#pragma once

#include <lz77/lpf_factorizer.hpp>

#include "sampled_lpf_factorizer.hpp"

#include "internal/benchmark.hpp"
#include "internal/io/block_coding.hpp"
#include "internal/sketch/topk_prefixes_misra_gries.hpp"

namespace zk {

class TopkLZ77 {
private:
    static constexpr uint64_t MAGIC =
        ((uint64_t)'T') << 56 |
        ((uint64_t)'O') << 48 |
        ((uint64_t)'P') << 40 |
        ((uint64_t)'K') << 32 |
        ((uint64_t)'F') << 24 |
        ((uint64_t)'A') << 16 |
        ((uint64_t)'C') << 8 |
        ((uint64_t)'T');

    static constexpr internal::TokenType TOK_TRIE_REF = 0;
    static constexpr internal::TokenType TOK_FACT_SRC = 1;
    static constexpr internal::TokenType TOK_FACT_LEN = 2;
    static constexpr internal::TokenType TOK_LITERAL = 3;
    static constexpr internal::TokenType TOK_FACT_REMAINDER = 4;

    static void setup_encoding(internal::BlockEncodingBase& enc, size_t const k, size_t const window_size) {
        enc.register_binary(k-1);                // TOK_TRIE_REF
        enc.register_binary(window_size-1);      // TOK_FACT_SRC
        enc.register_huffman();                  // TOK_FACT_LEN
        enc.register_binary(255, false);         // TOK_FACT_LITERAL
        enc.register_binary(window_size, false); // TOK_FACT_REMAINDER
    }
    
    static constexpr size_t MAX_LZ_REF_LEN = 255;
    static constexpr size_t MAX_WINDOW_SIZE = 2'147'483'647;

    using Index = uint32_t;
    using Node = Index;
    using Topk = internal::TopKPrefixesMisraGries<>;

    size_t k_;
    size_t window_size_;
    size_t max_freq_;
    size_t lz_sampling_;
    size_t block_size_;

    internal::Result result_;

    struct LZRef {
        Index pos, src, len;
        Index num_literals() const { return len; }
        Index end() const { return pos + len - 1; }
    } __attribute__((packed));

    enum WhatToDoNext { TRIE_REF, LZ_REF, LITERAL };

public:
    TopkLZ77(size_t const k, size_t const window_size, size_t const max_freq, size_t const lz_sampling, size_t const block_size)
        : k_(k),
          window_size_(std::min(window_size, MAX_WINDOW_SIZE)),
          max_freq_(max_freq),
          lz_sampling_(lz_sampling),
          block_size_(block_size) {
        
        if(window_size > MAX_WINDOW_SIZE) {
            std::cerr << "window size clamped to maximum value of " << MAX_WINDOW_SIZE << std::endl;
        }
    }

    TopkLZ77(TopkLZ77&&) = default;
    TopkLZ77& operator=(TopkLZ77&&) = default;
    TopkLZ77(TopkLZ77 const&) = default;
    TopkLZ77& operator=(TopkLZ77 const&) = default;

    template<iopp::STLInputStreamLike InputStream, iopp::BitSink Out>
    void compress(InputStream& in, Out out) {
        // init stats
        internal::MemoryTimePhase stats("topk-lz77");
        internal::TimePhase phase_read;
        internal::TimePhase phase_block_lz77;
        internal::TimePhase phase_process;
        result_ = {};

        size_t num_lz = 0;
        size_t num_trie = 0;
        size_t num_literal = 0;
        size_t trie_longest = 0;
        size_t lz_longest = 0;
        size_t total_trie_len = 0;
        size_t total_lz_len = 0;
        size_t num_relevant = 0;

        stats.start();

        // write header and initialize encoding
        out.write(MAGIC, 64);
        out.write(k_, 64);
        out.write(window_size_, 64);
        out.write(max_freq_, 64);

        internal::BlockEncoder enc(out, block_size_);
        setup_encoding(enc, k_, window_size_);

        // initialize top-k
        Topk topk(k_ - 1, max_freq_);

        // initialize factorizer
        std::vector<LZRef> lz_refs;

        // initialize buffers
        size_t block_offs = 0;
        auto block = std::make_unique<char[]>(window_size_);

        while(in) {
            // read next block
            Index block_num;
            {
                phase_read.resume();
                in.read(block.get(), window_size_);
                block_num = in.gcount();
                phase_read.pause();
            }

            if(block_num == 0) {
                continue;
            }

            // compute the LZ77 factorization of the block
            {
                phase_block_lz77.resume();
                lz_refs.clear();
                {
                    size_t i = 0;
                    auto emit_literal = [&](lz77::Factor){
                        i++;
                    };
                    auto emit_ref = [&](lz77::Factor f){
                        lz_refs.emplace_back(i, f.src, f.len);
                        i += f.num_literals();
                    };

                    if(lz_sampling_ > 0) {
                        SampledLPFFactorizer lpf(lz_sampling_, 16);
                        lpf.factorize(block.get(), block.get() + block_num, emit_literal, emit_ref);
                    } else {
                        lz77::LPFFactorizer lpf;
                        lpf.factorize(block.get(), block.get() + block_num, emit_literal, emit_ref);
                    }
                }
                phase_block_lz77.pause();
            }

            // parse and encode the block
            // at the beginning of each LZ77 factor, we attempt to find the longest possible string back in the top-k trie
            // if we find a string longer than the next LZ77 factor, we encode it using a trie reference and advance in the LZ77 factorization, potentially chopping
            // the factor that we reach into two fractions
            {
                phase_process.resume();

                auto topk_enter = [&](size_t const pos, size_t const len){
                    ++num_relevant;

                    typename Topk::StringState s = topk.empty_string();
                    Node node;
                    while(s.frequent && s.len < len && pos + s.len < block_num) {
                        node = s.node;
                        s = topk.extend(s, block[pos + s.len]);
                    }
                };

                Index z = 0; // the current or next LZ reference
                Index curpos = 0;
                while(curpos < block_num) {
                    auto const gpos = block_offs + curpos;

                    // find the longest string represented in the top-k trie starting at the current position
                    Node v;
                    Index dv = topk.find(block.get() + curpos, block_num - curpos, v);

                    // examine the current or next LZ reference
                    auto const& lz_ref = lz_refs[z];

                    // decide what to do
                    WhatToDoNext decision;

                    if(z >= lz_refs.size() || curpos < lz_ref.pos) {
                        // we have no current LZ reference that we can use
                        decision = (dv >= 1) ? TRIE_REF : LITERAL;
                    } else {
                        if(dv >= lz_ref.num_literals()) {
                            // trie reference is at least as good as LZ reference -- prefer it since the encoding is smaller
                            decision = TRIE_REF;
                        } else {
                            // LZ reference is longer than trie reference
                            decision = lz_ref.num_literals() > 1 ? LZ_REF : LITERAL;
                        }
                    }

                    if(decision == TRIE_REF) {
                        // encode a top-k trie reference
                        assert(v > 0);

                        enc.write_uint(TOK_FACT_LEN, 0);
                        enc.write_uint(TOK_TRIE_REF, v);
                        
                        ++num_trie;
                        trie_longest = std::max(trie_longest, (size_t)dv);
                        total_trie_len += dv;

                        // advance in LZ77 references
                        curpos += dv;
                        if(curpos < block_num) {
                            while(z < lz_refs.size() && curpos > lz_refs[z].end()) {
                                ++z;
                            }

                            if(z < lz_refs.size() && curpos > lz_refs[z].pos) {
                                // chop current LZ77 factor
                                auto const chop = curpos - lz_refs[z].pos;
                                lz_refs[z].pos += chop;
                                lz_refs[z].len -= chop;
                            }
                        }
                    } else if(decision == LZ_REF) {
                        // a real LZ77 reference
                        auto const fpos = curpos;
                        assert(fpos >= lz_ref.src);

                        if(lz_ref.len >= MAX_LZ_REF_LEN) {
                            // encode the maximum length, then encode the rest as a special token
                            enc.write_uint(TOK_FACT_LEN, MAX_LZ_REF_LEN);
                            enc.write_uint(TOK_FACT_REMAINDER, lz_ref.len - MAX_LZ_REF_LEN);
                        } else {
                            // simply encode the length
                            enc.write_uint(TOK_FACT_LEN, lz_ref.len);
                        }

                        // write source
                        enc.write_uint(TOK_FACT_SRC, lz_ref.src);

                        ++num_lz;
                        lz_longest = std::max(lz_longest, (size_t)lz_ref.len);
                        total_lz_len += lz_ref.len;

                        // enter
                        topk_enter(curpos, lz_ref.len);
                        curpos += lz_ref.len;

                        // advance to next LZ reference
                        ++z;
                    } else {
                        // encode a literal
                        enc.write_uint(TOK_FACT_LEN, 1);
                        enc.write_char(TOK_LITERAL, block[curpos]);

                        ++num_literal;

                        topk_enter(curpos, 1);
                        ++curpos;

                        if(curpos > lz_ref.end()) {
                            ++z;
                        }
                    }
                }
                phase_process.pause();
            }

            block_offs += block_num;
        }

        phase_process.resume();
        enc.flush();
        phase_process.pause();

        stats.stop();

        result_.add("z", num_lz + num_literal + num_trie);
        result_.add("phrases_ref", num_lz + num_trie);
        result_.add("phrases_literal", num_literal);
        result_.add("num_relevant", num_relevant);
        result_.add("phrases_longest", std::max(trie_longest, lz_longest));
        result_.add("phrases_longest_lz", lz_longest);
        result_.add("phrases_longest_trie", trie_longest);
        result_.add("phrases_ref_lz", num_lz);
        result_.add("phrases_ref_trie", num_trie);
        result_.add("phrases_avg_ref_len", std::round(100.0 * ((double)(total_lz_len + total_trie_len) / (double)(num_lz + num_trie))) / 100.0);
        result_.add("phrases_avg_ref_len_lz", std::round(100.0 * ((double)total_lz_len / (double)num_lz)) / 100.0);
        result_.add("phrases_avg_ref_len_trie", std::round(100.0 * ((double)total_trie_len / (double)num_trie)) / 100.0);
        result_.add("t", stats.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("t_read", phase_read.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("t_lz77", phase_block_lz77.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("t_process", phase_process.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("mem_peak", stats.get_metric<pm::MallocCounter::MemoryPeakMetric>());
    }

    template<iopp::BitSource In, iopp::STLOutputStreamLike OutputStream>
    static void decompress(In in, OutputStream& out) {
        // decode header
        uint64_t const magic = in.read(64);
        if(magic != MAGIC) {
            std::cerr << "wrong magic: 0x" << std::hex << magic << " (expected: 0x" << MAGIC << ")" << std::endl;
            std::abort();
        }

        auto const k = in.read(64);
        auto const window_size = in.read(64);
        auto const max_freq = in.read(64);

        // initialize decoding
        internal::BlockDecoder dec(in);
        setup_encoding(dec, k, window_size);
        Topk topk(k - 1, max_freq);
        
        auto block = std::make_unique<char[]>(window_size);
        auto block_offs = 0;
        size_t curpos = 0;

        while(in) {
            if(block_offs == 10485760) {
                block_offs = 10485760;
            }

            auto const len = dec.read_uint(TOK_FACT_LEN);
            size_t phrase_len;

            auto const gpos = block_offs + curpos;
            if(len == 0) {
                // a top-k trie reference
                auto const node = dec.read_uint(TOK_TRIE_REF);
                phrase_len = topk.get(node, block.get() + curpos);
            } else if(len == 1) {
                // a literal character
                auto const c = dec.read_char(TOK_LITERAL);
                block[curpos] = c;
                phrase_len = 1;
            } else {
                phrase_len = len;

                // a block-local LZ77 reference
                if(len == MAX_LZ_REF_LEN) {
                    // this factor may be even longer, decode remainder
                    phrase_len += dec.read_uint(TOK_FACT_REMAINDER);
                }

                auto const src = dec.read_uint(TOK_FACT_SRC);
                assert(curpos >= src);
                auto const srcpos = curpos - src;
                for(size_t i = 0; i < phrase_len; i++) {
                    block[curpos + i] = block[srcpos + i];
                }
            }

            // enter string into top-k structure
            {
                typename Topk::StringState s = topk.empty_string();
                Node node;
                while(s.frequent && s.len < phrase_len) {
                    assert(curpos + s.len < window_size);
                    node = s.node;
                    s = topk.extend(s, block[curpos + s.len]);
                }
            }

            // advance
            curpos += phrase_len;

            if(curpos >= window_size) {
                // emit and advance to new block
                out.write(block.get(), window_size);
                curpos = 0;
                block_offs += window_size;
            }
        }

        // emit final block
        out.write(block.get(), curpos);
    }

    auto&& consume_last_result() {
        return std::move(result_);
    }
};

}
