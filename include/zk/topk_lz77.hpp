#pragma once

#include <lz77/lpf_factorizer.hpp>

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

    using Index = uint32_t;
    using Node = Index;
    using Topk = internal::TopKPrefixesMisraGries<>;

    size_t threshold_;
    size_t k_;
    size_t window_size_;
    size_t max_freq_;
    size_t block_size_;

    internal::Result result_;

public:
    TopkLZ77(size_t const threshold, size_t const k, size_t const window_size, size_t const max_freq, size_t const block_size)
        : threshold_(threshold),
          k_(k),
          window_size_(window_size),
          max_freq_(max_freq),
          block_size_(block_size) {
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
        lz77::LPFFactorizer lpf;
        lpf.min_reference_length(threshold_);
        std::vector<lz77::Factor> factors;

        // initialize buffers
        auto block_offs = 0;
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
                factors.clear();
                lpf.factorize(block.get(), block.get() + block_num, std::back_inserter(factors));
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

                Index z = 0; // the current LZ77 factor
                Index curpos = 0;
                while(curpos < block_num) {
                    auto const gpos = block_offs + curpos;

                    // find the longest string represented in the top-k trie starting at the current position
                    Node v;
                    Index dv = topk.find(block.get() + curpos, block_num - curpos, v);

                    auto const& f = factors[z];
                    if(dv >= f.num_literals()) {
                        // encode a top-k trie reference
                        assert(v > 0);

                        enc.write_uint(TOK_FACT_LEN, 0);
                        enc.write_uint(TOK_TRIE_REF, v);
                        
                        ++num_trie;
                        trie_longest = std::max(trie_longest, (size_t)dv);
                        total_trie_len += dv;

                        // advance in LZ77 factorization
                        if(curpos + dv < block_num) {
                            auto d = dv;
                            while(d >= factors[z].num_literals()) {
                                d -= factors[z].num_literals();
                                ++z;
                            }

                            if(d > 0) {
                                // chop current LZ77 factor
                                assert(factors[z].len > d);
                                factors[z].len -= d;
                            }
                        }

                        // enter
                        topk_enter(curpos, dv);
                        curpos += dv;
                    } else {
                        // encode a LZ77 reference or a literal
                        if(f.is_literal() || f.num_literals() == 1) {
                            // a literal factor (possibly a reference of length one introduced due to chopping)
                            enc.write_uint(TOK_FACT_LEN, 1);
                            enc.write_char(TOK_LITERAL, block[curpos]);

                            ++num_literal;

                            topk_enter(curpos, 1);
                            ++curpos;
                        } else {
                            // a real LZ77 reference
                            auto const fpos = curpos;
                            assert(fpos >= f.src);

                            if(f.len >= MAX_LZ_REF_LEN) {
                                // encode the maximum length, then encode the rest as a special token
                                enc.write_uint(TOK_FACT_LEN, MAX_LZ_REF_LEN);
                                enc.write_uint(TOK_FACT_REMAINDER, f.len - MAX_LZ_REF_LEN);
                            } else {
                                // simply encode the length
                                enc.write_uint(TOK_FACT_LEN, f.len);
                            }

                            // write source
                            enc.write_uint(TOK_FACT_SRC, f.src);

                            ++num_lz;
                            lz_longest = std::max(lz_longest, (size_t)f.len);
                            total_lz_len += f.len;

                            // enter
                            topk_enter(curpos, f.len);
                            curpos += f.len;
                        }

                        // advance to next LZ77 factor
                        ++z;
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
