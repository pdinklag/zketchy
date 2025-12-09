#pragma once

#include <iopp/stream_input_iterator.hpp>
#include <iopp/stream_output_iterator.hpp>

#include "internal/benchmark.hpp"
#include "internal/io/block_coding.hpp"
#include "internal/sketch/topk_prefixes_misra_gries.hpp"

namespace zk {

class TopkLZ78 {
private:
    static constexpr uint64_t MAGIC =
        ((uint64_t)'T') << 56 |
        ((uint64_t)'O') << 48 |
        ((uint64_t)'P') << 40 |
        ((uint64_t)'K') << 32 |
        ((uint64_t)'L') << 24 |
        ((uint64_t)'Z') << 16 |
        ((uint64_t)'7') << 8 |
        ((uint64_t)'8');

    static constexpr internal::TokenType TOK_REF = 0;
    static constexpr internal::TokenType TOK_LITERAL = 1;

    static void setup_encoding(internal::BlockEncodingBase& enc, size_t const k) {
        enc.register_binary(k-1); // TOK_REF
        enc.register_huffman();   // TOK_LITERAL
    }
    
    static constexpr size_t MAX_LZ_REF_LEN = 255;
    static constexpr size_t MAX_WINDOW_SIZE = 2'147'483'647;

    using Index = uint32_t;
    using Node = Index;
    using Topk = internal::TopKPrefixesMisraGries<>;

    size_t k_;
    size_t max_freq_;
    size_t block_size_;

    internal::Result result_;

    struct LZRef {
        Index pos, src, len;
        Index end() const { return pos + len - 1; }
    } __attribute__((packed));

    enum WhatToDoNext { TRIE_REF, LZ_REF, LITERAL };

public:
    TopkLZ78(size_t const k, size_t const max_freq, size_t const block_size)
        : k_(k),
          max_freq_(max_freq),
          block_size_(block_size) {
    }

    TopkLZ78(TopkLZ78&&) = default;
    TopkLZ78& operator=(TopkLZ78&&) = default;
    TopkLZ78(TopkLZ78 const&) = default;
    TopkLZ78& operator=(TopkLZ78 const&) = default;

    template<iopp::STLInputStreamLike InputStream, iopp::BitSink Out>
    void compress(InputStream& in, Out out) {
        // init stats
        internal::MemoryTimePhase stats("topk-lz78");
        size_t longest = 0;
        size_t total_len = 0;
        result_ = {};

        stats.start();

        // write header and initialize encoding
        out.write(MAGIC, 64);
        out.write(k_, 64);
        out.write(max_freq_, 64);

        internal::BlockEncoder enc(out, block_size_);
        setup_encoding(enc, k_);

        // initialize top-k
        Topk topk(k_ + 1, max_freq_);

        // factorize
        using It = iopp::StreamInputIterator<InputStream>;
        auto const end = It::end(in);

        auto s = topk.empty_string();
        size_t z = 0;
        for(auto it = It(in); it != end; it++) {
            auto const c = *it;
            auto next = topk.extend(s, c);
            if(!next.frequent) {
                enc.write_uint(TOK_REF, s.node);
                enc.write_char(TOK_LITERAL, c);

                longest = std::max(longest, size_t(next.len));
                total_len += next.len;

                s = topk.empty_string();
                ++z;
            } else {
                s = next;
            }
        }

        // encode final phrase
        if(s.len > 0) {
            enc.write_uint(TOK_REF, s.node);
            ++z;
        }

        enc.flush();
        
        stats.stop();

        result_.add("z", z);
        result_.add("len_max", longest);
        result_.add("len_avg", std::round(100.0 * ((double)total_len / (double)z)) / 100.0);
        result_.add("t", stats.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("mem_peak", stats.get_metric<pm::MallocCounter::MemoryPeakMetric>());
    }

    template<iopp::BitSource In, iopp::STLOutputStreamLike OutputStream>
    static void decompress(In in, OutputStream& outs) {
        iopp::StreamOutputIterator out(outs);

        // decode header
        uint64_t const magic = in.read(64);
        if(magic != MAGIC) {
            std::cerr << "wrong magic: 0x" << std::hex << magic << " (expected: 0x" << MAGIC << ")" << std::endl;
            std::abort();
        }

        auto const k = in.read(64);
        auto const max_freq = in.read(64);

        // initialize decompression
        // - frequent substring 0 is reserved to indicate a literal character
        Topk topk(k + 1, max_freq);

        size_t n = 0;
        size_t z = 0;

        // initialize decoding
        internal::BlockDecoder dec(in);
        setup_encoding(dec, k);

        char* phrase = new char[k]; // phrases can be of length up to k...
        while(in) {
            // decode and handle phrase
            auto const x = dec.read_uint(TOK_REF);
            auto const phrase_len = topk.get(x, phrase);
            
            ++z;
            n += phrase_len;

            auto s = topk.empty_string();
            for(size_t i = 0; i < phrase_len; i++) {
                auto const c = phrase[i];
                s = topk.extend(s, c);
                *out++ = c;
            }

            // decode and handle literal
            if(in)
            {
                auto const literal = dec.read_char(TOK_LITERAL);
                topk.extend(s, literal);
                *out++ = literal;
                ++n;
            }
        }
    }

    auto&& consume_last_result() {
        return std::move(result_);
    }
};

}
