#pragma once

#include <bit>
#include <execution>
#include <tuple>
#include <vector>

#include <omp.h>

#include <fp/rk61.hpp>
#include <lz77/factor.hpp>

#include <ankerl/unordered_dense.h>

#include "internal/benchmark.hpp"
#include "internal/io/memory_input_stream.hpp"
#include "internal/io/overlapping_blocks.hpp"
#include "internal/io/vbyte_coding.hpp"
#include "internal/sketch/bloom_filter.hpp"
#include "internal/util/concurrent_map.hpp"
#include "internal/util/idiv_ceil.hpp"
#include "internal/util/si_iec_literals.hpp"

namespace zk {

class PSampLZ_FGGK15 {
private:
    static constexpr std::string MAGIC = "psamplz";

public:
    template<iopp::STLInputStreamLike InputStream>
    static std::string decompress(InputStream& in) {
        std::string s(MAGIC.length(), 0);

        in.read(s.data(), MAGIC.length());
        if(s != MAGIC) {
            return "";
        }
        s.clear();

        size_t i = 0;
        while(in.good()) {
            auto const num_literals = internal::decode_vbyte(in);
            if(!in.good()) break;

            for(size_t i = 0; i < num_literals; i++) {
                s.push_back(in.get());
            }
            i += num_literals;

            if(in.good()) {
                auto const lxp = internal::decode_vbyte(in);
                if(!in.good()) break;

                auto const src = internal::decode_vbyte(in);

                size_t const len = 1ULL << lxp;
                for(size_t i = 0; i < len; i++) {
                    s.push_back(s[src + i]);
                }
                i += len;
            }
        }

        return s;
    }

private:
    using RK = fp::RabinKarp61;
    using Fingerprint = RK::Fingerprint;
    using ConcurrentSampling = internal::ConcurrentMap<Fingerprint, size_t>;

    static constexpr Fingerprint rolling_fp_base_ = (1ULL << 16) - 39;

    struct LZRef2At {
        size_t src; // the (absolute) position to refer to (src < pos)
        size_t lxp; // the reference length exponent (len = 2^lxp)
        size_t pos; // the position at which to insert the reference

        size_t len() const { return 1ULL << lxp; }
        size_t end() const { return pos + len(); }
    };

    using BloomFilter = internal::BloomFilter<Fingerprint, 2>;

    size_t len_exp_min_;
    size_t len_exp_max_;
    size_t sampling_;
    size_t window_size_;
    size_t bloom_filter_scale_;
    
    internal::Result result_;

    template<typename InputStream>
    std::vector<LZRef2At> parse(InputStream& in, size_t const n) {
        if(n - 1 > UINT32_MAX) {
            return parse_impl<InputStream, uint64_t>(in, n);
        } else {
            return parse_impl<InputStream, uint32_t>(in, n);
        }
    }

    template<typename InputStream, std::unsigned_integral Index>
    std::vector<LZRef2At> parse_impl(InputStream& in, size_t const n) {
        result_ = {};
        internal::MemoryTimePhase overall, phase;

        size_t const window_size = std::min(window_size_, n);
        size_t const bloom_filter_size = bloom_filter_scale_ * (n / (1ULL << sampling_));
        size_t const num_threads = omp_get_max_threads();

        size_t const len_max = 1ULL << len_exp_max_;
        size_t const num_lens = len_exp_max_ - len_exp_min_ + 1;
        auto get_len = [&](size_t const l){ return 1ULL << (len_exp_min_ + l); };

        using BlockIndex = uint32_t;
        {
            size_t const max_num_blocks = internal::idiv_ceil(n, get_len(num_lens));
            if(max_num_blocks > UINT32_MAX) {
                std::cerr << "too many blocks" << std::endl;
                std::abort;
            }
        }

        overall.start();
        std::vector<LZRef2At> refs;
        {
            for(size_t l1 = num_lens; l1 > 0; l1--) {
                auto const l = l1 - 1;
                auto const len = get_len(l);

                std::cerr << "processing l=" << (len_exp_min_ + l) << " (len=" << len << ") ..." << std::endl;
                RK rk(rolling_fp_base_, len);
                
                // sample
                std::cerr << "\tsample ... "; std::cerr.flush();
                phase.start();

                // allocate block fingerprints
                size_t const num_blocks = internal::idiv_ceil(n, len);
                auto block_fp = std::make_unique<Fingerprint[]>(num_blocks);
                auto const bufsize = (window_size / len) * len;
                auto buffer = std::make_unique<char[]>(bufsize);

                // compute block fingerprints in parallel
                in.seekg(0, std::ios_base::beg);
                size_t block_offs = 0;
                while(in) {
                    in.read(buffer.get(), bufsize);
                    auto const read = in.gcount();
                    if(!read) break;

                    auto const num_blocks_read = read / len;

                    #pragma omp parallel for
                    for(size_t i = 0; i < num_blocks_read; i++) {
                        auto const block_beg = i * len;
                        auto const block_end = block_beg + len;

                        Fingerprint fp = 0;
                        for(size_t j = block_beg; j < block_end; j++) {
                            fp = rk.push(fp, buffer[j]);
                        }

                        block_fp[block_offs + i] = fp;
                    }

                    block_offs += num_blocks_read;
                }

                phase.stop();
                std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms" << std::endl;
                result_.add("t_sample_" + std::to_string(len_exp_min_ + l), phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
                
                // map
                ankerl::unordered_dense::map<Fingerprint, BlockIndex> fp_head;
                ankerl::unordered_dense::map<BlockIndex, BlockIndex> fp_next;

                std::cerr << "\tmap ... "; std::cerr.flush();
                phase.start();   
                for(size_t i = 0; i < num_blocks; i++) {
                    auto const fp = block_fp[i];
                    auto it = fp_head.find(fp);
                    if(it != fp_head.end()) {
                        auto const j = it->second;
                        fp_next[i] = j;
                        it->second = i;
                    } else {
                        fp_head.emplace(fp, i);
                    }
                }
                phase.stop();
                std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms (" << fp_head.size() << " distinct blocks out of " << num_blocks << " total)" << std::endl;
                result_.add("t_map_" + std::to_string(len_exp_min_ + l), phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());

                // parse
                std::cerr << "\tparse ... "; std::cerr.flush();
                phase.start();

                in.seekg(0, std::ios_base::beg);
                internal::OverlappingBlocks<InputStream> block(window_size, len);

                auto init_fingerprinting = [&](size_t const thread_num, size_t const num_per_thread){
                    // compute window part boundaries
                    char const* beg = block.begin() + thread_num * num_per_thread;
                    char const* end = std::min(beg + num_per_thread, block.end());

                    // initialize fingerprint with what comes right before the thread's block
                    Fingerprint fp = 0;
                    for(char const* p = beg - len; p < beg; p++) {
                        fp = rk.push(fp, *p);
                    }
                    return std::make_tuple(beg, end, fp);
                };

                auto roll_fingerprint = [&](Fingerprint& fp, char const* p){
                    fp = rk.roll(fp, *(p-len), *p);
                };

                struct BlockRef {
                    BlockIndex block;
                    Index src;
                } __attribute__((packed));

                std::unique_ptr<std::vector<BlockRef>> lrefs[num_threads];
                for(size_t x = 0; x < num_threads; x++) {
                    lrefs[x] = std::make_unique<std::vector<BlockRef>>();
                }

                block.init(in);
                do {
                    if(block.empty()) continue;

                    auto const num_per_thread = internal::idiv_ceil(block.size(), num_threads);
                    #pragma omp parallel
                    {
                        size_t const thread_num = omp_get_thread_num();
                        auto [beg, end, fp] = init_fingerprinting(thread_num, num_per_thread);
                        auto& local_block_refs = *lrefs[thread_num];

                        ssize_t i = block.offset() + (beg - block.begin()) - ssize_t(len) + 1;
                        char const* p = beg;
                        for(; p < beg + len; p++, i++) roll_fingerprint(fp, p);
                        for(; p < end; p++, i++) {
                            // compute new fingerprint
                            roll_fingerprint(fp, p);

                            // check if fingerprint matches any blocks
                            auto it = fp_head.find(fp);
                            if(it != fp_head.end()) {
                                // push potential refs
                                // TODO: check if even possible (i < block position)
                                auto block = it->second;
                                while(block != -1) {
                                    local_block_refs.push_back(BlockRef{block, Index(i)});

                                    auto it2 = fp_next.find(block);
                                    block = (it2 != fp_next.end()) ? it2->second : -1;
                                }
                            }
                        }
                    }

                    // advance
                    block.advance();
                } while(!block.last());

                phase.stop();
                std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms" << std::endl;
                result_.add("t_parse_" + std::to_string(len_exp_min_ + l), phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());

                // make refs
                std::cerr << "\tmake refs ... "; std::cout.flush();
                phase.start();

                auto block_src = std::make_unique<size_t[]>(num_blocks);
                for(size_t i = 0; i < num_blocks; i++) {
                    block_src[i] = SIZE_MAX;
                }

                for(size_t x = 0; x < num_threads; x++) {
                    for(auto& block_ref : *lrefs[x]) {
                        auto const i = block_ref.block;
                        size_t const pos = i * len;
                        auto const src = block_ref.src;
                        if(src < pos && src < block_src[i]) {
                            block_src[i] = src;
                        }
                    }
                }

                size_t num_refs = 0;
                size_t next_ref = 0; // the next relevant already existing reference
                size_t next_possible_ref_pos = 0;

                for(size_t i = 0; i < num_blocks; i++) {
                    size_t const block_pos = i * len;

                    // possibly advance next reference
                    while(next_ref < refs.size() && block_pos >= refs[next_ref].end()) {
                        next_possible_ref_pos = refs[next_ref].end();
                        ++next_ref;
                    }

                    if(block_pos > next_possible_ref_pos && (next_ref >= refs.size() || block_pos + len <= refs[next_ref].pos) && block_src[i] != SIZE_MAX) {
                        refs.push_back(LZRef2At{block_src[i], len_exp_min_ + l, block_pos});
                        next_possible_ref_pos = block_pos + len;
                        ++num_refs;
                    }
                }

                phase.stop();
                std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms (new refs: " << num_refs << ")" << std::endl;
                result_.add("t_mkrefs_" + std::to_string(len_exp_min_ + l), phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());

                // sort refs
                std::cout << "\tsort " << refs.size() << " refs ... "; std::cout.flush();
                phase.start();
                std::sort(std::execution::par_unseq, refs.begin(), refs.end(), [](LZRef2At const& a, LZRef2At const& b){ return a.pos <= b.pos; });
                phase.stop();
                std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms" << std::endl;
                result_.add("t_sort_" + std::to_string(len_exp_min_ + l), phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
            }
        }

        overall.stop();
        result_.add("p", num_threads);
        result_.add("sampling", sampling_);
        result_.add("window", window_size);
        result_.add("refs", refs.size());
        result_.add("t", overall.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());

        size_t const mem = overall.get_metric<pm::MallocCounter::MemoryPeakMetric>();
        double const exp_num_fingerprints = double(n) / double(1ULL << sampling_);
        result_.add("mem_peak", overall.get_metric<pm::MallocCounter::MemoryPeakMetric>());
        return refs;
    }

public:
    PSampLZ_FGGK15(size_t len_exp_min, size_t len_exp_max, size_t sampling, size_t bloom_filter_scale = 6, size_t window_size = 64_Mi)
        : len_exp_min_(len_exp_min),
          len_exp_max_(len_exp_max),
          sampling_(sampling),
          window_size_(window_size),
          bloom_filter_scale_(bloom_filter_scale) {
    }

    PSampLZ_FGGK15(PSampLZ_FGGK15&&) = default;
    PSampLZ_FGGK15& operator=(PSampLZ_FGGK15&&) = default;
    PSampLZ_FGGK15(PSampLZ_FGGK15 const&) = default;
    PSampLZ_FGGK15& operator=(PSampLZ_FGGK15 const&) = default;

    template<std::contiguous_iterator Input, std::output_iterator<lz77::Factor> Output>
    requires (sizeof(std::iter_value_t<Input>) == 1)
    void factorize(Input begin, Input const& end, Output out) {
        // parse
        std::string_view s(begin, end);
        size_t const n = s.size();
        internal::MemoryInputStream in(s.data(), n);
        auto refs = parse(in, n);

        // factorize
        size_t i = 0;
        for(auto it = refs.begin(); i < n && it != refs.end(); it++) {
            // create literal factors up to the next reference
            for(; i < it->pos; i++) {
                *out++ = lz77::Factor(s[i]);
            }

            // create reference
            *out++ = lz77::Factor(i - it->src, it->len());
            i += it->len();
        }

        // create literal factors for the remaining literals
        for(; i < n; i++) {
            *out++ = lz77::Factor(s[i]);
        }
    }

    template<iopp::STLInputStreamLike InputStream, iopp::STLOutputStreamLike OutputStream>
    requires requires(InputStream& subject, size_t const offs, std::ios_base::seekdir const dir){
        { subject.seekg(offs, dir) };
    }
    void compress(InputStream& in, size_t const n, OutputStream& out) {
        // parse
        auto refs = parse(in, n);

        // encode
        std::cerr << "encode ... "; std::cout.flush();
        internal::TimePhase phase;
        phase.start();

        size_t z = 0;
        size_t nout = 0;
        
        in.seekg(0, std::ios_base::beg);
        {
            // magic
            out.write(MAGIC.data(), MAGIC.length());
            nout += MAGIC.length();

            size_t i = 0;

            auto copy_literals = [&](size_t const num_literals){
                nout += internal::encode_vbyte(out, num_literals);
                
                for(size_t j = 0; j < num_literals; j++) {
                    out.put(in.get());
                }

                z += num_literals;
                nout += num_literals;
                i += num_literals;
            };
            
            for(auto it = refs.begin(); it != refs.end(); it++) {
                // copy literals up to next reference
                copy_literals(it->pos - i);

                // encode reference
                ++z;
                nout += internal::encode_vbyte(out, it->lxp);
                nout += internal::encode_vbyte(out, it->src);
                
                // skip literals replaced by reference
                i += it->len();
                in.seekg(it->len(), std::ios_base::cur);
            }

            // copy remaining literals
            if(i < n) {
                copy_literals(n - i);
            }
        }

        phase.stop();
        std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms" << std::endl;
        result_.add("t_encode", phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());

        // result
        result_.add("z", z);
        result_.add("nout", nout);
    }

    auto&& consume_last_result() {
        return std::move(result_);
    }
};

}
