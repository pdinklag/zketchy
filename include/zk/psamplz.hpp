#pragma once

#include <bit>
#include <execution>
#include <tuple>
#include <vector>

#include <omp.h>

#include <fp/rk61.hpp>
#include <iopp/util/overlapping_blocks.hpp>
#include <lz77/factor.hpp>

#include "internal/benchmark.hpp"
#include "internal/io/memory_input_stream.hpp"
#include "internal/io/vbyte_coding.hpp"
#include "internal/sketch/bloom_filter.hpp"
#include "internal/util/concurrent_map.hpp"
#include "internal/util/idiv_ceil.hpp"
#include "internal/util/si_iec_literals.hpp"

namespace zk {

class PSampLZ {
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

    struct UpdateLeftmost {
        using mapped_type = size_t;

        mapped_type operator()(mapped_type& lhs, const mapped_type& rhs) const { return lhs = std::min(lhs, rhs); }
    };

    static constexpr Fingerprint rolling_fp_base_ = (1ULL << 16) - 39;

    struct CHTLookupStats {
        size_t bloom_lookups;
        size_t map_lookups;
        size_t false_positives;

        CHTLookupStats() : bloom_lookups(0), map_lookups(0), false_positives(0) {
        }

        CHTLookupStats(CHTLookupStats const&) = default;
        CHTLookupStats& operator=(CHTLookupStats const&) = default;

        CHTLookupStats& operator+=(CHTLookupStats const& x) {
            bloom_lookups += x.bloom_lookups;
            map_lookups += x.map_lookups;
            false_positives += x.false_positives;
            return *this;
        }

        CHTLookupStats operator+(CHTLookupStats const& x) {
            CHTLookupStats sum = *this;
            sum += x;
            return sum;
        }
    };

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
    internal::MemoryTimePhase stats_;

    template<typename InputStream>
    std::vector<LZRef2At> parse(InputStream& in, size_t const n) {
        result_ = {};
        CHTLookupStats cht_total_stats;

        stats_ = internal::MemoryTimePhase("psamplz");
        stats_.start();

        size_t const window_size = std::min(window_size_, n);
        size_t const bloom_filter_size = bloom_filter_scale_ * (n / (1ULL << sampling_));
        size_t const num_threads = omp_get_max_threads();

        size_t const len_max = 1ULL << len_exp_max_;
        size_t const num_lens = len_exp_max_ - len_exp_min_ + 1;
        auto get_len = [&](size_t const l){ return 1ULL << (len_exp_min_ + l); };

        // sample and parse for each length
        std::vector<LZRef2At> refs;
        {
            internal::TimePhase phase_sample_and_parse("sample and parse");
            phase_sample_and_parse.start();
            iopp::OverlappingBlocks<InputStream> block(window_size, len_max);

            for(size_t l1 = num_lens; l1 > 0; l1--) {
                auto const l = l1 - 1;
                auto const len = get_len(l);

                internal::TimePhase phase_process_len("len=" + std::to_string(len));
                phase_process_len.data()["len"] = len;
                phase_process_len.start();

                // initialize fingerprinting
                RK rk(rolling_fp_base_, len);

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

                // sample
                internal::TimePhase phase_sample("sample");
                phase_sample.start();
                
                ConcurrentSampling growt_map(n / (1ULL << sampling_));

                std::unique_ptr<BloomFilter> lbloom[num_threads];
                for(size_t thread_num = 0; thread_num < num_threads; thread_num++) {
                    lbloom[thread_num] = std::make_unique<BloomFilter>(bloom_filter_size);
                }

                {
                    auto const s = (1ULL << sampling_) - 1;

                    in.seekg(0, std::ios_base::beg);
                    block.init(in);
                    do {
                        auto const num_per_thread = internal::idiv_ceil(block.size(), num_threads);

                        #pragma omp parallel
                        {
                            // process block
                            size_t const thread_num = omp_get_thread_num();
                            auto& bloom = *lbloom[thread_num];
                            auto growt_handle = growt_map.get_handle();

                            auto [beg, end, fp] = init_fingerprinting(thread_num, num_per_thread);
                            ssize_t i = block.offset() + (beg - block.begin()) - ssize_t(len) + 1;

                            char const* p = beg;
                            for(; i < 0; p++, i++) roll_fingerprint(fp, p);
                            for(; p < end; p++, i++) {
                                // compute new fingerprint
                                roll_fingerprint(fp, p);

                                // sample
                                if((fp & s) == 0) {
                                    growt_handle.insert_or_update(fp, i, UpdateLeftmost(), i);
                                    bloom.emplace(fp);
                                }
                            }
                        }
                    } while(block.advance());
                }

                phase_sample.data()["num_fingerprints_approx"] = growt_map.get_handle().element_count_approx();
                phase_sample.stop();
                phase_process_len.append_child(phase_sample);

                // merge bloom filters
                internal::TimePhase phase_merge_bloom("merge bloom filters");
                phase_merge_bloom.start();

                BloomFilter bloom(bloom_filter_size);
                for(size_t thread_num = 0; thread_num < num_threads; thread_num++) {
                    bloom |= *lbloom[thread_num];
                    lbloom[thread_num].reset();
                }

                phase_merge_bloom.stop();
                phase_process_len.append_child(phase_merge_bloom);

                // parse
                internal::TimePhase phase_parse("parse");
                phase_parse.start();

                std::unique_ptr<std::vector<LZRef2At>> lrefs[num_threads];

                std::unique_ptr<CHTLookupStats> lcht_stats[num_threads];
                for(size_t x = 0; x < num_threads; x++) {
                    lrefs[x] = std::make_unique<std::vector<LZRef2At>>();
                    lcht_stats[x] = std::make_unique<CHTLookupStats>();
                }

                {
                    in.seekg(0, std::ios_base::beg);
                    block.init(in);
                    do {
                        auto const num_per_thread = internal::idiv_ceil(block.size(), num_threads);
                        #pragma omp parallel
                        {
                            size_t const thread_num = omp_get_thread_num();
                            auto [beg, end, fp] = init_fingerprinting(thread_num, num_per_thread);
                            auto& new_lrefs = *lrefs[thread_num];
                            auto& cht_stats = *lcht_stats[thread_num];

                            // process portion
                            auto growt_handle = growt_map.get_handle();

                            size_t next_ref = 0; // the next relevant reference
                            size_t next_possible_ref_pos = 0;
                            // TODO: find initial value of next_ref via binary search 

                            ssize_t i = block.offset() + (beg - block.begin()) - ssize_t(len) + 1;

                            char const* p = beg;
                            for(; p < beg + len; p++, i++) roll_fingerprint(fp, p);
                            for(; p < end; p++, i++) {
                                // compute new fingerprint
                                roll_fingerprint(fp, p);

                                // possibly advance next reference
                                while(next_ref < refs.size() && i >= refs[next_ref].end()) {
                                    next_possible_ref_pos = refs[next_ref].end();
                                    ++next_ref;
                                }

                                // test whether we are good to introduce a reference here
                                if(i >= next_possible_ref_pos && (next_ref >= refs.size() || i + len <= refs[next_ref].pos))
                                {
                                    ++cht_stats.bloom_lookups;
                                    if(!bloom.lookup(fp)) continue;

                                    ++cht_stats.map_lookups;
                                    auto it = growt_handle.find(fp);
                                    if(it != growt_handle.end()) {
                                        auto const j = (*it).second;
                                        if(i > ssize_t(j)) {
                                            new_lrefs.push_back(LZRef2At{j, len_exp_min_ + l, size_t(i)});
                                            next_possible_ref_pos = i + len;
                                        }
                                    } else {
                                        ++cht_stats.false_positives;
                                    }
                                }
                            }
                        }
                    } while(block.advance());
                }

                phase_parse.stop();

                // process lookup stats
                {
                    CHTLookupStats cht_round_stats;
                    for(size_t thread_num = 0; thread_num < num_threads; thread_num++) {
                        cht_round_stats += *lcht_stats[thread_num];
                    }
                    
                    phase_parse.data()["bloom_lookups"] = cht_round_stats.bloom_lookups;
                    phase_parse.data()["map_lookups"] = cht_round_stats.map_lookups;
                    phase_parse.data()["false_positives"] = cht_round_stats.false_positives;
                    phase_parse.data()["bloom_fpr"] = double(cht_round_stats.false_positives) / double(cht_round_stats.bloom_lookups);
                    phase_parse.data()["bloom_load_factor"] = double(bloom.count_set_bits()) / double(bloom.size());
                    cht_total_stats += cht_round_stats;
                }
                
                phase_process_len.append_child(phase_parse);

                // merge refs and lrefs, then sort
                internal::TimePhase phase_merge_and_sort("merge and sort");
                phase_merge_and_sort.start();

                for(size_t thread_num = 0; thread_num < num_threads; thread_num++) {
                    for(auto ref : *lrefs[thread_num]) {
                        refs.push_back(ref);
                    }
                }
                std::sort(std::execution::par_unseq, refs.begin(), refs.end(), [](LZRef2At const& a, LZRef2At const& b){ return a.pos <= b.pos; });

                phase_merge_and_sort.stop();
                phase_process_len.append_child(phase_merge_and_sort);

                phase_process_len.stop();
                phase_sample_and_parse.append_child(phase_process_len);
            }

            phase_sample_and_parse.stop();
            stats_.append_child(phase_sample_and_parse);
        }

        stats_.stop();
        result_.add("p", num_threads);
        result_.add("sampling", sampling_);
        result_.add("window", window_size);
        result_.add("bloom_size", bloom_filter_size);
        result_.add("bloom_fpr", double(cht_total_stats.false_positives) / double(cht_total_stats.bloom_lookups));
        result_.add("refs", refs.size());
        result_.add("t", stats_.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());

        size_t const mem = stats_.get_metric<pm::MallocCounter::MemoryPeakMetric>();
        double const exp_num_fingerprints = double(n) / double(1ULL << sampling_);
        result_.add("mem_peak", stats_.get_metric<pm::MallocCounter::MemoryPeakMetric>());
        result_.add("mem_avg_per_fp", double(mem - BloomFilter::memory_usage(bloom_filter_size) - window_size_) / double(exp_num_fingerprints));
        return refs;
    }

public:
    PSampLZ(size_t len_exp_min, size_t len_exp_max, size_t sampling, size_t bloom_filter_scale = 6, size_t window_size = 64_Mi)
        : len_exp_min_(len_exp_min),
          len_exp_max_(len_exp_max),
          sampling_(sampling),
          window_size_(window_size),
          bloom_filter_scale_(bloom_filter_scale) {
    }

    PSampLZ(PSampLZ&&) = default;
    PSampLZ& operator=(PSampLZ&&) = default;
    PSampLZ(PSampLZ const&) = default;
    PSampLZ& operator=(PSampLZ const&) = default;

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
        auto buffer = std::make_unique<char[]>(window_size_);
        size_t z = 0;
        size_t nout = 0;
        {
            internal::TimePhase phase_encode("encode");
            phase_encode.start();

            in.seekg(0, std::ios_base::beg);
            {
                // magic
                out.write(MAGIC.data(), MAGIC.length());
                nout += MAGIC.length();

                size_t i = 0;

                auto copy_literals = [&](size_t const num_literals){
                    nout += internal::encode_vbyte(out, num_literals);
                    
                    size_t num_written = 0;
                    while(num_written < num_literals) {
                        in.read(buffer.get(), std::min(num_literals, window_size_));
                        auto const w = in.gcount();
                        out.write(buffer.get(), w);
                        num_written += w;
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

            phase_encode.stop();
            stats_.append_child(phase_encode);
        }

        stats_.stop();

        // result
        result_.add("z", z);
        result_.add("nout", nout);
    }

    auto&& consume_last_result() {
        return std::move(result_);
    }

    auto&& consume_last_stats() {
        return std::move(stats_);
    }
};

}
