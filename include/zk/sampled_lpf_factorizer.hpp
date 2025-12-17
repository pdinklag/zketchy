#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <functional>
#include <iterator>
#include <memory>
#include <numeric>
#include <vector>

#include <omp.h>

#include <fp/rk31.hpp>
#include <fp/rk61.hpp>
#include <ankerl/unordered_dense.h>

#include <lz77/emit_function.hpp>
#include <libsais.h>
#include <libsais64.h>

#include "internal/benchmark.hpp"
#include "internal/util/idiv_ceil.hpp"
#include "internal/util/si_iec_literals.hpp"

namespace zk {

class SampledLPFFactorizer {
private:
    #ifdef _ZK_SAMPLED_LPF_DEBUG
    static constexpr bool debug_ = true;
    #else
    static constexpr bool debug_ = false;
    #endif

    using RK = fp::RabinKarp31;
    using RK64 = fp::RabinKarp61;
    using Fingerprint = RK::Fingerprint;
    using Fingerprint64 = RK64::Fingerprint;

    static constexpr Fingerprint rolling_fp_base_ = (1ULL << 16) - 39;
    static constexpr size_t MAX_SIZE_32BIT = 1ULL << 31;

    size_t sampling_;
    size_t fp_window_;

    using MIndex = uint32_t; // nb: we generally assume that we won't ever have more than 4G metacharacters...
    using MLength = uint32_t;

    template<bool require_64bit>
    void factorize(std::string_view const& t, lz77::EmitFunction emit_literal, lz77::EmitFunction emit_reference) {
        using Index = std::conditional_t<require_64bit, uint64_t, uint32_t>;
        
        struct Metachar {
            Index occ;
            MLength len;
            Fingerprint64 fp;
        } __attribute__((packed));

        // prefix free parsing
        Index const n = t.size();
        auto const num_threads = omp_get_max_threads();
        Index const num_per_thread = internal::idiv_ceil(n, num_threads);

        internal::MemoryTimePhase phase;

        size_t gap_total = 0, gap_num = 0;

        if constexpr(debug_) {
            std::cout << "parallel pre-parse ... ";
            std::cout.flush();
            phase.start();
        }

        std::unique_ptr<std::vector<Metachar>> lpre_parse[num_threads];
        {
            RK rk_trigger(rolling_fp_base_, fp_window_);
            RK64 rk_meta(rolling_fp_base_);

            size_t const s = (1ULL << sampling_) - 1;
            size_t const min_metachar_len = s / 2;

            #pragma omp parallel
            {
                Index const thread_num = omp_get_thread_num();
                Index const beg = thread_num * num_per_thread;
                Index const end = std::min(beg + num_per_thread, n);

                lpre_parse[thread_num] = std::make_unique<std::vector<Metachar>>();
                auto& local_pre_parse = *lpre_parse[thread_num];

                Fingerprint fp_trigger = 0;
                Fingerprint64 fp_meta = 0;

                if(thread_num > 0) {
                    // consider previous characters for consistent triggering
                    for(Index j = beg - fp_window_; j < beg; j++) {
                        fp_trigger = rk_trigger.push(fp_trigger, t[j]);
                    }
                }

                Index last = beg;
                Index i = beg;

                for(; i < end && i < fp_window_; i++) {
                    if(thread_num == 0) {
                        // initialize
                        fp_trigger = rk_trigger.push(fp_trigger, t[i]);
                    } else {
                        // roll
                        fp_trigger = rk_trigger.roll(fp_trigger, t[i - fp_window_], t[i]);
                    }
                    fp_meta = rk_meta.push(fp_meta, t[i]);
                }

                for(; i < end; i++) {
                    if(i - last >= min_metachar_len && (fp_trigger & s) == 0) {
                        local_pre_parse.push_back(Metachar{ last, MLength(i - last), fp_meta });
                        last = i;
                        fp_meta = 0;
                    }

                    fp_trigger = rk_trigger.roll(fp_trigger, t[i - fp_window_], t[i]);
                    fp_meta = rk_meta.push(fp_meta, t[i]);
                }

                if(last < i) {
                    local_pre_parse.push_back(Metachar{ last, MLength(i - last), fp_meta });
                }

                local_pre_parse.shrink_to_fit();
            }
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;

            size_t pre_parse_size = 0;
            for(size_t thread_num = 0; thread_num < num_threads; thread_num++) {
                pre_parse_size += lpre_parse[thread_num]->size();
            }
            std::cout << "\tparse size: " << pre_parse_size << std::endl;
        }

        if constexpr(debug_) {
            std::cout << "compute distinct metacharacters ... ";
            std::cout.flush();
            phase.start();
        }

        std::vector<Metachar> meta;
        std::vector<MIndex> parse;
        {
            ankerl::unordered_dense::map<Fingerprint64, MIndex> meta_fps;
            for(size_t thread_num = 0; thread_num < num_threads; thread_num++) {
                for(auto& x : *lpre_parse[thread_num]) {
                    auto const fp = x.fp;
                    auto it = meta_fps.find(fp);
                    if(it != meta_fps.end()) {
                        parse.push_back(it->second);
                    } else {
                        parse.push_back(MIndex(meta.size()));

                        meta_fps.emplace(fp, meta.size());
                        meta.push_back(x);
                    }
                }
                lpre_parse[thread_num].reset();
            }
        }

        if(meta.size() >= 4_Gi) std::abort(); // if this happens, you wouldn't want to wait for the result anyway

        meta.shrink_to_fit();
        parse.shrink_to_fit();

        auto const m = parse.size();
        auto const sigma = meta.size();

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
            std::cout << "\tdistinct metacharacters: " << sigma << std::endl;
        }

        // sort meta characters
        if constexpr(debug_) {
            std::cout << "sort metacharacters ... ";
            std::cout.flush();
            phase.start();
        }

        auto meta_sorted = std::make_unique<MIndex[]>(sigma);
        {
            for(size_t i = 0; i < sigma; i++) {
                meta_sorted[i] = i;
            }

            std::sort(std::execution::par_unseq, meta_sorted.get(), meta_sorted.get() + sigma, [&](MIndex const a, MIndex const b){
                return std::string_view(t.data() + meta[a].occ, meta[a].len).compare(std::string_view(t.data() + meta[b].occ, meta[b].len)) < 0;
            });
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
        }

        // parse text
        if constexpr(debug_) {
            std::cout << "compute parsing ... ";
            std::cout.flush();
            phase.start();
        }

        auto parse_beg = std::make_unique<Index[]>(m);
        {
            // compute starting positions of phrases
            {
                size_t i = 0;
                for(size_t j = 0; j < m; j++) {
                    parse_beg[j] = i;
                    i += meta[parse[j]].len;
                }
            }

            // inverse metacharacter order
            auto meta_sorted_inv = std::make_unique<MIndex[]>(sigma);
            for(size_t i = 0; i < sigma; i++) {
                meta_sorted_inv[meta_sorted[i]] = i;
            }

            // rewrite parsing
            for(size_t j = 0; j < m; j++) {
                parse[j] = meta_sorted_inv[parse[j]];
            }
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
        }

        // compute suffix array of parsing
        if constexpr(debug_) {
            std::cout << "compute suffix array ... ";
            std::cout.flush();
            phase.start();
        }

        auto const sa_extra_space = 6 * sigma; // recommended for libsais
        auto sa = std::make_unique<Index[]>(m + sa_extra_space);
        if constexpr(require_64bit) {
            libsais64_long((int64_t*)parse.data(), (int64_t*)sa.get(), m, meta.size(), sa_extra_space);
        } else {
            libsais_int((int32_t*)parse.data(), (int32_t*)sa.get(), m, meta.size(), sa_extra_space);
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
        }

        // compute inverse
        if constexpr(debug_) {
            std::cout << "compute inverse suffix array ... ";
            std::cout.flush();
            phase.start();
        }

        auto isa = std::make_unique<Index[]>(m);
        {
            for(size_t i = 0; i < m; i++) {
                isa[sa[i]] = i;
            }
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
        }

        // factorize
        auto get_meta = [&](MIndex const i){ return meta[meta_sorted[i]]; };

        if constexpr(debug_) {
            std::cout << "factorizing ... ";
            std::cout.flush();
            phase.start();
        }

        internal::TimePhase phase_gather, phase_emit;
        {
            struct Ref {
                Index beg, src, len;
                Index end() const { return beg + len - 1; }
            } __attribute__((packed));
            auto const num_threads = omp_get_max_threads();
            std::unique_ptr<std::vector<Ref>> lrefs[num_threads];
            for(size_t x = 0; x < num_threads; x++) {
                lrefs[x] = std::make_unique<std::vector<Ref>>();
            }

            phase_gather.start();

            #pragma omp parallel for
            for(size_t j = 0; j < m; j++) {
                auto const thread_num = omp_get_thread_num();

                // get SA position for parse suffix j
                size_t const cur_pos = isa[j];

                // longest common extension
                auto lce = [&](size_t const a, size_t const b, size_t& matched_meta, size_t& rext){
                    size_t l = 0;

                    // first compare meta characters
                    matched_meta = 0;
                    while(a + matched_meta < m && b + matched_meta < m && parse[a + matched_meta] == parse[b + matched_meta]) {
                        l += get_meta(parse[a + matched_meta]).len;
                        ++matched_meta;
                    }

                    // once we have a mismatch, extend ...
                    rext = 0;

                    if(matched_meta >= 1) {
                        // ... to the right
                        {
                            size_t const x = parse_beg[a] + l;
                            size_t const y = parse_beg[b] + l;
                            while(x + rext < n && y + rext < n && t[x + rext] == t[y + rext]) {
                                ++rext;
                            }
                        }
                    }
                    return l + rext;
                };

                // compute PSV and NSV as well as longest common prefixes
                ssize_t psv_pos = (ssize_t)cur_pos - 1;
                while (psv_pos >= 0 && sa[psv_pos] > j) --psv_pos;

                size_t psv_matched_meta, psv_rext;
                size_t const psv_lcp = psv_pos >= 0 ? lce(j, (size_t)sa[psv_pos], psv_matched_meta, psv_rext) : 0;

                size_t nsv_pos = cur_pos + 1;
                while(nsv_pos < m && sa[nsv_pos] > j) ++nsv_pos;

                size_t nsv_matched_meta, nsv_rext;
                size_t const nsv_lcp = nsv_pos < m ? lce(j, (size_t)sa[nsv_pos], nsv_matched_meta, nsv_rext) : 0;

                if(psv_matched_meta >= 1 || nsv_matched_meta >= 1) {
                    // select maximum
                    size_t dst, src, lcp, matched_meta, rext;
                    if(psv_lcp > nsv_lcp) {
                        lcp = psv_lcp;
                        dst = parse_beg[j];
                        src = dst - parse_beg[sa[psv_pos]];
                        matched_meta = psv_matched_meta;
                        rext = psv_rext;
                    } else {
                        lcp = nsv_lcp;
                        dst = parse_beg[j];
                        src = dst - parse_beg[sa[nsv_pos]];
                        matched_meta = nsv_matched_meta;
                        rext = nsv_rext;
                    }

                    if(lcp > 1)[[likely]] { // nb: may occur as a rare bordercase
                        // emit reference
                        lrefs[thread_num]->emplace_back(dst, src, lcp);

                        j += matched_meta - 1;
                        if(rext > 0) {
                            // we have encoded characters from the following meta characters
                            ++j;
                            while(j < m && rext >= get_meta(parse[j]).len) {
                                // TODO: we may have skipped additional metacharacters... but HOW???
                                rext -= get_meta(parse[j]).len;
                                ++j;
                            }
                        }
                    }
                }
            }

            for(size_t x = 0; x < num_threads; x++) {
                lrefs[x]->shrink_to_fit();
            }
            phase_gather.stop();

            // emit
            phase_emit.start();
            {
                size_t cur_gap = 0;

                size_t i = 0;


                size_t x = 0;
                size_t j = 0;
                auto has_next_ref = [&](){ return x < num_threads && j < lrefs[x]->size(); };
                auto next_ref = [&](){ return (*lrefs[x])[j]; };
                auto advance_ref = [&](){
                    ++j;
                    if(x < num_threads && j >= lrefs[x]->size()) {
                        ++x;
                        j = 0;
                    }
                };
                
                while(i < n) {
                    while(has_next_ref() && i > next_ref().end()) {
                        advance_ref();
                    }

                    if(has_next_ref() && i >= next_ref().beg) {
                        auto const& ref = next_ref();
                        auto const d = i - ref.beg;
                        emit_reference(lz77::Factor(ref.src, ref.len - d));
                        i = ref.end() + 1;
                        advance_ref();

                        if(cur_gap > 0) ++gap_num;
                        cur_gap = 0;
                    } else {
                        emit_literal(t[i++]);

                        ++gap_total;
                        ++cur_gap;
                    }
                }
                if(cur_gap > 0) ++gap_num;
            }
            phase_emit.stop();
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
            std::cout << "\tt_gather=" << (size_t)phase_gather.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>()
                      << ", t_emit=" << (size_t)phase_emit.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>()
                      << std::endl;

            double const avg_gap_len = double(gap_total) / double(gap_num);
            std::cout << "\taverage gap length: " << avg_gap_len << " (of " << gap_num << " gaps with total length " << gap_total << ")" << std::endl;
        }
    }

public:
    SampledLPFFactorizer(size_t sampling, size_t fp_window)
        : sampling_(sampling), fp_window_(fp_window) {
    }

    template<std::contiguous_iterator Input>
    requires (sizeof(std::iter_value_t<Input>) == 1)
    void factorize(Input begin, Input const& end, lz77::EmitFunction emit_literal, lz77::EmitFunction emit_reference) {
        std::string_view const t(begin, end);
        size_t const n = t.size();

        // FIXME: whether or not we need 64 bits should depend on the size of the PARSING, not the input
        if(n < MAX_SIZE_32BIT) {
            factorize<false>(t, emit_literal, emit_reference);
        } else {
            factorize<true>(t, emit_literal, emit_reference);
        }
    }

    template<std::contiguous_iterator Input, std::output_iterator<lz77::Factor> Output>
    requires (sizeof(std::iter_value_t<Input>) == 1)
    void factorize(Input begin, Input const& end, Output out) {
        auto emit = [&](lz77::Factor f){ *out++ = f; };
        factorize(begin, end, emit, emit);
    }
};

}
