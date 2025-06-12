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

#include <fp/rk31.hpp>
#include <fp/rk61.hpp>
#include <ankerl/unordered_dense.h>

#include <lz77/emit_function.hpp>
#include <libsais.h>
#include <libsais64.h>

#include "internal/benchmark.hpp"
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

        RK rk_trigger(rolling_fp_base_, fp_window_);
        RK64 rk_meta(rolling_fp_base_);
        Fingerprint fp_trigger = 0;
        Fingerprint64 fp_meta = 0;
        Fingerprint64 fp_short = 0;

        internal::MemoryTimePhase phase;

        size_t gap_total = 0, gap_num = 0;

        if constexpr(debug_) {
            std::cout << "compute metacharacters ... ";
            std::cout.flush();
            phase.start();
        }

        std::vector<Metachar> meta;
        std::vector<MIndex> parse;
        {
            size_t const s = (1ULL << sampling_) - 1;

            Index beg = 0;
            Index i = 0;

            ankerl::unordered_dense::map<Fingerprint64, MIndex> meta_fps;
            auto on_trigger = [&](){
                if(i - beg > std::numeric_limits<MLength>::max()) std::abort();

                Metachar x { beg, MLength(i - beg), fp_meta };

                auto it = meta_fps.find(fp_meta);
                if(it != meta_fps.end()) {
                    parse.push_back(it->second);
                } else {
                    meta_fps.emplace(fp_meta, meta.size());
                    parse.push_back(MIndex(meta.size()));
                    meta.push_back(x);
                }

                beg = i;
                fp_meta = 0;
            };

            for(; i < n && i < fp_window_; i++) {
                fp_trigger = rk_trigger.push(fp_trigger, t[i]);
                fp_meta = rk_meta.push(fp_meta, t[i]);
            }

            for(; i < n; i++) {
                if((fp_trigger & s) == 0) {
                    on_trigger();
                }

                fp_trigger = rk_trigger.roll(fp_trigger, t[i - fp_window_], t[i]);
                fp_meta = rk_meta.push(fp_meta, t[i]);
            }

            if(beg < i) {
                on_trigger();
            }
        }

        if(meta.size() >= 4_Gi) std::abort(); // if this happens, you wouldn't want to wait for the result anyway

        meta.shrink_to_fit();
        parse.shrink_to_fit();

        auto const m = parse.size();
        auto const sigma = meta.size();

        if constexpr(debug_) {
            phase.stop();
            std::cout << " found " << sigma << " distinct meta characters, parsing size: " << m
                << " (" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
            std::cout << "\tsizeof(meta)=" << sigma * sizeof(Metachar) << ", sizeof(pre_parse)=" << m * sizeof(Fingerprint64) << std::endl;
        }

        // sort meta characters
        if constexpr(debug_) {
            std::cout << "sort meta characters ... ";
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

        internal::TimePhase phase_lce, phase_nsv_psv, phase_emit;
        
        {
            for(size_t j = 0; j < m; j++) {
                // keep track of how many characters from the current meta character we have already encoded
                size_t joffs;

                // get SA position for parse suffix j
                size_t const cur_pos = isa[j];

                // longest common extension
                auto lce = [&](size_t const a, size_t const b, size_t& matched_meta, size_t& ext){
                    phase_lce.resume();
                    size_t l = 0;

                    // first compare meta characters
                    matched_meta = 0;
                    while(a + matched_meta < m && b + matched_meta < m && parse[a + matched_meta] == parse[b + matched_meta]) {
                        l += get_meta(parse[a + matched_meta]).len;
                        ++matched_meta;
                    }

                    // once we have a mismatch, extend by comparing remaining characters
                    size_t const x = parse_beg[a] + l;
                    size_t const y = parse_beg[b] + l;
                    ext = 0;
                    while(x + ext < n && y + ext < n && t[x + ext] == t[y + ext]) {
                        ++ext;
                    }

                    phase_lce.pause();
                    return l + ext;
                };

                // compute PSV and NSV as well as longest common prefixes
                phase_nsv_psv.resume();
                ssize_t psv_pos = (ssize_t)cur_pos - 1;
                while (psv_pos >= 0 && sa[psv_pos] > j) --psv_pos;
                phase_nsv_psv.pause();

                size_t psv_matched_meta, psv_ext;
                size_t const psv_lcp = psv_pos >= 0 ? lce(j, (size_t)sa[psv_pos], psv_matched_meta, psv_ext) : 0;

                phase_nsv_psv.resume();
                size_t nsv_pos = cur_pos + 1;
                while(nsv_pos < m && sa[nsv_pos] > j) ++nsv_pos;
                phase_nsv_psv.pause();

                size_t nsv_matched_meta, nsv_ext;
                size_t const nsv_lcp = nsv_pos < m ? lce(j, (size_t)sa[nsv_pos], nsv_matched_meta, nsv_ext) : 0;

                // select maximum
                size_t lcp, src, matched_meta, ext;
                if(psv_lcp > nsv_lcp) {
                    lcp = psv_lcp;
                    src = parse_beg[j] - parse_beg[sa[psv_pos]];
                    matched_meta = psv_matched_meta;
                    ext = psv_ext;
                } else {
                    lcp = nsv_lcp;
                    src = parse_beg[j] - parse_beg[sa[nsv_pos]];
                    matched_meta = nsv_matched_meta;
                    ext = nsv_ext;
                }

                // emit reference
                phase_emit.resume();
                if(lcp >= 2) {
                    emit_reference(lz77::Factor(src, lcp));

                    j += matched_meta - 1;
                    if(ext > 0) {
                        // we have encoded characters from the following meta characters
                        ++j;
                        while(j < m && ext >= get_meta(parse[j]).len) {
                            // TODO: we may have skipped additional metacharacters... but HOW???
                            ext -= get_meta(parse[j]).len;
                            ++j;
                        }
                        joffs = ext;
                    } else {
                        // we have fully encoded all meta characters with previous occurrence
                        joffs = get_meta(parse[j]).len;
                    }
                } else {
                    joffs = 0; // we have not encoded anything
                }

                // emit any remaining characters from current meta character as literals
                if(j < m) {
                    size_t const gap_len = get_meta(parse[j]).len - joffs;
                    gap_total += gap_len;
                    ++gap_num;

                    for(; joffs < get_meta(parse[j]).len; joffs++) {
                        // TODO: keep table for short repetitions?
                        emit_literal(lz77::Factor(t[get_meta(parse[j]).occ + joffs]));
                    }
                }
                phase_emit.pause();
            }
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;

            std::cout << "\tt_lce=" << (size_t)phase_lce.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>()
                      << ", t_nsv_psv=" << (size_t)phase_nsv_psv.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>()
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
