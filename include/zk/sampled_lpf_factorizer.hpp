#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <iterator>
#include <memory>
#include <numeric>
#include <vector>

#include <fp/rk31.hpp>

#include <lz77/factor.hpp>

#include <pm/stopwatch.hpp> // ONLY FOR DEVELOPMENT

namespace zk {

class SampledLPFFactorizer {
private:
    using RK = fp::RabinKarp31;
    using Fingerprint = RK::Fingerprint;

    static constexpr Fingerprint rolling_fp_base_ = (1ULL << 16) - 39;
    static constexpr size_t MAX_SIZE_32BIT = 1ULL << 31;

    static size_t lce(std::string_view const& t, size_t const i, size_t const j) {
        auto const n = t.length();

        size_t l = 0;
        while(i + l < n && j + l < n && t[i + l] == t[j + l]) ++l;

        return l;
    }

    size_t sampling_;
    size_t fp_window_;
    size_t min_ref_len_;

    template<bool require_64bit, std::output_iterator<lz77::Factor> Output>
    void factorize(std::string_view const& t, Output& out) {
        using Index = std::conditional_t<require_64bit, uint64_t, uint32_t>;

        // sample
        Index const n = t.size();

        RK rk(rolling_fp_base_, fp_window_);
        Fingerprint fp = 0;

        pm::Stopwatch sw;

        std::cout << "sample ... ";
        std::cout.flush();

        sw.start();
        std::vector<Index> samples;
        {
            size_t const s = (1ULL << sampling_) - 1;

            size_t i = 0;
            for(; i < n && i < fp_window_; i++) {
                fp = rk.push(fp, t[i]);
            }

            for(; i < n; i++) {
                if((fp & s) == 0) {
                    samples.push_back(i - fp_window_);
                }
                fp = rk.roll(fp, t[i - fp_window_], t[i]);
            }

            if((fp & s) == 0) {
                samples.push_back(i - fp_window_);
            }
        }
        sw.stop();
        std::cout << samples.size() << " (" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;

        samples.shrink_to_fit();
        size_t const m = samples.size();

        // construct suffix array of samples
        std::cout << "suffix array ... ";
        std::cout.flush();

        auto sa = std::make_unique<Index[]>(m);
        std::iota(sa.get(), sa.get() + m, 0);
        sw.start();
        {
            auto compare = [&](Index const a, Index const b) {
                size_t i = samples[a];
                size_t j = samples[b];
                assert(i != j);
                
                for(; i < n && j < n && t[i] == t[j]; i++, j++);

                if(i == n) return true;
                else if(j == n) return false;
                else return t[i] < t[j];
            };
            std::sort(std::execution::par_unseq, sa.get(), sa.get() + m, compare); // FIXME: way too slow -- gsaca may help?
        }
        sw.stop();
        std::cout << "(" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;

        // construct inverse
        std::cout << "inverse suffix array ... ";
        std::cout.flush();
        sw.start();

        auto isa = std::make_unique<Index[]>(m);
        for(size_t i = 0; i < m; i++) isa[sa[i]] = i;

        sw.stop();
        std::cout << "(" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;

        // factorize
        std::cout << "inverse suffix array ... ";
        std::cout.flush();
        sw.start();

        // TODO: implement

        sw.stop();
        std::cout << "(" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;
    }

public:
    SampledLPFFactorizer(size_t sampling, size_t fp_window)
        : sampling_(sampling), fp_window_(fp_window), min_ref_len_(2) {
    }

    template<std::contiguous_iterator Input, std::output_iterator<lz77::Factor> Output>
    requires (sizeof(std::iter_value_t<Input>) == 1)
    void factorize(Input begin, Input const& end, Output out) {
        std::string_view const t(begin, end);
        size_t const n = t.size();

        if(n < MAX_SIZE_32BIT) {
            factorize<false>(t, out);
        } else {
            factorize<true>(t, out);
        }
    }

    /**
     * \brief Reports the minimum length of a referencing factor
     * 
     * If a referencing factor is shorter than this length, a literal factor is emitted instead
     * 
     * \return the minimum reference length
     */
    size_t min_reference_length() const { return min_ref_len_; }

    /**
     * \brief Sets the minimum length of a referencing factor
     * 
     * If a referencing factor is shorter than this length, a literal factor is emitted instead
     * 
     * \param min_ref_len the minimum reference length
     */
    void min_reference_length(size_t min_ref_len) { min_ref_len_ = min_ref_len; }
};

}
