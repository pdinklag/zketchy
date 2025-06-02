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
#include <fp/rk61.hpp>
#include <ankerl/unordered_dense.h>

#include <lz77/factor.hpp>
#include <libsais.h>
#include <libsais64.h>

#include <pm/stopwatch.hpp> // ONLY FOR DEVELOPMENT

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

    struct Metachar {
        std::string_view s;
        Fingerprint64 fp;
        size_t rank;

        size_t size() const {
            return s.size();
        }
    };

    template<bool require_64bit, std::output_iterator<lz77::Factor> Output>
    void factorize(std::string_view const& t, Output& out) {
        using Index = std::conditional_t<require_64bit, uint64_t, uint32_t>;

        // prefix free parsing
        Index const n = t.size();

        RK rk(rolling_fp_base_, fp_window_);
        RK64 rk64(rolling_fp_base_);
        Fingerprint fp = 0;
        Fingerprint64 fp64 = 0;

        pm::Stopwatch sw;

        size_t gap_total = 0, gap_num = 0;

        if constexpr(debug_) {
            std::cout << "compute metacharacters ... ";
            std::cout.flush();
            sw.start();
        }

        std::vector<Metachar> meta;
        std::vector<Fingerprint64> pre_parse;
        {
            size_t const s = (1ULL << sampling_) - 1;

            size_t beg = 0;
            size_t i = 0;

            ankerl::unordered_dense::set<Fingerprint64> meta_fps;
            auto on_trigger = [&](){
                Metachar x { t.substr(beg, i - beg), fp64, 0 };

                if(!meta_fps.contains(fp64)) {
                    meta.push_back(x);
                    meta_fps.emplace(fp64);
                }
                pre_parse.push_back(fp64);

                beg = i;
                fp64 = 0;
            };

            for(; i < n && i < fp_window_; i++) {
                fp = rk.push(fp, t[i]);
                fp64 = rk64.push(fp64, t[i]);
            }

            for(; i < n; i++) {
                if((fp & s) == 0) {
                    on_trigger();
                }

                fp = rk.roll(fp, t[i - fp_window_], t[i]);
                fp64 = rk64.push(fp64, t[i]);
            }

            if(beg < i) {
                on_trigger();
            }
        }

        if constexpr(debug_) {
            sw.stop();
            std::cout << " found " << meta.size() << " distinct meta characters, parsing size: " << pre_parse.size() << " (" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;
        }

        // sort meta characters
        if constexpr(debug_) {
            std::cout << "sort meta characters ... ";
            std::cout.flush();
            sw.start();
        }
        
        {
            std::sort(std::execution::par_unseq, meta.begin(), meta.end(), [&](Metachar const& a, Metachar const& b){
                return a.s.compare(b.s) < 0;
            });

            for(size_t i = 0; i < meta.size(); i++) {
                meta[i].rank = i;
            }
        }

        if constexpr(debug_) {
            sw.stop();
            std::cout << "(" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;
        }

        // parse text
        if constexpr(debug_) {
            std::cout << "compute parsing ... ";
            std::cout.flush();
            sw.start();
        }

        std::vector<Index> parse;
        std::vector<Index> parse_beg;
        {
            // map fingerprints to ranks
            ankerl::unordered_dense::map<Fingerprint64, Index> metachar_ranks;
            for(auto const& x : meta) {
                metachar_ranks.emplace(x.fp, x.rank);
            }

            // parse
            parse.reserve(pre_parse.size());
            for(auto fp64 : pre_parse) {
                parse.push_back(metachar_ranks.find(fp64)->second);
            }

            // compute starting positions of phrases
            parse_beg.reserve(pre_parse.size());
            {
                size_t i = 0;
                for(auto x : parse) {
                    parse_beg.push_back(i);
                    i += meta[x].size();
                }
            }

            // discard no longer needed stuff
            pre_parse = {};
        }

        if constexpr(debug_) {
            sw.stop();
            std::cout << "(" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;
        }

        // compute suffix array of parsing
        if constexpr(debug_) {
            std::cout << "compute suffix array ... ";
            std::cout.flush();
            sw.start();
        }

        auto const m = parse.size();

        auto const sa_extra_space = 6 * meta.size(); // recommended for libsais
        auto sa = std::make_unique<Index[]>(m + sa_extra_space);
        if constexpr(require_64bit) {
            // TODO: libsais64_int ???
            std::abort();
        } else {
            libsais_int((int32_t*)parse.data(), (int32_t*)sa.get(), m, meta.size(), sa_extra_space);
        }

        if constexpr(debug_) {
            sw.stop();
            std::cout << "(" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;
        }

        // compute inverse
        if constexpr(debug_) {
            std::cout << "compute inverse suffix array ... ";
            std::cout.flush();
            sw.start();
        }

        auto isa = std::make_unique<Index[]>(m);
        {
            for(size_t i = 0; i < m; i++) {
                isa[sa[i]] = i;
            }
        }

        if constexpr(debug_) {
            sw.stop();
            std::cout << "(" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;
        }

        // factorize
        if constexpr(debug_) {
            std::cout << "factorizing ... ";
            std::cout.flush();
            sw.start();
        }

        {
            for(size_t j = 0; j < m; j++) {
                // keep track of how many characters from the current meta character we have already encoded
                size_t joffs;

                // get SA position for parse suffix j
                size_t const cur_pos = isa[j];

                // longest common extension
                auto lce = [&](size_t const a, size_t const b, size_t& matched_meta, size_t& ext){
                    size_t l = 0;

                    // first compare meta characters
                    matched_meta = 0;
                    while(a + matched_meta < m && b + matched_meta < m && parse[a + matched_meta] == parse[b + matched_meta]) {
                        l += meta[parse[a + matched_meta]].size();
                        ++matched_meta;
                    }

                    // once we have a mismatch, extend by comparing remaining characters
                    size_t const x = parse_beg[a] + l;
                    size_t const y = parse_beg[b] + l;
                    ext = 0;
                    while(x + ext < n && y + ext < n && t[x + ext] == t[y + ext]) {
                        ++ext;
                    }

                    return l + ext;
                };

                // compute PSV and NSV as well as longest common prefixes
                ssize_t psv_pos = (ssize_t)cur_pos - 1;
                while (psv_pos >= 0 && sa[psv_pos] > j) --psv_pos;
                size_t psv_matched_meta, psv_ext;
                size_t const psv_lcp = psv_pos >= 0 ? lce(j, (size_t)sa[psv_pos], psv_matched_meta, psv_ext) : 0;

                size_t nsv_pos = cur_pos + 1;
                while(nsv_pos < m && sa[nsv_pos] > j) ++nsv_pos;
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
                if(lcp >= 2) {
                    *out++ = lz77::Factor(src, lcp);

                    j += matched_meta - 1;
                    if(ext > 0) {
                        // we have encoded characters from the following meta characters
                        ++j;
                        while(j < m && ext >= meta[parse[j]].size()) {
                            // TODO: we may have skipped additional metacharacters... but HOW???
                            ext -= meta[parse[j]].size();
                            ++j;
                        }
                        joffs = ext;
                    } else {
                        // we have fully encoded all meta characters with previous occurrence
                        joffs = meta[parse[j]].size();
                    }
                } else {
                    joffs = 0; // we have not encoded anything
                }

                // emit any remaining characters from current meta character as literals
                if(j < m) {
                    size_t const gap_len = meta[parse[j]].size() - joffs;
                    gap_total += gap_len;
                    ++gap_num;

                    for(; joffs < meta[parse[j]].size(); joffs++) {
                        // TODO: keep table for short repetitions?
                        *out++ = lz77::Factor(meta[parse[j]].s[joffs]);
                    }
                }
            }
        }

        if constexpr(debug_) {
            sw.stop();
            std::cout << "(" << (size_t)sw.elapsed_time_millis() << "ms)" << std::endl;

            double const avg_gap_len = double(gap_total) / double(gap_num);
            std::cout << "average gap length: " << avg_gap_len << " (of " << gap_num << " gaps)" << std::endl;
        }
    }

public:
    SampledLPFFactorizer(size_t sampling, size_t fp_window)
        : sampling_(sampling), fp_window_(fp_window) {
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
};

}
