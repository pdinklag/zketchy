#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <execution>
#include <functional>
#include <iterator>
#include <limits>
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
#include "internal/io/memory_input_stream.hpp"
#include "internal/io/overlapping_blocks.hpp"
#include "internal/util/idiv_ceil.hpp"
#include "internal/util/si_iec_literals.hpp"

namespace zk {

template<std::unsigned_integral Index>
class SampledLPFFactorizer {
private:
    static constexpr bool debug_ = true;

    using RK = fp::RabinKarp31;
    using RK64 = fp::RabinKarp61;
    using Fingerprint = RK::Fingerprint;
    using Fingerprint64 = RK64::Fingerprint;

    using MIndex = uint32_t; // nb: we generally assume that we won't ever have more than 4G metacharacters...
    using MLength = uint32_t;

    struct Metachar {
        Index occ;
        MLength len;
        Fingerprint64 fp;
    } __attribute__((packed));

    static constexpr Fingerprint rolling_fp_base_ = (1ULL << 16) - 39;

    size_t sampling_;
    size_t fp_window_;

    template<typename InputStream, bool has_text_access>
    void factorize(InputStream& in, size_t const block_size, std::string_view const& t, size_t const n, lz77::EmitFunction emit_literal, lz77::EmitFunction emit_reference) {
        auto const num_threads = omp_get_max_threads();
        size_t const s = (1ULL << sampling_) - 1;

        internal::MemoryTimePhase phase;

        // parsing
        std::vector<Metachar> meta;
        std::unique_ptr<char[]> meta_buf;
        std::unique_ptr<Index[]> meta_ptr;
        std::vector<Index> parsing;

        {
            if constexpr(debug_) {
                std::cout << "parallel pre-parse ... ";
                std::cout.flush();
                phase.start();
            }

            RK rk_trigger(rolling_fp_base_, fp_window_);
            RK64 rk_meta(rolling_fp_base_);

            auto push_trigger = [&](Fingerprint& fp, char const* p){ fp = rk_trigger.push(fp, *p); };
            auto roll_trigger = [&](Fingerprint& fp, char const* p){ fp = rk_trigger.roll(fp, *(p-fp_window_), *p); };
            auto push_meta = [&](Fingerprint64& fp, char const* p){ fp = rk_meta.push(fp, *p); };

            std::vector<std::unique_ptr<std::vector<Metachar>>> pre_parsing;

            internal::OverlappingBlocks<InputStream> block;
            bool done;
            if constexpr(has_text_access) {
                // we only scan one "block" - the whole input
                done = true;
            } else {
                // we will process the input blockwise
                block = internal::OverlappingBlocks(in, block_size, fp_window_);
            }

            // when the last thread reaches the end of a block (not the last) and has not yet found a trigger string,
            // it leaves this delta for the first thread processing the next block
            Index prev_block_delta;
            do {
                size_t const pre_parsing_offs = pre_parsing.size();
                Index const num_per_thread = internal::idiv_ceil(has_text_access ? n : block.size(), num_threads);

                // add p partial parsings for each block
                for(size_t thread_num = 0; thread_num < num_threads; thread_num++) {
                    pre_parsing.emplace_back(std::make_unique<std::vector<Metachar>>());
                }

                #pragma omp parallel
                {
                    Index const thread_num = omp_get_thread_num();

                    auto& local_pre_parsing = *pre_parsing[pre_parsing_offs + thread_num];
                    local_pre_parsing.reserve(size_t(1.2 * double(has_text_access ? n : block.size()) / double(s * num_threads))); // exaggerate a little bit to account for standard deviation

                    char const* t_beg, *t_end;
                    if constexpr(has_text_access) {
                        t_beg = t.data();
                        t_end = t_beg + n;
                    } else {
                        t_beg = block.begin();
                        t_end = block.end();
                    }

                    char const* beg = t_beg + thread_num * num_per_thread;
                    char const* end = std::min(beg + num_per_thread, t_end);

                    Fingerprint fp_trigger = 0;
                    Fingerprint64 fp_meta = 0;

                    if(thread_num > 0 || !block.first()) {
                        // consider previous characters for consistent triggering
                        // (when scanning blockwise, the first thread must do this too unless this is the first block)
                        for(char const* p = beg - fp_window_; p < beg; p++) {
                            push_trigger(fp_trigger, p);
                        }
                    }

                    char const* last = beg;
                    if(thread_num == 0 && !block.first()) {
                        // if this is not the first block, the first thread must use whatever memo the last thread left from the previous block
                        last -= prev_block_delta;
                    }

                    char const* p = beg;

                    // the first thread must fingerprint the initial window if this is the first block
                    if(thread_num == 0 && block.first()) {
                        for(; p < end && p < beg + fp_window_; p++) {
                            push_trigger(fp_trigger, p);
                            push_meta(fp_meta, p);
                        }
                    }

                    for(; last < end && p < t_end; p++) {
                        if((fp_trigger & s) == 0) {
                            local_pre_parsing.push_back(Metachar{ Index(block.offset() + (last - t_beg)), MLength(p - last), fp_meta });
                            last = p - fp_window_;

                            fp_meta = 0;
                            for(char const* q = p - fp_window_; q < p; q++) {
                                push_meta(fp_meta, q);
                            }
                        }

                        roll_trigger(fp_trigger, p);
                        push_meta(fp_meta, p);
                    }

                    // the last thread must either introduce the final metacharacter, or leave information for the first thread regarding the next block
                    if(thread_num == num_threads - 1 && last < t_end) {
                        if(block.last()) {
                            // we are in the last block -- introduce final metacharacter
                            // nb: this is safe if not scanning blockwise; block.last() will then always return true
                            local_pre_parsing.push_back(Metachar{ Index(block.offset() + (last - t_beg)), MLength(t_end - last), fp_meta });
                        } else {
                            // leave a memo for the first thread processing the next block
                            prev_block_delta = t_end - last;
                        }
                    }
                }

                if constexpr(!has_text_access) {
                    // advance to next block
                    done = !block.advance();
                }
            } while(!done);

            size_t pre_parsing_length = 0;
            for(size_t i = 0; i < pre_parsing.size(); i++) {
                pre_parsing_length += pre_parsing[i]->size();
            }

            if constexpr(debug_) {
                phase.stop();
                std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
                std::cout << "\tpreliminary parsing length: " << pre_parsing_length << std::endl;
            }

            if constexpr(debug_) {
                std::cout << "compute distinct metacharacters ... ";
                std::cout.flush();
                phase.start();
            }

            std::vector<Metachar> pre_meta;
            size_t meta_len_total = 0;
            parsing.reserve(pre_parsing_length);
            {
                ankerl::unordered_dense::map<Fingerprint64, MIndex> meta_fps;
                size_t pos = 0;
                size_t skipped = 0;
                for(size_t i = 0; i < pre_parsing.size(); i++) {
                    for(size_t j = 0; j < pre_parsing[i]->size(); j++) {
                        auto const& x = (*pre_parsing[i])[j];
                        if(x.occ < pos) {
                            // a thread may have parsed beyond its block boundary to the next trigger string
                            ++skipped;
                            continue;
                        };
                        if(x.occ > pos)[[unlikely]] std::abort();

                        auto const fp = x.fp;
                        auto it = meta_fps.find(fp);
                        if(it != meta_fps.end()) {
                            parsing.push_back(it->second);
                        } else {
                            parsing.push_back(Index(pre_meta.size()));

                            meta_fps.emplace(fp, pre_meta.size());
                            pre_meta.push_back(x);

                            meta_len_total += x.len;
                        }

                        pos += x.len - fp_window_;
                    }
                    pre_parsing[i].reset();
                }
                std::cout << "(skipped=" << skipped << ", final pos=" << pos << ") ";
            }

            if(pre_meta.size() >= 4_Gi) std::abort(); // if this happens, you wouldn't want to wait for the result anyway
            pre_meta.shrink_to_fit();

            auto const m = parsing.size();
            auto const sigma = pre_meta.size();

            if constexpr(debug_) {
                phase.stop();

                std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
                std::cout << "\tdistinct metacharacters: " << sigma << std::endl;
                std::cout << "\taverage length: " << double(meta_len_total) / double(sigma) << " (total: " << meta_len_total << ")" << std::endl;
            }

            // load metacharacters
            if constexpr(!has_text_access) {
                if constexpr(debug_) {
                    std::cout << "load metacharacters ... ";
                    std::cout.flush();
                    phase.start();
                }

                size_t const meta_bufsize = meta_len_total + fp_window_ * (sigma-1); // nb: account for the fact that we prepend the previous trigger string (except for the first metacharacter)
                meta_buf = std::make_unique<char[]>(meta_bufsize);
                meta_ptr = std::make_unique<Index[]>(sigma);

                char* p = meta_buf.get();

                // load first metacharacter
                in.seekg(0, std::ios_base::beg);
                in.read(p, pre_meta.front().len);
                meta_ptr[0] = 0;
                p += pre_meta.front().len;

                // load remaining metacharacters
                for(size_t i = 1; i < sigma; i++) {
                    meta_ptr[i] = Index(p - meta_buf.get());
                    in.seekg(pre_meta[i].occ, std::ios_base::beg);

                    auto const len = pre_meta[i].len;
                    in.read(p, len);
                    p += len;
                }

                if constexpr(debug_) {
                    phase.stop();
                    std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
                    std::cout << "\tsize=" << meta_bufsize << std::endl;
                }
            }

            // sort meta characters
            if constexpr(debug_) {
                std::cout << "sort metacharacters ... ";
                std::cout.flush();
                phase.start();
            }

            auto meta_order = std::make_unique<MIndex[]>(sigma);
            {
                std::iota(meta_order.get(), meta_order.get() + sigma, 0);

                if constexpr(has_text_access) {
                    // sort by accessing the text
                    std::sort(std::execution::par_unseq, meta_order.get(), meta_order.get() + sigma, [&](MIndex const a, MIndex const b){
                        auto const la = pre_meta[a].len;
                        size_t pa = pre_meta[a].occ;

                        auto const lb = pre_meta[a].len;
                        size_t pb = pre_meta[b].occ;

                        return std::string_view(t.data() + pa, la).compare(std::string_view(t.data() + pb, lb)) < 0;
                    });
                } else {
                    // sort by accessing the buffer
                    std::sort(std::execution::par_unseq, meta_order.get(), meta_order.get() + sigma, [&](MIndex const a, MIndex const b){
                        auto const la = pre_meta[a].len;
                        char const* pa = meta_buf.get() + meta_ptr[a];

                        auto const lb = pre_meta[b].len;
                        char const* pb = meta_buf.get() + meta_ptr[b];

                        return std::string_view(pa, la).compare(std::string_view(pb, lb)) < 0;
                    });
                }
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

            {
                // inverse metacharacter order
                auto meta_sorted_inv = std::make_unique<MIndex[]>(sigma);
                for(size_t i = 0; i < sigma; i++) {
                    meta_sorted_inv[meta_order[i]] = i;
                }

                // rewrite parsing
                for(size_t j = 0; j < m; j++) {
                    parsing[j] = meta_sorted_inv[parsing[j]];
                }
            }

            // rewrite meta
            // if we are streaming, we also need to rewrite the pointers into the dictionary
            meta.reserve(sigma);

            std::unique_ptr<Index[]> new_meta_ptr;
            if constexpr(!has_text_access) {
                new_meta_ptr = std::make_unique<Index[]>(sigma);
            }
            
            for(size_t i = 0; i < sigma; i++) {
                meta.push_back(pre_meta[meta_order[i]]);
                if constexpr(!has_text_access) new_meta_ptr[i] = meta_ptr[meta_order[i]];
            }

            if constexpr(!has_text_access) {
                meta_ptr = std::move(new_meta_ptr);
            }

            if constexpr(debug_) {
                phase.stop();
                std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
                std::cout << "\tparsing length: " << parsing.size() << std::endl;
            }
        }

        auto const m = parsing.size();
        auto const sigma = meta.size();

        // compute starting positions of phrases
        auto parsing_beg = std::make_unique<Index[]>(m);
        {
            size_t i = 0;
            for(size_t j = 0; j < m; j++) {
                parsing_beg[j] = i;
                i += meta[parsing[j]].len - fp_window_; // nb: account for overlap
            }
        }

        // compute suffix array of parsing
        if constexpr(debug_) {
            std::cout << "compute suffix array ... ";
            std::cout.flush();
            phase.start();
        }

        auto const sa_extra_space = 6 * sigma; // recommended for libsais
        auto sa = std::make_unique<Index[]>(m + sa_extra_space);
        if constexpr(std::numeric_limits<Index>::digits > 32) {
            libsais64_long_omp((int64_t*)parsing.data(), (int64_t*)sa.get(), m, sigma, sa_extra_space, num_threads);
        } else {
            libsais_int_omp((int32_t*)parsing.data(), (int32_t*)sa.get(), m, sigma, sa_extra_space, num_threads);
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
            #pragma omp parallel for
            for(size_t i = 0; i < m; i++) {
                isa[sa[i]] = i;
            }
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
        }

        // factorize
        struct Ref {
            Index beg, src, len;
            Index end() const { return beg + len - 1; }
        } __attribute__((packed));

        std::unique_ptr<std::vector<Ref>> lrefs[num_threads];
        for(size_t x = 0; x < num_threads; x++) {
            lrefs[x] = std::make_unique<std::vector<Ref>>();
        }

        if constexpr(debug_) {
            std::cout << "factorize ... ";
            std::cout.flush();
            phase.start();
        }

        #pragma omp parallel for
        for(size_t j = 0; j < m; j++) {
            auto const thread_num = omp_get_thread_num();

            // get SA position for parse suffix j
            size_t const cur_pos = isa[j];

            // longest common extension
            auto lce = [&](size_t const meta_dst, size_t const meta_src, size_t& matched_meta, size_t& rext, size_t& lext){
                size_t l = 0;

                // first compare meta characters
                matched_meta = 0;
                while(meta_dst + matched_meta < m && parsing[meta_dst + matched_meta] == parsing[meta_src + matched_meta]) {
                    l += meta[parsing[meta_dst + matched_meta]].len;
                    ++matched_meta;
                }

                // once we have a mismatch, extend ...
                rext = 0;
                lext = 0;

                if(matched_meta >= 1) {
                    l -= (matched_meta - 1) * fp_window_; // we have to account for matched_meta-1 overlaps in the phrase

                    if constexpr(has_text_access) {
                        // ... by accessing the text ...
                        // ... to the right
                        {
                            size_t const x = parsing_beg[meta_dst] + l;
                            size_t const y = parsing_beg[meta_src] + l;
                            while(x + rext < n && t[x + rext] == t[y + rext]) {
                                ++rext;
                            }
                        }
                        // ... and to the left
                        {
                            size_t x = parsing_beg[meta_dst];
                            size_t y = parsing_beg[meta_src];
                            while(y > 0 && t[x-1] == t[y-1]) {
                                --x;
                                --y;
                                ++lext;
                            }
                        }
                    } else {
                        // ... by accessing the buffer ...
                        auto get_meta_str = [&](size_t const i, char const*& begin, char const*& end){
                            auto const x = parsing[i];
                            begin = meta_buf.get() + meta_ptr[x];
                            end = meta_buf.get() + meta_ptr[x] + meta[x].len;
                        };

                        // ... to the right
                        {
                            size_t x = meta_dst + matched_meta;
                            size_t y = meta_src + matched_meta;
                            if(x < m)[[likely]] {
                                char const *px, *xend, *py, *yend;
                                get_meta_str(x, px, xend); px += fp_window_;
                                get_meta_str(y, py, yend); py += fp_window_;
                                while(*px == *py) {
                                    ++rext;

                                    if(++px >= xend) {
                                        if(++x >= m) break;
                                        get_meta_str(x, px, xend); px += fp_window_;
                                    }

                                    if(++py >= yend) {
                                        get_meta_str(++y, py, yend); py += fp_window_;
                                    }
                                }
                            }
                        }
                        // ... and to the left
                        if(meta_src > 0)[[likely]] {
                            size_t x = meta_dst - 1;
                            size_t y = meta_src - 1;
                            char const *px, *xbeg, *py, *ybeg;
                            get_meta_str(x, xbeg, px); px -= fp_window_ + 1;
                            get_meta_str(y, ybeg, py); py -= fp_window_ + 1;
                            while(*px == *py) {
                                ++lext;

                                if(--py < ybeg) {
                                    if(y-- == 0) break; // nb: important to post-decrement here
                                    get_meta_str(y, ybeg, py); py -= fp_window_ + 1;
                                }

                                if(--px < xbeg) {
                                    get_meta_str(--x, xbeg, px); px -= fp_window_ + 1;
                                }
                            }
                        }
                    }
                }
                return l + rext + lext;
            };

            // compute PSV and NSV as well as longest common prefixes
            ssize_t psv_pos = (ssize_t)cur_pos - 1;
            while (psv_pos >= 0 && sa[psv_pos] > j) --psv_pos;

            size_t psv_matched_meta = 0, psv_rext, psv_lext;
            size_t const psv_lcp = psv_pos >= 0 ? lce(j, (size_t)sa[psv_pos], psv_matched_meta, psv_rext, psv_lext) : 0;

            size_t nsv_pos = cur_pos + 1;
            while(nsv_pos < m && sa[nsv_pos] > j) ++nsv_pos;

            size_t nsv_matched_meta = 0, nsv_rext, nsv_lext;
            size_t const nsv_lcp = nsv_pos < m ? lce(j, (size_t)sa[nsv_pos], nsv_matched_meta, nsv_rext, nsv_lext) : 0;

            if(psv_matched_meta >= 1 || nsv_matched_meta >= 1) {
                // select maximum
                size_t dst, src, lcp, matched_meta, rext;
                if(psv_lcp > nsv_lcp) {
                    lcp = psv_lcp;
                    dst = parsing_beg[j] - psv_lext;
                    src = dst - (parsing_beg[sa[psv_pos]] - psv_lext);
                    matched_meta = psv_matched_meta;
                    rext = psv_rext;
                } else {
                    lcp = nsv_lcp;
                    dst = parsing_beg[j] - nsv_lext;
                    src = dst - (parsing_beg[sa[nsv_pos]] - nsv_lext);
                    matched_meta = nsv_matched_meta;
                    rext = nsv_rext;
                }

                if(lcp > 1)[[likely]] { // nb: may occur as a rare bordercase
                    // emit reference
                    lrefs[thread_num]->emplace_back(dst, src, lcp);

                    j += matched_meta - 1; // nb: -1 because the for loop will already advance by one
                    if(rext > 0) {
                        // we have encoded characters from the following meta characters
                        ++j;
                        while(j < m && rext >= meta[parsing[j]].len) {
                            // TODO: we may have skipped additional metacharacters... but HOW??? (happens very rarely)
                            rext -= meta[parsing[j]].len;
                            ++j;
                        }
                    }
                }
            }
        }

        size_t max_refs = 0;
        for(size_t x = 0; x < num_threads; x++) {
            max_refs += lrefs[x]->size();
            lrefs[x]->shrink_to_fit();
        }

        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
            std::cout << "\tpreliminary refs: " << max_refs << std::endl;
        }

        if constexpr(debug_) {
            std::cout << "emit ... ";
            std::cout.flush();
            phase.start();
        }

        // emit
        size_t gap_total = 0, gap_num = 0;
        size_t copy_total = 0, copy_num = 0;
        {
            size_t cur_gap = 0;
            size_t cur_gap_begin = 0;
            std::string gap_buffer; // nb: only for streaming
            auto emit_current_gap = [&](){
                if(cur_gap > 0) {
                    if constexpr(has_text_access) {
                        for(size_t c = 0; c < cur_gap; c++) {
                            emit_literal(t[cur_gap_begin + c]);
                        }
                    } else {
                        in.seekg(cur_gap_begin, std::ios_base::beg);

                        gap_buffer.resize_and_overwrite(cur_gap, [&](char* data, size_t num){
                            in.read(data, num);
                            return in.gcount();
                        });

                        for(auto c : gap_buffer) {
                            emit_literal(c);
                        }
                    }

                    ++gap_num;
                    cur_gap = 0;
                }
            };

            size_t x = 0;
            size_t i = 0;
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

            if constexpr(!has_text_access) {
                gap_buffer.reserve(2 * s); // the expected gap length would be s, reserve a little more
            }

            while(i < n) {
                while(has_next_ref() && i > next_ref().end()) {
                    advance_ref();
                }

                if(has_next_ref() && i >= next_ref().beg) {
                    // emit gap up until reference
                    emit_current_gap();

                    // emit reference
                    auto const& ref = next_ref();
                    auto const d = i - ref.beg;
                    emit_reference(lz77::Factor(ref.src, ref.len - d));

                    copy_total += ref.len - d;
                    ++copy_num;

                    i = ref.end() + 1;
                    cur_gap_begin = i;

                    // advance
                    advance_ref();
                } else {
                    ++gap_total;
                    ++cur_gap;
                    ++i;
                }
            }
            emit_current_gap();
        }
        
        if constexpr(debug_) {
            phase.stop();
            std::cout << "(" << (size_t)phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;

            double const avg_gap_len = double(gap_total) / double(gap_num);
            double const avg_copy_len = double(copy_total) / double(copy_num);
            std::cout << std::endl;
            std::cout << "average gap length: " << avg_gap_len << " (of " << gap_num << " gaps with total length " << gap_total << ")" << std::endl;
            std::cout << "average copy length: " << avg_copy_len << " (of " << copy_num << " copy phrases with total length " << copy_total << ")" << std::endl;
        }
    }

public:
    SampledLPFFactorizer(size_t sampling, size_t fp_window)
        : sampling_(sampling), fp_window_(fp_window) {
    }

    template<typename InputStream>
    void factorize(InputStream& in, size_t const n, size_t const block_size, lz77::EmitFunction emit_literal, lz77::EmitFunction emit_reference) {
        factorize<InputStream, false>(in, block_size, std::string_view(), n, emit_literal, emit_reference);
    }

    template<std::contiguous_iterator Input>
    requires (sizeof(std::iter_value_t<Input>) == 1)
    void factorize(Input begin, Input const& end, lz77::EmitFunction emit_literal, lz77::EmitFunction emit_reference) {
        std::string_view const t(begin, end);
        size_t const n = t.size();

        struct NoStream {
            // define types for overlapping blocks to compile
            using char_type = char;
            using int_type = int;
        };
        NoStream no_stream;
        factorize<NoStream, true>(no_stream, 0, t, n, emit_literal, emit_reference);
    }

    template<std::contiguous_iterator Input, std::output_iterator<lz77::Factor> Output>
    requires (sizeof(std::iter_value_t<Input>) == 1)
    void factorize(Input begin, Input const& end, Output out) {
        auto emit = [&](lz77::Factor f){ *out++ = f; };
        factorize(begin, end, emit, emit);
    }
};

}
