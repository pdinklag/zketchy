#pragma once

#include <cassert>
#include <libsais.h>
#include <string>
#include <tuple>

#include "binary_rank.hpp"

namespace zk::internal {

class FMIndex {
private:
    static constexpr size_t SIGMA = 256;

    size_t n_;

    size_t c_array_[SIGMA+1];
    uint64_t effective_alphabet_[SIGMA / 64];
    uint8_t effective_mapping_[SIGMA];

    std::unique_ptr<int32_t[]> sa_;
    size_t wt_levels_;
    std::unique_ptr<std::unique_ptr<uint64_t[]>[]> wt_data_;
    std::unique_ptr<BinaryRank[]> wt_rank_;

    bool try_map(char const c, uint8_t& out) const {
        uint8_t const x = c;

        auto effective_alphabet = word_packing::bit_accessor(effective_alphabet_);
        if(effective_alphabet[x]) {
            out = effective_mapping_[x];
            return true;
        } else {
            return false;
        }
    }

    size_t wt_rank(uint8_t const x, size_t const pos) const {
        ssize_t i = pos;
        size_t vl = 0;
        size_t vr = n_;
        for(size_t l = 0; l < wt_levels_; l++)  {
            auto const& level_rank = wt_rank_[l];
            auto const b = (x >> (wt_levels_ - 1 - l)) & 1;
            if(b) {
                // 1-bit, navigate right
                auto const v_offset0 = vl > 0 ? level_rank.rank0(vl - 1) : 0;
                auto const v_offset1 = vl - v_offset0; // nb: needed later, but we modify vl so we compute it now
                auto const v_left_child_size = level_rank.rank0(vr - 1) - v_offset0;
                vl += v_left_child_size;
                assert(vl <= vr);
                
                auto const i_rank1 = level_rank.rank1(i);
                if(i_rank1 == 0) return 0;
                assert(i_rank1 >= v_offset1);
                assert(i_rank1 <= i + 1);
                auto const i_offs_in_v = i_rank1 - 1 - v_offset1;
                i = vl + i_offs_in_v;
                assert(i <= vr);
            } else {
                // 0-bit, navigate left
                auto const v_offset0 = vl > 0 ? level_rank.rank0(vl - 1) : 0;
                auto const v_left_child_size = level_rank.rank0(vr - 1) - v_offset0;
                vr = vl + v_left_child_size;
                assert(vr <= n_);
                assert(vl <= vr);

                auto const i_rank0 = level_rank.rank0(i);
                if(i_rank0 == 0) return 0;
                assert(i_rank0 >= v_offset0);
                assert(i_rank0 <= i + 1);
                auto const i_offs_in_v = i_rank0 - 1 - v_offset0;
                i = vl + i_offs_in_v;
                assert(i < vr);
            }
        }
        return i + 1 - vl;
    }

public:
    FMIndex(std::string_view const& s) : n_(s.length() + 1) {
        // compute suffix array of reverse string, then the BWT
        sa_ = std::make_unique<int32_t[]>(n_);
        std::string bwtr;

        {
            // reverse
            std::string rs;
            rs.reserve(n_);
            std::reverse_copy(s.begin(), s.end(), std::back_inserter(rs));
            rs.push_back(0);
            libsais((uint8_t const*)rs.data(), sa_.get(), n_, 0, nullptr);

            // compute BWT
            bwtr.reserve(n_);
            for(size_t i = 0; i < n_; i++) {
                bwtr.push_back(sa_[i] == 0 ? rs[n_-1] : rs[sa_[i]-1]);
            }
        }

        // compute histogram
        size_t occ[SIGMA];
        for(size_t x = 0; x < SIGMA; x++) {
            occ[x] = 0;
        }

        for(auto c : bwtr) {
            ++occ[uint8_t(c)];
        }

        if(occ[0] != 1) std::abort(); // expect exactly one dollar

        // compute C array, and effective alphabet
        size_t sigma = 0;
        auto effective_alphabet = word_packing::bit_accessor(effective_alphabet_);
        size_t hist[SIGMA];

        for(size_t x = 0; x < SIGMA; x++) {
            if(occ[x] > 0) {
                hist[sigma] = occ[x];

                effective_alphabet[x] = 1;
                effective_mapping_[x] = sigma++;
            } else {
                effective_alphabet[x] = 0;
            }
        }

        // compute C array
        c_array_[0] = 0;
        for(size_t x = 1; x <= sigma; x++) {
            c_array_[x] = c_array_[x-1] + hist[x-1];
        }

        // compute wavelet tree bottom-up
        wt_levels_ = std::bit_width(sigma - 1);
        wt_data_ = std::make_unique<std::unique_ptr<uint64_t[]>[]>(wt_levels_);

        size_t start[SIGMA];
        for(size_t l1 = wt_levels_; l1 > 0; l1--) {
            size_t const l = l1 - 1;
            size_t const level_nodes = (1ULL << l) - 1;

            hist[0] += hist[1];
            start[0] = 0;
            for(size_t v = 1; v < level_nodes; v++) {
                hist[v] = hist[2 * v] + hist[2 * v + 1];
                start[v] = start[v-1] + hist[v - 1];
            }

            auto bits = word_packing::bit_alloc(wt_data_[l], n_);
            auto const rsh = wt_levels_ - l;
            for(size_t i = 0; i < n_; i++) {
                auto const x = effective_mapping_[uint8_t(bwtr[i])];
                auto const bit_prefix = x >> rsh;
                auto const pos = start[bit_prefix]++;
                bits[pos] = (x >> (rsh - 1)) & 1;
            }
        }

        // compute rank data structures
        wt_rank_ = std::make_unique<BinaryRank[]>(wt_levels_);
        for(size_t l = 0; l < wt_levels_; l++) {
            wt_rank_[l] = BinaryRank(wt_data_[l].get(), n_);
        }
    }

    // find the longest possible prefix of pattern p that still occurs in the input
    // returns the position of one occurrence, the number of total occurrences and the matched length
    std::tuple<size_t, size_t, size_t> find(std::string_view const& p) const {
        // backward search for p
        uint8_t x;
        if(try_map(p[0], x)) {
            size_t d = 1;
            auto l = c_array_[x];
            auto r = c_array_[x + 1] - 1;
            assert(r < n_);
            assert(l <= r);

            for(size_t i = 1; i < p.length(); i++) {
                if(try_map(p[i], x)) {
                    // try and refine the suffix array interval
                    auto const target = c_array_[x];

                    auto const off_l = (l > 0) ? wt_rank(x, l - 1) : 0;
                    auto const new_l = target + off_l;
                    assert(new_l < n_);

                    auto const off_r = wt_rank(x, r);
                    auto const new_r = target + off_r - 1;
                    assert(new_r <= n_);

                    if(new_l <= new_r) {
                        // matched another character
                        l = new_l;
                        r = new_r;
                        ++d;
                    } else {
                        // done matching
                        break;
                    }
                } else {
                    // next character does not exist in input, we're done matching
                    break;
                }
            }

            // reconstruct the position of an occurrence from [l, r]
            auto const pos = (n_ - 1) - sa_[l] - d; // -1 for the sentinel
            return {pos, r - l + 1, d};
        } else {
            // no occurrence of even the first character of p
            return {-1, 0, 0};
        }
    }
};

}
