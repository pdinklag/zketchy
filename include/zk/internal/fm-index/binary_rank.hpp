#pragma once

#include <bit>

#include <word_packing.hpp>

#include "../util/idiv_ceil.hpp"

namespace zk::internal {

class BinaryRank {
public:
    static constexpr uint8_t popcount_prefix(const uint64_t v, const uint8_t x) {
        return std::popcount(v & (UINT64_MAX >> (~x & 63ULL))); // 63 - x
    }

private:
    static constexpr size_t SUPERBLOCK_WIDTH = 12;
    static constexpr size_t SUPERBLOCK_SIZE = 1ULL << SUPERBLOCK_WIDTH;
    static constexpr size_t BLOCKS_PER_SUPERBLOCK = SUPERBLOCK_SIZE >> 6ULL;

    uint64_t const* data_;
    std::unique_ptr<uint64_t[]> block_ranks_;
    std::unique_ptr<uint64_t[]> supblock_ranks;

public:
    BinaryRank() : data_(nullptr) {
    }

    BinaryRank(uint64_t const* data, size_t const n) : data_(data) {
        size_t const num_blocks = idiv_ceil(n, 64);
        auto block_ranks = word_packing::alloc(block_ranks_, idiv_ceil(n, 64), SUPERBLOCK_WIDTH);
        supblock_ranks = std::make_unique<uint64_t[]>(idiv_ceil(n, SUPERBLOCK_SIZE));

        // construct
        {
            size_t rank_bv = 0; // 1-bits in whole BV
            size_t rank_sb = 0; // 1-bits in current superblock
            size_t cur_sb = 0;  // current superblock

            for(size_t i = 0; i < num_blocks; i++) {
                if(i % BLOCKS_PER_SUPERBLOCK == 0) {
                    // we reached a new superblock
                    supblock_ranks[cur_sb++] = rank_bv;
                    rank_sb = 0;
                }
                
                block_ranks[i] = rank_sb;

                const auto rank_b = std::popcount(data[i]);
                rank_sb += rank_b;
                rank_bv += rank_b;
            }
        }
    }

    BinaryRank(BinaryRank&&) = default;
    BinaryRank& operator=(BinaryRank&&) = default;

    BinaryRank(BinaryRank const&) = delete;
    BinaryRank& operator=(BinaryRank const&) = delete;

    inline size_t rank1(size_t const x) const {
        auto block_ranks = word_packing::accessor(block_ranks_.get(), SUPERBLOCK_WIDTH);

        const size_t r_sb = supblock_ranks[x / SUPERBLOCK_SIZE];
        const size_t j   = x / 64;
        const size_t r_b = block_ranks[j];
        
        return r_sb + r_b + popcount_prefix(data_[j], x % 64);
    }

    inline size_t rank0(size_t const x) const {
        return x + 1 - rank1(x);
    }
};

}
