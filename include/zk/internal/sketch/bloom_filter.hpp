#pragma once

#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>

#include <word_packing.hpp>

namespace zk::internal {

template<std::unsigned_integral T, size_t num_hashes_ = 2>
class BloomFilter {
private:
    static constexpr uintmax_t shuf_[] = {
        0x319ED645ULL,
        0xD645319EULL,
        0x19ED6453ULL,
        0x645319EDULL,
        0x9ED64531ULL,
        0xED645319ULL,
        0x45319ED6ULL,
        0x5319ED64ULL,
    };
    static constexpr uintmax_t prime_[] = {
        (1ULL << 32) - 17,
        (1ULL << 32) - 65,
        (1ULL << 32) - 99,
        (1ULL << 32) - 107,
        (1ULL << 32) - 135,
        (1ULL << 32) - 153,
        (1ULL << 32) - 185,
        (1ULL << 32) - 209,
    };
    static_assert(num_hashes_ > 0);
    static_assert(num_hashes_ < 8);

    size_t size_;
    uintmax_t mask_;
    std::unique_ptr<uintmax_t[]> bits_;

    inline uintmax_t hash(T const x, size_t const i) const {
        return ((x * shuf_[i]) % prime_[i]) & mask_;
    }

    inline size_t num_packs_required() const {
        return word_packing::num_packs_required<uintmax_t>(size_, 1);
    }

public:
    BloomFilter() : mask_(0) {
    }

    BloomFilter(size_t const size) : size_(std::bit_ceil(size)), mask_(size_ - 1) {
        auto const num_packs = num_packs_required();
        bits_ = std::make_unique<uintmax_t[]>(num_packs);
    }

    BloomFilter(BloomFilter&&) = default;
    BloomFilter& operator=(BloomFilter&&) = default;

    BloomFilter(BloomFilter const&) = delete;
    BloomFilter& operator=(BloomFilter const&) = delete;

    inline void emplace(T const x) {
        for(size_t i = 0; i < num_hashes_; i++) {
            auto const h_i = hash(x, i);
            
            auto bits = word_packing::bit_accessor(bits_.get());
            bits[h_i] = 1;
        }
    }

    inline bool lookup(T const x) const {
        for(size_t i = 0; i < num_hashes_; i++) {
            auto bits = word_packing::bit_accessor(bits_.get());

            auto const h_i = hash(x, i);
            if(!bits[h_i]) return false;
        }
        return true;
    }

    inline BloomFilter& operator|=(BloomFilter const& other) {
        assert(size_ == other.size_);

        auto const n = num_packs_required();
        for(size_t i = 0; i < n; i++) {
            bits_[i] |= other.bits_[i];
        }
        return *this;
    }

    inline size_t count_set_bits() const {
        auto const n = num_packs_required();

        size_t count = 0;
        for(size_t i = 0; i < n; i++) {
            count += std::popcount(bits_[i]);
        }
        return count;
    }

    size_t size() const {
        return size_;
    }
};

}
