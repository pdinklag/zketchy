#pragma once

#include <algorithm>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>

namespace zk::internal {

// mantains an array of trie edges
template<std::integral Character = char, std::unsigned_integral NodeIndex = uint32_t>
class TrieEdgeArray {
private:
    static constexpr size_t capacity_for(size_t const n) {
        return (n > 0) ? std::bit_ceil(n) : 0;
    }

    using UCharacter = std::make_unsigned_t<Character>;
    using BitPack = uintmax_t;
    using Size = UCharacter;

    static constexpr size_t bits_per_pack_ = 8 * sizeof(BitPack);
    static constexpr size_t bit_pack_mask_ = bits_per_pack_ - 1;
    static constexpr size_t sigma_ = ((size_t)std::numeric_limits<UCharacter>::max() + 1);
    static_assert((sigma_ % bits_per_pack_) == 0);
    static constexpr size_t num_bit_packs_ = sigma_ / bits_per_pack_;

    struct ExternalArray {
        BitPack ind[num_bit_packs_];
        NodeIndex* links;

        #ifndef NDEBUG
        // only used for debugging
        size_t size() const {
            size_t s = 0;
            for(size_t i = 0; i < num_bit_packs_; i++) {
                s += std::popcount(ind[i]);
            }
            return s;
        }
        #endif

        inline void clear() {
            for(size_t i = 0; i < num_bit_packs_; i++) {
                ind[i] = 0;
            }
        }

        inline void set(UCharacter const i) {
            size_t const b = i / bits_per_pack_;
            size_t const j = i % bits_per_pack_;
            ind[b] |= (1ULL << j);
        }

        inline void unset(UCharacter const i) {
            size_t const b = i / bits_per_pack_;
            size_t const j = i % bits_per_pack_;
            ind[b] &= ~(1ULL << j);
        }

        inline bool get(UCharacter const i) const {
            size_t const b = i / bits_per_pack_;
            size_t const j = i % bits_per_pack_;
            return (ind[b] & (1ULL << j)) != 0;
        }
        
        inline size_t rank(UCharacter const i) const {
            assert(get(i));

            size_t r = 0;
            size_t const b = i / bits_per_pack_;
            size_t const j = i % bits_per_pack_;
            for(size_t i = 0; i < b; i++) {
                r += std::popcount(ind[i]);
            }

            BitPack const mask = std::numeric_limits<BitPack>::max() >> (std::numeric_limits<BitPack>::digits - 1 - j);
            return r + std::popcount(ind[b] & mask) - 1;
        }

        inline UCharacter select(size_t k) const {
            for(size_t i = 0; i < num_bit_packs_; i++) {
                auto x = ind[i];
                size_t rsh = 0;

                while(x) {
                    auto const j = std::countr_zero(x);
                    if(k == 0) {
                        // this is the bit we were looking for, reconstruct the character
                        return UCharacter(i * bits_per_pack_ + rsh + j);
                    } else {
                        // continue
                        x >>= j;
                        x >>= 1; // nb: j may be 63, and shifting by j+1 = 64 would be undefined
                        rsh += j + 1;
                        --k;
                    }
                }
            }
            assert(false);
            return 0;
        }
    } __attribute__((packed));

public:
    static constexpr size_t inline_size_ = sizeof(ExternalArray) / (sizeof(NodeIndex) + sizeof(Character));
    static constexpr size_t inline_align_ = sizeof(ExternalArray) - inline_size_ * (sizeof(NodeIndex) + sizeof(Character));

private:
    struct InlineArray {
        Character labels[inline_size_];
        NodeIndex links[inline_size_];
    } __attribute__((packed));

    Size size_;
    union {
        ExternalArray ext;
        InlineArray   inl;
    } data_;

    inline size_t find(Character const label) const {
        if(is_inline()) {
            NodeIndex found = 0;
            while(found < size_ && data_.inl.labels[found] != label) ++found;
            return found;
        } else {
            return data_.ext.get(label) ? data_.ext.rank(label) : size_;
        }
    }

public:
    TrieEdgeArray() : size_(0) {
        data_.inl.labels[0] = 0; // emptiness convention
        data_.ext.links = nullptr;
    }

    ~TrieEdgeArray() {
        if(!is_inline()) {
            delete[] data_.ext.links;
        }
    }

    TrieEdgeArray(TrieEdgeArray const& other) {
        *this = other;
    }

    TrieEdgeArray& operator=(TrieEdgeArray const& other) {
        clear();

        size_ = other.size_;
        data_ = other.data_;

        if(!other.is_inline()) {
            // deep copy of links
            data_.ext.links = new NodeIndex[capacity_for(size_)];
            for(size_t i = 0; i < size_; i++) {
                data_.ext.links[i] = other.data_.ext.links[i];
            }
        }
        return *this;
    }

    TrieEdgeArray(TrieEdgeArray&& other) { size_ = 0; *this = std::move(other); }
    TrieEdgeArray& operator=(TrieEdgeArray&& other) {
        // deallocate children
        clear();

        // copy data
        size_ = other.size_;
        data_ = other.data_;

        // invalidate other
        other.data_.ext.links = nullptr;
        other.size_ = 0;

        // done
        return *this;
    }

    inline void clear() {
        if(!is_inline()) {
            delete[] data_.ext.links;
            data_.ext.links = nullptr;
        }
        size_ = 0;
        data_.inl.labels[0] = 0; // emptiness convention
    }

    inline bool is_leaf() const {
        return size_ == 0 && data_.inl.labels[0] == 0;
    }

    inline bool is_inline() const {
        return is_leaf() || size_ <= inline_size_;
    }

    inline size_t size() const {
        if(size_ != 0)[[likely]] {
            return size_;
        } else {
            return data_.inl.labels[0] == 0 ? 0 : sigma_;
        }
    }

    size_t allocated_extra_memory() const {
        return is_inline() ? 0 : capacity_for(size_) * sizeof(NodeIndex);
    }

    inline NodeIndex operator[](size_t const i) const {
        if(is_inline()) {
            return data_.inl.links[i];
        } else {
            return data_.ext.links[i];
        }
    }

    inline Character label(size_t const i) const {
        if(is_inline()) {
            return data_.inl.labels[i];
        } else {
            return (Character)data_.ext.select(i);
        }
    }

    inline void sort() {
        if(is_inline()) {
            std::pair<UCharacter, NodeIndex> c[size_];
            for(size_t i = 0; i < size_; i++) {
                auto const label = (UCharacter)data_.inl.labels[i];
                auto const link = data_.inl.links[i];
                c[i] = std::make_pair(label, link);
            }
            std::sort(c, c + size_, [](auto const& a, auto const& b){ return a.first < b.first; });

            for(size_t i = 0; i < size_; i++) {
                data_.inl.labels[i] = (Character)c[i].first;
                data_.inl.links[i] = c[i].second;
            }
        } else {
            // nothing to do, children in external array are always sorted
        }
    }

    inline bool contains(NodeIndex const what) const {
        NodeIndex discard;
        return find(what, discard);
    }

    inline bool find(NodeIndex const what, NodeIndex& child_index) const {
        if(is_inline()) {
            for(NodeIndex i = 0; i < size_; i++) {
                if(data_.inl.links[i] == what) {
                    child_index = i;
                    return true;
                }
            }
        } else {
            for(NodeIndex i = 0; i < size_; i++) {
                if(data_.ext.links[i] == what) {
                    child_index = i;
                    return true;
                }
            }
        }
        return false;
    }

    void insert(Character const label, NodeIndex const link) {
        // possibly allocate slots
        // nb: because we are only keeping track of the size and assume the capacity to always be its hyperceil,
        //     we may re-allocate here even though it would not be necessary
        //     however, if there were many removals, this may also (unknowingly) actually shrink the capacity
        //     in any event: it's not a bug, it's a feature!
        if(size_ == inline_size_ || (!is_inline() && size_ == capacity_for(size_))) {
            auto* new_links = new NodeIndex[capacity_for(size_ + 1)];
            if(is_inline()) {
                assert(size_ == inline_size_);

                // becoming large
                Character old_labels[inline_size_];

                // we need to make sure now that the children are transferred in the order of their labels
                std::pair<UCharacter, NodeIndex> old_child_info[inline_size_];
                for(size_t j = 0; j < size_; j++) {
                    auto const child_j = data_.inl.links[j];
                    old_child_info[j] = { data_.inl.labels[j], child_j };
                }
                std::sort(old_child_info, old_child_info + size_, [](auto const& a, auto const& b){ return a.first < b.first; });

                // clear all bits first
                data_.ext.clear();

                // set bits that need to be set and write children
                for(size_t j = 0; j < size_; j++) {
                    data_.ext.set(old_child_info[j].first);
                    new_links[j] = old_child_info[j].second;
                }
            } else {
                // staying large
                std::copy(data_.ext.links, data_.ext.links + size_, new_links);
                delete[] data_.ext.links;
            }
            data_.ext.links = new_links;
        }
        
        // insert
        ++size_;
        if(is_inline()) {
            auto const i = size_ - 1;
            data_.inl.labels[i] = label;
            data_.inl.links[i] = link;
        } else {
            data_.ext.set(label);
            assert(size_ == data_.ext.size());

            auto const i = data_.ext.rank(label);
            for(size_t j = size_ - 1; j > i; j--) {
                data_.ext.links[j] = data_.ext.links[j - 1];
            }
            data_.ext.links[i] = link;
        }

        assert(contains(link));
    }

    void remove(Character const label) {
        auto const i = find(label);
        if(i < size_) {
            // remove from link array if necessary
            if(is_inline()) {
                // swap with last child if necessary
                if(size_ > 1) {
                    NodeIndex const last = size_ - 1;
                    data_.inl.labels[i] = data_.inl.labels[last];
                    data_.inl.links[i] = data_.inl.links[last];
                }
            } else {
                // remove
                data_.ext.unset(label);
                assert(data_.ext.size() == size_ - 1);

                if(size_ > 1) {
                    for(size_t j = i; j + 1 < size_; j++) {
                        data_.ext.links[j] = data_.ext.links[j+1];
                    }
                }
            }

            // possibly convert into inline node
            auto const was_inline = is_inline();
            --size_; // nb: must be done before moving data

            if(!was_inline && is_inline()) {
                assert(size_ == inline_size_);
                Character new_labels[inline_size_];
                NodeIndex new_links[inline_size_];

                {
                    size_t j = 0;
                    for(size_t c = 0; c <= std::numeric_limits<UCharacter>::max(); c++) {
                        if(data_.ext.get(c)) {
                            new_labels[j] = (Character)c;
                            new_links[j] = data_.ext.links[j];
                            ++j;
                        }
                    }
                    assert(j == inline_size_);
                }

                delete[] data_.ext.links;

                for(size_t j = 0; j < size_; j++) {
                    data_.inl.labels[j] = new_labels[j];
                    data_.inl.links[j] = new_links[j];
                }
            }

            if(size_ == 0) {
                data_.inl.labels[0] = 0; // emptiness convention
            }
        } else {
            // "this kid is not my son"
            assert(false);
        }
    }

    inline bool try_get(Character const label, NodeIndex& out_link) const {
        if(is_inline()) {
            for(NodeIndex i = 0; i < size_; i++) {
                if(data_.inl.labels[i] == label) {
                    out_link = data_.inl.links[i];
                    return true;
                }
            }
            return false;
        } else {
            if(data_.ext.get(label)) {
                out_link = data_.ext.links[data_.ext.rank(label)];
                return true;
            } else {
                return false;
            }
        }
    }

} __attribute__((packed));

}
