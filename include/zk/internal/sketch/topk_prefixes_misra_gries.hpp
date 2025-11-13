#pragma once

#include <zk/internal/trie/trie.hpp>
#include <zk/internal/trie/basic_trie_node.hpp>
#include <zk/internal/sketch/space_saving.hpp>

namespace zk::internal {

template<std::unsigned_integral TrieNodeIndex = uint32_t>
class TopKPrefixesMisraGries {
private:
    static constexpr bool gather_stats_ = true;

    struct NodeData;
    static constexpr auto NIL = SpaceSaving<NodeData>::NIL;

    struct NodeData : public BasicTrieNode<TrieNodeIndex> {
        using Character = BasicTrieNode<TrieNodeIndex>::Character;
        using Index = BasicTrieNode<TrieNodeIndex>::Index;

    private:
        TrieNodeIndex freq_; // the current frequency
        TrieNodeIndex prev_; // the previous node in frequency order
        TrieNodeIndex next_; // the next node in frequency order

    public:
        NodeData() {
        }

        NodeData(Index v, Character c) : BasicTrieNode<TrieNodeIndex>(v, c), freq_(0), prev_(NIL), next_(NIL) {
        }

        // SpaceSavingItem
        inline TrieNodeIndex freq() const { return freq_; }
        inline TrieNodeIndex prev() const { return prev_; }
        inline TrieNodeIndex next() const { return next_; }
        
        inline bool is_linked() const { return this->is_leaf(); }

        inline void freq(TrieNodeIndex const f) { freq_ = f; }
        inline void prev(TrieNodeIndex const x) { prev_ = x; }
        inline void next(TrieNodeIndex const x) { next_ = x; }
    } __attribute__((packed));

    using TrieNodeDepth = TrieNodeIndex;

    size_t k_;

    Trie<NodeData> trie_;
    SpaceSaving<NodeData> space_saving_;

    inline bool insert(TrieNodeIndex const parent, char const label, TrieNodeIndex& out_node) {
        TrieNodeIndex v;
        if(space_saving_.get_garbage(v)) {
            // recycle something from the garbage
            assert(v < k_);
            assert(v != 0);

            auto& vdata = trie_.node(v);
            assert(vdata.is_leaf());
            assert(vdata.freq() <= space_saving_.threshold());

            // extract from trie
            auto const old_parent = trie_.extract(v);

            // old parent may have become a leaf
            if(trie_.is_valid_nonroot(old_parent) && trie_.is_leaf(old_parent)) {
                space_saving_.link(old_parent);
            }

            // new parent can no longer be a leaf
            if(trie_.is_valid_nonroot(parent) && trie_.is_leaf(parent)) {
                space_saving_.unlink(parent);
            }

            // insert into trie with new parent
            trie_.insert_child(v, parent, label);

            // now simply increment
            space_saving_.increment(v);

            out_node = v;
            return true;
        } else {
            // there is no garbage to recycle
            out_node = trie_.root();
            return false;
        }
    }

public:
    inline TopKPrefixesMisraGries() : k_(0) {
    }

    inline TopKPrefixesMisraGries(size_t const k, size_t const sketch_columns, size_t const fp_window_size = 8)
        : trie_(k),
          k_(k),
          space_saving_(trie_.nodes(), 1, k_ - 1, sketch_columns - 1) {
        
        // initialize all k nodes as orphans in trie
        trie_.fill();

        // make all of them garbage (except the root)
        space_saving_.init_garbage();
    }

    TopKPrefixesMisraGries(TopKPrefixesMisraGries&&) = default;
    TopKPrefixesMisraGries& operator=(TopKPrefixesMisraGries&&) = default;

    TopKPrefixesMisraGries(TopKPrefixesMisraGries const& other) {
        *this = other;
    }

    TopKPrefixesMisraGries& operator=(TopKPrefixesMisraGries const& other) {
        k_ = other.k_;
        trie_ = other.trie_;
        space_saving_ = other.space_saving_;
        space_saving_.set_items(trie_.nodes());
        return *this;
    }

    struct StringState {
        TrieNodeIndex len;         // length of the string
        TrieNodeIndex node;        // the string's node in the trie filter
        bool          frequent;    // whether or not the string is frequent
    };

    // returns a string state for the empty string to start with
    inline StringState empty_string() const {
        StringState s;
        s.len = 0;
        s.node = trie_.root();
        s.frequent = true;
        return s;
    }

    // extends a string to the right by a new character
    inline StringState extend(StringState const& s, char const c) {
        auto const i = s.len;

        // try and find extension in trie
        StringState ext;
        ext.len = i + 1;

        bool const edge_exists = s.frequent && trie_.try_get_child(s.node, c, ext.node);
        if(edge_exists) {
            // the current prefix is frequent, increment
            space_saving_.increment(ext.node);

            // done
            ext.frequent = true;
        } else {
            // the current prefix is non-frequent

            // attempt to insert it
            if(!insert(s.node, c, ext.node)) {
                // that failed, decrement everything else in turn
                space_saving_.decrement_all();
            }

            // we dropped out of the trie, so no extension can be frequent (not even if the current prefix was inserted or swapped in)
            ext.frequent = false;
        }

        // advance
        return ext;
    }

    // read the string with the given index into the buffer
    TrieNodeDepth get(TrieNodeIndex const index, char* buffer) const {
        return trie_.spell(index, buffer);
    }

    size_t freq(TrieNodeIndex const index) const {
        return trie_.node(index).freq();
    }

    // try to find the string in the trie and report its depth and node
    TrieNodeDepth find(char const* s, size_t const max_len, TrieNodeIndex& out_node) const {
        auto v = trie_.root();
        TrieNodeDepth dv = 0;
        while(dv < max_len) {
            TrieNodeIndex u;
            if(trie_.try_get_child(v, s[dv], u)) {
                v = u;
                ++dv;
            } else {
                break;
            }
        }

        out_node = v;
        return dv;
    }
};

}
