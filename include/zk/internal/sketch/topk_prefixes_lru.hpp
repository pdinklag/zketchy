#pragma once

#include <map>

#include <zk/internal/trie/trie.hpp>
#include <zk/internal/trie/basic_trie_node.hpp>
#include <zk/internal/sketch/space_saving.hpp>

namespace zk::internal {

template<std::unsigned_integral TrieNodeIndex = uint32_t>
class TopKPrefixesLRU {
private:
    struct NodeData : public BasicTrieNode<TrieNodeIndex> {
        using Character = BasicTrieNode<TrieNodeIndex>::Character;
        using Index = BasicTrieNode<TrieNodeIndex>::Index;

    private:
        size_t time_; // the insert time

    public:
        NodeData() : time_(-1) {
        }

        NodeData(Index v, Character c) : BasicTrieNode<TrieNodeIndex>(v, c) {
        }

        // SpaceSavingItem
        inline size_t time() const { return time_; }
        inline void time(size_t time) { time_ = time; }
    } __attribute__((packed));

    using TrieNodeDepth = TrieNodeIndex;
    struct PQEntry {
        TrieNodeIndex node;
        size_t time;
    } __attribute__((packed));

    struct TimeCompare {
        bool operator()(PQEntry const& a, PQEntry const& b) const {
            return a.time > b.time;
        }
    };

    size_t k_;

    Trie<NodeData> trie_;
    std::map<size_t, TrieNodeIndex> map_; // time -> leaf that was inserted at that time
    size_t time_;

public:
    inline TopKPrefixesLRU() : k_(0), time_(0) {
    }

    inline TopKPrefixesLRU(size_t const k)
        : trie_(k),
          k_(k),
          time_(0) {
    }

    TopKPrefixesLRU(TopKPrefixesLRU&&) = default;
    TopKPrefixesLRU& operator=(TopKPrefixesLRU&&) = default;

    TopKPrefixesLRU(TopKPrefixesLRU const& other) = delete;
    TopKPrefixesLRU& operator=(TopKPrefixesLRU const& other) = delete;

    struct StringState {
        TrieNodeIndex len;      // length of the string
        TrieNodeIndex node;     // the string's node in the trie filter
        bool          frequent; // whether or not the string is frequent
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
            // done
            ext.frequent = true;
        } else {
            // the current prefix is non-frequent

            TrieNodeIndex v;
            if(trie_.size() == k_) {
                // extract LRU leaf from trie
                auto lru = map_.begin();
                if(lru->second == s.node) ++lru; // this may happen!

                v = lru->second;
                map_.erase(lru->first);
                
                auto const old_parent = trie_.extract(v);
                if(trie_.is_valid_nonroot(old_parent) && trie_.is_leaf(old_parent)) {
                    // old parent became a leaf, insert it into the map
                    map_.emplace(trie_.node(old_parent).time(), old_parent);
                }
            } else {
                v = trie_.new_node();
            }

            // insert new leaf
            auto& node = trie_.insert_child(v, s.node, c);
            node.time(time_);

            map_.emplace(time_, v);
            ++time_;

            // emovre parent, which is no longer a leaf
            map_.erase(trie_.node(s.node).time());

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
