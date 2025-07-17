#pragma once

#include "trie_edge_array.hpp"

namespace zk::internal {

template<std::unsigned_integral NodeIndex = uint32_t>
struct BasicTrieNode {
    static constexpr NodeIndex NIL = std::numeric_limits<NodeIndex>::max(); // nb: used to denote orphans and is only ever used if orphans are allowed

    using Character = char;
    using Index = NodeIndex;

    using ChildArray = TrieEdgeArray<Character, Index>;

    ChildArray children;
    Character inlabel;
    NodeIndex parent;

    BasicTrieNode(NodeIndex const _parent, Character const _inlabel) : parent(_parent), inlabel(_inlabel) {
    }
    
    BasicTrieNode() : BasicTrieNode(0, 0) {
    }

    inline size_t size() const {
        return children.size();
    }

    inline bool is_leaf() const {
        return size() == 0;
    }
} __attribute__((packed));

}
