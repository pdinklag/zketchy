/**
 * internal/trie/basic_trie_node.hpp
 * part of pdinklag/zketchy
 * 
 * MIT License
 * 
 * Copyright (c) Patrick Dinklage
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _ZKETCHY_INTERNAL_TRIE_BASIC_TRIE_NODE_HPP
#define _ZKETCHY_INTERNAL_TRIE_BASIC_TRIE_NODE_HPP

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

#endif
