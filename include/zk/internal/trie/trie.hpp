#pragma once

#include <algorithm>
#include <bit>
#include <cassert>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <string>

namespace zk::internal {

template<typename Node>
concept TrieNode =
    requires {
        typename Node::Character;
        typename Node::Index;
    }
    && std::constructible_from<Node> // default
    && std::constructible_from<Node, typename Node::Index, typename Node::Character> // parent / label
    && requires(Node node) {
        { node.children };
        { node.inlabel };
        { node.parent };
    } && requires(Node const& node) {
        { node.size() } -> std::unsigned_integral;
        { node.is_leaf() } -> std::same_as<bool>;
    };

template<TrieNode Node>
class Trie {
private:
    using Character = typename Node::Character;
    using NodeIndex = typename Node::Index;

public:
    static constexpr NodeIndex ROOT = 0;
    static constexpr NodeIndex NIL = std::numeric_limits<NodeIndex>::max(); // nb: used to denote orphans and is only ever used if orphans are allowed

    static constexpr bool is_nil(NodeIndex const x) { return x == NIL; }
    static constexpr bool is_valid(NodeIndex const x) { return !is_nil(x); }
    static constexpr bool is_root(NodeIndex const x) { return x == ROOT; }
    static constexpr bool is_root_or_nil(NodeIndex const x) { return is_root(x) || is_nil(x); }
    static constexpr bool is_valid_nonroot(NodeIndex x) { return !is_root_or_nil(x); }

private:
    NodeIndex capacity_;
    NodeIndex size_;

    std::unique_ptr<Node[]> nodes_;

    #ifndef NDEBUG
    bool is_child_of(NodeIndex const node, NodeIndex const parent) const {
        return nodes_[parent].children.contains(node);
    }
    #endif

public:
    Trie() : capacity_(0), size_(0) {
    }

    Trie(NodeIndex const capacity) : capacity_(capacity), size_(1) {
        nodes_ = std::make_unique<Node[]>(capacity_);
        for(NodeIndex i = 0; i < capacity_; i++) {
            nodes_[i] = Node(NIL, {});
        }
    }

    Trie(Trie&&) = default;
    Trie& operator=(Trie&&) = default;

    Trie(Trie const& other) {
        *this = other;
    }

    Trie& operator=(Trie const& other) {
        capacity_ = other.capacity_;
        size_ = other.size_;
        nodes_ = std::make_unique<Node[]>(capacity_);
        for(NodeIndex i = 0; i < capacity_; i++) {
            nodes_[i] = other.nodes_[i];
        }
        return *this;
    }

    void fill() {
        size_ = capacity_;
    }

    Node& insert_child(NodeIndex const node, NodeIndex const parent, Character const label) {
        assert(is_valid_nonroot(node));
        assert(is_valid(parent));
        assert(node < capacity_);
        assert(node != parent);
        assert(!is_child_of(node, parent));

        NodeIndex discard;
        assert(!try_get_child(parent, label, discard));

        nodes_[parent].children.insert(label, node);
        nodes_[node].parent = parent;
        nodes_[node].inlabel = label;
        nodes_[node].children.clear();

        assert(is_leaf(node));

        return nodes_[node];
    }

    // extract node from trie and return parent
    NodeIndex extract(NodeIndex const node) {
        assert(!is_root(node));
        assert(is_leaf(node)); // cannot extract an inner node

        auto& v = nodes_[node];

        // orphan node
        auto const parent = v.parent;
        if(is_valid(parent)) {
            auto const label = v.inlabel;
            assert(is_child_of(node, parent));
            nodes_[parent].children.remove(label);

            NodeIndex discard;
            assert(!try_get_child(parent, label, discard));
            assert(!is_child_of(node, parent));
        }
        v.parent = NIL;

        return parent;
    }

    inline NodeIndex new_node() {
        assert(size_ < capacity_);
        return size_++;
    }

    inline  bool try_get_child(NodeIndex const node, Character const label, NodeIndex& out_child) const {
        return nodes_[node].children.try_get(label, out_child);
    }

    NodeIndex child_count(NodeIndex const node) const {
        return nodes_[node].children.size();
    }

    NodeIndex index_in_parent(NodeIndex const node) const {
        auto const parent = nodes_[node].parent;
        auto const label = nodes_[node].inlabel;
        return nodes_[parent].children.find(label);
    }

    NodeIndex depth(NodeIndex const node) const {
        NodeIndex d = 0;
        auto v = node;
        while(is_valid_nonroot(v)) {
            ++d;
            v = parent(v);
        }
        return d;
    }

    bool is_leaf(NodeIndex const node) const {
        return nodes_[node].is_leaf();
    }

    NodeIndex root() const {
        return ROOT;
    }

    NodeIndex parent(NodeIndex const node) const {
        return nodes_[node].parent;
    }

    bool full() const { 
        return size_ == capacity_;
    }

    NodeIndex size() const {
        return size_;
    }

    Node& node(NodeIndex const v) {
        return nodes_[v];
    }

    Node const& node(NodeIndex const v) const {
        return nodes_[v];
    }

    Node* nodes() {
        return nodes_.get();
    }

    auto const& children_of(NodeIndex const v) const {
        return nodes_[v].children;
    }

    size_t spell_reverse(NodeIndex const node, Character* buffer) const {
        size_t d = 0;
        auto v = node;
        while(is_valid_nonroot(v)) {
            buffer[d++] = nodes_[v].inlabel;
            v = nodes_[v].parent;
        }
        return d;
    }

    size_t spell(NodeIndex const node, Character* buffer) const {
        auto const d = spell_reverse(node, buffer);
        std::reverse(buffer, buffer +  d);
        return d;
    }

    size_t mem_size() const {
        auto sz = sizeof(Trie) + capacity_ * sizeof(Node);
        for(size_t v = 0; v < capacity_; v++) {
            sz += nodes_[v].children.allocated_extra_memory();
        }
        return sz;
    }
};

}
