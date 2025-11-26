#include <cstdint>
#include <concepts>
#include <functional>
#include <memory>
#include <type_traits>
#include <vector>

namespace zk::internal {

// a naive (!) suffix tree implementation for conceptual uses, not succinctness or performance
template<std::unsigned_integral Index = uint32_t>
class SuffixTree {
public:
    struct Node {
        Index depth;
        Index parent;
        std::vector<Index> children;

        Node() : depth(0), parent(0) {
        }
    };

private:
    using SignedIndex = std::make_signed_t<Index>;

    size_t num_leaves_;
    std::unique_ptr<Node[]> nodes_;
    size_t num_internal_;

    size_t leaf(size_t const i) const { return 1 + i; }
    size_t internal(size_t const i) const { return num_leaves_ + 1 + i; }

    size_t new_internal_node() {
        return internal(num_internal_++);
    }

public:
    SuffixTree(SignedIndex const* sa, SignedIndex* const lcp, size_t const n) : num_leaves_(n), num_internal_(0), nodes_(std::make_unique<Node[]>(2 * n)) {
        auto v = root();
        for(size_t i = 0; i < n; i++) {
            auto const x = leaf(sa[i]);
            auto const d = lcp[i];

            // std::cout << "next leaf: x=" << x << ", d=" << d << std::endl;
            while(nodes_[v].depth > d) {
                v = nodes_[v].parent;
            }

            // std::cout << "\t-> navigate up to v=" << v << ", d(v)=" << depth(v) << std::endl;

            Index y;
            if(nodes_[v].depth == d) {
                // simply insert a new leaf as a child of v
                y = v;
            } else {
                // split edge
                auto const w = nodes_[v].children.back();
                y = new_internal_node();
                // std::cout << "\t-> split edge (" << v << ", " << w << ") using new node " << y << " at depth " << d << std::endl;
                
                nodes_[y].children.push_back(w);
                nodes_[y].parent = v;
                nodes_[y].depth = d;
                nodes_[w].parent = y;
                nodes_[v].children.back() = y;
            }

            // insert a new leaf
            // std::cout << "\t-> insert leaf as child of " << y << std::endl;
            nodes_[x].parent = y;
            nodes_[x].depth = n - sa[i];
            nodes_[y].children.push_back(x);

            // advance
            v = x;
        }
    }

    Index root() const { return 0; }
    int32_t depth(Index const v) const { return nodes_[v].depth; }
    int32_t parent(Index const v) const { return nodes_[v].parent; }
    bool is_leaf(Index const v) const { return nodes_[v].children.empty(); }
    std::vector<Index> const& children(Index const v) const { return nodes_[v].children; }

    using VisitFunc = std::function<void(size_t const)>;

    void dfs(Index const v, VisitFunc visit) const {
        for(auto const x : nodes_[v].children) {
            dfs(x, visit);
        }
        visit(v);
    }

    size_t freq(Index const v) const {
        if(is_leaf(v)) {
            return 1;
        } else {
            size_t f = 0;
            for(auto const x : nodes_[v].children) {
                f += freq(x);
            }
            return f;
        }
    }

    size_t occ(Index const v) const {
        if(is_leaf(v)) {
            return v - 1;
        } else {
            return occ(nodes_[v].children[0]);
        }
    }

    size_t size() const { return 1 + num_leaves_ + num_internal_; }
};

}