#pragma once

#include <concepts>
#include <memory>

namespace zk::internal {
 
template<std::unsigned_integral T>
class UnionFind {
    std::unique_ptr<T[]> root_;
    std::unique_ptr<size_t[]> size_;

public:
    UnionFind(size_t const n) : root_(std::make_unique<T[]>(n)), size_(std::make_unique<size_t[]>(n)) {
        for(size_t i = 0; i < n; i++) {
            root_[i] = T(i);
            size_[i] = 1;
        }
    }

    T find(T const i) {
        if(root_[i] == i) return i;
        return root_[i] = find(root_[i]); // path compression
    }

    void unite(T const i, T const j){
        auto const root_i = find(i);
        auto const root_j = find(j);
        if (root_i != root_j) {
            if (size_[root_i] < size_[root_j]) {
                root_[root_i] = root_j;
            } else if (size_[root_i] > size_[root_j]) {
                root_[root_j] = root_i;
            } else {
                root_[root_j] = root_i;
                size_[root_i]++;
            }
        }
    }
};

}
