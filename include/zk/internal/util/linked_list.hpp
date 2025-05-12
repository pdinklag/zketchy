#pragma once

#include <algorithm>
#include <cassert>
#include <concepts>
#include <functional>
#include <iostream>
#include <memory>

namespace zk::internal {

template<typename T>
concept LinkedListItem =
    requires { typename T::Index; } &&
    requires(T const& item) {
        { item.prev() } -> std::convertible_to<typename T::Index>;
        { item.next() } -> std::convertible_to<typename T::Index>;
    } &&
    requires(T& item, typename T::Index const x) {
        { item.prev(x) };
        { item.next(x) };
    };

template<LinkedListItem T>
class LinkedList {
private:
    using Index = typename T::Index;

    static constexpr Index NIL = -1;

    Index head_;

public:
    LinkedList() : head_(NIL) {
    }

    LinkedList(LinkedList&&) = default;
    LinkedList& operator=(LinkedList&&) = default;

    LinkedList(LinkedList const& other) = default;
    LinkedList& operator=(LinkedList const& other) = default;

    inline void push_front(T* items, Index const i) {
        auto& item = items[i];

        if(head_ != NIL) {
            auto& head = items[head_];
            head.prev(i);
        }

        item.prev(NIL);
        item.next(head_);

        head_ = i;
    }

    inline void pop_front(T* items) {
        erase(head_, items);
    }

    inline void erase(T* items, Index const i) {
        assert(!empty());

        auto const& item_to_delete = items[i];
        auto const iprev = item_to_delete.prev();
        auto const inext = item_to_delete.next();

        if(iprev != NIL) items[iprev].next(inext);
        if(inext != NIL) items[inext].prev(iprev);

        if(head_ == i) {
            head_ = inext;
        }
    }

    inline void append(T* items, LinkedList<T> const& other) {
        if(empty()) {
            // trivial
            *this = other;
        } else {
            // find tail
            Index last;
            for(auto x = head_; x != NIL; last = x, x = items[x].next());

            // link tail to head of other
            auto& tail = items[last];
            tail.next(other.front());

            if(other.front() != NIL) {
                auto& link = items[other.front()];
                link.prev(last);
            }
        }
    }

    inline void clear() {
        head_ = NIL; // well ...
    }

    inline bool contains(T const* items, Index const i) const {
        for(auto cur = head_; cur != NIL; cur = items[cur].next()) {
            if(cur == i) return true;
        }
        return false;
    }

    inline Index front() const { return head_; }

    inline size_t size(T const* items) const { 
        size_t sz = 0;
        for(auto cur = head_; cur != NIL; cur = items[cur].next()) ++sz;
        return sz;
    }

    inline bool empty() const { return head_ == NIL; }
} __attribute__((packed));

}
