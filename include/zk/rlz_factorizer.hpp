#pragma once

#include <memory>
#include <string>
#include <tuple>

#include <lz77/emit_function.hpp>
#include <libsais.h>

namespace zk {

class RLZFactorizer {
private:
    static constexpr size_t MAX_REF_LEN = (1ULL << 31) - 1;

    std::string_view ref_;
    std::unique_ptr<int32_t[]> sa_;

    uint8_t access(size_t const q, size_t const d) {
        auto const i = sa_[q] + d;
        return i < ref_.length() ? ref_[i] : 0; // TODO: is this good?
    }

    // nb: assumes that min <= x <= max
    size_t lmost_or_succ(size_t const l, size_t const r, uint8_t const x, size_t const d) {
        ssize_t a = ssize_t(l)- 1;
        size_t b = r;
        while(a + 1 < b) {
            auto const q = (a + b) / 2;
            if(access(q, d) < x) {
                a = q;
            } else {
                b = q;
            }
        }
        return b;
    }

    // nb: assumes that min <= x <= max
    size_t rmost_or_pred(size_t const l, size_t const r, uint8_t const x, size_t const d) {
        size_t a = l;
        size_t b = r + 1;
        while(a + 1 < b) {
            auto const q = (a + b) / 2;
            if(access(q, d) <= x) {
                a = q;
            } else {
                b = q;
            }
        }
        return a;
    }

    std::pair<size_t, size_t> step(size_t const l, size_t const r, uint8_t const x, size_t const d) {
        auto const min = access(l, d);
        if(x < min)[[unlikely]] return {r+1, l};
        auto const max = access(r, d);
        if(x > max)[[unlikely]] return {r+1, l};

        auto const ll = lmost_or_succ(l, r, x, d);
        auto const rr = rmost_or_pred(l, r, x, d);
        return {ll, rr};
    }

public:
    RLZFactorizer(std::string_view const& r) : ref_(r) {
        auto const n = ref_.length();
        if(n > MAX_REF_LEN)[[unlikely]] {
            std::abort();
        }

        // compute index for reference
        sa_ = std::make_unique<int32_t[]>(n);
        libsais((uint8_t const*)ref_.data(), sa_.get(), n, 0, nullptr);
    }

    template<std::input_iterator Input>
    requires (sizeof(std::iter_value_t<Input>) == 1)
    void factorize(Input it, Input const& end, lz77::EmitFunction emit_literal, lz77::EmitFunction emit_reference) {
        while(it != end) {
            // match in reference
            size_t d = 0;
            size_t l = 0;
            size_t r = ref_.length() - 1;
            char last_literal;

            while(l <= r && it != end) {
                auto [new_l, new_r] = step(l, r, *it, d);
                if(new_l <= new_r) {
                    l = new_l;
                    r = new_r;
                    ++d;
                    last_literal = *it++;
                } else {
                    break;
                }
            }

            // we matched d characters
            if(d > 1) {
                emit_reference(lz77::Factor(sa_[l], d));
            } else if(d == 1) {
                emit_literal(lz77::Factor(last_literal));
            } else {
                // d == 0
                emit_literal(lz77::Factor(*it));
                ++it; // nb: make sure we actually advance
            }
        }
    }
};

}
