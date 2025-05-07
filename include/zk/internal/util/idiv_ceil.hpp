#pragma once

#include <cstdint>

namespace zk::internal {

constexpr uintmax_t idiv_ceil(uintmax_t const a, uintmax_t const b) {
    return ((a + b) - 1ULL) / b;
}

}
