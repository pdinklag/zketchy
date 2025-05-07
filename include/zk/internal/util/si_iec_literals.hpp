#pragma once

#include <cstddef>

/// \brief Custom literal for SI unit "kilo", multiplying a size by <tt>10^3</tt>.
constexpr size_t operator"" _K(unsigned long long s) { return s * 1'000ULL; }

/// \brief Custom literal for SI unit "mega", multiplying a size by <tt>10^6</tt>.
constexpr size_t operator"" _M(unsigned long long s) { return s * 1'000'000ULL; }

/// \brief Custom literal for SI unit "giga", multiplying a size by <tt>10^9</tt>.
constexpr size_t operator"" _G(unsigned long long s) { return s * 1'000'000'000ULL; }

/// \brief Custom literal for IEC unit "kibi", multiplying a size by <tt>2^10</tt>.
constexpr size_t operator"" _Ki(unsigned long long s) { return s << 10ULL; }

/// \brief Custom literal for IEC unit "mebi", multiplying a size by <tt>2^20</tt>.
constexpr size_t operator"" _Mi(unsigned long long s) { return s << 20ULL; }

/// \brief Custom literal for IEC unit "gibi", multiplying a size by <tt>2^30</tt>.
constexpr size_t operator"" _Gi(unsigned long long s) { return s << 30ULL; }
