/**
 * internal/util/linked_list.hpp
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

#ifndef _ZKETCHY_INTERNAL_UTIL_SI_IEC_LITERALS_HPP
#define _ZKETCHY_INTERNAL_UTIL_SI_IEC_LITERALS_HPP

#include <cstddef>

/// \brief Custom literal for SI unit "kilo", multiplying a size by <tt>10^3</tt>.
constexpr size_t operator""_K(unsigned long long s) { return s * 1'000ULL; }

/// \brief Custom literal for SI unit "mega", multiplying a size by <tt>10^6</tt>.
constexpr size_t operator""_M(unsigned long long s) { return s * 1'000'000ULL; }

/// \brief Custom literal for SI unit "giga", multiplying a size by <tt>10^9</tt>.
constexpr size_t operator""_G(unsigned long long s) { return s * 1'000'000'000ULL; }

/// \brief Custom literal for IEC unit "kibi", multiplying a size by <tt>2^10</tt>.
constexpr size_t operator""_Ki(unsigned long long s) { return s << 10ULL; }

/// \brief Custom literal for IEC unit "mebi", multiplying a size by <tt>2^20</tt>.
constexpr size_t operator""_Mi(unsigned long long s) { return s << 20ULL; }

/// \brief Custom literal for IEC unit "gibi", multiplying a size by <tt>2^30</tt>.
constexpr size_t operator""_Gi(unsigned long long s) { return s << 30ULL; }

#endif
