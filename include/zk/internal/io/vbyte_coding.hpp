/**
 * internal/io/vbyte_coding.hpp
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

#ifndef _ZKETCHY_INTERNAL_IO_VBYTE_CODING_HPP
#define _ZKETCHY_INTERNAL_IO_VBYTE_CODING_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include <iopp/concepts.hpp>

namespace zk::internal {

template<typename _Char>
struct VByteTraits {
    using Char = _Char;
    using UChar = std::make_unsigned_t<Char>;

    static constexpr uintmax_t mask = std::numeric_limits<UChar>::max() >> 1;
    static constexpr uintmax_t rsh = std::numeric_limits<UChar>::digits - 1;
};

template<iopp::STLOutputStreamLike OutputStream>
inline size_t encode_vbyte(OutputStream& out, uintmax_t x) {
    using Traits = VByteTraits<typename OutputStream::char_type>;
    
    size_t written = 0;
    do {        
        auto byte = typename Traits::UChar(x & Traits::mask);

        x >>= Traits::rsh;
        if(x) byte |= ~Traits::mask;

        out.put(byte);
        ++written;
    } while(x);

    return written;
}

template<iopp::STLInputStreamLike InputStream>
inline uintmax_t decode_vbyte(InputStream& in) {
    using Traits = VByteTraits<typename InputStream::char_type>;

    uintmax_t x = 0;
    size_t lsh = 0;
    bool has_next;
    do {
        auto const byte = typename Traits::UChar(in.get());
        if(!in.good()) break;

        x |= (byte & Traits::mask) << lsh;
        lsh += Traits::rsh;

        has_next = (byte >> Traits::rsh) != 0;
    } while(has_next);

    return x;
}

}

#endif
