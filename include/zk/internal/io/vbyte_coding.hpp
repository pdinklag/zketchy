#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include <iopp/concepts.hpp>

namespace zk::internal {

template<iopp::STLOutputStreamLike OutputStream>
inline size_t encode_vbyte(OutputStream& out, uintmax_t x) {
    using Char = typename OutputStream::char_type;
    using UChar = std::make_unsigned_t<Char>;
    
    static constexpr uintmax_t mask = std::numeric_limits<UChar>::max() >> 1;
    static constexpr uintmax_t rsh = std::numeric_limits<UChar>::digits - 1;

    size_t written = 0;
    do {        
        UChar byte = UChar(x & mask);

        x >>= rsh;
        if(x) byte |= ~mask;

        out.put(Char(byte));
        ++written;
    } while(x);

    return written;
}

}
