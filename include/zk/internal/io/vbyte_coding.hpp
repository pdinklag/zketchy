#pragma once

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
