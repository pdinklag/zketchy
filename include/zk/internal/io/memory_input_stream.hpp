#pragma once

#include <algorithm>
#include <iostream>

namespace zk::internal {

class MemoryInputStream {
public:
    using char_type = char;

private:
    char const* data_;
    size_t const size_;

    size_t pos_;
    size_t gcount_;

    inline void reset() {
        pos_ = 0;
        gcount_ = 0;
    }

public:
    MemoryInputStream(char const* data, size_t const size) : data_(data), size_(size) {
        reset();
    }

    inline MemoryInputStream& seekg(ssize_t off, std::ios_base::seekdir dir) {
        size_t new_pos;
        switch(dir) {
            case std::ios::beg:
                new_pos = off;
                break;
            
            case std::ios::cur:
                new_pos = pos_ + off;
                break;
            
            case std::ios::end:
                new_pos = size_ + off;
                break;
        }

        pos_ = new_pos;
        gcount_ = 0;
        return *this;
    }

    inline MemoryInputStream& read(char* outp, size_t const num) {
        gcount_ = std::min(num, size_ - pos_);
        for(size_t i = 0; i < gcount_; i++) {
            *outp++ = data_[pos_++];
        }
        return *this;
    }

    inline size_t gcount() const { return gcount_; }
};

}
