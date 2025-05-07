#pragma once

#include <cstddef>
#include <memory>
#include <iopp/concepts.hpp>

namespace zk::internal {

// utility to process an input blockwise such that the blocks overlap by a given margin
// this allows accessing "negative" positions in each block up to the given overlap
template<iopp::STLInputStreamLike InputStream>
class OverlappingBlocks {
private:
    using Char = typename InputStream::char_type;

    size_t block_size_;
    size_t overlap_;

    std::unique_ptr<Char[]> buffer_;
    size_t cur_size_;
    size_t cur_offs_;
    Char* cur_begin_;

    InputStream* stream_;

    void read_next() {
        stream_->read(cur_begin_, block_size_);
        cur_size_ = stream_->gcount();
    }

public:
    OverlappingBlocks() : block_size_(0), overlap_(0), buffer_(), cur_begin_(nullptr), cur_size_(0), cur_offs_(0), stream_(nullptr) {
    }

    OverlappingBlocks(size_t const block_size, size_t const overlap) : block_size_(block_size),
        overlap_(overlap),
        buffer_(std::make_unique<Char[]>(block_size_ + overlap_)),
        cur_begin_(buffer_.get() + overlap_),
        cur_size_(block_size_ + overlap_),
        cur_offs_(0),
        stream_(nullptr) {
    }

    OverlappingBlocks(InputStream& stream, size_t const block_size, size_t const overlap) : OverlappingBlocks(block_size, overlap) {
        init(stream);
    }

    OverlappingBlocks(OverlappingBlocks&&) = default;
    OverlappingBlocks& operator=(OverlappingBlocks&&) = default;

    OverlappingBlocks(OverlappingBlocks const&) = delete;
    OverlappingBlocks& operator=(OverlappingBlocks const&) = delete;

    // initialize the given input stream
    void init(InputStream& stream) {
        stream_ = &stream;

        // initialize the initial overlap region
        for(size_t i = 0; i < overlap_; i++) {
            buffer_[i] = 0;
        }

        // read the first block
        cur_offs_ = 0;
        read_next();
    }

    // advance to the next block
    void advance() {
        // slide overlap, then read next
        for(ssize_t i = 0; i < ssize_t(overlap_); i++) {
            buffer_[i] = cur_begin_[ssize_t(cur_size_) - ssize_t(overlap_) + i];
        }

        cur_offs_ += cur_size_;
        read_next();
    }

    inline Char operator[](ssize_t const i) const{ return cur_begin_[i]; }

    inline Char const* begin() const { return cur_begin_; }
    inline Char const* end() const { return cur_begin_ + cur_size_; }
    inline size_t size() const { return cur_size_; }
    inline size_t offset() const { return cur_offs_; }

    inline bool empty() const { return cur_size_ == 0; }
    inline bool last() const { return cur_size_ < block_size_; }

    Char* buffer() const { return buffer_.get(); }
};

}
