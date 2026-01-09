#pragma once

#include <cstddef>
#include <memory>
#include <iopp/concepts.hpp>

namespace zk::internal {

// utility to process an input blockwise such that the blocks overlap by a given margin
// this allows accessing "negative" positions in each block up to the given overlap
template<typename InputStream>
class OverlappingBlocks {
private:
    using Char = typename InputStream::char_type;
    static constexpr auto EOF_TOKEN = std::char_traits<Char>::eof();

    size_t block_size_;
    size_t overlap_;

    std::unique_ptr<Char[]> buffer_;
    InputStream::int_type probe_;
    size_t cur_size_;
    size_t cur_offs_;
    Char* cur_begin_;

    InputStream* stream_;

    void read_next() {
        if(first()) {
            stream_->read(cur_begin_, block_size_);
            cur_size_ = stream_->gcount();
        } else {
            *cur_begin_ = Char(probe_);
            stream_->read(cur_begin_ + 1, block_size_ - 1);
            cur_size_ = stream_->gcount() + 1;
        }
        probe_ = stream_->get();
    }

public:
    OverlappingBlocks() : block_size_(0), overlap_(0), buffer_(), probe_(EOF_TOKEN), cur_begin_(nullptr), cur_size_(0), cur_offs_(0), stream_(nullptr) {
    }

    OverlappingBlocks(size_t const block_size, size_t const overlap) : block_size_(block_size),
        overlap_(overlap),
        buffer_(std::make_unique<Char[]>(block_size_ + overlap_)),
        probe_(EOF_TOKEN),
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
    bool advance() {
        if(last()) return false;

        // slide overlap, then read next
        for(ssize_t i = 0; i < ssize_t(overlap_); i++) {
            buffer_[i] = cur_begin_[ssize_t(cur_size_) - ssize_t(overlap_) + i];
        }

        cur_offs_ += cur_size_;
        read_next();
        return cur_size_ > 0;
    }

    inline Char operator[](ssize_t const i) const{ return cur_begin_[i]; }

    inline Char const* begin() const { return cur_begin_; }
    inline Char const* end() const { return cur_begin_ + cur_size_; }
    inline size_t size() const { return cur_size_; }
    inline size_t offset() const { return cur_offs_; }

    inline bool empty() const { return cur_size_ == 0; }
    inline bool first() const { return cur_offs_ == 0; }
    inline bool last() const { return probe_ == EOF_TOKEN; }

    Char* buffer() const { return buffer_.get(); }
};

}
