#include <algorithm>
#include <bit>
#include <memory>
#include <vector>

#include <ankerl/unordered_dense.h>
#include <fp/rk31.hpp>
#include <iopp/concepts.hpp>

#include "internal/io/overlapping_blocks.hpp"
#include "internal/util/si_iec_literals.hpp"

namespace zk {

class MinimizerSampling {
public:
    using Fingerprint = fp::RabinKarp31::Fingerprint;
    using TextIndex = size_t;
    using SampleIndex = uint32_t;
    using MinimizerIndex = uint16_t;

    struct Sample {
        TextIndex index;
        SampleIndex len;

        MinimizerIndex max_minimizers;
        MinimizerIndex num_minimizers;
        Fingerprint* minimizers;

        Sample() : index(0), len(), max_minimizers(0), num_minimizers(0), minimizers(nullptr) {
        }

        Sample(TextIndex const _index, MinimizerIndex const _max_minimizers) : index(_index), len(0), max_minimizers(_max_minimizers), num_minimizers(0), minimizers(new Fingerprint[max_minimizers]) {
        }

        Sample(Sample&& other) {
            *this = std::move(other);
        }

        Sample& operator=(Sample&& other) {
            index = other.index;
            len = other.len;
            minimizers = other.minimizers;
            max_minimizers = other.max_minimizers;
            num_minimizers = other.num_minimizers;

            other.minimizers = nullptr; // nb: make sure it's not deleted
            return *this;
        }

        Sample(Sample const&) = delete;
        Sample& operator=(Sample const&) = delete;

        ~Sample() {
            if(minimizers) {
                delete[] minimizers;
                minimizers = nullptr;
            }
        }

        void insert(Fingerprint const fp) {
            if(num_minimizers == max_minimizers && fp > minimizers[num_minimizers])[[unlikely]] return;

            MinimizerIndex rank = 0;
            while(rank < num_minimizers && minimizers[rank] < fp) {
                ++rank;
            }

            if(num_minimizers < max_minimizers) {
                ++num_minimizers;
            }

            for(size_t i = num_minimizers; i > rank; i--) {
                minimizers[i] = minimizers[i-1];
            }
            minimizers[rank] = fp;
        }

        double similarity(Sample const& other) const {
            size_t cut_size = 0;
            size_t i = 0, j = 0;
            while(i < num_minimizers && j < other.num_minimizers) {
                if(minimizers[i] < other.minimizers[j]) {
                    ++i;
                } else if(minimizers[i] > other.minimizers[j]) {
                    ++j;
                } else {
                    ++cut_size;
                    ++i;
                    ++j;
                }
            }

            auto const union_size = num_minimizers + other.num_minimizers - cut_size;
            return double(cut_size) / double(union_size);
        }

        Fingerprint hash() const {
            Fingerprint h = 17;
            for(size_t i = 0; i < num_minimizers; i++) {
                h = 31 * h + minimizers[i];
            }
            return h;
        }
    } __attribute__((packed));

    struct Cluster {
        std::vector<SampleIndex> sample_indices;

        double similarity(std::vector<Sample> const& samples, Cluster const& other) const {
            return samples[sample_indices.front()].similarity(samples[other.sample_indices.front()]);
        }
    };

    static std::vector<Cluster> compute_clusters(std::vector<Sample> const& samples) {
        std::vector<Cluster> clusters;
        auto const num_samples = samples.size();

        ankerl::unordered_dense::map<uint32_t, uint32_t> map;
        for(uint32_t i = 0; i < num_samples; i++) {
            auto const h = samples[i].hash();

            Cluster* cluster;
            auto it = map.find(h);
            if(it != map.end()) {
                cluster = &clusters[it->second];
            } else {
                map.emplace(h, clusters.size());

                clusters.emplace_back();
                cluster = &clusters.back();
            }
            cluster->sample_indices.push_back(i);
            }
        return clusters;
    }

private:
    static constexpr Fingerprint fp_base_ = 257;

    size_t sampling_;
    size_t fp_window_;
    MinimizerIndex max_minimizers_;

    SampleIndex min_sample_len_;
    Fingerprint sample_mask_;
    SampleIndex max_sample_len_;

public:
    MinimizerSampling(size_t const sampling, size_t const fp_window, MinimizerIndex const max_minimizers_per_sample)
        : sampling_(sampling),
          min_sample_len_(std::bit_floor(sampling / 3)),
          sample_mask_(2 * min_sample_len_ - 1),
          max_sample_len_(9 * min_sample_len_),
          fp_window_(fp_window),
          max_minimizers_(max_minimizers_per_sample) {
    }

    template<iopp::STLInputStreamLike InputStream>
    std::vector<Sample> sample(InputStream& in, size_t const buffer_size = 64_Mi) {
        std::vector<Sample> samples;

        // initialize I/O
        zk::internal::OverlappingBlocks block(in, buffer_size, fp_window_);
        if(block.size() < fp_window_) std::abort();
        char const* p = block.begin();

        // initialize fingerprinting
        fp::RabinKarp31 rk(fp_base_, fp_window_);
        Fingerprint fp = 0;

        size_t i = 0;
        for(; i < fp_window_; i++) {
            fp = rk.push(fp, *p++);
        }

        // initialize samples
        samples.emplace_back(0, max_minimizers_);
        samples.back().insert(fp);

        // process blocks
        size_t next_allowed = min_sample_len_;
        size_t next_mandatory = max_sample_len_;

        do {
            if(block.empty()) continue;

            while(p < block.end()) {
                // read next character
                fp = rk.roll(fp, *(p - fp_window_), *p);

                if(((fp & sample_mask_) == 0 && i >= next_allowed) || i >= next_mandatory) {
                    // finalize previous sample
                    size_t const new_beg = i - fp_window_ + 1;
                    samples.back().len = new_beg - samples.back().index;
                    samples.emplace_back(new_beg, max_minimizers_);

                    next_allowed = new_beg + min_sample_len_;
                    next_mandatory = new_beg + max_sample_len_;
                }
                samples.back().insert(fp);

                // advance character
                ++i;
                ++p;
            }

            // advance block
            block.advance();
            p = block.begin();
        } while(!block.empty());

        samples.back().len = i - samples.back().index;
        return samples;
    }

    size_t min_sample_len() const { return min_sample_len_; }
    size_t max_sample_len() const { return max_sample_len_; }
};

}
