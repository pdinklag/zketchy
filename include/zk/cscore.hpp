#pragma once

#include "internal/benchmark.hpp"
#include "internal/sketch/min_hash.hpp"
#include "internal/util/idiv_ceil.hpp"
#include "internal/util/si_iec_literals.hpp"
#include "minimizer_sampling.hpp"

#include <iopp/concepts.hpp>
#include <fp/rk31.hpp>

#include <omp.h>
#include <memory>
#include <vector>

namespace zk {

class CScore {
private:
    using RK = fp::RabinKarp31;
    using Fingerprint = RK::Fingerprint;
    using Sig = zk::internal::MinHashSignature<32>;

    static constexpr Fingerprint rolling_fp_base_ = (1ULL << 16) - 39;

    inline static double all_to_all_similarity(std::vector<Sig> const& signatures) {
        size_t const num_threads = omp_get_max_threads();
        size_t const num_blocks = signatures.size();
        size_t const num_per_thread = zk::internal::idiv_ceil(num_blocks, num_threads);
    
        std::unique_ptr<double> thread_score[num_threads];
        for(size_t x = 0; x < num_threads; x++) {
            thread_score[x] = std::make_unique<double>(0.0);
        }
    
        for(size_t i = 0; i < num_blocks; i++) {
            auto const& sig_i = signatures[i];
    
            #pragma omp parallel
            {
                size_t const x = omp_get_thread_num();
                double& score = *thread_score[x];
    
                size_t const beg = x * num_per_thread;
                size_t const end = (x < num_threads - 1) ? (beg + num_per_thread) : num_blocks;
    
                for(size_t j = beg; j < end; j++) {
                    if(i != j) {
                        score += sig_i.similarity(signatures[j]);
                    }
                }
            }
        }
    
        double score = 0.0;
        for(size_t x = 0; x < num_threads; x++) {
            score += *thread_score[x];
        }
        score /= double(num_blocks * (num_blocks-1));
        return score;
    }

    size_t sampling_;
    size_t len_;
    size_t max_minimizers_;
    size_t buffer_size_;
    double skip_probability_;

    internal::Result result_;

public:
    CScore(size_t const sampling, size_t const len, size_t const max_minimizers, size_t const buffer_size = 64_Mi, double skip_probability = 0.0)
        : sampling_(sampling),
          len_(len),
          max_minimizers_(max_minimizers),
          buffer_size_(buffer_size),
          skip_probability_(skip_probability) {
    }

    CScore(CScore&&) = default;
    CScore& operator=(CScore&&) = default;
    CScore(CScore const&) = default;
    CScore& operator=(CScore const&) = default;

    template<iopp::STLInputStreamLike InputStream>
    double compute(InputStream& in) {
        internal::MemoryTimePhase phase;
        
        phase.start();
        MinimizerSampling s(sampling_, len_, max_minimizers_);
        auto samples = s.sample(in, buffer_size_, skip_probability_);
        auto clusters = zk::MinimizerSampling::compute_clusters(samples);
        phase.stop();

        result_ = internal::Result();
        result_.add("t", phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("mem", phase.get_metric<pm::MallocCounter::MemoryPeakMetric>());
        result_.add("num_samples", samples.size());
        result_.add("num_clusters", clusters.size());

        return double(clusters.size()) / double(samples.size());
    }

    internal::Result&& consume_last_result() {
        return std::move(result_);
    }
};

}
