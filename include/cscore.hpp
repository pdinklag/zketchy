#pragma once

#include "internal/benchmark.hpp"
#include "internal/hashing/min_hash.hpp"
#include "internal/util/idiv_ceil.hpp"
#include "internal/util/si_iec_literals.hpp"

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

    size_t pattern_len_;
    size_t sampling_exp_;
    size_t block_size_;
    size_t buffer_size_;

    internal::Result result_;

public:
    CScore(size_t const pattern_len, size_t const sampling_exp, size_t const block_size = 64_Ki, size_t const buffer_size = 32_Mi)
        : pattern_len_(pattern_len),
          sampling_exp_(sampling_exp),
          block_size_(block_size),
          buffer_size_(buffer_size) {
    }

    CScore(CScore&&) = default;
    CScore& operator=(CScore&&) = default;
    CScore(CScore const&) = default;
    CScore& operator=(CScore const&) = default;

    template<iopp::STLInputStreamLike InputStream>
    double compute(InputStream& in) {
        result_ = {};
        result_.add("len", pattern_len_);
        result_.add("s", sampling_exp_);
        result_.add("block_size", block_size_);
        result_.add("buffer_size", buffer_size_);

        internal::MemoryTimePhase phase_all;
        internal::TimePhase phase_signatures, phase_similarity;
        phase_all.start();
        phase_signatures.start();
    
        RK rk(rolling_fp_base_, pattern_len_);
        size_t const num_threads = omp_get_max_threads();
    
        std::vector<Sig> block_signatures;
    
        size_t t_signatures;
        size_t t_gather;
        {
            Fingerprint const s = (1ULL << sampling_exp_) - 1;
    
            auto buffer = std::make_unique<char[]>(buffer_size_);
    
            std::unique_ptr<std::vector<Sig>> thread_sig[num_threads];
            for(size_t x = 0; x < num_threads; x++) {
                thread_sig[x] = std::make_unique<std::vector<Sig>>();
            }
    
            while(in.good()) {
                in.read(buffer.get(), buffer_size_);
    
                auto const num_read = in.gcount();
                auto const num_blocks = num_read / block_size_;
                if(num_blocks <= 0) continue;
    
                #pragma omp parallel
                {
                    size_t const x = omp_get_thread_num();
                    auto& local_sigs = *thread_sig[x];
    
                    #pragma omp for
                    for(size_t j = 0; j < num_blocks; j++) {
                        auto const* block = buffer.get() + j * block_size_;
    
                        Sig sig;
                        Fingerprint fp = 0;
    
                        for(size_t i = 0; i < pattern_len_; i++) {
                            fp = rk.push(fp, block[i]);
                        }
                        for(size_t i = pattern_len_; i < block_size_; i++) {
                            if((fp & s) == 0) sig.update(fp);
                            fp = rk.roll(fp, block[i - pattern_len_], block[i]);
                        }
                        if((fp & s) == 0) sig.update(fp);
                        
                        local_sigs.emplace_back(std::move(sig));
                    }
                }
            }
    
            // gather
            for(size_t x = 0; x < num_threads; x++) {
                for(auto const& sig : *thread_sig[x]) {
                    block_signatures.push_back(sig);
                }
            }
        }
        phase_signatures.stop();
    
        phase_similarity.start();
        double const score = all_to_all_similarity(block_signatures);
        phase_similarity.stop();
    
        phase_all.stop();
    
        result_.add("score", score);
        result_.add("t", phase_all.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("t_sig", phase_signatures.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("t_score", phase_similarity.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>());
        result_.add("mem_peak", phase_all.get_metric<pm::MallocCounter::MemoryPeakMetric>());
        return score;
    }

    internal::Result&& consume_last_result() {
        return std::move(result_);
    }
};

}
