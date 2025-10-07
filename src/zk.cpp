#include <zk/cscore.hpp>
#include <zk/psamplz.hpp>
#include <zk/topk_lz77.hpp>

#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

#include <cmdline/program.hpp>

#include <zk/internal/benchmark.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

class ZK : public cmdline::Program {
private:
    std::string filename;
    size_t memory = 0;
    double psamplz_window_ratio = 0.1;
    double cscore_window_ratio = 0.1;
    double topk_lz77_window_ratio = 0.75;
    bool decompress = false;

public:
    ZK() : cmdline::Program("zketchy compression utility", "Compresses the input") {
        required_arg("file", filename, "The input file.");
        option('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
        option('m', "memory", memory, "The memory limit.");
        option('c', "cscore-window-ratio", cscore_window_ratio, "The ratio of memory to use for the window in cscore.");
        option('y', "psamplz-window-ratio", psamplz_window_ratio, "The ratio of memory to use for the window in psamplz.");
        option('x', "topk-lz77-window-ratio", topk_lz77_window_ratio, "The ratio of memory to use for the window in top-k-lz77.");
    }

    virtual int main() override {
        auto const n = std::filesystem::file_size(filename);
        if(decompress) {

        } else {
            iopp::FileInputStream in(filename);

            // determine parameters respecting the given memory limit
            if(memory == 0) {
                memory = n; // by default, never use more memory than the input size
            }
            std::cout << "zk on " << filename << " with memory constrained to " << memory << " bytes" << std::endl;

            auto precompress = false;
            auto strong_lz77 = false;
            size_t len_exp_min, len_exp_max;

            // probe compressibility score
            double score;

            {
                static constexpr size_t cscore_bytes_per_w = 1;
                static constexpr size_t cscore_num_minimizers = 8;
                static constexpr size_t cscore_len = 8;
                static constexpr size_t cscore_bytes_per_sample = 150;
                static constexpr size_t cscore_window_max = 64_Mi;
                static constexpr size_t cscore_sampling_min = 16_Ki;

                double const mem_window = cscore_window_ratio * double(memory);
                size_t const w = std::min(cscore_window_max, size_t(cscore_bytes_per_w * mem_window));
                size_t const mem_samples = memory - w * cscore_bytes_per_w;
                double const max_samples = double(mem_samples) / double(cscore_bytes_per_sample);
                size_t const cscore_sampling = std::max(cscore_sampling_min, size_t(double(n) / max_samples));

                std::cout << "cscore (s=" << cscore_sampling << ", w=" << w << ") ... "; std::cout.flush();

                zk::internal::MemoryTimePhase phase;
                phase.start();
                score = zk::CScore(cscore_sampling, cscore_len, cscore_num_minimizers, w).compute(in);
                phase.stop();

                std::cout << "time=" << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << ", mem=" << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ", score=" << score << std::endl;
            }

            if(score < 0.25) {
                precompress = true;
                if(score < 0.025) {
                    len_exp_min = 10;
                    len_exp_max = 12;
                } else {
                    len_exp_min = 12;
                    len_exp_max = 12;
                }
            } else if(score > 0.75) {
                strong_lz77 = true;
            }

            // psamplz
            in.seekg(0, std::ios_base::beg);

            auto const tmp_filename = precompress ? filename + ".zk_tmp" : filename;
            if(precompress) {
                static constexpr size_t psamplz_bytes_per_fp = 72;
                static constexpr size_t psamplz_window_max = 64_Mi;
                static constexpr size_t psamplz_bytes_per_w = 1;
                static constexpr size_t psamplz_bloom_scale = 6UL;

                double const mem_window = psamplz_window_ratio * double(memory);
                size_t const w = std::min(psamplz_window_max, size_t(psamplz_bytes_per_w * mem_window));
                size_t const mem_fp = memory - w * psamplz_bytes_per_w;
                size_t const psamplz_sampling = size_t(std::ceil(std::log2(size_t(double(psamplz_bytes_per_fp * n) / double(mem_fp)))));
                
                std::cout << "psamplz (s=" << psamplz_sampling << ", w=" << w << ", min=" << len_exp_min << ", max=" << len_exp_max << ") ... "; std::cout.flush();

                zk::internal::MemoryTimePhase phase;
                phase.start();
                {
                    iopp::FileOutputStream out_tmp(tmp_filename);
                    zk::PSampLZ psamplz(len_exp_min, len_exp_max, psamplz_sampling, psamplz_bloom_scale, w);
                    psamplz.compress(in, n, out_tmp);
                }
                phase.stop();
                auto const nout = std::filesystem::file_size(tmp_filename);
                auto const cratio = 100.0 * double(nout) / double(n);
                std::cout << "time=" << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << ", mem=" << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ", nout=" << nout << " (" << cratio << "%)" << std::endl;
            }

            // topk-lz77
            {
                static constexpr size_t topk_lz_sampling_weak = 4;
                static constexpr size_t topk_lz_sampling_strong = 0;
                static constexpr size_t topk_bytes_per_k = 59;
                static constexpr size_t topk_bytes_per_w_weak = 3;
                static constexpr size_t topk_bytes_per_w_strong = 12;
                static constexpr size_t topk_window_max = 2_Gi - 1; // nb: never go beyond 31-bit

                double const mem_window = topk_lz77_window_ratio * double(memory);
                size_t const sampling = strong_lz77 ? topk_lz_sampling_strong : topk_lz_sampling_weak;
                size_t const topk_bytes_per_w = strong_lz77 ? topk_bytes_per_w_strong : topk_bytes_per_w_weak;
                size_t const w = std::min(topk_window_max, size_t(mem_window / topk_bytes_per_w));
                size_t const k = size_t(double(memory - w * topk_bytes_per_w) / double(topk_bytes_per_k));

                std::cout << "topk-lz77 (k=" << k << ", w=" << w << ", s=" << sampling << ") ... "; std::cout.flush();

                zk::internal::MemoryTimePhase phase;
                auto const out_filename = filename + ".zk";
                phase.start();
                {
                    in = iopp::FileInputStream(tmp_filename);
                    iopp::FileOutputStream fout(out_filename);
                    
                    zk::TopkLZ77 topk_lz77(k, w, 1_Ki, sampling, 1, 32_Ki);
                    topk_lz77.compress(in, iopp::bitwise_output_to(fout));
                }
                phase.stop();

                auto const nout = std::filesystem::file_size(out_filename);
                auto const cratio = 100.0 * double(nout) / double(n);
                std::cout << "time=" << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << ", mem=" << phase.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ", nout=" << nout << " (" << cratio << "%)" << std::endl;
            }

            // clean up
            if(precompress) {
                std::filesystem::remove(tmp_filename);
            }
        }
        return 0;
    }
};

int main(int argc, char** argv) {
    return ZK().run(argc, argv);
}
