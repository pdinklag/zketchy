#include <cmath>

#include <alz/approximate_lz77.hpp>
#include <zk/topk_lz77.hpp>

#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

#include <cmdline/program.hpp>

#include <zk/internal/benchmark.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <zk/internal/io/vbyte_coding.hpp>

class ZK : public cmdline::Program {
private:
    std::string filename;
    size_t memory = 0;
    double alz_threshold = 0.3;
    double alz_block_ratio = 0.05;
    double topk_lz77_window_ratio = 0.75;
    bool disable_precompression = false;
    bool decompress = false;

public:
    ZK() : cmdline::Program("zketchy compression utility", "Compresses the input") {
        required_arg("file", filename, "The input file.");
        option('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
        option('m', "memory", memory, "The memory limit.");
        option('a', "alz-threshold", alz_threshold, "Do precompression only if metacharacter cardinality falls below this ratio.");
        option('x', "alz-block-ratio", alz_block_ratio, "The ratio of memory to use for the blocks in alz.");
        option('y', "topk-lz77-window-ratio", topk_lz77_window_ratio, "The ratio of memory to use for the window in topk-lz77.");
        option("disable-precompression", disable_precompression, "Completely disable precompression.");
    }

    virtual int main() override {
        auto const n = std::filesystem::file_size(filename);
        if(decompress) {

        } else {
            // determine parameters respecting the given memory limit
            if(memory == 0) {
                memory = n; // by default, never use more memory than the input size
            }
            std::cout << "zk on " << filename << " with memory constrained to " << memory << " bytes" << std::endl;

            iopp::FileInputStream in(filename);

            // alz precompression
            bool precompress = false;
            auto tmp_filename = filename + ".zk_tmp";
            if(!disable_precompression) {
                constexpr size_t alz_min_sampling = 6;
                constexpr size_t alz_max_sampling = 10;
                constexpr size_t const alz_fp_window = 10;
                constexpr size_t const alz_block_size_max = 16 * 1024 * 1024;
                size_t const alz_memory_headroom = size_t(0.8 * memory);

                bool const use_64bit = (n >= 4_Gi);
                size_t const sizeof_index = use_64bit ? 8 : 4;
                size_t const sizeof_metachar = 12 + sizeof_index;
                size_t const alz_block_size = std::min(alz_block_size_max, size_t(alz_block_ratio * memory));
                
                // first coarse estimation of the required sampling to be able to do anything
                size_t const alz_initial_sampling = std::max(double(alz_min_sampling), std::ceil(1.0 + std::log2(double(n * sizeof_metachar) / double(memory - alz_block_size))));
                size_t alz_sampling = alz_initial_sampling;
                if(alz_sampling <= alz_max_sampling) {
                    std::cout << "alz<" << (use_64bit ? "64" : "32") << "> -s " << alz_sampling << " -l " << alz_fp_window << " -w " << alz_block_size << " ..." << std::endl;

                    zk::internal::MemoryTimePhase phase;

                    phase.start();
                    {                
                        iopp::FileOutputStream out(tmp_filename);

                        auto estimate_memory = [&](size_t const parsing, size_t const card){
                            return alz_block_size + parsing * (sizeof_metachar + 2 * sizeof_index) + (1ULL << (alz_sampling + 1)) * card;
                        };

                        auto pre_parse_callback = [&](size_t parsing, size_t card){
                            size_t est_mem = estimate_memory(parsing, card);
                            double card_ratio = double(card) / double(parsing);
                            double est_mem_ratio = double(est_mem) / double(memory);
                            std::cout << "\tparsing=" << parsing;
                            std::cout << ", est_card=" << card << " (" << 100.0 * card_ratio << "%)";
                            std::cout << " -> est_mem=" << est_mem << " (" << 100.0 * est_mem_ratio << "% of limit)";
                            std::cout << std::endl;

                            if(card_ratio <= alz_threshold) {
                                precompress = true;

                                while(estimate_memory(parsing*2, card/2) < alz_memory_headroom && alz_sampling > alz_min_sampling) {
                                    --alz_sampling;
                                    card /= 2;
                                    parsing *= 2;
                                    est_mem = estimate_memory(parsing, card);
                                    precompress = false;
                                }

                                while(est_mem > alz_memory_headroom && alz_sampling <= alz_max_sampling) {
                                    ++alz_sampling;
                                    card *= 2;
                                    parsing /= 2;
                                    est_mem = estimate_memory(parsing, card);
                                    precompress = false;
                                }

                                if(alz_sampling == alz_initial_sampling - 1) {
                                    // don't waste time on shoving off a sub-percent in precompression
                                    std::cout << "\tnot retrying with -s " << alz_sampling << " by heuristic" << std::endl;
                                    alz_sampling = alz_initial_sampling;
                                    precompress = true;
                                }

                                if(!precompress) {
                                    if(alz_sampling > alz_max_sampling) {
                                        std::cout << "\tmemory limit too low" << std::endl;
                                    } else {
                                        double est_mem_ratio = double(est_mem) / double(memory);
                                        std::cout << "\tretrying with -s " << alz_sampling << " -> est_mem=" << est_mem << " (" << 100.0 * est_mem_ratio << "% of limit) ..." << std::endl;
                                    }
                                }
                            } else {
                                std::cout << "\tscore too low" << std::endl;
                                precompress = false;
                            }
                            return precompress;
                        };

                        auto emit_literal = [&](lz77::Factor f){
                            zk::internal::encode_vbyte(out, 0);
                            out.put(f.literal());
                        };

                        auto emit_copy = [&](lz77::Factor f){
                            zk::internal::encode_vbyte(out, f.len);
                            zk::internal::encode_vbyte(out, f.src);
                        };

                        auto run_alz = [&](auto alz, bool callback){
                            if(callback) {
                                alz.pre_parse_callback = pre_parse_callback;
                            }
                            alz.factorize(in, n, alz_block_size, emit_literal, emit_copy);
                        };
                        
                        // attempt #1
                        if(use_64bit) {
                            // 64-bit
                            run_alz(alz::ApproximateLZ77<uint64_t>(alz_sampling, alz_fp_window), true);
                        } else {
                            // 32-bit
                            run_alz(alz::ApproximateLZ77<uint32_t>(alz_sampling, alz_fp_window), true);
                        }

                        if(!precompress && alz_sampling != alz_initial_sampling && alz_sampling <= alz_max_sampling) {
                            // attempt #2
                            in.seekg(0, std::ios::beg);
                            if(use_64bit) {
                                // 64-bit
                                run_alz(alz::ApproximateLZ77<uint64_t>(alz_sampling, alz_fp_window), false);
                            } else {
                                // 32-bit
                                run_alz(alz::ApproximateLZ77<uint32_t>(alz_sampling, alz_fp_window), false);
                            }
                            precompress = true;
                        }
                    }
                    phase.stop();

                    auto const time = phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>();
                    auto const mem = phase.get_metric<pm::MallocCounter::MemoryPeakMetric>();
                    auto const mratio = double(mem) / double(memory);
                    std::cout << "\ttime=" << time << ", mem=" << mem << " (" << 100.0 * mratio << " % of limit)";
                    if(precompress) {
                        auto const nout = std::filesystem::file_size(tmp_filename);
                        auto const cratio = double(nout) / double(n);
                        std::cout << ", nout=" << nout << " (" << 100.0 * cratio << "% of input)";
                    }
                    std::cout << std::endl;

                    if(!precompress) {
                        std::filesystem::remove(tmp_filename);
                    }
                } else {
                    precompress = false;
                    std::cout << "skipping precompression -- memory limit too low (sampling rate would have to be " << alz_sampling << ")" << std::endl;
                }
            }

            if(!precompress) {
                tmp_filename = filename;
            }
            
            // topk-lz77
            {
                static constexpr size_t topk_lz_sampling_weak = 4;
                static constexpr size_t topk_lz_sampling_strong = 0;
                static constexpr size_t topk_bytes_per_k = 58;
                static constexpr double topk_bytes_per_w_weak = 3.75;
                static constexpr double topk_bytes_per_w_strong = 12;
                static constexpr size_t topk_window_max = 2_Gi - 1; // nb: never go beyond 31-bit

                auto const tmp_n = std::filesystem::file_size(tmp_filename);
                double const mem_window = topk_lz77_window_ratio * double(memory);
                size_t const sampling = topk_lz_sampling_strong;
                double const topk_bytes_per_w = topk_bytes_per_w_strong;
                size_t const w = std::min(tmp_n, std::min(topk_window_max, size_t(mem_window / topk_bytes_per_w)));
                size_t const k = size_t(double(memory - w * topk_bytes_per_w) / double(topk_bytes_per_k));

                std::cout << "topk-lz77 (k=" << k << ", w=" << w << ", s=" << sampling << ") ... " << std::endl;

                zk::internal::MemoryTimePhase phase;
                auto const out_filename = filename + ".zk";
                phase.start();
                size_t z;
                {
                    in = iopp::FileInputStream(tmp_filename);
                    iopp::FileOutputStream fout(out_filename);
                    
                    zk::TopkLZ77 topk_lz77(k, w, 1_Ki, sampling, 1, 32_Ki);
                    z = topk_lz77.compress(in, iopp::bitwise_output_to(fout));
                }
                phase.stop();

                auto const nout = std::filesystem::file_size(out_filename);
                auto const cratio = double(nout) / double(n);

                auto const time = phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>();
                auto const mem = phase.get_metric<pm::MallocCounter::MemoryPeakMetric>();
                auto const mratio = double(mem) / double(memory);
                std::cout << "\ttime=" << time << ", mem=" << mem << " (" << 100.0 * mratio << " % of limit)" << ", z=" << z << ", nout=" << nout << " (" << 100.0 * cratio << "% of input)" << std::endl;
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
