#include <bit>

#include <zk/internal/io/overlapping_blocks.hpp>
#include <zk/internal/util/complete_graph.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <code/binary.hpp>
#include <fp/rk31.hpp>
#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>
#include <oocmd.hpp>
#include <pm/result.hpp>

class RearrangeTool : public oocmd::ConfigObject {
private:
    using Fingerprint = fp::RabinKarp31::Fingerprint;
    static constexpr Fingerprint fp_base_ = 257;

    size_t sampling = 16_Ki;
    size_t len = 8;
    size_t max_minimizers = 64;
    size_t buffer_size = 64_Mi;
    size_t prefix = SIZE_MAX;
    bool verbose = false;
    std::string output_filename;

    struct Sample {
        size_t index;
        uint32_t len;

        std::unique_ptr<Fingerprint[]> minimizers;
        uint32_t num_minimizers;

        Sample(size_t const _index, size_t const max_minimizers) : index(_index), len(0), minimizers(std::make_unique<Fingerprint[]>(max_minimizers)), num_minimizers(0) {
        }

        void insert(Fingerprint const fp, size_t const max) {
            if(num_minimizers == max && fp > minimizers[num_minimizers])[[unlikely]] return;

            size_t rank = 0;
            while(rank < num_minimizers && minimizers[rank] < fp) {
                ++rank;
            }

            if(num_minimizers < max) {
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

            size_t const union_size = num_minimizers + other.num_minimizers - cut_size;
            return double(cut_size) / double(union_size);
        }
    };

    struct Edge {
        uint32_t target;
        float cost;
    } __attribute__((packed));

public:
    RearrangeTool() : oocmd::ConfigObject("Rearrange", "Rearrange blocks of the input according to their similarity") {
        param('s', "sample", sampling, "The sampling rate.");
        param('l', "len", len, "The pattern length.");
        param('m', "minimizers", max_minimizers, "The maximum number of minimizers per sample.");
        param('w', "buffer_size", buffer_size, "The buffer size.");
        param('p', "prefix", prefix, "Process only this prefix of the input file.");
        param('v', "verbose", verbose, "Print a lot of info.");
        param('o', "out", output_filename, "The output filename.");
    }

    int run(oocmd::Application const& app) {
        if(!app.args().empty()) {
            auto const& filename = app.args()[0];
            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            iopp::FileInputStream fis(filename, 0, n);

            // sample
            std::cout << "scan ..." << std::endl;

            size_t const min_sample_len = std::bit_floor(sampling / 3); // rate * (1/3)
            size_t const max_sample_len = 9 * min_sample_len; // rate * 3

            if(verbose) {
                std::cout << "\tmin_sample_len=" << min_sample_len << std::endl;
                std::cout << "\tmax_sample_len=" << max_sample_len << std::endl;
            }

            std::vector<Sample> samples;
            {
                Fingerprint const s = 2 * min_sample_len - 1; // rate * (2/3)

                // initialize I/O
                zk::internal::OverlappingBlocks block(fis, buffer_size, len);
                if(block.size() < len) std::abort();
                char const* p = block.begin();

                // initialize fingerprinting
                fp::RabinKarp31 rk(fp_base_, len);
                Fingerprint fp = 0;

                size_t i = 0;
                for(; i < len; i++) {
                    fp = rk.push(fp, *p++);
                }

                // initialize samples
                samples.emplace_back(0, max_minimizers);
                samples.back().insert(fp, max_minimizers);

                // process blocks
                size_t next_allowed = min_sample_len;
                size_t next_mandatory = max_sample_len;

                do {
                    if(block.empty()) continue;

                    while(p < block.end()) {
                        // read next character
                        fp = rk.roll(fp, *(p - len), *p);

                        if(((fp & s) == 0 && i >= next_allowed) || i >= next_mandatory) {
                            // finalize previous sample
                            size_t const new_beg = i - len + 1;
                            samples.back().len = new_beg - samples.back().index;
                            samples.emplace_back(new_beg, max_minimizers);

                            next_allowed = new_beg + min_sample_len;
                            next_mandatory = new_beg + max_sample_len;
                        }
                        samples.back().insert(fp, max_minimizers);

                        // advance character
                        ++i;
                        ++p;
                    }

                    // advance block
                    block.advance();
                    p = block.begin();
                } while(!block.last());

                samples.back().len = n - samples.back().index;
            }

            if(verbose) {
                std::cout << "Samples:" << std::endl;
                for(auto& x : samples) {
                    std::cout << "\tindex=" << x.index << ", len=" << x.len << ", minimizers=" << x.num_minimizers << std::endl;
                }
            }

            auto const num_samples = samples.size();

            // construct graph
            std::cout << "construct graph (num_samples=" << num_samples << ") ..." << std::endl;
            zk::internal::CompleteGraph g(num_samples);
            for(size_t i = 0; i < num_samples; i++) {
                g.dist(i, i) = 0.0f;
                for(size_t j = i+1; j < num_samples; j++) {
                    float const sim = float(samples[i].similarity(samples[j]));
                    float const d = 1.0f - sim;
                    g.dist(j, i) = d;
                    g.dist(i, j) = d;
                }
            }

            // nb: similarity is NOT a metric
            // compute APSP -- SLOW!
            /*
            g.all_pairs_shortest_paths();
            
            // convert similary to cost
            for(size_t i = 0; i < num_samples; i++) {
                for(size_t j = 0; j < num_samples; j++) {
                    g.dist(i, j) = 1.0f - g.dist(i, j);
                }
            }
            */

            // compute TSP tour via Christofides
            std::cout << "compute TSP tour ..." << std::endl;
            auto const tour = g.tsp_approx();            
            if(verbose) {
                float total_cost = 0;
                std::cout << std::endl << "tour:" << std::endl;
                for(size_t i = 1; i < tour.size(); i++) {
                    auto const cost = g.dist(tour[i-1], tour[i]);
                    total_cost += cost;

                    std::cout << "\t" << tour[i-1] << " -> " << tour[i] << " (cost " << cost << ")" << std::endl;
                }
                std::cout << "-> " << tour.size() << " blocks, total_cost=" << total_cost << std::endl;
            }

            // encode
            std::cout << "encode ..." << std::endl;
            {
                if(output_filename.empty()) {
                    output_filename = filename + ".rearr";

                    // encode tour
                    iopp::FileOutputStream fos(output_filename);
                    
                    auto write_uint32 = [&](uint32_t const x){
                        char* p = (char*)&x;
                        for(size_t i = 0; i < sizeof(uint32_t); i++) {
                            fos.put(p[i]);
                        }
                    };

                    write_uint32(tour.size());
                    for(auto v : tour) {
                        write_uint32(v);
                        write_uint32(samples[v].len);
                    }

                    // rearrange
                    auto buffer = std::make_unique<char[]>(max_sample_len);
                    for(auto v : tour) {
                        auto const offs = samples[v].index;
                        fis.seekg(offs, std::ios::beg);
                        fis.read(buffer.get(), samples[v].len);
                        fos.write(buffer.get(), samples[v].len);
                    }
                    // done
                }
            }

            pm::Result result;
            result.add("algo", "rearrange");
            result.add("file", std::filesystem::path(filename).filename().string());
            result.add("n", n);
            result.sort();
            result.print();
            return 0;
        } else {
            app.print_usage(*this);
            return -1;
        }
    }
};

int main(int argc, char** argv) {
    RearrangeTool app;
    return oocmd::Application::run(app, argc, argv);
}
