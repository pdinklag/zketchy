#include <zk/cscore.hpp>
#include <zk/internal/util/complete_graph.hpp>

#include <oocmd.hpp>
#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>
#include <code/binary.hpp>

class RearrangeTool : public oocmd::ConfigObject {
private:
    uint64_t sampling = 8;
    uint64_t len = 8;
    size_t block_size = 64_Ki;
    size_t buffer_size = 32_Mi;
    size_t prefix = SIZE_MAX;
    bool verbose = false;
    std::string output_filename;

    struct Edge {
        uint32_t target;
        float cost;
    } __attribute__((packed));

public:
    RearrangeTool() : oocmd::ConfigObject("Rearrange", "Rearrange blocks of the input according to their similarity") {
        param('s', "sample", sampling, "The sampling rate (2^value).");
        param('l', "len", len, "The pattern length.");
        param('b', "block_size", block_size, "The block size.");
        param('w', "buffer_size", buffer_size, "The buffer size.");
        param('p', "prefix", prefix, "Process only this prefix of the input file.");
        param('v', "verbose", verbose, "Print a lot of info.");
        param('o', "out", output_filename, "The output filename.");
    }

    int run(oocmd::Application const& app) {
        if((buffer_size % block_size) != 0) {
            std::cerr << "buffer size must be a multiple of block size" << std::endl;
            return -1;
        }

        if(!app.args().empty()) {
            auto const& filename = app.args()[0];
            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            if((n + block_size - 1) / block_size > UINT32_MAX) {
                std::abort();
            }

            iopp::FileInputStream fis(filename, 0, n);

            // compute block signatures
            std::cout << "compute block signatures ..." << std::endl;
            zk::CScore cscore(len, sampling, block_size, buffer_size);
            auto sig = cscore.compute_block_signatures(fis);
            auto const num_blocks = sig.size();

            // construct graph
            std::cout << "construct graph (num_blocks=" << num_blocks << ") ..." << std::endl;
            zk::internal::CompleteGraph g(num_blocks);
            for(size_t i = 0; i < num_blocks; i++) {
                g.dist(i, i) = 0.0f;
                for(size_t j = i+1; j < num_blocks; j++) {
                    float const sim = float(sig[i].similarity(sig[j]));
                    g.dist(j, i) = sim;
                    g.dist(i, j) = sim;
                }
            }

            // nb: similarity is NOT a metric
            // compute APSP -- SLOW!
            /*
            g.all_pairs_shortest_paths();
            */
            
            // convert similary to cost
            for(size_t i = 0; i < num_blocks; i++) {
                for(size_t j = 0; j < num_blocks; j++) {
                    g.dist(i, j) = 1.0f - g.dist(i, j);
                }
            }

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
                    for(auto v : tour) write_uint32(v);

                    // rearrange
                    auto buffer = std::make_unique<char[]>(block_size);
                    for(auto v : tour) {
                        auto const offs = v * block_size;
                        fis.seekg(offs, std::ios::beg);
                        fis.read(buffer.get(), block_size);
                        fos.write(buffer.get(), block_size);
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
