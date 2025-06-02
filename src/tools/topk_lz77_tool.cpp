#include <zk/topk_lz77.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <oocmd.hpp>
#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

class TopkLZ77Tool : public oocmd::ConfigObject {
private:
    std::string output_filename;
    bool decompress = false;

    uint64_t block_size = 32_Ki; // best value according to many many experiments
    uint64_t prefix = UINTMAX_MAX;
    uint64_t k = 1_Mi;
    uint64_t max_freq = 1_Ki;
    uint64_t lz_sampling = 4;
    uint64_t window_size = 1_Mi;

public:
    TopkLZ77Tool() : oocmd::ConfigObject("Top-k LZ77", "Best of both worlds approach to blockwise LZ77 and top-k LZ78.") {
        param('o', "out", output_filename, "The output filename.");
        param('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
        param('b', "block-size", block_size, "The block size for encoding.");
        param('p', "prefix", prefix, "The prefix of the input file to consider.");
        param('k', "num-frequent", k, "The number of frequent substrings to maintain.");
        param('c', "max-freq", max_freq, "The maximum frequency of a frequent pattern.");
        param('s', "lz-sampling", lz_sampling, "The LZ77 sampling rate (2^value, 0 for exact).");
        param('w', "window", window_size, "The window size.");
    }

    int run(oocmd::Application const& app) {
        if(!app.args().empty()) {
            auto const& filename = app.args()[0];
            if(decompress) {
                if(output_filename.empty()) {
                    output_filename = filename + ".dec";
                }

                iopp::FileInputStream in(filename);
                iopp::FileOutputStream out(output_filename);
                zk::TopkLZ77::decompress(iopp::bitwise_input_from(in), out);
            } else {
                size_t const n = std::min(std::filesystem::file_size(filename), prefix);

                if(output_filename.empty()) {
                    output_filename = filename + ".topklz77";
                }

                zk::TopkLZ77 topk_lz77(k, window_size, max_freq, lz_sampling, block_size);
                {
                    iopp::FileInputStream in(filename, 0, n);
                    iopp::FileOutputStream out(output_filename);
                    topk_lz77.compress(in, iopp::bitwise_output_to(out));
                }

                auto result = topk_lz77.consume_last_result();
                result.add("algo", "topk_lz77");
                result.add("file", std::filesystem::path(filename).filename().string());
                result.add("n", n);
                result.add("nout", std::filesystem::file_size(output_filename));
                result.sort();
                result.print();
            }
            return 0;
        } else {
            app.print_usage(*this);
            return -1;
        }
    }
};

int main(int argc, char** argv) {
    TopkLZ77Tool app;
    return oocmd::Application::run(app, argc, argv);
}
