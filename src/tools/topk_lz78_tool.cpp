#include <zk/topk_lz78.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <cmdline/program.hpp>
#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

class TopkLZ78Tool : public cmdline::Program {
private:
    std::string filename;
    std::string output_filename;
    bool decompress = false;

    uint64_t block_size = 32_Ki; // best value according to many many experiments
    uint64_t prefix = UINTMAX_MAX;
    uint64_t k = 1_Mi;
    uint64_t max_freq = 1_Ki;

public:
    TopkLZ78Tool() : cmdline::Program("Top-k LZ78", "Top-k LZ78.") {
        required_arg("file", filename, "The input file.");
        option('o', "out", output_filename, "The output filename.");
        option('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
        option('b', "block-size", block_size, "The block size for encoding.");
        option('p', "prefix", prefix, "The prefix of the input file to consider.");
        option('k', "num-frequent", k, "The number of frequent substrings to maintain.");
        option('c', "max-freq", max_freq, "The maximum frequency of a frequent pattern.");
    }

    virtual int main() override {
        if(decompress) {
            if(output_filename.empty()) {
                output_filename = filename + ".dec";
            }

            iopp::FileInputStream in(filename);
            iopp::FileOutputStream out(output_filename);
            zk::TopkLZ78::decompress(iopp::bitwise_input_from(in), out);
        } else {
            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            if(output_filename.empty()) {
                output_filename = filename + ".topklz78";
            }

            zk::TopkLZ78 topk_lz78(k, max_freq, block_size);
            {
                iopp::FileInputStream in(filename, 0, n);
                iopp::FileOutputStream out(output_filename);
                topk_lz78.compress(in, iopp::bitwise_output_to(out));
            }

            auto result = topk_lz78.consume_last_result();
            result.add("algo", "topk_lz78");
            result.add("file", std::filesystem::path(filename).filename().string());
            result.add("n", n);
            result.add("k", k);
            result.add("fmax", max_freq);
            result.add("nout", std::filesystem::file_size(output_filename));
            result.sort();
            result.print();
        }
        return 0;
    }
};

int main(int argc, char** argv) {
    return TopkLZ78Tool().run(argc, argv);
}
