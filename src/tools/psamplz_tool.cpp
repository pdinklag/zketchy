#include <zk/psamplz.hpp>
#include <zk/psamplz_fggk15.hpp>

#include <cmdline/program.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

class PSampLZTool : public cmdline::Program {
private:
    std::string filename;
    std::string output_filename;
    size_t window_size = 64_Mi;
    size_t len_exp_min = 10;
    size_t len_exp_max = 16;
    uint64_t sampling = 4;
    size_t bloom_filter_scale = 6;
    size_t prefix = SIZE_MAX;
    bool decompress = false;
    bool alt = false;

    template<typename Compressor>
    int main_using() {
        if(decompress) {
            iopp::FileInputStream in(filename, 0);
            auto s = Compressor::decompress(in);

            if(s.empty()) {
                std::cerr << "ill-formed input" << std::endl;
                return -1;
            } else {
                if(output_filename.empty()) {
                    output_filename = filename + ".dec";
                }

                iopp::FileOutputStream out(output_filename);
                out.write(s.data(), s.length());
            }
        } else {
            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            if(output_filename.empty()) {
                output_filename = filename + ".psamplz";
            }

            Compressor psamplz(len_exp_min, len_exp_max, sampling, bloom_filter_scale, window_size);
            {
                iopp::FileInputStream in(filename, 0, n);
                iopp::FileOutputStream out(output_filename);
                psamplz.compress(in, n, out);
            }

            auto result = psamplz.consume_last_result();
            result.add("algo", "psamplz");
            result.add("file", std::filesystem::path(filename).filename().string());
            result.add("n", n);
            result.sort();
            result.print();
        }
        return 0;
    }

public:
    PSampLZTool() : cmdline::Program("Compressibility score", "Quickly estimate how compressible the input is") {
        required_arg("file", filename, "The input file.");
        option('d', "decompress", decompress, "Decompress the input.");
        option('o', "out", output_filename, "The output file.");
        option("min", len_exp_min, "The minimum pattern length (2^value).");
        option("max", len_exp_max, "The maximum pattern length (2^value).");
        option('s', "sample", sampling, "The sampling rate (2^value).");
        option('w', "window", window_size, "The window size.");
        option('f', "bloom-scale", bloom_filter_scale, "The scale of the bloom filter - the number of bits will be this times the number of expected samples (0 to disable).");
        option('p', "prefix", prefix, "Process only this prefix of the input file.");
        option("alt", alt, "Use alternative approximation.");
    }

    virtual int main() override {
        return alt ? main_using<zk::PSampLZ_FGGK15>() : main_using<zk::PSampLZ>();
    }
};

int main(int argc, char** argv) {
    return PSampLZTool().run(argc, argv);
}
