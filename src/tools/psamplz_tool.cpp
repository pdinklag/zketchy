#include <zk/psamplz.hpp>

#include <oocmd.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

class PSampLZTool : public oocmd::ConfigObject {
private:
    std::string output_filename;
    size_t window_size = 64_Mi;
    size_t len_exp_min = 10;
    size_t len_exp_max = 16;
    uint64_t sampling = 4;
    size_t bloom_filter_scale = 6;
    size_t prefix = SIZE_MAX;
    bool decompress = false;

public:
    PSampLZTool() : oocmd::ConfigObject("Compressibility score", "Quickly estimate how compressible the input is") {
        param('d', "decompress", decompress, "Decompress the input.");
        param('o', "out", output_filename, "The output file.");
        param("min", len_exp_min, "The minimum pattern length (2^value).");
        param("max", len_exp_max, "The maximum pattern length (2^value).");
        param('s', "sample", sampling, "The sampling rate (2^value).");
        param('w', "window", window_size, "The window size.");
        param('f', "bloom-scale", bloom_filter_scale, "The scale of the bloom filter - the number of bits will be this times the number of expected samples (0 to disable).");
        param('p', "prefix", prefix, "Process only this prefix of the input file.");
    }

    int run(oocmd::Application const& app) {
        if(!app.args().empty()) {
            auto const& filename = app.args()[0];
            if(decompress) {
                iopp::FileInputStream in(filename, 0);
                auto s = zk::PSampLZ::decompress(in);

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

                zk::PSampLZ psamplz(len_exp_min, len_exp_max, sampling, bloom_filter_scale, window_size);
                {
                    iopp::FileInputStream in(filename, 0, n);
                    iopp::FileOutputStream out(output_filename);
                    psamplz.compress(in, n, out);
                }

                auto stats = psamplz.consume_last_stats();
                std::cout << stats.gather_data().dump(4) << std::endl;

                auto result = psamplz.consume_last_result();
                result.add("algo", "psamplz");
                result.add("file", std::filesystem::path(filename).filename().string());
                result.add("n", n);
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
    PSampLZTool app;
    return oocmd::Application::run(app, argc, argv);
}
