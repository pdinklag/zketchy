#include <cscore.hpp>
#include <oocmd.hpp>

#include <internal/util/si_iec_literals.hpp>
#include <iopp/file_input_stream.hpp>

#include <iostream>

class CScoreTool : public oocmd::ConfigObject {
private:
    uint64_t sampling = 8;
    uint64_t len = 8;
    size_t block_size = 64_Ki;
    size_t buffer_size = 32_Mi;
    size_t prefix = SIZE_MAX;

public:
    CScoreTool() : oocmd::ConfigObject("Compressibility score", "Quickly estimate how compressible the input is") {
        param('s', "sample", sampling, "The sampling rate (2^value).");
        param('l', "len", len, "The pattern length.");
        param('b', "block_size", block_size, "The block size.");
        param('w', "buffer_size", buffer_size, "The buffer size.");
        param('p', "prefix", prefix, "Process only this prefix of the input file.");
    }

    int run(oocmd::Application const& app) {
        if((buffer_size % block_size) != 0) {
            std::cerr << "buffer size must be a multiple of block size" << std::endl;
            return -1;
        }

        if(!app.args().empty()) {
            auto const& filename = app.args()[0];
            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            iopp::FileInputStream fis(filename, 0, n);

            double const score = zk::cscore(fis, len, sampling, block_size, buffer_size);
            std::cout << score << std::endl;
            return 0;
        } else {
            app.print_usage(*this);
            return -1;
        }
    }
};

int main(int argc, char** argv) {
    CScoreTool app;
    return oocmd::Application::run(app, argc, argv);
}
