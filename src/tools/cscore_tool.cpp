#include <zk/cscore.hpp>

#include <cmdline/program.hpp>
#include <iopp/file_input_stream.hpp>

class CScoreTool : public cmdline::Program {
private:
    std::string filename;
    uint64_t sampling = 8;
    uint64_t len = 8;
    size_t block_size = 64_Ki;
    size_t buffer_size = 32_Mi;
    size_t prefix = SIZE_MAX;

public:
    CScoreTool() : cmdline::Program("Compressibility score", "Quickly estimate how compressible the input is") {
        required_arg("file", filename, "The input file.");
        option('s', "sample", sampling, "The sampling rate (2^value).");
        option('l', "len", len, "The pattern length.");
        option('b', "block_size", block_size, "The block size.");
        option('w', "buffer_size", buffer_size, "The buffer size.");
        option('p', "prefix", prefix, "Process only this prefix of the input file.");
    }

    virtual int main() override {
        if((buffer_size % block_size) != 0) {
            std::cerr << "buffer size must be a multiple of block size" << std::endl;
            return -1;
        }

        size_t const n = std::min(std::filesystem::file_size(filename), prefix);

        iopp::FileInputStream fis(filename, 0, n);

        zk::CScore cscore(len, sampling, block_size, buffer_size);
        double const score = cscore.compute(fis);

        auto result = cscore.consume_last_result();
        result.add("algo", "cscore");
        result.add("file", std::filesystem::path(filename).filename().string());
        result.add("n", n);
        result.sort();
        result.print();

        std::cout << score << std::endl;
        return 0;
    }
};

int main(int argc, char** argv) {
    return CScoreTool().run(argc, argv);
}
