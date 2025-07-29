#include <iostream>

#include <cmdline/program.hpp>

#include <iopp/file_input_stream.hpp>
#include <iopp/load_file.hpp>

#include <zk/rlz_factorizer.hpp>

class RLZTool : public cmdline::Program {
private:
    std::string ref_filename;
    std::string text_filename;

    static void progress(size_t const i, size_t const n, size_t& s) {
        size_t const step = (size_t)(0.05 * double(n));
        if(i >= s * step) {
            double const p = 100.0 * double(i) / double(n);
            std::cout << "\t" << p << "%" << std::endl;

            s = i / step + 1;
        }
    }

public:
    RLZTool() : cmdline::Program("RLZ", "Relative Lempel-Ziv") {
        required_arg("ref", ref_filename);
        required_arg("text", text_filename);
    }

    int main() {
        std::cout << "building index ..."; std::cout.flush();
        auto r = iopp::load_file_str(ref_filename);
        zk::RLZFactorizer rlz(r);

        std::cout << std::endl;
        std::cout << "compressing ..." << std::endl;

        size_t const n = std::filesystem::file_size(text_filename);
        iopp::FileInputStream in(text_filename);
        size_t i = 0;
        size_t s = 1;
        size_t z = 0;
        size_t z_literal = 0;
        size_t z_ref = 0;
        size_t ref_len_sum = 0;

        rlz.factorize(in.begin(), in.end(), [&](lz77::Factor f){
            ++i; progress(i, n, s);
            ++z;
            ++z_literal;
        },
        [&](lz77::Factor f){
            i += f.num_literals(); progress(i, n, s);
            ++z;
            ++z_ref;
            ref_len_sum += f.len;
        });

        double const ref_len_avg = double(ref_len_sum) / double(z_ref);
        std::cout << "-> z=" << z << " (z_literal=" << z_literal << ", z_ref=" << z_ref << ", ref_len_avg=" << ref_len_avg << ")" << std::endl;
        return 0;
    }
};

int main(int argc, char** argv) {
    RLZTool rlz;
    rlz.run([&](){ return rlz.main(); }, argc, argv);
}
