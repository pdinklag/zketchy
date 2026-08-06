#include <filesystem>

#include <cmdline/program.hpp>
#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>
#include <iopp/load_file.hpp>
#include <libsais64.h>
#include <omp.h>

#include <zk/internal/benchmark.hpp>

class RepTool : public cmdline::Program {
private:
    std::string filename;

public:
    RepTool() : cmdline::Program("Delta", "Compute the delta compressibility measure") {
        required_arg("file", filename, "The input file.");
    }

    virtual int main() override {
        size_t const n = std::filesystem::file_size(filename);
        std::cout << "file=" << filename << ", n=" << n; std::cout.flush();

        // load file
        auto s = iopp::load_file_str(filename, n);

        // alphabet and H0 entropy
        std::cout << ", sigma="; std::cout.flush();
        size_t sigma = 0;
        double h0 = 0;
        {
            size_t hist[256];
            for(size_t c = 0; c < 256; c++) hist[c] = 0;

            for(size_t i = 0; i < n; i++) ++hist[uint8_t(s[i])];

            for(size_t c = 0; c < 256; c++) {
                auto const nc = hist[c];
                if(nc) {
                    ++sigma;
                    h0 += (double(nc) / double(n)) * std::log2(double(n) / double(nc));
                }
            }
        }
        std::cout << sigma << ", h0=" << h0; std::cout.flush();

        // compute suffix array
        auto sa = std::make_unique<uint64_t[]>(n);
        #ifdef LIBSAIS_OPENMP
        libsais64_omp((uint8_t const*)s.data(), (int64_t*)sa.get(), n, 0, nullptr, omp_get_max_threads());
        #else
        libsais64((uint8_t const*)s.data(), (int64_t*)sa.get(), n, 0, nullptr);
        #endif

        // compute # of BWT runs
        std::cout << ", r="; std::cout.flush();
        size_t r = 0;
        {
            auto bwt = [&](size_t const i){
                auto const j = sa[i];
                return j > 0 ? s[j-1] : s[n-1];
            };
            
            uint8_t last = bwt(0);
            for(size_t i = 1; i < n; i++) {
                auto const c = bwt(i);
                if(last != 0 && c != last) ++r;
                last = c;
            }
        }
        std::cout << r; std::cout.flush();

        // compute phi
        auto work = std::make_unique<uint64_t[]>(n);
        for(size_t i = 1, prev = sa[0]; i < n; i++) {
            work[sa[i]] = prev;
            prev = sa[i];
        }
        work[sa[0]] = sa[n-1];

        // write suffix array to disk
        std::string sa_filename = filename + ".sa";
        {
            iopp::FileOutputStream f_sa(sa_filename);
            for(size_t i = 1; i < n; i++) {
                f_sa.write((char const*)&sa[i], sizeof(uint64_t));
            }
        }

        // discard suffix array
        sa.reset();

        // compute PLCP in-place
        for(size_t i = 0, l = 0; i < n - 1; ++i) {
            const size_t phi_i = work[i];
            while(s[i+l] == s[phi_i+l]) ++l;
            work[i] = l;
            if(l) --l;
        }

        // compute LCP streaming the suffix array
        auto lcp = std::make_unique<uint64_t[]>(n);
        {
            iopp::FileInputStream f_sa(sa_filename);
            lcp[0] = 0;
            uint64_t sa_i;
            for(size_t i = 1; i < n; i++) {
                f_sa.read((char*)&sa_i, sizeof(uint64_t));
                lcp[i] = work[sa_i];
            }
        }

        // discard work array and SA
        work.reset();
        std::filesystem::remove(sa_filename);

        // compute delta
        std::cout << ", delta="; std::cout.flush();
        double delta = 0;
        {
            auto dk = std::make_unique<uint32_t[]>(n);
            for(size_t i = 1; i < n; i++) {
                dk[lcp[i]+1]++;
                if(dk[lcp[i+1]] == UINT32_MAX) {
                    std::cerr << "int overflow" << std::endl;
                    return -1;
                }
            }

            double x = dk[1];
            delta = x;
            for(size_t k = 2; k < n; k++) {
                x = x + dk[k] - 1;
                delta = std::max(delta, x / k);
            }
        }
        std::cout << std::fixed << delta << std::endl;

        return 0;
    }
};

int main(int argc, char** argv) {
    return RepTool().run(argc, argv);
}
