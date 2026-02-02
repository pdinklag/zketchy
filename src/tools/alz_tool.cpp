#include <cmdline/program.hpp>
#include <iopp/load_file.hpp>
#include <iopp/file_output_stream.hpp>

#define _ZK_SAMPLED_LPF_DEBUG
#include <zk/approximate_lz77.hpp>
#include <zk/internal/io/vbyte_coding.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <zk/internal/benchmark.hpp>

class ALZTool : public cmdline::Program {
private:
    static constexpr size_t MAX_SIZE_32BIT = 1ULL << 31 - 1;

    static constexpr char const* MAGIC = "ALZ";
    static constexpr size_t MAGIC_LEN = 3;

    static constexpr char MODE_VBYTE = 'V';

    std::string filename;
    std::string output_filename;
    bool decompress = false;
    
    size_t prefix = SIZE_MAX;
    uint64_t block_size = 32_Ki; // best value according to many many experiments

    uint64_t sampling = 4;
    uint64_t fp_window = 16;
    uint64_t window = 0;

    std::string load_input(size_t const n) {
        zk::internal::MemoryTimePhase t;
        std::cout << "load file " << filename << " (n=" << n << ") ... "; std::cout.flush();
        t.start();
        auto s = iopp::load_file_str(filename, n);
        t.stop();
        std::cout << "(" << (size_t)t.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms)" << std::endl;
        return s;
    }

    class VByteEmitter {
    private:
        iopp::FileOutputStream* out_;
        size_t z_;
    
    public:
        VByteEmitter(iopp::FileOutputStream& out) : out_(&out), z_(0) {
        }

        void emit_literal(char const c) {
            zk::internal::encode_vbyte(*out_, 0);
            out_->put(c);
            ++z_;
        }

        void emit_copy(uintmax_t const src, uintmax_t const len) {
            zk::internal::encode_vbyte(*out_, len);
            zk::internal::encode_vbyte(*out_, src);
            ++z_;
        }

        size_t num_phrases() const {
            return z_;
        }
    };

    template<typename Factorizer, typename Emitter>
    size_t factorize_vbyte(Factorizer& factorizer, Emitter&& emitter, size_t const n) {
        if(window > 0) {
            iopp::FileInputStream in(filename);
            factorizer.factorize(in, n, window, emitter);
        } else {
            auto s = load_input(n);
            factorizer.factorize(s.begin(), s.end(), emitter);
        }
        return emitter.num_phrases();
    }

    template<std::unsigned_integral Index>
    size_t compress(size_t const n) {
        zk::ApproximateLZ77<Index> alz(sampling, fp_window);

        iopp::FileOutputStream fout(output_filename);
        fout.write(MAGIC, MAGIC_LEN);
        fout.put(0); // reserved for flags
        zk::internal::encode_vbyte(fout, n);
        return factorize_vbyte(alz, VByteEmitter(fout), n);
    }

public:
    ALZTool() : cmdline::Program("alz", "Compute and encode an approximate LZ77 factorization") {
        required_arg("file", filename, "The input file.");
        option('s', "sampling", sampling, "The sampling rate (2^value).");
        option('l', "len", fp_window, "The fingerprint window size.");
        option('b', "block-size", block_size, "The block size for encoding.");
        option('p', "prefix", prefix, "Process only this prefix of the input file.");
        option('w', "window", window, "The input window size; leave at 0 to load entire input into RAM.");
        option('o', "out", output_filename, "The output filename.");
        option('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
    }

    virtual int main() override {
        if(decompress) {
            if(output_filename.empty()) {
                output_filename = filename + ".dec";
            }

            std::string s;
            {
                iopp::FileInputStream fin(filename);

                // check magic
                {
                    char magic[MAGIC_LEN];
                    fin.read(magic, MAGIC_LEN);
                    for(size_t i = 0; i < MAGIC_LEN; i++) {
                        if(magic[i] != MAGIC[i]) {
                            std::cerr << "wrong magic" << std::endl;
                            std::abort();
                        }
                    }
                }

                fin.get(); // skip currently unused flags byte

                auto const n = zk::internal::decode_vbyte(fin);
                while(s.length() < n) {
                    auto const len = zk::internal::decode_vbyte(fin);
                    if(len == 0) {
                        s.push_back(fin.get());
                    } else {
                        auto const src = s.length() - zk::internal::decode_vbyte(fin);
                        for(size_t i = 0; i < len; i++) {
                            s.push_back(s[src + i]);
                        }
                    }
                }
            }

            iopp::FileOutputStream fout(output_filename);
            fout.write(s.data(), s.length());
        } else {
            zk::internal::MemoryTimePhase t;

            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            if(window > n) {
                window = 0; // load whole file
            }

            if(output_filename.empty()) {
                output_filename = filename + ".alz";
            }

            size_t z = 0;
            {
                t.start();
                if(n <= MAX_SIZE_32BIT) {
                    z = compress<uint32_t>(n);
                } else {
                    z = compress<uint64_t>(n);
                }
                t.stop();
            }
            
            if constexpr(zk::internal::do_benchmark) {
                auto const nout = std::filesystem::file_size(output_filename);
                std::cout << "n=" << n << ", z=" << z << ", nout=" << nout << ", t=" << t.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << ", m=" << t.get_metric<pm::MallocCounter::MemoryPeakMetric>() << std::endl;
            }
        }
        return 0;
    }
};

int main(int argc, char** argv) {
    return ALZTool().run(argc, argv);
}