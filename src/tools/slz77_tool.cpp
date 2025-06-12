#include <oocmd.hpp>
#include <iopp/bitwise_io.hpp>
#include <iopp/load_file.hpp>
#include <iopp/file_output_stream.hpp>

#define _ZK_SAMPLED_LPF_DEBUG
#include <zk/sampled_lpf_factorizer.hpp>
#include <zk/internal/io/block_coding.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <zk/internal/benchmark.hpp>

class SLZ77Tool : public oocmd::ConfigObject {
private:
    static constexpr uint64_t MAGIC =
        ((uint64_t)'S') << 24 |
        ((uint64_t)'Z') << 16 |
        ((uint64_t)'7') << 8 |
        ((uint64_t)'7');

    static constexpr zk::internal::TokenType TOK_LITERAL = 0;
    static constexpr zk::internal::TokenType TOK_REF_LEN = 1;
    static constexpr zk::internal::TokenType TOK_REF_SRC = 2;

    static void setup_encoding(zk::internal::BlockEncodingBase& enc, size_t const n) {
        enc.register_binary(255, false); // TOK_FACT_LITERAL
        enc.register_huffman();          // TOK_FACT_LEN
        enc.register_binary(n, false);   // TOK_REF_SRC
    }

    std::string output_filename;
    bool decompress = false;
    
    size_t prefix = SIZE_MAX;
    uint64_t block_size = 32_Ki; // best value according to many many experiments

    uint64_t sampling = 8;
    uint64_t fp_window = 16;

public:
    SLZ77Tool() : oocmd::ConfigObject("LZ77", "Compute and encode the exact LZ77 factorization") {
        param('s', "sampling", sampling, "The sampling rate (2^value).");
        param('l', "len", fp_window, "The fingerprint window size.");
        param('b', "block-size", block_size, "The block size for encoding.");
        param('p', "prefix", prefix, "Process only this prefix of the input file.");
        param('o', "out", output_filename, "The output filename.");
        param('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
    }

    int run(oocmd::Application const& app) {
        if(!app.args().empty()) {
            auto const& filename = app.args()[0];

            if(decompress) {
                if(output_filename.empty()) {
                    output_filename = filename + ".dec";
                }

                std::string s;
                {
                    iopp::FileInputStream fin(filename);
                    auto in = iopp::bitwise_input_from(fin);

                    uint64_t const magic = in.read(32);
                    if(magic != MAGIC) {
                        std::cerr << "wrong magic: 0x" << std::hex << magic << " (expected: 0x" << MAGIC << ")" << std::endl;
                        std::abort();
                    }

                    auto const n = in.read(64);
                    zk::internal::BlockDecoder dec(in);
                    setup_encoding(dec, n);

                    while(in) {
                        auto const len = dec.read_uint(TOK_REF_LEN);
                        if(len == 0) {
                            s.push_back(dec.read_char(TOK_LITERAL));
                        } else {
                            auto const src = s.length() - dec.read_uint(TOK_REF_SRC);
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

                if(output_filename.empty()) {
                    output_filename = filename + ".slz77";
                }

                size_t z = 0;
                t.start();
                {
                    auto s = iopp::load_file_str(filename, n);

                    iopp::FileOutputStream fout(output_filename);
                    auto out = iopp::bitwise_output_to(fout);

                    out.write(MAGIC, 32);
                    out.write(n, 64);

                    zk::internal::BlockEncoder enc(out, block_size);
                    setup_encoding(enc, n);

                    zk::SampledLPFFactorizer lz77(sampling, fp_window);
                    lz77.factorize(s.begin(), s.end(),
                        [&](lz77::Factor literal){
                            enc.write_uint(TOK_REF_LEN, 0);
                            enc.write_char(TOK_LITERAL, literal.literal());
                            ++z;
                        },
                        [&](lz77::Factor ref){
                            enc.write_uint(TOK_REF_LEN, ref.len);
                            enc.write_uint(TOK_REF_SRC, ref.src);
                            ++z;
                        });

                    enc.flush();
                }
                t.stop();

                if constexpr(zk::internal::do_benchmark) {                
                    std::cout << "n=" << n << ", z=" << z << ", t=" << t.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << ", m=" << t.get_metric<pm::MallocCounter::MemoryPeakMetric>() << std::endl;
                }
            }
            return 0;
        } else {
            app.print_usage(*this);
            return -1;
        }
    }
};

int main(int argc, char** argv) {
    SLZ77Tool app;
    return oocmd::Application::run(app, argc, argv);
}
