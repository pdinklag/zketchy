#include <cmdline/program.hpp>
#include <iopp/bitwise_io.hpp>
#include <iopp/load_file.hpp>
#include <iopp/file_output_stream.hpp>

#define _ZK_SAMPLED_LPF_DEBUG
#include <zk/sampled_lpf_factorizer.hpp>
#include <zk/internal/io/block_coding.hpp>
#include <zk/internal/io/vbyte_coding.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <zk/internal/benchmark.hpp>

class SLZ77Tool : public cmdline::Program {
private:
    static constexpr char const* MAGIC = "SZ77";
    static constexpr size_t MAGIC_LEN = 4;

    static constexpr zk::internal::TokenType TOK_LITERAL = 0;
    static constexpr zk::internal::TokenType TOK_REF_LEN = 1;
    static constexpr zk::internal::TokenType TOK_REF_SRC = 2;

    static void setup_encoding(zk::internal::BlockEncodingBase& enc, size_t const n) {
        enc.register_binary(255, false); // TOK_FACT_LITERAL
        enc.register_huffman();          // TOK_FACT_LEN
        enc.register_binary(n, false);   // TOK_REF_SRC
    }

    std::string filename;
    std::string output_filename;
    bool decompress = false;
    
    size_t prefix = SIZE_MAX;
    uint64_t block_size = 32_Ki; // best value according to many many experiments

    uint64_t sampling = 4;
    uint64_t fp_window = 16;
    uint64_t window = 0;

    bool count_only = false;
    bool vbyte_coding = false;

    size_t factorize(zk::SampledLPFFactorizer& factorizer, std::string const& s) {
        size_t z = 0;
        factorizer.factorize(s.begin(), s.end(),
            [&](lz77::Factor literal){
                ++z;
            },
            [&](lz77::Factor ref){
                ++z;
            });
        return z;
    }

    size_t factorize_vbyte(zk::SampledLPFFactorizer& factorizer, std::string const& s, iopp::FileOutputStream& out) {
        size_t z = 0;
        if(window > 0) {
            std::cout << ">>> DISCLAIMER: testing streaming version -- text is still loaded fully in RAM <<<" << std::endl;
            factorizer.factorize_debug_streaming(s.begin(), s.end(), window,
                [&](lz77::Factor literal){
                    zk::internal::encode_vbyte(out, 0);
                    out.put(literal.literal());
                    ++z;
                },
                [&](lz77::Factor ref){
                    zk::internal::encode_vbyte(out, ref.len);
                    zk::internal::encode_vbyte(out, ref.src);
                    ++z;
                });
        } else {
            factorizer.factorize(s.begin(), s.end(),
                [&](lz77::Factor literal){
                    zk::internal::encode_vbyte(out, 0);
                    out.put(literal.literal());
                    ++z;
                },
                [&](lz77::Factor ref){
                    zk::internal::encode_vbyte(out, ref.len);
                    zk::internal::encode_vbyte(out, ref.src);
                    ++z;
                });
        }
        return z;
    }

    size_t factorize_block_enc(zk::SampledLPFFactorizer& factorizer, std::string const& s, auto& enc) {
        size_t z = 0;
        if(window > 0) {
            std::cout << ">>> DISCLAIMER: testing streaming version -- text is still loaded fully in RAM <<<" << std::endl;
            factorizer.factorize_debug_streaming(s.begin(), s.end(), window,
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
        } else {
            factorizer.factorize(s.begin(), s.end(),
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
        }
        return z;
    }

public:
    SLZ77Tool() : cmdline::Program("LZ77", "Compute and encode the exact LZ77 factorization") {
        required_arg("file", filename, "The input file.");
        option('s', "sampling", sampling, "The sampling rate (2^value).");
        option('l', "len", fp_window, "The fingerprint window size.");
        option('b', "block-size", block_size, "The block size for encoding.");
        option('p', "prefix", prefix, "Process only this prefix of the input file.");
        option('w', "window", window, "The input window size; leave at 0 to load entire input into RAM.");
        option('o', "out", output_filename, "The output filename.");
        option('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
        option("vbyte", vbyte_coding, "Use V-Byte encoding of phrases -- better suits the output for pipelining to other compressors.");
        option("count", count_only, "Only count the factors, don't actually write to the output file.");
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

                char const mode = fin.get();
                if(mode == 'B') {
                    auto in = iopp::bitwise_input_from(fin);

                    zk::internal::BlockDecoder dec(in);
                    auto const n = in.read(64);
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
                } else {
                    std::cerr << "decompression not implemented for mode '" << mode << "'" << std::endl;
                    std::abort();
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
            
            {
                std::cout << "load file " << filename << " (n=" << n << ") ... "; std::cout.flush();
                t.start();
                auto s = iopp::load_file_str(filename, n);
                t.stop();
                std::cout << "(" << (size_t)t.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms, peak mem " << t.get_metric<pm::MallocCounter::MemoryPeakMetric>() << ")" << std::endl;
                
                zk::SampledLPFFactorizer lz77(sampling, fp_window);
                t.start();
                if(!count_only) {
                    iopp::FileOutputStream fout(output_filename);
                    fout.write(MAGIC, MAGIC_LEN);

                    if(vbyte_coding) {
                        fout.put('V');
                        zk::internal::encode_vbyte(fout, n);
                        z = factorize_vbyte(lz77, s, fout);
                    } else {
                        fout.put('B');
                        auto out = iopp::bitwise_output_to(fout);

                        zk::internal::BlockEncoder enc(out, block_size);

                        out.write(n, 64);
                        setup_encoding(enc, n);

                        z = factorize_block_enc(lz77, s, enc);

                        enc.flush();
                    }
                } else {
                    z = factorize(lz77, s);
                }
                t.stop();
            }
            
            if constexpr(zk::internal::do_benchmark) {                
                std::cout << "n=" << n << ", z=" << z << ", t=" << t.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << ", m=" << t.get_metric<pm::MallocCounter::MemoryPeakMetric>() << std::endl;
            }
        }
        return 0;
    }
};

int main(int argc, char** argv) {
    return SLZ77Tool().run(argc, argv);
}