#include <cmdline/program.hpp>
#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>
#include <iopp/stream_input_iterator.hpp>
#include <iopp/stream_output_iterator.hpp>

#include <zk/internal/sketch/topk_prefixes_lru.hpp>
#include <zk/internal/util/si_iec_literals.hpp>
#include <zk/internal/benchmark.hpp>
#include <zk/internal/io/block_coding.hpp>


class LZ78LRUTool : public cmdline::Program {
private:
    static constexpr zk::internal::TokenType TOK_REF = 0;
    static constexpr zk::internal::TokenType TOK_LITERAL = 1;

    static void setup_encoding(zk::internal::BlockEncodingBase& enc, size_t const k) {
        enc.register_binary(k-1); // TOK_REF
        enc.register_huffman();   // TOK_LITERAL
    }

    std::string filename;
    std::string output_filename;
    bool decompress = false;

    uint64_t block_size = 32_Ki; // best value according to many many experiments
    uint64_t prefix = UINTMAX_MAX;
    uint64_t k = 1_Mi;
    uint64_t max_freq = 1_Ki;

public:
    LZ78LRUTool() : cmdline::Program("Top-k LZ78 LRU", "Top-k LZ78 using the LRU heuristic.") {
        required_arg("file", filename, "The input file.");
        option('o', "out", output_filename, "The output filename.");
        option('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
        option('b', "block-size", block_size, "The block size for encoding.");
        option('p', "prefix", prefix, "The prefix of the input file to consider.");
        option('k', "num-frequent", k, "The number of frequent substrings to maintain.");
    }

    virtual int main() override {
        if(decompress) {
            // not implemented
        } else {
            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            if(output_filename.empty()) {
                output_filename = filename + ".lz78lru";
            }

            size_t z = 0;
            {
                iopp::FileInputStream in(filename, 0, n);
                iopp::FileOutputStream fout(output_filename);
                auto out = iopp::bitwise_output_to(fout);

                zk::internal::BlockEncoder enc(out, block_size);
                setup_encoding(enc, k);

                // initialize top-k
                zk::internal::TopKPrefixesLRU<> topk(k);

                // factorize
                using It = iopp::StreamInputIterator<iopp::FileInputStream>;
                auto const end = It::end(in);

                auto s = topk.empty_string();
                for(auto it = It(in); it != end; it++) {
                    auto const c = *it;
                    auto next = topk.extend(s, c);
                    if(!next.frequent) {
                        enc.write_uint(TOK_REF, s.node);
                        enc.write_char(TOK_LITERAL, c);

                        s = topk.empty_string();
                        ++z;
                    } else {
                        s = next;
                    }
                }

                // encode final phrase
                if(s.len > 0) {
                    enc.write_uint(TOK_REF, s.node);
                    ++z;
                }

                enc.flush();
            }

            zk::internal::Result result;
            result.add("algo", "topk_lz78");
            result.add("file", std::filesystem::path(filename).filename().string());
            result.add("n", n);
            result.add("nout", std::filesystem::file_size(output_filename));
            result.add("z", z);
            result.sort();
            result.print();
        }
        return 0;
    }
};

int main(int argc, char** argv) {
    return LZ78LRUTool().run(argc, argv);
}
