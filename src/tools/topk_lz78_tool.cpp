/**
 * tools/topk_lz78_tool.cpp
 * part of pdinklag/zketchy
 * 
 * MIT License
 * 
 * Copyright (c) Patrick Dinklage
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <zk/topk_lz78.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <cmdline/program.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

class TopkLZ78Tool : public cmdline::Program {
private:
    std::string filename;
    std::string output_filename;
    bool decompress = false;

    uint64_t block_size = 32_Ki; // best value according to many many experiments
    uint64_t prefix = UINTMAX_MAX;
    uint64_t k = 1_Mi;
    uint64_t max_freq = 1_Ki;

public:
    TopkLZ78Tool() : cmdline::Program("Top-k LZ78", "Top-k LZ78.") {
        required_arg("file", filename, "The input file.");
        option('o', "out", output_filename, "The output filename.");
        option('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
        option('b', "block-size", block_size, "The block size for encoding.");
        option('p', "prefix", prefix, "The prefix of the input file to consider.");
        option('k', "num-frequent", k, "The number of frequent substrings to maintain.");
        option('c', "max-freq", max_freq, "The maximum frequency of a frequent pattern.");
    }

    virtual int main() override {
        if(decompress) {
            if(output_filename.empty()) {
                output_filename = filename + ".dec";
            }

            iopp::FileInputStream in(filename);
            iopp::FileOutputStream out(output_filename);
            zk::TopkLZ78::decompress(in, out);
        } else {
            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            if(output_filename.empty()) {
                output_filename = filename + ".topklz78";
            }

            zk::TopkLZ78 topk_lz78(k, max_freq, block_size);
            {
                iopp::FileInputStream in(filename, 0, n);
                iopp::FileOutputStream out(output_filename);
                topk_lz78.compress(in, out);
            }

            auto result = topk_lz78.consume_last_result();
            result.add("algo", "topk_lz78");
            result.add("file", std::filesystem::path(filename).filename().string());
            result.add("n", n);
            result.add("k", k);
            result.add("fmax", max_freq);
            result.add("nout", std::filesystem::file_size(output_filename));
            result.sort();
            result.print();
        }
        return 0;
    }
};

int main(int argc, char** argv) {
    return TopkLZ78Tool().run(argc, argv);
}
