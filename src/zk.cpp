#include <zk/cscore.hpp>
#include <zk/psamplz.hpp>
#include <zk/topk_lz77.hpp>

#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

#include <oocmd.hpp>

#include <zk/internal/util/si_iec_literals.hpp>

class ZK : public oocmd::ConfigObject {
private:
    bool decompress = false;

public:
    ZK() : oocmd::ConfigObject("zketchy compression utility", "Compresses the input") {
        param('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
    }

    int run(oocmd::Application const& app) {
        if(!app.args().empty()) {
            auto const& filename = app.args()[0];
            auto const n = std::filesystem::file_size(filename);
            if(decompress) {

            } else {
                iopp::FileInputStream in(filename);

                auto precompress = false;
                auto exact_lz77 = false;
                size_t len_exp_min, len_exp_max;

                // probe compressibility scores
                double score = zk::CScore(256, 8).compute(in);
                std::cout << "score for len=256: " << score << std::endl;
                in.seekg(0, std::ios_base::beg);
                if(score > 0.003) {
                    precompress = true;

                    score = zk::CScore(4096, 8).compute(in);
                    std::cout << "score for len=4096: " << score << std::endl;
                    if(score > 0.01) {
                        len_exp_min = 12;
                        len_exp_max = 12;
                    } else {
                        len_exp_min = 8;
                        len_exp_max = 12;
                    }
                } else {
                    score = zk::CScore(8, 8).compute(in);
                    std::cout << "score for len=8: " << score << std::endl;
                    if(score < 0.03) {
                        exact_lz77 = true;
                    }
                }
                
                // psamplz
                in.seekg(0, std::ios_base::beg);

                auto const tmp_filename = precompress ? filename + ".zk_tmp" : filename;
                if(precompress) {
                    std::cout << "psamplz ..." << std::endl;
                    std::cout.flush();

                    iopp::FileOutputStream out_tmp(tmp_filename);
                    zk::PSampLZ(len_exp_min, len_exp_max, 10).compress(in, n, out_tmp);
                }

                // topk-lz77
                std::cout << "topk-lz77 ..." << std::endl;
                std::cout.flush();
                in = iopp::FileInputStream(tmp_filename);
                {
                    auto const out_filename = filename + ".zk";
                    iopp::FileOutputStream fout(out_filename);
                    zk::TopkLZ77(16_Mi, 2_Gi-1, 1_Ki, exact_lz77 ? 0 : 4, 32_Ki).compress(in, iopp::bitwise_output_to(fout));
                }

                // clean up
                if(precompress) {
                    std::filesystem::remove(tmp_filename);
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
    ZK app;
    return oocmd::Application::run(app, argc, argv);
}
