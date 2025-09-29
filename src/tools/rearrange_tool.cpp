#include <zk/minimizer_sampling.hpp>
#include <zk/internal/benchmark.hpp>
#include <zk/internal/util/complete_graph.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <iopp/bitwise_io.hpp>
#include <iopp/file_input_stream.hpp>
#include <iopp/file_output_stream.hpp>

#include <cmdline/program.hpp>

class RearrangeTool : public cmdline::Program {
private:
    using Fingerprint = fp::RabinKarp31::Fingerprint;
    static constexpr Fingerprint fp_base_ = 257;

    bool do_decode = false;
    size_t sampling = 16_Ki;
    size_t len = 8;
    size_t max_minimizers = 64;
    size_t buffer_size = 64_Mi;
    double tsp_threshold = 0.75;
    size_t prefix = SIZE_MAX;
    bool verbose = false;

    std::string filename;
    std::string output_filename;

    struct Sample {
        size_t index;
        uint32_t len;

        std::unique_ptr<Fingerprint[]> minimizers;
        uint32_t num_minimizers;

        Sample(size_t const _index, size_t const max_minimizers) : index(_index), len(0), minimizers(std::make_unique<Fingerprint[]>(max_minimizers)), num_minimizers(0) {
        }

        void insert(Fingerprint const fp, size_t const max) {
            if(num_minimizers == max && fp > minimizers[num_minimizers])[[unlikely]] return;

            size_t rank = 0;
            while(rank < num_minimizers && minimizers[rank] < fp) {
                ++rank;
            }

            if(num_minimizers < max) {
                ++num_minimizers;
            }

            for(size_t i = num_minimizers; i > rank; i--) {
                minimizers[i] = minimizers[i-1];
            }
            minimizers[rank] = fp;
        }

        double similarity(Sample const& other) const {
            size_t cut_size = 0;
            size_t i = 0, j = 0;
            while(i < num_minimizers && j < other.num_minimizers) {
                if(minimizers[i] < other.minimizers[j]) {
                    ++i;
                } else if(minimizers[i] > other.minimizers[j]) {
                    ++j;
                } else {
                    ++cut_size;
                    ++i;
                    ++j;
                }
            }

            auto const union_size = num_minimizers + other.num_minimizers - cut_size;
            return double(cut_size) / double(union_size);
        }

        Fingerprint hash() const {
            Fingerprint h = 17;
            for(size_t i = 0; i < num_minimizers; i++) {
                h = 31 * h + minimizers[i];
            }
            return h;
        }
    };

    struct Edge {
        uint32_t target;
        float cost;
    } __attribute__((packed));

public:
    RearrangeTool() : cmdline::Program("Rearrange", "Rearrange blocks of the input according to their similarity") {
        required_arg("file", filename, "The input file.");
        option('d', "decode", do_decode, "Decode the input.");
        option('s', "sample", sampling, "The sampling rate.");
        option('l', "len", len, "The pattern length.");
        option('m', "minimizers", max_minimizers, "The maximum number of minimizers per sample.");
        option('w', "buffer-size", buffer_size, "The buffer size.");
        option('p', "prefix", prefix, "Process only this prefix of the input file.");
        option('t', "tsp-threshold", tsp_threshold, "Use an arbitrary tour over TSP if the cluster ratio goes above this value.");
        option('v', "verbose", verbose, "Print a lot of info.");
        option('o', "out", output_filename, "The output filename.");
    }

    void encode() {
        pm::Result result;
        zk::internal::TimePhase phase;

        size_t const n = std::min(std::filesystem::file_size(filename), prefix);
        iopp::FileInputStream fis(filename, 0, n);

        // sample
        zk::MinimizerSampling s(sampling, len, max_minimizers);
        result.add("min_sample_len", s.min_sample_len());
        result.add("max_sample_len", s.max_sample_len());
        result.add("max_minimizers", max_minimizers);

        std::cout << "scan ... "; std::cout.flush();

        phase.start();
        auto samples = s.sample(fis, buffer_size);
        auto const num_samples = samples.size();
        phase.stop();
        std::cout << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms -> " << num_samples << " samples" << std::endl;

        if(verbose) {
            std::cout << "Samples:" << std::endl;
            for(auto& x : samples) {
                std::cout << "\tindex=" << x.index << ", len=" << x.len << ", minimizers=" << x.num_minimizers << ", hash=" << x.hash() << std::endl;
            }
        }

        // cluster
        std::cout << "cluster ... "; std::cout.flush();

        phase.start();
        auto clusters = zk::MinimizerSampling::compute_clusters(samples);
        auto const num_clusters = clusters.size();
        auto const cluster_ratio = double(num_clusters) / double(num_samples);
        phase.stop();
        std::cout << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms -> " << num_clusters << " clusters (ratio=" << cluster_ratio << ")" << std::endl;

        result.add("num_samples", num_samples);
        result.add("num_clusters", num_clusters);
        result.add("cluster_ratio", cluster_ratio);

        std::vector<uint32_t> tour;

        if(cluster_ratio > tsp_threshold) {
            std::cout << "continuing with arbitrary tour due to a large amount of clusters" << std::endl;
            tour.reserve(num_clusters);
            for(uint32_t i = 0; i < num_clusters; i++) {
                tour.push_back(i);
            }
        } else {
            // construct graph
            std::cout << "construct graph ... "; std::cout.flush();

            phase.start();
            zk::internal::CompleteGraph g(num_clusters);
            for(size_t i = 0; i < num_clusters; i++) {
                g.dist(i, i) = 0.0f;
                for(size_t j = i+1; j < num_clusters; j++) {
                    float const sim = float(clusters[i].similarity(samples, clusters[j]));
                    float const d = 1.0f - sim;
                    g.dist(j, i) = d;
                    g.dist(i, j) = d;
                }
            }
            phase.stop();
            std::cout << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms" << std::endl;

            // nb: similarity is NOT a metric
            // compute APSP -- SLOW!
            /*
            g.all_pairs_shortest_paths();
            
            // convert similary to cost
            for(size_t i = 0; i < num_samples; i++) {
                for(size_t j = 0; j < num_samples; j++) {
                    g.dist(i, j) = 1.0f - g.dist(i, j);
                }
            }
            */

            // compute TSP tour via Christofides        
            std::cout << "compute TSP tour ... "; std::cout.flush();
            phase.start();
            tour = g.tsp_approx();            
            phase.stop();
            std::cout << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms" << std::endl;

            if(verbose) {
                float total_cost = 0;
                std::cout << std::endl << "tour:" << std::endl;
                for(size_t i = 1; i < tour.size(); i++) {
                    auto const cost = g.dist(tour[i-1], tour[i]);
                    total_cost += cost;

                    std::cout << "\t" << tour[i-1] << " -> " << tour[i] << " (cost " << cost << ")" << std::endl;
                }
                std::cout << "-> " << tour.size() << " blocks, total_cost=" << total_cost << std::endl;
            }
        }

        // encode
        std::cout << "encode ... "; std::cout.flush();
        phase.start();
        {
            if(output_filename.empty()) {
                output_filename = filename + ".re";
            }

            // encode tour
            iopp::FileOutputStream fos(output_filename);
            
            auto write_uint32 = [&](uint32_t const x){
                char* p = (char*)&x;
                for(size_t i = 0; i < sizeof(uint32_t); i++) {
                    fos.put(p[i]);
                }
            };

            {
                write_uint32(num_samples);

                size_t count = 0;
                for(auto v : tour) {
                    for(auto i : clusters[v].sample_indices) {
                        ++count;
                        write_uint32(i);
                        write_uint32(samples[i].len);
                    }
                }
                assert(count == num_samples);
            }

            // rearrange
            auto buffer = std::make_unique<char[]>(s.max_sample_len());
            for(auto v : tour) {
                for(auto i : clusters[v].sample_indices) {
                    auto const offs = samples[i].index;
                    fis.seekg(offs, std::ios::beg);
                    fis.read(buffer.get(), samples[i].len);
                    fos.write(buffer.get(), samples[i].len);
                }
            }
        }
        phase.stop();
        std::cout << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << " ms" << std::endl;

        result.print();
    }

    void decode() {
        // decode tour
        iopp::FileInputStream fis(filename);

        auto read_uint32 = [&](){
            uint32_t x;
            char* p = (char*)&x;
            for(size_t i = 0; i < sizeof(uint32_t); i++) {
                *p++ = fis.get();
            }
            return x;
        };

        struct DecSample {
            uint32_t rank;
            uint32_t len;
            size_t index;
        } __attribute__((packed));

        auto const num_samples = read_uint32();
        std::cout << "decoding " << num_samples << " samples" << std::endl;

        std::vector<DecSample> samples;
        samples.reserve(num_samples);

        size_t index = 0;
        uint32_t longest = 0;
        for(size_t i = 0; i < num_samples; i++) {
            auto const rank = read_uint32();
            auto const len = read_uint32();
            
            samples.emplace_back(DecSample{rank, len, index});

            longest = std::max(longest, len);
            index += len;
        }

        size_t const offset = fis.tellg();

        // compute original order
        std::sort(samples.begin(), samples.end(), [](DecSample const& a, DecSample const& b){
            return a.rank < b.rank;
        });

        // unarrange
        if(output_filename.empty()) {
            output_filename = filename + ".dec";
        }

        {
            iopp::FileOutputStream fos(output_filename);
            auto buffer = std::make_unique<char[]>(longest);

            for(auto const& x : samples) {
                if(verbose) {
                    std::cout << "\tindex=" << x.index << ", rank=" << x.rank << ", len=" << x.len << std::endl;
                }

                fis.seekg(offset + x.index, std::ios::beg);
                fis.read(buffer.get(), x.len);
                fos.write(buffer.get(), x.len);
            }
        }

        
    }

    virtual int main() override {
        if(do_decode) {
            decode();
        } else {
            encode();
        }
        return 0;
    }
};

int main(int argc, char** argv) {
    return RearrangeTool().run(argc, argv);
}
