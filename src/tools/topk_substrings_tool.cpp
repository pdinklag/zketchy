#include <fstream>
#include <libsais.h>
#include <memory>
#include <stack>

#include <cmdline/program.hpp>
#include <fp/rk61.hpp>
#include <iopp/load_file.hpp>
#include <iopp/file_output_stream.hpp>

#include <zk/internal/benchmark.hpp>
#include <zk/internal/sketch/topk_prefixes_misra_gries.hpp>
#include <zk/internal/util/suffix_tree.hpp>

class StringFrequencyTable {
public:
    struct Entry {
        uint64_t fp;
        size_t len;
        size_t freq;
    } __attribute__((packed));

private:
    std::vector<Entry> table_;

public:
    StringFrequencyTable() {
    }

    StringFrequencyTable(zk::internal::SuffixTree const& st, std::string_view const& text, size_t const k) {
        // PASS 1 - ignoring implcit nodes
        // compute node frequencies and occurrences
        struct NodeFreq {
            uint32_t v;
            uint32_t freq;
            uint32_t occ;

            NodeFreq() : v(0), freq(0), occ(0) {
            }

            NodeFreq(uint32_t const node) : v(node), freq(0), occ(0) {
            }
        } __attribute__((packed));

        std::vector<NodeFreq> nodes_with_freq;
        nodes_with_freq.reserve(st.size());
        for(uint32_t v = 0; v < st.size(); v++) {
            nodes_with_freq.emplace_back(v);
        }

        // compute frequencies and occurrences
        st.dfs(st.root(), [&](size_t const v) {
            if(v == st.root()) return;

            if(st.is_leaf(v)) {
                nodes_with_freq[v].freq = 1;
                nodes_with_freq[v].occ = v - 1;
            } else {
                size_t f = 0;
                for(auto u : st.children(v)) {
                    f += nodes_with_freq[u].freq;
                }
                nodes_with_freq[v].freq = f;
                nodes_with_freq[v].occ = nodes_with_freq[st.children(v)[0]].occ;
            }
        });

        // sort descending by frequency
        std::sort(nodes_with_freq.begin(), nodes_with_freq.end(),
            [&](NodeFreq const& a, NodeFreq const& b){
                return a.freq > b.freq;
            });

        // keep only the top k (so far)
        if(nodes_with_freq.size() > k) {
            nodes_with_freq.resize(k);
        }
        nodes_with_freq.shrink_to_fit();

        // PASS 2 - account for implicit nodes
        struct SubstringFreq {
            uint32_t pos;
            uint32_t len;
            uint32_t freq;
        } __attribute__((packed));

        std::vector<SubstringFreq> substrings_with_freq;
        substrings_with_freq.reserve(k);
        {
            for(auto const& x : nodes_with_freq) {
                uint32_t const d = st.depth(x.v);
                uint32_t const parent_d = st.depth(st.parent(x.v));

                // 
                for(uint32_t j = 1; j <= d - parent_d; j++) {
                    substrings_with_freq.push_back(SubstringFreq{ x.occ, parent_d + j, x.freq });
                }
            }
        }

        // keep only the top k
        if(substrings_with_freq.size() > k) {
            substrings_with_freq.resize(k);
        }

        // translate
        fp::RabinKarp61 rk(257);

        table_.reserve(substrings_with_freq.size());
        for(auto x : substrings_with_freq) {
            auto const s = text.substr(x.pos, x.len);

            uint64_t fp = 0;
            for(auto const c : s) {
                fp = rk.push(fp, c);
            }

            table_.push_back(Entry{fp, s.length(), x.freq});
        }
    }
    
    StringFrequencyTable(zk::internal::TopKPrefixesMisraGries<> const& topk, size_t const k) {
        // extract strings and produce entries
        char buffer[1024 * 1024]; // 1 meg should suffice...
        fp::RabinKarp61 rk(257);

        table_.reserve(k);
        for(size_t v = 1; v <= k; v++) {
            auto const len = topk.get(v, buffer);
            uint64_t fp = 0;
            for(size_t i = 0; i < len; i++) {
                fp = rk.push(fp, buffer[i]);
            }

            table_.push_back(Entry{fp, len, topk.freq(v)});
        }

        // sort descending by frequency
        std::sort(table_.begin(), table_.end(),
            [&](Entry const& a, Entry const& b){
                return a.freq > b.freq;
            });
    }

    std::vector<Entry> const& table() const { return table_; }

    void sort_by_fingerprint() {
        // sort ascending by fingerprint
        std::sort(table_.begin(), table_.end(),
            [&](Entry const& a, Entry const& b){
                return a.fp < b.fp;
            });
    }

    template<typename Out>
    void serialize(Out& out) {
        for(auto const& x : table_) {
            out << x.fp << "," << std::dec << x.len << "," << x.freq << std::endl;
        }
    }
};

class TopkSubstringsTool : public cmdline::Program {
private:
    std::string filename;
    size_t prefix = SIZE_MAX;
    bool exact = false;

    std::string outfilename;

    size_t k = SIZE_MAX;
    size_t max_freq = 1024 * 1024;

    StringFrequencyTable compute_exact(std::string_view const& t) {
        auto const n = t.length();
        zk::internal::MemoryTimePhase phase;

        std::cerr << "construct SA ... "; std::cerr.flush();
        phase.start();

        auto sa = std::make_unique<int32_t[]>(n);
        libsais((uint8_t const*)t.data(), sa.get(), n, 0, nullptr);

        phase.stop();
        std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms" << std::endl;

        std::cerr << "construct LCP ... "; std::cerr.flush();
        auto lcp = std::make_unique<int32_t[]>(n);
        {
            auto plcp = std::make_unique<int32_t[]>(n);
            libsais_plcp((uint8_t const*)t.data(), sa.get(), plcp.get(), n);
            libsais_lcp(plcp.get(), sa.get(), lcp.get(), n);
        }

        phase.stop();
        std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms" << std::endl;

        std::cerr << "compute suffix tree ... "; std::cerr.flush();
        phase.start();

        zk::internal::SuffixTree st(sa.get(), lcp.get(), n);
        
        phase.stop();
        std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms" << std::endl;

        std::cerr << "compute top-k substrings ... "; std::cerr.flush();
        StringFrequencyTable freqs(st, t, k);
        std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms" << std::endl;

        return freqs;
    }

    StringFrequencyTable compute_approx(std::string_view const& t) {
        auto const n = t.length();
        zk::internal::MemoryTimePhase phase;

        std::cerr << "process ... "; std::cerr.flush();
        phase.start();

        zk::internal::TopKPrefixesMisraGries<> topk(k+1, max_freq);
        for(size_t i = 0; i < n; i++) {
            auto s = topk.empty_string();
            size_t j = i;
            while(s.frequent && j < n) {
                s = topk.extend(s, t[j++]);
            }
        }

        phase.stop();
        std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms" << std::endl;

        std::cerr << "compute top-k substrings ... "; std::cerr.flush();
        phase.start();

        StringFrequencyTable freqs(topk, k);

        phase.stop();
        std::cerr << phase.get_metric<pm::Stopwatch::ElapsedTimeMillisMetric>() << "ms (renormalizations: " << topk.num_renormalizations() << ")" << std::endl;

        return freqs;
    }

public:
    TopkSubstringsTool() : cmdline::Program("topk-substrings", "Compute the top-k substrings") {
        required_arg("file", filename, "The input file.");
        option('p', "prefix", prefix, "Process only this prefix of the input file.");
        option('k', "count", k, "The number of substrings to extract");
        option("fmax", max_freq, "The maximum allowed frequency for the approximation.");
        option('e', "exact", exact, "Compute the exact top-k substrings using the suffix tree (memory heavy!).");
        option('o', "out", outfilename, "The output file; stdout if empty");
    }

    virtual int main() override {
        std::cerr << "load " << filename << " ... "; std::cerr.flush();

        auto t = iopp::load_file_str(filename, prefix);
        t.push_back(0);
        t.shrink_to_fit();

        std::cerr << "n=" << t.length() << std::endl;
        
        auto freqs = exact ? compute_exact(t) : compute_approx(t);
        if(outfilename.empty()) {
            freqs.serialize(std::cout);
        } else {
            std::ofstream f(outfilename);
            freqs.serialize(f);
            f.close();
        }

        return 0;
    }
};

int main(int argc, char** argv) {
    return TopkSubstringsTool().run(argc, argv);
}
