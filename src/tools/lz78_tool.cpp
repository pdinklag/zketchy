#include <cmdline/program.hpp>

#include <code/binary.hpp>

#include <iopp/bitwise_io.hpp>
#include <iopp/load_file.hpp>
#include <iopp/file_output_stream.hpp>

#include <zk/internal/io/block_coding.hpp>
#include <zk/internal/trie/trie_edge_array.hpp>
#include <zk/internal/util/si_iec_literals.hpp>

#include <zk/internal/benchmark.hpp>

#include <stack>

class LZ78Tool : public cmdline::Program {
private:
    static constexpr uint64_t MAGIC =
        ((uint64_t)'L') << 24 |
        ((uint64_t)'Z') << 16 |
        ((uint64_t)'7') << 8 |
        ((uint64_t)'8');

    using Node = uint64_t;

    std::string filename;
    std::string output_filename;
    bool decompress = false;
    
    size_t prefix = SIZE_MAX;

    class Trie {
    public:
        using NodeData = zk::internal::TrieEdgeArray<char, Node>;

    private:
        static constexpr size_t block_size_ = 1ULL << 20;
        static constexpr size_t block_mask_ = block_size_ - 1;

        std::vector<std::unique_ptr<NodeData[]>> blocks_;
        size_t size_;
        size_t capacity_;

        NodeData& node(Node const i) {
            auto const block = i / block_size_;
            return blocks_[block][i & block_mask_];
        }

        NodeData const& node(Node const i) const {
            auto const block = i / block_size_;
            return blocks_[block][i & block_mask_];
        }

        Node insert_child(Node const parent, char const c) {
            auto const v = size_++;
            if(v == capacity_) {
                blocks_.emplace_back(alloc_block());
                capacity_ += block_size_;
            }

            node(v) = NodeData();
            node(parent).insert(c, v);
            return v;
        }

        std::unique_ptr<NodeData[]> alloc_block() {
            return std::make_unique<NodeData[]>(block_size_);
        }

        template<typename Trie>
        void construct(Trie const& other, size_t const other_v, Node const v) {
            auto const& children = other.children_of(other_v);
            for(size_t i = 0; i < children.size(); i++) {
                Node child;
                follow_edge(v, children.label(i), child);
                construct(other, children[i], child);
            }
        }

    public:
        Trie() {
            clear();
        }

        Node root() const { return 0; }

        size_t size() const { return size_; }

        bool follow_edge(Node const v, char const c, Node& out_node) {
            auto const found = node(v).try_get(c, out_node);
            if(!found) {
                out_node = insert_child(v, c);
            }
            return found;
        }

        bool try_get_child(Node const v, char const c, Node& out_node) const {
            return node(v).try_get(c, out_node);
        }

        void clear() {
            size_ = 1;
            capacity_ = block_size_;
            blocks_.clear();
            blocks_.emplace_back(alloc_block());
            node(0) = NodeData();
        }
    };

public:
    LZ78Tool() : cmdline::Program("LZ78", "Compute and encode the exact LZ78 factorization") {
        required_arg("file", filename, "The input file.");
        option('p', "prefix", prefix, "Process only this prefix of the input file.");
        option('o', "out", output_filename, "The output filename.");
        option('d', "decompress", decompress, "Decompress the input file rather than compressing it.");
    }

    virtual int main() override {
        if(decompress) {
            if(output_filename.empty()) {
                output_filename = filename + ".dec";
            }

            std::string s;

            // TODO

            iopp::FileOutputStream fout(output_filename);
            fout.write(s.data(), s.length());
        } else {
            zk::internal::MemoryTimePhase t;
            size_t const n = std::min(std::filesystem::file_size(filename), prefix);

            if(output_filename.empty()) {
                output_filename = filename + ".lz78";
            }

            Trie trie;

            iopp::FileInputStream fin(filename, 0, n);
            iopp::FileOutputStream fout(output_filename);
            auto sink = iopp::bitwise_output_to(fout);

            size_t z = 0;
            size_t longest = 0;

            auto current = trie.root();
            size_t current_len = 0;
            for(auto it = fin.begin(); it != fin.end(); it++) {
                auto const c = *it;
                if(fin.good()) {
                    ++current_len;

                    Node v;
                    if(trie.follow_edge(current, c, v)) {
                        current = v;
                    } else {
                        // encode
                        code::Binary::encode(sink, current, std::bit_width(z));
                        code::Binary::encode(sink, c, code::Universe::of<unsigned char>());

                        ++z;
                        current = trie.root();

                        longest = std::max(longest, current_len);
                        current_len = 0;
                    }
                }
            }

            if(current != trie.root()) {
                longest = std::max(longest, current_len);
                code::Binary::encode(sink, current, std::bit_width(z));
                ++z;
            }

            std::cout << "z=" << z << ", longest=" << longest << std::endl;
        }
        return 0;
    }
};

int main(int argc, char** argv) {
    return LZ78Tool().run(argc, argv);
}
