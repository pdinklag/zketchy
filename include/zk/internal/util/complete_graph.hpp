#pragma once

#include <algorithm>
#include <limits>
#include <list>
#include <random>
#include <stack>
#include <unordered_map>
#include <vector>

#include "union_find.hpp"

namespace zk::internal {

class CompleteGraph {
public:
    using Node = uint32_t;

private:
    size_t num_vertices_;
    std::unique_ptr<float[]> dist_;

    struct Edge {
        Node u, v;
        float weight;
    } __attribute__((packed));

    // Kruskal
    std::vector<Edge> minimum_spanning_tree() {
        std::vector<bool> visited(num_vertices_, false);

        std::vector<Edge> edges;
        edges.reserve((num_vertices_ * (num_vertices_-1)) / 2);
        for(Node u = 0; u < num_vertices_; u++) {
            for(Node v = u + 1; v < num_vertices_; v++) {
                edges.push_back(Edge{u, v, dist(u,v)});
            }
        }

        std::sort(edges.begin(), edges.end(), [](Edge const& a, Edge const& b){ return a.weight < b.weight; });

        std::vector<Edge> mst;
        UnionFind<Node> uf(num_vertices_);
        for(auto const& e : edges) {
            if(uf.find(e.u) != uf.find(e.v)) {
                mst.push_back(e);
                uf.unite(e.u, e.v);
                if(mst.size() == num_vertices_-1) break;
            }
        }
        mst.shrink_to_fit();
        return mst;
    }

    // vertices with odd degree
    std::vector<Node> find_odd_vertices(std::vector<Edge> const& edges) {
        auto degree = std::make_unique<Node[]>(num_vertices_);
        for(const auto& edge : edges) {
            degree[edge.u]++;
            degree[edge.v]++;
        }

        std::vector<Node> odd_vertices;
        for(Node u = 0; u < num_vertices_; u++) {
            if((degree[u] % 2) != 0) {
                odd_vertices.push_back(u);
            }
        }
        odd_vertices.shrink_to_fit();
        return odd_vertices;
    }

    // greedy heuristic -- shuffles odd_vertices
    std::vector<Edge> minimum_weight_matching(std::vector<Edge> const& mst, std::vector<Node>& odd_vertices) {
        std::vector<Edge> matching;
        if(odd_vertices.empty())[[unlikely]] return matching;

        matching.reserve(mst.size());
        for(auto& e : mst) {
            matching.push_back(e);
        }

        // shuffled copy of odd vertices
        auto const num_odd = odd_vertices.size();
        std::shuffle(odd_vertices.begin(), odd_vertices.end(), std::mt19937());
        std::vector<bool> visited(num_odd, false);

        for(size_t i = 0; i < num_odd; i++) {
            if(visited[i]) continue;

            auto const v = odd_vertices[i];
            auto min_distance = std::numeric_limits<float>::infinity();
            auto closest_u_index = SIZE_MAX;

            // find the closest unmatched odd vertex occurring after v in cur_odd
            for(size_t j = i+1; j < num_odd; j++) {
                if(!visited[j]) {
                    auto const u = odd_vertices[j];
                    auto const d = dist(v, u);
                    if(d < min_distance) {
                        min_distance = d;
                        closest_u_index = j;
                    }
                }
            }

            if(closest_u_index < SIZE_MAX) {
                auto const u = odd_vertices[closest_u_index];
                matching.push_back(Edge{v, u, min_distance});

                visited[i] = true;
                visited[closest_u_index] = true;
            } else {
                std::abort();
            }
        }
        return matching;
    }

    // Hierholzer
    std::vector<Node> eulerian_tour(std::vector<Edge> const& matching) {
        std::vector<Node> tour;
        if(matching.empty())[[unlikely]] return tour;

        // build adjacency list representation of the multigraph (MST + matching)
        struct Adjacency {
            Node neighbour;
            size_t edge_index;
        } __attribute__((packed));

        size_t const m = matching.size();
        std::vector<std::list<Adjacency>> adj(num_vertices_);
        
        for (size_t i = 0; i < m; i++) {
            auto const& e = matching[i];
            adj[e.u].push_back(Adjacency{e.v, i});
            adj[e.v].push_back(Adjacency{e.u, i});
        }
        std::vector<bool> used(m, false);

        std::stack<Node> path;
        auto u = matching[0].u;
        path.push(u);
        
        while(!path.empty()) {
            u = path.top();
            bool found_edge = false;

            // find an unused edge from the current node
            for(auto x : adj[u]) {
                auto const v = x.neighbour;
                auto const i = x.edge_index;

                if(!used[i]) {
                    used[i] = true;
                    path.push(v);
                    found_edge = true;
                    break;
                }
            }
        
            // if no unused edge was found from current_node, backtrack
            if(!found_edge) {
                tour.push_back(path.top());
                path.pop();
            }
        }

        // revert tour
        std::reverse(tour.begin(), tour.end());
        return tour;
    }

    std::vector<Node> hamiltonian_cycle(std::vector<Node> const& tour) {
        std::vector<Node> cycle;
        std::vector<bool> visited(tour.size(), false);

        for(size_t i = 0; i < tour.size(); i++) {
            auto const v = tour[i];
            if(!visited[v]) {
                cycle.push_back(v);
                visited[v] = true;
            }
        }
        return cycle;
    }

public:
    CompleteGraph() : num_vertices_(0) {
    }

    CompleteGraph(size_t const n) : num_vertices_(n), dist_(std::make_unique<float[]>(n * n)) {
    }

    CompleteGraph(CompleteGraph&&) = default;
    CompleteGraph& operator=(CompleteGraph&&) = default;

    CompleteGraph(CompleteGraph const&) = delete;
    CompleteGraph& operator=(CompleteGraph const&) = delete;

    // Christofides
    std::vector<Node> tsp_approx() {
        auto const mst = minimum_spanning_tree();
        auto odd = find_odd_vertices(mst);
        auto const matching = minimum_weight_matching(mst, odd);
        auto const tour = eulerian_tour(matching);
        return hamiltonian_cycle(tour);
    }

    // Floyd-Warshall
    void all_pairs_shortest_paths() {
        for(Node u = 0; u < num_vertices_; u++) {
            for(Node v = 0; v < num_vertices_; v++) {
                for(Node w = 0; w < num_vertices_; w++) {
                    auto const x = dist(u, v);
                    auto const y = dist(v, w);
                    auto& z = dist(u, w);
                    if(x + y < z) {
                        z = x + y;
                    }
                }
            }
        }
    }

    float& dist(Node const u, Node const v) {
        return dist_[u * num_vertices_ + v];
    }

    float dist(Node const u, Node const v) const {
        return dist_[u * num_vertices_ + v];
    }

    size_t num_vertices() const {
        return num_vertices_;
    }
};

}
