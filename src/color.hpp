#pragma once
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include "CKT.hpp"
#include "bsim4v82/bsim4v82.hpp"

class BSIM4Coloring {
private:
    std::vector<std::vector<size_t>> color_groups;
    size_t num_colors;
    
    // Extract all nodes used by a BSIM4 instance
    // Note: Excludes node 0 (ground) since grounded nodes don't create matrix conflicts
    std::unordered_set<int> getInstanceNodes(const BSIM4& instance) {
        std::unordered_set<int> nodes;

        // Add all non-ground nodes that this instance uses
        // Node 0 is excluded because grounded nodes don't appear in MNA matrix equations
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::D_NODE]
            && instance.bsim4v82Instance.BSIM4dNode != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4dNode);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::G_NODE_EXT]
            && instance.bsim4v82Instance.BSIM4gNodeExt != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4gNodeExt);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::S_NODE]
            && instance.bsim4v82Instance.BSIM4sNode != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4sNode);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::B_NODE]
            && instance.bsim4v82Instance.BSIM4bNode != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4bNode);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::D_NODE_PRIME]
            && instance.bsim4v82Instance.BSIM4dNodePrime != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4dNodePrime);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::G_NODE_PRIME]
            && instance.bsim4v82Instance.BSIM4gNodePrime != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4gNodePrime);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::G_NODE_MID]
            && instance.bsim4v82Instance.BSIM4gNodeMid != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4gNodeMid);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::S_NODE_PRIME]
            && instance.bsim4v82Instance.BSIM4sNodePrime != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4sNodePrime);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::B_NODE_PRIME]
            && instance.bsim4v82Instance.BSIM4bNodePrime != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4bNodePrime);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::DB_NODE]
            && instance.bsim4v82Instance.BSIM4dbNode != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4dbNode);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::SB_NODE]
            && instance.bsim4v82Instance.BSIM4sbNode != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4sbNode);
        if (instance.bsim4v82Instance.BSIM4nodeValid[bsim4::BSIM4V82::NodeType::Q_NODE]
            && instance.bsim4v82Instance.BSIM4qNode != 0)
            nodes.insert(instance.bsim4v82Instance.BSIM4qNode);

        return nodes;
    }
    
    // Check if two node sets intersect
    bool nodesIntersect(const std::unordered_set<int>& nodes1, 
                       const std::unordered_set<int>& nodes2) {
        // Use smaller set for iteration (optimization)
        if (nodes1.size() > nodes2.size()) {
            return nodesIntersect(nodes2, nodes1);
        }
        
        for (int node : nodes1) {
            if (nodes2.count(node) > 0) {
                return true;
            }
        }
        return false;
    }

public:
    void computeColoring(const std::vector<BSIM4>& instances) {
        const size_t n = instances.size();
        
        // Step 1: Extract node sets for all instances
        std::vector<std::unordered_set<int>> instance_nodes(n);
        for (size_t i = 0; i < n; ++i) {
            instance_nodes[i] = getInstanceNodes(instances[i]);
        }
        
        // Step 2: Build conflict graph (adjacency list)
        std::vector<std::unordered_set<size_t>> adj(n);
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = i + 1; j < n; ++j) {
                if (nodesIntersect(instance_nodes[i], instance_nodes[j])) {
                    adj[i].insert(j);
                    adj[j].insert(i);
                }
            }
        }
        
        // Step 3: Greedy coloring with degree-based ordering
        std::vector<std::pair<size_t, size_t>> degree_order;
        for (size_t i = 0; i < n; ++i) {
            degree_order.push_back({adj[i].size(), i});
        }
        // Sort by degree (descending) - color high-degree nodes first
        std::sort(degree_order.begin(), degree_order.end(), 
                  std::greater<std::pair<size_t, size_t>>());
        
        std::vector<int> colors(n, -1);
        num_colors = 0;
        
        for (const auto& [degree, i] : degree_order) {
            // Find the smallest color not used by neighbors
            std::unordered_set<int> neighbor_colors;
            for (size_t neighbor : adj[i]) {
                if (colors[neighbor] != -1) {
                    neighbor_colors.insert(colors[neighbor]);
                }
            }
            
            int color = 0;
            while (neighbor_colors.count(color) > 0) {
                color++;
            }
            
            colors[i] = color;
            num_colors = std::max(num_colors, static_cast<size_t>(color + 1));
        }
        
        // Step 4: Group instances by color
        color_groups.clear();
        color_groups.resize(num_colors);
        for (size_t i = 0; i < n; ++i) {
            color_groups[colors[i]].push_back(i);
        }
    }
    
    const std::vector<std::vector<size_t>>& getColorGroups() const {
        return color_groups;
    }
    
    size_t getNumColors() const {
        return num_colors;
    }
};