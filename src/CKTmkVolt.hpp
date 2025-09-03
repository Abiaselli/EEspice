#pragma once
#include <unordered_map>
#include <string>
#include <stdexcept>
#include <iostream>

// Insert function specifically for internal_nodes and return a new id for internal node
int CKTmkVolt(std::unordered_map<std::string, int>& map_internal_nodes,
                         const std::string& key) 
{
    int new_id = static_cast<int>(map_internal_nodes.size()) + 1;
    auto result = map_internal_nodes.insert({key, new_id});
    if (!result.second) {
        std::cerr << "Error: Duplicate internal node detected for " << key << std::endl;
        throw std::runtime_error("Duplicate internal node: " + key);
    }
    return result.first->second;
}