#pragma once
#include <unordered_map>
#include <string>
#include <stdexcept>
#include <iostream>
#include "CKT.hpp"

// Insert function specifically for internal_nodes and return a new id for internal node
int CKTmkVolt(CKTcircuit &ckt, const std::string& key)
{
    int external_size = static_cast<int>(ckt.map.map_nodes.size());
    int internal_size = static_cast<int>(ckt.map.map_internal_nodes.size());
    int new_id = external_size + internal_size + 1;
    auto result = ckt.map.map_internal_nodes.insert({key, new_id});
    if (!result.second) {
        std::cerr << "Error: Duplicate internal node detected for " << key << std::endl;
        throw std::runtime_error("Duplicate internal node: " + key);
    }
    return result.first->second;
}