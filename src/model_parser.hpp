#pragma once
#include <iostream>
#include <string>
#include <strings.h>
#include <sstream>
#include <fstream>
#include <map>
#include <variant>
#include <algorithm>
#include <vector>
#include <cctype>     // for std::tolower

#include "models.hpp"
#include "map.hpp"
#include "bsim4v82/bsim4v82paser.hpp"

void NMOSParamLV1::setFromMap(const std::map<std::string, std::string>& kvMap) {
    if (kvMap.count("vto")) vt0 = std::stod(kvMap.at("vto"));
    if (kvMap.count("kp")) kp = std::stod(kvMap.at("kp"));
    if (kvMap.count("gamma")) gamma = std::stod(kvMap.at("gamma"));
    if (kvMap.count("phi")) phi = std::stod(kvMap.at("phi"));
    if (kvMap.count("lambda")) LAMBDA = std::stod(kvMap.at("lambda"));
    if (kvMap.count("cgso")) CGSO = std::stod(kvMap.at("cgso"));
    if (kvMap.count("cgdo")) CGDO = std::stod(kvMap.at("cgdo"));
    if (kvMap.count("cgbo")) CGBO = std::stod(kvMap.at("cgbo"));
    if (kvMap.count("cbd")) CBD = std::stod(kvMap.at("cbd"));
    if (kvMap.count("cbs")) CBS = std::stod(kvMap.at("cbs"));
    if (kvMap.count("rd")) RD = std::stod(kvMap.at("rd"));
    if (kvMap.count("rs")) RS = std::stod(kvMap.at("rs"));
    if (kvMap.count("rg")) RG = std::stod(kvMap.at("rg"));
}

void PMOSParamLV1::setFromMap(const std::map<std::string, std::string>& kvMap) {
    if (kvMap.count("vto")) vt0 = std::stod(kvMap.at("vto"));
    if (kvMap.count("kp")) kp = std::stod(kvMap.at("kp"));
    if (kvMap.count("gamma")) gamma = std::stod(kvMap.at("gamma"));
    if (kvMap.count("phi")) phi = std::stod(kvMap.at("phi"));
    if (kvMap.count("lambda")) LAMBDA = std::stod(kvMap.at("lambda"));
    if (kvMap.count("cgso")) CGSO = std::stod(kvMap.at("cgso"));
    if (kvMap.count("cgdo")) CGDO = std::stod(kvMap.at("cgdo"));
    if (kvMap.count("cgbo")) CGBO = std::stod(kvMap.at("cgbo"));
    if (kvMap.count("cbd")) CBD = std::stod(kvMap.at("cbd"));
    if (kvMap.count("cbs")) CBS = std::stod(kvMap.at("cbs"));
    if (kvMap.count("rd")) RD = std::stod(kvMap.at("rd"));
    if (kvMap.count("rs")) RS = std::stod(kvMap.at("rs"));
    if (kvMap.count("rg")) RG = std::stod(kvMap.at("rg"));
}

void parseModel(std::istringstream &iss, const std::string &line, Modelmap &modmap) {
    // Example formats:
    // .model MyNMOS nmos (level=1 vto=0.7 kp=2e-5 ...)  - with parentheses
    // .model MyNMOS nmos level=1 vto=0.7 kp=2e-5 ...    - without parentheses
    
    std::string modelName;
    int level = 1;
    std::string version;
    
    iss >> modelName; // e.g., "MyNMOS" or "N90"
    
    // Get the rest of the line (type and parameters)
    std::string typeAndParams;
    std::getline(iss, typeAndParams);
    typeAndParams.erase(0, typeAndParams.find_first_not_of(" \t")); // Trim leading spaces
    
    // Extract model type and parameters
    std::string modelType;
    std::string paramsStr;
    
    // Check if parameters are enclosed in parentheses
    size_t openParen = typeAndParams.find('(');
    size_t closeParen = typeAndParams.find(')');
    
    if (openParen != std::string::npos && closeParen != std::string::npos && closeParen > openParen) {
        // Format with parentheses
        modelType = typeAndParams.substr(0, openParen);
        modelType.erase(std::remove_if(modelType.begin(), modelType.end(), ::isspace), modelType.end());
        paramsStr = typeAndParams.substr(openParen + 1, closeParen - openParen - 1);
    } else {
        // Format without parentheses
        // First token is the model type, rest are parameters
        std::istringstream typeStream(typeAndParams);
        typeStream >> modelType;
        
        // Get the remaining string as parameters
        std::streampos pos = typeStream.tellg();
        if (pos != std::streampos(-1)) {
            paramsStr = typeAndParams.substr(pos);
        } else {
            paramsStr = ""; // No parameters after model type
        }
    }
    
    // Parse parameters into key-value map
    std::map<std::string, std::string> kvMap;
    
    // Replace '=' with ' = ' to ensure proper tokenization
    // This handles both "level=14" and "level = 14" formats
    std::string processedParams;
    for (size_t i = 0; i < paramsStr.length(); i++) {
        if (paramsStr[i] == '=') {
            // Add spaces around '=' if not already present
            if (i > 0 && paramsStr[i-1] != ' ') {
                processedParams += ' ';
            }
            processedParams += '=';
            if (i < paramsStr.length() - 1 && paramsStr[i+1] != ' ') {
                processedParams += ' ';
            }
        } else {
            processedParams += paramsStr[i];
        }
    }
    
    std::istringstream paramsStream(processedParams);
    std::string token;
    std::string currentKey;
    bool expectingValue = false;
    
    while (paramsStream >> token) {
        if (token == "=") {
            expectingValue = true;
        } else if (expectingValue) {
            // This is a value for the current key
            if (!currentKey.empty()) {
                // Convert key to lowercase for consistency
                std::transform(currentKey.begin(), currentKey.end(), currentKey.begin(), 
                            [](unsigned char c){ return std::tolower(c); });
                kvMap[currentKey] = token;
                currentKey.clear();
            }
            expectingValue = false;
        } else {
            // parsed by this block "level=14" in one go.
            // Check if this token contains '='
            size_t eqPos = token.find('=');
            if (eqPos != std::string::npos) {
                // Token is in "key=value" format
                std::string key = token.substr(0, eqPos);
                std::string value = token.substr(eqPos + 1);
                // Convert key to lowercase
                std::transform(key.begin(), key.end(), key.begin(), 
                            [](unsigned char c){ return std::tolower(c); });
                kvMap[key] = value;
            } else {
                // This is a key
                currentKey = token;
            }
        }
    }
    
    // Handle NMOS models
    if (strcasecmp(modelType.c_str(), "nmos") == 0) {
        // Set level if specified
        if (kvMap.count("level")) {
            level = std::stoi(kvMap.at("level"));
        }
        
        // Set version if specified
        if (kvMap.count("version")) {
            version = kvMap.at("version");
        }
        
        // Assign parameters based on level
        if (level == 1) {
            NMOSModel model;
            model.level = level;
            model.version = version;
            NMOSParamLV1 param;
            param.setFromMap(kvMap);
            model.params = std::move(param);
            modmap.nmosModels[modelName] = model;
        }
        else if (level == 14) {
            std::shared_ptr<bsim4::BSIM4model> bsim4Model = bsim4::paserBSIM4Model("nmos", modelName, kvMap);
            modmap.bsim4Models[modelName] = bsim4Model;
        }
        else {
            throw ParsingException("Unsupported NMOS level: " + std::to_string(level) + " in line: " + line, "parseModel");
        }
    }
    // Handle PMOS models
    else if (strcasecmp(modelType.c_str(), "pmos") == 0) {
        // Set level if specified
        if (kvMap.count("level")) {
            level = std::stoi(kvMap.at("level"));
        }
        
        // Set version if specified
        if (kvMap.count("version")) {
            version = kvMap.at("version");
        }
        
        // Assign parameters based on level
        if (level == 1) {
            PMOSModel model;
            model.level = level;
            model.version = version;
            PMOSParamLV1 param;
            param.setFromMap(kvMap);
            model.params = std::move(param);
            modmap.pmosModels[modelName] = model;
        }
        else if (level == 14) {
            std::shared_ptr<bsim4::BSIM4model> bsim4Model = bsim4::paserBSIM4Model("pmos", modelName, kvMap);
            modmap.bsim4Models[modelName] = bsim4Model;
        }
        else {
            throw ParsingException("Unsupported PMOS level: " + std::to_string(level) + " in line: " + line, "parseModel");
        }
    }
    else {
        throw ParsingException("Unknown model type: " + modelType + " in line: " + line, "parseModel");
    }
}