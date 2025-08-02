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

#include "models.hpp"
#include "map.hpp"
#include "bsim4v82/bsim4v82paser.hpp"

void NMOSParamLV1::setFromMap(const std::map<std::string, std::string>& kvMap) {
    if (kvMap.count("vto")) vt0 = std::stod(kvMap.at("vto"));
    if (kvMap.count("kp")) kp = std::stod(kvMap.at("kp"));
    if (kvMap.count("gamma")) gamma = std::stod(kvMap.at("gamma"));
    if (kvMap.count("phi")) phi = std::stod(kvMap.at("phi"));
    if (kvMap.count("lambda")) LAMBDA = std::stod(kvMap.at("lambda"));
    if (kvMap.count("Cgso")) CGSO = std::stod(kvMap.at("Cgso"));
    if (kvMap.count("Cgdo")) CGDO = std::stod(kvMap.at("Cgdo"));
    if (kvMap.count("Cgbo")) CGBO = std::stod(kvMap.at("Cgbo"));
    if (kvMap.count("Cbd")) CBD = std::stod(kvMap.at("Cbd"));
    if (kvMap.count("Cbs")) CBS = std::stod(kvMap.at("Cbs"));
    if (kvMap.count("RD")) RD = std::stod(kvMap.at("RD"));
    if (kvMap.count("RS")) RS = std::stod(kvMap.at("RS"));
    if (kvMap.count("RG")) RG = std::stod(kvMap.at("RG"));
}

void PMOSParamLV1::setFromMap(const std::map<std::string, std::string>& kvMap) {
    if (kvMap.count("vto")) vt0 = std::stod(kvMap.at("vto"));
    if (kvMap.count("kp")) kp = std::stod(kvMap.at("kp"));
    if (kvMap.count("gamma")) gamma = std::stod(kvMap.at("gamma"));
    if (kvMap.count("phi")) phi = std::stod(kvMap.at("phi"));
    if (kvMap.count("lambda")) LAMBDA = std::stod(kvMap.at("lambda"));
    if (kvMap.count("Cgso")) CGSO = std::stod(kvMap.at("Cgso"));
    if (kvMap.count("Cgdo")) CGDO = std::stod(kvMap.at("Cgdo"));
    if (kvMap.count("Cgbo")) CGBO = std::stod(kvMap.at("Cgbo"));
    if (kvMap.count("Cbd")) CBD = std::stod(kvMap.at("Cbd"));
    if (kvMap.count("Cbs")) CBS = std::stod(kvMap.at("Cbs"));
    if (kvMap.count("RD")) RD = std::stod(kvMap.at("RD"));
    if (kvMap.count("RS")) RS = std::stod(kvMap.at("RS"));
    if (kvMap.count("RG")) RG = std::stod(kvMap.at("RG"));
}

void parseModel(std::istringstream &iss, const std::string &line, Modelmap &modmap) {
    // Example: .model MyNMOS nmos (level=1 vto=0.7 kp=2e-5 ...)
    std::string modelName;
    int level = 1;
    std::string version;
    iss >> modelName; // e.g., "MyNMOS"

    // Get the rest of the line (type and parameters)
    std::string typeAndParams;
    std::getline(iss, typeAndParams);
    typeAndParams.erase(0, typeAndParams.find_first_not_of(" \t")); // Trim leading spaces

    // Extract model type and parameters
    size_t openParen = typeAndParams.find('(');
    size_t closeParen = typeAndParams.find(')');
    if (openParen == std::string::npos || closeParen == std::string::npos || closeParen < openParen) {
        throw ParsingException("Invalid .model line (missing parentheses): " + line, "parseModel");
    }

    std::string modelType = typeAndParams.substr(0, openParen);
    modelType.erase(std::remove_if(modelType.begin(), modelType.end(), ::isspace), modelType.end()); // Remove whitespace

    std::string paramsStr = typeAndParams.substr(openParen + 1, closeParen - openParen - 1);
    std::istringstream paramsStream(paramsStr);
    std::string paramToken;

    // Collect all parameters into a map
    std::map<std::string, std::string> kvMap;
    while (paramsStream >> paramToken) {
        size_t eqPos = paramToken.find('=');
        if (eqPos == std::string::npos) continue; // Skip malformed tokens
        std::string key = paramToken.substr(0, eqPos);
        std::string valueStr = paramToken.substr(eqPos + 1);
        kvMap[key] = valueStr;
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
            NMOSModel model; // Default level is 1
            model.level = level;
            model.version = version;
            NMOSParamLV1 param;
            param.setFromMap(kvMap);
            model.params = std::move(param);
            modmap.nmosModels[modelName] = model;
        }
        // Add more levels as needed, e.g.:
        // else if (model.level == 2) {
        //     NMOSParamLV2 param;
        //     param.setFromMap(kvMap);
        //     model.params = param;
        // }
        else if(level == 14){
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
            PMOSModel model; // Default level is 1
            model.level = level;
            model.version = version;
            PMOSParamLV1 param;
            param.setFromMap(kvMap);
            model.params = std::move(param);
            modmap.pmosModels[modelName] = model;
        }
        // Add more levels as needed, e.g.:
        // else if (model.level == 2) {
        //     PMOSParamLV2 param;
        //     param.setFromMap(kvMap);
        //     model.params = param;
        // }
        else if(level == 14){
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