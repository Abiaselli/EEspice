#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <map>

#include "CKT.hpp"
#include "global.hpp"

/*
    This new helper file will contain the implementation details for the batch simulation, 
    such as the logic for applying a specific set of parameters to a circuit instance.
*/

namespace batch {

// A structure to hold information about a single parameter being swept.
struct BatchParam {
    std::string device_id_str;
    char device_type;
    std::string param_name; // e.g., "value", "W", "L"
    std::vector<double> values;
};

// A map to represent a single circuit configuration, e.g., {"V1.value": 1.0, "M1.W": 10e-6}
using CircuitConfig = std::map<std::string, double>;

// Function to find all parameters that have batch values specified.
std::vector<BatchParam> find_batch_params(const CircuitElements& elements) {
    std::vector<BatchParam> params;

    for (const auto& vol : elements.voltageSources) {
        if (!vol.batchValues.empty()) {
            params.push_back({vol.id_str, 'V', "value", vol.batchValues});
        }
    }
    for (const auto& res : elements.resistors) {
        if (!res.batchValues.empty()) {
            params.push_back({res.id_str, 'R', "value", res.batchValues});
        }
    }
    for (const auto& cap : elements.capacitors) {
        if (!cap.batchValues.empty()) {
            params.push_back({cap.id_str, 'C', "value", cap.batchValues});
        }
    }
    for (const auto& cur : elements.currentSources) {
        if (!cur.batchValues.empty()) {
            params.push_back({cur.id_str, 'I', "value", cur.batchValues});
        }
    }
    for (const auto& vccs : elements.vccs) {
        if (!vccs.batchValues.empty()) {
            params.push_back({vccs.id_str, 'G', "value", vccs.batchValues});
        }
    }
    for (const auto& nmos : elements.nmos) {
        if (!nmos.batchW.empty()) {
            params.push_back({nmos.id_str, 'M', "W", nmos.batchW});
        }
        if (!nmos.batchL.empty()) {
            params.push_back({nmos.id_str, 'M', "L", nmos.batchL});
        }
    }
    for (const auto& pmos : elements.pmos) {
        if (!pmos.batchW.empty()) {
            params.push_back({pmos.id_str, 'M', "W", pmos.batchW});
        }
        if (!pmos.batchL.empty()) {
            params.push_back({pmos.id_str, 'M', "L", pmos.batchL});
        }
    }
    for (const auto& bsim4 : elements.bsim4) {
        if (!bsim4.batchW.empty()) {
            params.push_back({bsim4.id_str, 'M', "W", bsim4.batchW});
        }
        if (!bsim4.batchL.empty()) {
            params.push_back({bsim4.id_str, 'M', "L", bsim4.batchL});
        }
    }
    // Diodes have two potential batch parameters
    for (const auto& diode : elements.diodes) {
        if (!diode.batchIs.empty()) {
            params.push_back({diode.id_str, 'D', "Is", diode.batchIs});
        }
        if (!diode.batchVT.empty()) {
            params.push_back({diode.id_str, 'D', "VT", diode.batchVT});
        }
    }

    return params;
}

// Compute the Cartesian-product size up front
std::size_t estimate_num_configs(const std::vector<BatchParam>& params)
{
    using limit = std::numeric_limits<std::size_t>;
    std::size_t product = 1;

    for (const BatchParam& p : params) {
        std::size_t n = p.values.size();          // how many values for this param?
        if (n == 0) continue;                     // empty vector ⇒ no influence

        // --- overflow guard -------------------------------------------------
        if (n > limit::max() / product) {         // would overflow size_t?
            throw std::overflow_error(
                "Sweep is too large to materialise (would overflow size_t)");
        }
        // --------------------------------------------------------------------
        product *= n;
    }
    return product;                               // 0 → no parameters
}

// Recursive helper function to generate all combinations of parameter values.
void generate_configs_recursive(
    const std::vector<BatchParam>& all_params,
    std::vector<CircuitConfig>& all_configs,
    CircuitConfig current_config,
    size_t param_index)
{
    // Base case: If we have processed all parameters, this configuration is complete.
    if (param_index == all_params.size()) {
        all_configs.push_back(current_config);
        return;
    }

    const auto& param = all_params[param_index];
    std::string config_key = param.device_id_str + "." + param.param_name;

    // Recursive step: For the current parameter, iterate through all its possible values.
    for (double value : param.values) {
        current_config[config_key] = value;
        generate_configs_recursive(all_params, all_configs, current_config, param_index + 1);
    }
}

// Generates a list of all possible circuit configurations to be simulated.
std::vector<CircuitConfig> generate_all_configs(const std::vector<BatchParam>& batch_params) {
    std::vector<CircuitConfig> all_configs;
    if (!batch_params.empty()) {
        const std::size_t n_cfg = estimate_num_configs(batch_params);
        all_configs.reserve(n_cfg);  // Reserve space for all configurations
        generate_configs_recursive(batch_params, all_configs, {}, 0);
    }
    return all_configs;
}

// Applies a specific configuration to a circuit by modifying its element parameters.
void apply_circuit_config(CircuitElements& CKTelements, const CircuitConfig& config) {
    for (const auto& [key, value] : config) {
        size_t dot_pos = key.find('.');
        if (dot_pos == std::string::npos) continue;

        std::string device_id = key.substr(0, dot_pos);
        std::string param_name = key.substr(dot_pos + 1);
        // BUGFIX: Convert the first character to uppercase to handle case-insensitivity.
        char device_type = std::toupper(device_id[0]);

        bool found = false;
        switch(device_type) {
            case 'V':
                for (auto& dev : CKTelements.voltageSources) {
                    if (dev.id_str == device_id) { dev.value = value; found = true; break; }
                }
                break;
            case 'I':
                for (auto& dev : CKTelements.currentSources) {
                    if (dev.id_str == device_id) { dev.value = value; found = true; break; }
                }
                break;
            case 'R':
                for (auto& dev : CKTelements.resistors) {
                    if (dev.id_str == device_id) { dev.value = value; found = true; break; }
                }
                break;
            case 'C':
                for (auto& dev : CKTelements.capacitors) {
                    if (dev.id_str == device_id) { dev.value = value; found = true; break; }
                }
                break;
            case 'G':
                 for (auto& dev : CKTelements.vccs) {
                    if (dev.id_str == device_id) { dev.value = value; found = true; break; }
                }
                break;
            case 'D':
                for (auto& dev : CKTelements.diodes) {
                    if (dev.id_str == device_id) {
                        if (param_name == "Is") dev.Is = value;
                        else if (param_name == "VT") dev.VT = value;
                        found = true; break;
                    }
                }
                break;
            case 'M':
                for (auto& dev : CKTelements.nmos) {
                    if (dev.id_str == device_id) {
                        if (param_name == "W") dev.W = value;
                        else if (param_name == "L") dev.L = value;
                        found = true; break;
                    }
                }
                if (found) break;
                for (auto& dev : CKTelements.pmos) {
                    if (dev.id_str == device_id) {
                        if (param_name == "W") dev.W = value;
                        else if (param_name == "L") dev.L = value;
                        found = true; break;
                    }
                }
                if (found) break;
                for (auto& dev : CKTelements.bsim4) {
                    if (dev.id_str == device_id) {
                        if (param_name == "W") dev.W = value;
                        else if (param_name == "L") dev.L = value;
                        found = true; break;
                    }
                }
                break;
        }
        if (!found) {
             std::cerr << "Warning: Device with ID " << device_id << " not found for batch update." << std::endl;
        }
    }
}

} // namespace batch
