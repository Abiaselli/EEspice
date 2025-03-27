#pragma once
#include <iostream>
#include <string>
#include <unordered_map>
#include "bsim4v82.hpp"


namespace bsim4{

// Parameter map for BSIM4
const std::unordered_map<std::string, int> bsim4ParamMap = {
{"w", BSIM4_W},
{"l", BSIM4_L},
{"version", BSIM4_MOD_VERSION}
};

} // namespace bsim4