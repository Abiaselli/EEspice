#pragma once
#include <iostream>
#include "map.hpp"
#include "bsim4v82/bsim4v82setup.hpp"
#include "bsim4v82/bsim4v82temp.hpp"
#include "sim_variables.hpp"

// 1. setup the golbal parameters for the model (modelSetup)
// 2. setup the temperature dependent parameters for the model (modelTemp)
// 3. check the model (modelTemp->modelCheck)

void modelSetup(Modelmap &modmap, double ckttemp){
    for(auto &model : modmap.bsim4Models){
        // Setup the model using the default temperature
        bsim4::modelSetup(*model.second, nomTemp);
        // Setup the model using the actual temperature
        bsim4::modelTemp(*model.second, ckttemp);
    }
}