#pragma once
#include <iostream>
#include "map.hpp"
#include "bsim4v82/bsim4v82setup.hpp"
#include "sim_variables.hpp"

void modelSetup(Modelmap &modmap){
    for(auto &model : modmap.bsim4Models){
        bsim4::modelSetup(*model.second, nomTemp);
    }
    
}