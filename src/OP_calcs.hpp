#pragma once
#include <armadillo>
#include <string>
#include <vector>

struct MosfetOpData
{
    std::string name;
    std::string model;
    double id{}, vgs{}, vds{}, vbs{}, vth{}, vdsat{};
    double gm{}, gds{}, gmb{}, cbd{}, cbs{};
};

struct OPResult{
    arma::vec solution;
    std::vector<MosfetOpData> mosfet_data;
};
