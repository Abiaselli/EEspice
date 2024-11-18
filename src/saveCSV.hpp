#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <variant>

#include "device.hpp"
#include "matrix.hpp"
#include "CKT.hpp"
#include "Transient.hpp"

void save_csv(const CKTcircuit &ckt, std::vector<Transient> &vec_trans)
{
    std::ofstream file("final_solution.csv");

    // Write the header
    file << "Time, Time Step";
    for (size_t j = 0; j < ckt.external_nodes; ++j)
    {
        file << ", Voltage " << (j + 1);
    }
    for (size_t z = ckt.no_of_V_sources; z > 0; --z)
    {
        file << ", Current " << (ckt.no_of_V_sources - z + 1);
    }
    // file << ", Current 2";
    file << std::endl;

    // Write the data
    for (size_t i = 0; i < vec_trans.size(); ++i)
    {
        file << std::scientific << std::setprecision(20);                // Set precision to 20 decimal places
        file << vec_trans.at(i).time_trans << ", " << vec_trans.at(i).h; // Time and Timestep

        for (size_t j = 0; j < ckt.external_nodes; ++j)
        {
            file << ", " << vec_trans.at(i).solution(j, 0); // Voltages
        }

        for (size_t z = ckt.no_of_V_sources; z > 0; --z)
        {
            file << ", " << vec_trans.at(i).solution(ckt.cktdematrix->n_rows - z, 0); // Currents
        }
        // file << ", " << vec_trans.at(i).C_current(0,0);
        file << std::endl;
    }

    file.close();
}