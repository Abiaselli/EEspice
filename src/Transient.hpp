#pragma once

#include "Transient_calcs.hpp"
#include "Transient_multi_solvers.hpp"
#include "Transient_single_solvers.hpp"
#include "CKTsetup.hpp"
#include "OP.hpp"

Transient Fixed_TimeStep(CKTcircuit &ckt, TransientSimulator &trans_sim, const Modelmap &modmap){
    Transient trans;
    trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation
    trans.h = trans_sim.trans_config.init_h; // Initial time step
    trans.time_trans = trans_sim.vec_trans.back().time_trans;

    std::cout << "Fixed time step" << std::endl;
    arma::vec solution = NewtonRaphson_system(ckt, trans.h, 1, trans.time_trans, trans_sim.vec_trans.back().solution, modmap);
    // solution.print("The transient analysis of the circuit is: ");
    ARMA_PRINT(solution, "The transient analysis of the circuit is: ");

    trans.time_trans = trans_sim.vec_trans.back().time_trans + trans.h;
    trans.solution = solution;

    trans.trans_count = trans_sim.vec_trans.back().trans_count + 1;
    return trans;
}

Transient Varibale_TimeStep(CKTcircuit &ckt, TransientSimulator &trans_sim, const Modelmap &modmap){
    Transient trans;
    trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation

    if (trans_sim.vec_trans.size() == 1)
    {
        // The first step of the transient simulation
        trans.h = trans_sim.trans_config.init_h; // Initial time step
        trans.time_trans = trans_sim.trans_config.t_start + trans.h;
        NIcomCof(ckt, trans.h);

        trans.solution = NewtonRaphson_system(ckt, trans.h, 1, trans.time_trans, trans_sim.vec_trans.back().solution, modmap);
        // solution.print("The transient analysis of the circuit is: ");
        ARMA_PRINT(trans.solution, "The transient analysis of the circuit is: ");

        trans.CapState = get_cap_state(ckt, trans.solution, trans.h, trans_sim.vec_trans);

        trans.trans_count += 1;
        trans.next_h = trans.h;
    }

    else
    {
        // Just to prevent single_next_h from directly using time_trans as the time_trans of the previous timestep (not recommended).
        trans.time_trans = trans_sim.vec_trans.back().time_trans; 

        single_timestep single_h = single_next_h(trans, ckt, trans_sim, modmap);
        
        // The current time of the transient simulation. Use to compare with the breakpoints!
        double current_time = single_h.h + trans_sim.vec_trans.back().time_trans;

        /*  If the time difference between the previous simulation time point and the breakpoint is less than 10*hmin,
            it will be counted as 1 time point and no additional breakpoint simulation will be performed.*/
        if (trans_sim.breakpoints.empty() == false && trans_sim.breakpoints.front() - trans_sim.vec_trans.back().time_trans < 10 * trans_sim.trans_config.h_MIN)
        {
            trans_sim.breakpoints.pop_front();
            if (trans_sim.breakpoints.empty())
            {
                trans_sim.trans_end = true;
            }
        }

        if (trans_sim.breakpoints.empty() == false && current_time > trans_sim.breakpoints.front())
        {
            Transient breakpoints_trans;
            breakpoints_trans.mode = 1; // 0 to do OP analysis, 1 to do transient simulation

            std::cout << "The time step is too large, back to breakpoint" << std::endl;
            breakpoints_trans.time_trans = trans_sim.breakpoints.front();
            breakpoints_trans.h = breakpoints_trans.time_trans - trans_sim.vec_trans.back().time_trans;

            breakpoints_trans.solution = NewtonRaphson_system(ckt, breakpoints_trans.h, 1, breakpoints_trans.time_trans, trans_sim.vec_trans.back().solution, modmap);

            breakpoints_trans.CapState = get_cap_state(ckt, breakpoints_trans.solution, breakpoints_trans.h, trans_sim.vec_trans);

            breakpoints_trans.trans_count = trans_sim.vec_trans.back().trans_count + 1;
            breakpoints_trans.next_h = breakpoints_trans.h;
            
            std::cout << "time step: " << breakpoints_trans.h << std::endl;
            std::cout << "time_trans: " << breakpoints_trans.time_trans << std::endl;

            trans_sim.breakpoints.pop_front();
            trans = std::move(breakpoints_trans);
        }

        else
        {
            if (trans_sim.breakpoints.empty() == false && current_time == trans_sim.breakpoints.front())
            {
                trans_sim.breakpoints.pop_front();
                if (trans_sim.breakpoints.empty())
                {
                    trans_sim.trans_end = true;
                }
            }

            // trans has already been updated in the multi_next_h function
            trans.h = single_h.h;
            trans.time_trans = trans_sim.vec_trans.back().time_trans + trans.h;
            trans.next_h = single_h.next_h;                          // Sometimes temp_h(next_h) sometimes is 2 * single_h.h!
            trans.solution = std::move(single_h.solution);
            trans.CapState = std::move(single_h.CapState);
            trans.trans_count = trans_sim.vec_trans.back().trans_count + 1;
        }
    }

    return trans;
}

std::vector<Transient> Transient_ops(CKTcircuit &ckt, TransientSimulator &trans_sim, const Modelmap modmap)
{
    ScopedTimer analysisTimer(ckt.sim_stats.simTime.analysis_time); // Time the analysis
    Transient trans_op;

    // LHS and RHS are stored in ckt.cktmatrix with baseline reset mechanism
    // No need to copy to trans_op - OperatingPointAnalysis uses ckt.cktmatrix directly
    trans_op.h = 0;

    // OPERATING POINT ANALYSIS SYSTEM
    trans_op.mode = 0;                          // 0 to do OP analysis, 1 to do transient simulation

    // OP analysis used as initial condition for next evaluation
    ckt.spiceCompatible.setFlagsTranOP();
    // Benchmarking for OP analysis
    auto tstart_op = std::chrono::high_resolution_clock::now();

    trans_op.solution = OperatingPointAnalysis(ckt, modmap, trans_sim.trans_config.non_linear);
    trans_op.next_h = trans_sim.trans_config.init_h;

    if (trans_sim.trans_config.timestep_control)
    {   
        // Get the capacitance state for the OP analysis if the time step control is on
        trans_op.CapState = get_cap_state(ckt, trans_op.solution, trans_op.h, trans_sim.vec_trans);
    }
    auto tstop_op = std::chrono::high_resolution_clock::now();


    if (debugMode == false && batchMode == false)
    {   
        printOperatingPointWithNames(trans_op.solution, ckt.map);
    }
    ARMA_PRINT(trans_op.solution, "The OP analysis (include internal nodes) of the circuit is: ");

    {
        ScopedTimer t(ckt.sim_stats.simTime.update_device_time);
        updateDeviceState(ckt);
    }
    {
        ScopedTimer t(ckt.sim_stats.simTime.history_update_time);
        history_trans_update(trans_op, trans_sim);
    }

    /*-----------------------------------------------------------*/
    //                Transient simulation start
    /*-----------------------------------------------------------*/
    auto tstart_trans = std::chrono::high_resolution_clock::now();
    std::cout << "transient simulation start" << std::endl;

    if (ckt.pulse_num != 0)
    {
        trans_sim.breakpoints = get_breakpoints(ckt, trans_sim); // Get the time of breakpoints for the transient simulation
        std::cout << "The breakpoints are: ";
        for (double num : trans_sim.breakpoints)
        {
            std::cout << num << " ";
        }
        std::cout << std::endl;
    }

    do
    {   ckt.spiceCompatible.setFlagsTR();
        // std::cout << "-------------------------------------------------------" << std::endl;
        // std::cout << "Next step" << std::endl;

        // Fixed time step
        if (trans_sim.trans_config.timestep_control == false)
        {
           Transient trans = Fixed_TimeStep(ckt, trans_sim, modmap);
           {
               ScopedTimer t(ckt.sim_stats.simTime.history_update_time);
               history_trans_update(trans, trans_sim);
           }
           {
               ScopedTimer t(ckt.sim_stats.simTime.update_device_time);
               updateDeviceState(ckt); // Copy state0 to state1
           }
        }
        // time step control: varibale time step
        else{
            Transient trans = Varibale_TimeStep(ckt, trans_sim, modmap);
            {
                ScopedTimer t(ckt.sim_stats.simTime.history_update_time);
                history_trans_update(trans, trans_sim);
            }
            {
                ScopedTimer t(ckt.sim_stats.simTime.update_device_time);
                updateDeviceState(ckt); // Copy state0 to state1
            }
        }

    } while (trans_sim.vec_trans.back().time_trans < trans_sim.trans_config.t_end && trans_sim.trans_end == false);

    ckt.sim_stats.num_data_points = static_cast<int>(trans_sim.vec_trans.size());

    return trans_sim.vec_trans;
}
