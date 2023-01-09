#include <tcf.hpp>
#include <iostream>

// Solvation time correlation function 
std::tuple<std::vector<double>, std::vector<double>, double>
time_corr_function(Eigen::Ref<Eigen::RowVectorXd> energy_fluc, 
    Eigen::Ref<Eigen::RowVectorXd> time, 
    int max_tau, int start_t0, int start_tau, int delta_tau){
    //
    std::vector<double> tau_values; // Vector of lag times
    std::vector<double> tcf_values; // Vector of the unnormalized TCF
    double current_tau, current_tcf;  
    Eigen::RowVectorXd energy_fluc2(time.cols()); // Squared energy fluctuations
    int n_origins = max_tau - start_t0; // Number of time origins 
    // Time origins to average over, for each lag time tau
    Eigen::RowVectorXd tcf_origins(n_origins);  
    int t0, t, counter; 
    double mean_t0; // For the first time origin, mean over the entire range 

    // for(auto x : time) std::cout << x <<"\n";
    // for (auto it = time.begin(); it < time.end(); std::advance(it, 2)){
    //     std::cout << *it <<"\n";
    // }

    // Calculation of the TCF at time=0
    tau_values.push_back( time(start_t0) ); // update the lag time
    // Square of the energy fluctuations
    energy_fluc2 = energy_fluc.array().square();
    // Mean of the squared energy fluctuations upto and including max_tau 
    current_tcf = energy_fluc2(Eigen::seq(start_t0, max_tau)).mean(); 
    tcf_values.push_back( current_tcf ); // Update the TCF vector
    // Mean over the entire range (will return later)
    mean_t0 = energy_fluc2(Eigen::seq(start_t0, Eigen::last)).mean();

    // Calculation of the TCF at different lag times,
    // starting from start_tau 
    for (int i_tau = start_tau; i_tau <= max_tau; i_tau+=delta_tau)
    {
        current_tau = time(i_tau); // Current lag time 
        // Get the TCF for this lag time over the desired time origins
        counter = 0; // initialize 
        // Loop through time origins 
        for (auto& i_origin : tcf_origins)
        {
            t0 = counter + start_t0; // current time origin
            t = t0 + i_tau; // time t = t0 + lag time 

            // Calculate TCF for t0 and tau
            i_origin = energy_fluc(t0) * energy_fluc(t);

            ++counter; // Increment counter 
        } // end of loop through time origins at i_tau

        // The current TCF is the average over all time origins
        current_tcf = tcf_origins.mean();

        // Update the TCF and time vectors
        tau_values.push_back(current_tau); 
        tcf_values.push_back(current_tcf);
    } // end of loop through lag times 

    // Return the lag times, TCF and the mean over the entire time series
    // for C(0). The value in the TCF vector is uptil max_tau 
    return std::make_tuple(tau_values, tcf_values, mean_t0);

}
