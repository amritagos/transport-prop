#include <tcf.hpp>
#include <iostream>

// Solvation time correlation function 
void test(Eigen::Ref<Eigen::RowVectorXd> energy_fluc, 
    Eigen::Ref<Eigen::RowVectorXd> time, 
    int max_tau, int start_t0, int start_tau, int delta_tau){
    //
    std::vector<double> tau_values; // Vector of lag times
    std::vector<double> tcf_values; // Vector of the unnormalized TCF
    double current_tau, current_tcf;  
    Eigen::RowVectorXd energy_fluc2(time.cols()); // Squared energy fluctuations 

    // for(auto x : time) std::cout << x <<"\n";
    // for (auto it = time.begin(); it < time.end(); std::advance(it, 2)){
    //     std::cout << *it <<"\n";
    // }

    for (auto it = time.begin()+start_t0; it < time.begin()+max_tau+1; std::advance(it, delta_tau)){
        std::cout << *it <<"\n";
    }

    // Calculation of the TCF at time=0
    tau_values.push_back( time(start_t0) ); // update the lag time
    // Square of the energy fluctuations
    energy_fluc2 = energy_fluc.array().square();
    // Mean of the squared energy fluctuations upto and including max_tau 
    current_tcf = energy_fluc2(Eigen::seq(start_t0, max_tau)).mean(); 
    tcf_values.push_back( current_tcf ); // Update the TCF vector  

    int j=0;
    std::cout<<"Number of data elements is "<<time.cols()<<"\n";
}
