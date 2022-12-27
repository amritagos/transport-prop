#include <tcf.hpp>
#include <iostream>

// Solvation time correlation function 
void test(Eigen::Ref<Eigen::RowVectorXd> energy_fluc, 
    Eigen::Ref<Eigen::RowVectorXd> time, 
    int max_tau){
    //
    // for(auto x : time) std::cout << x <<"\n";
    int j=0;
    std::cout<<"Number of data elements is "<<time.cols()<<"\n";
}
