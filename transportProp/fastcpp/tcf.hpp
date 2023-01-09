#ifndef TCF_H
#define TCF_H

#include <Eigen/Dense>
#include <tuple>

std::tuple<std::vector<double>, std::vector<double>, double>
time_corr_function(Eigen::Ref<Eigen::RowVectorXd> energy_fluc, 
    Eigen::Ref<Eigen::RowVectorXd> time,
    int max_tau, int start_t0, int start_tau, int delta_tau);

#endif /* TCF_H */
