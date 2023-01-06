#ifndef TCF_H
#define TCF_H

#include <Eigen/Dense>

void test(Eigen::Ref<Eigen::RowVectorXd> energy_fluc, 
    Eigen::Ref<Eigen::RowVectorXd> time,
    int max_tau, int start_t0, int start_tau, int delta_tau);

#endif /* TCF_H */
