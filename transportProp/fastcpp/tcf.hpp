#ifndef TCF_H
#define TCF_H

#include <Eigen/Dense>

void test(Eigen::Ref<Eigen::RowVectorXd> energy_fluc, 
    Eigen::Ref<Eigen::RowVectorXd> time,
    int max_tau);

#endif /* TCF_H */
