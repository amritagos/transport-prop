#ifndef GENERIC_H
#define GENERIC_H

#include <Eigen/Dense>
#include <math.h>

const double PI = atan(1.0)*4;

// Given the box lengths, finds the periodic distances
// between a particle i and particle j
double periodic_distance(Eigen::RowVectorXd r_i, 
    Eigen::RowVectorXd r_j,
    Eigen::RowVectorXd box);

#endif /* GENERIC_H */