#ifndef RDF_H
#define RDF_H

#include <Eigen/Dense>
#include <tuple>

Eigen::RowVectorXd calc_rdf(Eigen::MatrixXd pos_A, 
    Eigen::MatrixXd pos_B,
    Eigen::RowVectorXd box,
    double binwidth, int nbin, double cutoff);

#endif /* RDF_H */