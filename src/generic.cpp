#include <generic.hpp>

// Finds the periodic distance between particles i and j 
double periodic_distance(Eigen::RowVectorXd r_i, 
    Eigen::RowVectorXd r_j,
    Eigen::RowVectorXd box){
    //
    int dim = box.cols(); // Number of dimensions 
    Eigen::VectorXd dr(dim); // Difference in the distance 
    double r2 = 0.0; // Squared absolute distance

    // Get x1-x2 for each dimension
    dr = r_i - r_j; 

    // Squaring and normalizing by the box lengths  
    for (int k = 0; k < dim; ++k)
    {
        dr(k) -= box(k) * std::round(dr(k) / box(k));
        r2 += pow(dr(k), 2.0);
    } // end of getting the relative distance

    return sqrt(r2);
}