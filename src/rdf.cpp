#include <rdf.hpp>
#include <iostream>
#include <generic.hpp>

//Radial distribution function for a single time step 
// Eigen::Ref<Eigen::MatrixXd>
Eigen::RowVectorXd calc_rdf(Eigen::MatrixXd pos_A, 
    Eigen::MatrixXd pos_B,
    Eigen::RowVectorXd box,
    double binwidth, int nbin, double cutoff){
    //
    Eigen::RowVectorXd hist = Eigen::RowVectorXd::Zero(nbin); // Number of unnormalized counts per timestep; histogram (nbin)
    Eigen::RowVectorXd rdf_val = Eigen::RowVectorXd::Zero(nbin); // RDF normalized by the bin volume and solvent density 
    int nop_A = pos_A.rows(); // Number of the central particle A 
    int nop_B = pos_B.rows(); // Number of the solvent particles B  
    int dim = box.cols();
    double r_ij; // Unwrapped distance between particles i and j
    int i_bin; // Current bin in which the particle falls  
    double total_volume = 1; // Simulation box volume
    double rho_B; // Number density of B in the simulation box  
    double bin_volume; // Bin volume of a particular ith bin 
    //
    // Loop through the central particles A
    for (int i = 0; i < nop_A; ++i)
    {
        // Loop through all molecules of type B 
        for (int j = 0; j < nop_B; ++j)
        {
            // Calculate the distance of jatom from iatom 
            r_ij = periodic_distance(pos_A.row(i), pos_B.row(j), box);
            // If the distance is within the cutoff, update the 
            // histogram count 
            if (r_ij <= cutoff)
            {
                // Find the bin in which the particle count falls
                i_bin = int(r_ij/binwidth); 
                // Update the histogram
                hist(i_bin) += 1;
            } // end of check for r_ij within the cutoff 
        } // end of loop through B particles 
    } // end of loop through particles A 

    // Normalize the histogram counts by the local density of B
    // and the number of particles of A.
    // 
    // Get the simulation box volume 
    for (int k = 0; k < dim; ++k)
     {
         total_volume *= box(k);
     } // end of getting the simulation box volume 
     // 
     // Getting the volume density of B
     rho_B = nop_B/total_volume;

     // Loop through all the bins 
     for (int i = 0; i < nbin; ++i)
     {
        // Bin volume of the ith bin 
        bin_volume = 4.0*PI*( pow(binwidth,3.0) )*(pow(i+1,3)-pow(i,3))/3.0;
        // Update the RDF
        rdf_val(i) = hist(i)/(rho_B*nop_A*bin_volume);
     } // end of loop through all bins 

    return rdf_val;
} // end of rdf function