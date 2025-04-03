#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <Eigen/Dense>

#include "tcf.hpp"
#include "rdf.hpp"

PYBIND11_MODULE(fastcpp, m) {
    m.doc() = "Module for correlation function loops"; // module docstring

    m.def("time_corr_function", &time_corr_function, "Given an energy of energy difference fluctuations, calculates the time correlation function.");
    m.def("calc_rdf", &calc_rdf, "Calculate the radial distribution function for a single time step.");
}
