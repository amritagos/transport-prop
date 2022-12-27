#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

#include "tcf.hpp"

PYBIND11_MODULE(fastcpp, m) {
    m.doc() = "Module for correlation function loops"; // module docstring

    m.def("test", &test, "Given an energy of energy difference fluctuations, calculates the time correlation function.");
}
