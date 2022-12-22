#include <pybind11/pybind11.h>

#include "add.h"

PYBIND11_MODULE(fastcpp, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
}
