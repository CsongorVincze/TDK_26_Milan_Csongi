#pragma once
#include <vector>

// The dynamical state of the system on the radial grid.
struct State {
    std::vector<double> phi;   // field
    std::vector<double> pi;    // d(phi)/dt (conjugate momentum)
};
