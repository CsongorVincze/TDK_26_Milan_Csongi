#include "initial_conditions.hpp"
#include "params.hpp"
#include <cmath>

void set_initial_conditions(State& s) {
    s.phi.assign(params::N, 0.0);
    s.pi.assign(params::N, 0.0);

    for (int i = 0; i < params::N; ++i) {
        double r = params::radius(i);
        double w = params::width0;
        s.phi[i] = params::amp0 * std::exp(-0.5 * r * r / (w * w));
    }

    // The outer point is a fixed Dirichlet boundary in the solver.
    s.phi[params::N - 1] = 0.0;
}
