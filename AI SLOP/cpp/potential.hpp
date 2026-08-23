#pragma once
#include "params.hpp"

// Attractive phi^4 potential from the simulation note:
//   V(phi) = 1/2 phi^2 - 1/4 phi^4.
// The quadratic term supplies the mass scale and the negative quartic term
// provides the self-focusing nonlinearity.
namespace physics {

inline double V_prime(double f) {
    return params::mass * params::mass * f - params::lambda * f * f * f;
}

inline double V(double f) {          // potential itself, for energy diagnostic
    return 0.5 * params::mass * params::mass * f * f
         - 0.25 * params::lambda * f * f * f * f;
}

} // namespace physics
