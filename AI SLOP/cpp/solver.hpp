#pragma once
#include "fields.hpp"

// Right-hand side of the semi-discrete system:
//   d(phi_i)/dt = pi_i
//   d(pi_i)/dt  = laplacian(phi)_i - V'(phi_i) - gamma_i*pi_i
// on the stretched radial grid described in params.hpp.
void rhs(const State& y, State& dydt);

// One classic 4th-order Runge-Kutta step of the whole state, size dt.
void rk4_step(State& y);

// E = integral of [pi^2/2 + phi_r^2/2 + V(phi)] * 4*pi*r^2 dr.
// The diagnostic is expected to decrease after radiation reaches the sponge.
double total_energy(const State& y);
