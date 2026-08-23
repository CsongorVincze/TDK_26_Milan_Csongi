#include "solver.hpp"
#include "params.hpp"
#include "potential.hpp"

namespace {

State k1, k2, k3, k4, ys;   // RK4 stage scratch buffers

void rhs_impl(const State& y, State& dy) {
    const int N = params::N;
    dy.phi.resize(N);
    dy.pi.resize(N);

    for (int i = 0; i < N; ++i) {
        if (i == N - 1) {                     // outer boundary: hold phi ~ 0
            dy.phi[i] = 0.0;
            dy.pi[i] = 0.0;
            continue;
        }

        dy.phi[i] = y.pi[i];

        if (i == 0) {
            // Regularity at the origin gives laplacian(phi) = 3 phi_rr.
            // The ghost-point relation phi[-1] = phi[1] gives this stencil
            // directly in the computational coordinate xi.
            const double rp0 = params::radius_prime(0);
            const double laplacian = 6.0 * (y.phi[1] - y.phi[0])
                                   / (rp0 * rp0 * params::dxi * params::dxi);
            dy.pi[i] = laplacian
                       - physics::V_prime(y.phi[i]);
            continue;
        }

        // Chain rule for the stretched grid:
        // phi_rr = phi_xixi/r'^2 - r'' phi_xi/r'^3,
        // (2/r) phi_r = 2 phi_xi/(r r').
        const double phi_xi = (y.phi[i + 1] - y.phi[i - 1])
                            / (2.0 * params::dxi);
        const double phi_xixi = (y.phi[i + 1] - 2.0 * y.phi[i] + y.phi[i - 1])
                              / (params::dxi * params::dxi);
        const double r = params::radius(i);
        const double rp = params::radius_prime(i);
        const double rpp = params::radius_second(i);
        const double laplacian = phi_xixi / (rp * rp)
                               + (2.0 / (r * rp) - rpp / (rp * rp * rp)) * phi_xi;

        dy.pi[i] = laplacian - physics::V_prime(y.phi[i])
                   - params::sponge(r) * y.pi[i];
    }
}

} // namespace

void rhs(const State& y, State& dydt) {
    rhs_impl(y, dydt);
}

void rk4_step(State& y) {
    const int    N  = params::N;
    const double dt = params::dt;

    ys.phi.resize(N);
    ys.pi.resize(N);

    rhs(y, k1);
    for (int i = 0; i < N; ++i) {
        ys.phi[i] = y.phi[i] + 0.5 * dt * k1.phi[i];
        ys.pi[i]  = y.pi[i]  + 0.5 * dt * k1.pi[i];
    }

    rhs(ys, k2);
    for (int i = 0; i < N; ++i) {
        ys.phi[i] = y.phi[i] + 0.5 * dt * k2.phi[i];
        ys.pi[i]  = y.pi[i]  + 0.5 * dt * k2.pi[i];
    }

    rhs(ys, k3);
    for (int i = 0; i < N; ++i) {
        ys.phi[i] = y.phi[i] + dt * k3.phi[i];
        ys.pi[i]  = y.pi[i]  + dt * k3.pi[i];
    }

    rhs(ys, k4);
    for (int i = 0; i < N; ++i) {
        y.phi[i] += dt / 6.0 * (k1.phi[i] + 2.0 * k2.phi[i]
                               + 2.0 * k3.phi[i] + k4.phi[i]);
        y.pi[i]  += dt / 6.0 * (k1.pi[i] + 2.0 * k2.pi[i]
                               + 2.0 * k3.pi[i] + k4.pi[i]);
    }
}

double total_energy(const State& y) {
    const int N = params::N;
    double E = 0.0;
    auto energy_density = [&](int i) {
        const double r = params::radius(i);
        double phi_r = 0.0;
        if (i == N - 1) {
            const double dr = params::radius(i) - params::radius(i - 1);
            phi_r = (y.phi[i] - y.phi[i - 1]) / dr;
        } else if (i > 0) {
            phi_r = (y.phi[i + 1] - y.phi[i - 1])
                  / (2.0 * params::dxi * params::radius_prime(i));
        }
        const double rho = 0.5 * y.pi[i] * y.pi[i]
                         + 0.5 * phi_r * phi_r
                         + physics::V(y.phi[i]);
        return 4.0 * params::pi * r * r * rho;
    };

    for (int i = 0; i < N - 1; ++i) {
        const double dr = params::radius(i + 1) - params::radius(i);
        E += 0.5 * (energy_density(i) + energy_density(i + 1)) * dr;
    }
    return E;
}
