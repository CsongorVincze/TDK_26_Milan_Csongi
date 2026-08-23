// Driver: sets up the Gaussian, evolves it with RK4, writes results.
//
// Output files:
//   output/data/snapshots.dat   lines "t r phi", blank line between snapshots
//   output/data/diagnostics.dat lines "t phi(0,t) E(t)"
//   output/data/run_info.txt    exact parameters and run date
#include "params.hpp"
#include "fields.hpp"
#include "initial_conditions.hpp"
#include "solver.hpp"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace {

void write_run_info(std::ostream& out) {
    const std::time_t now = std::time(nullptr);
    const std::tm* local = std::localtime(&now);

    if (local)
        out << "date=" << std::put_time(local, "%Y-%m-%d") << '\n';
    else
        out << "date=unknown\n";

    out << std::setprecision(17)
        << "N=" << params::N << '\n'
        << "R_max=" << params::R_max << '\n'
        << "alpha=" << params::alpha << '\n'
        << "cfl=" << params::cfl << '\n'
        << "t_end=" << params::t_end << '\n'
        << "out_every=" << params::out_every << '\n'
        << "amp0=" << params::amp0 << '\n'
        << "width0=" << params::width0 << '\n'
        << "initial_pi=0\n"
        << "mass=" << params::mass << '\n'
        << "lambda=" << params::lambda << '\n'
        << "sponge_start_fraction=" << params::sponge_start_fraction << '\n'
        << "sponge_strength=" << params::sponge_strength << '\n'
        << "sponge_power=" << params::sponge_power << '\n'
        << "dxi=" << params::dxi << '\n'
        << "dr_min=" << params::dr_min << '\n'
        << "dt=" << params::dt << '\n';
}

} // namespace

int main() {
    State y;
    set_initial_conditions(y);

    std::ofstream snap("output/data/snapshots.dat");
    std::ofstream diag("output/data/diagnostics.dat");
    std::ofstream info("output/data/run_info.txt");
    if (!snap || !diag || !info) {
        std::cerr << "failed to open output files\n";
        return 1;
    }
    write_run_info(info);

    const int n_steps = static_cast<int>(std::floor(params::t_end / params::dt));
    int milestone = std::max(1, n_steps / 10);

    for (int step = 0; step <= n_steps; ++step) {
        double t = step * params::dt;

        if (step % params::out_every == 0 || step == n_steps) {
            for (int i = 0; i < params::N; ++i)
                snap << t << ' ' << params::radius(i) << ' ' << y.phi[i] << '\n';
            snap << '\n';
            diag << t << ' ' << y.phi[0] << ' ' << total_energy(y) << '\n';
        }

        if (step < n_steps)
            rk4_step(y);

        if (step % milestone == 0)
            std::cout << "t = " << t << '\n';
    }

    std::cout << "done: " << n_steps << " steps\n";
    return 0;
}
