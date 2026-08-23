#pragma once
#include <cmath>

// All tunable simulation parameters live here.
namespace params {

// --- spatial grid ---
inline constexpr int    N     = 400;           // number of radial grid points
inline constexpr double R_max = 30.0;          // outer boundary radius
inline constexpr double alpha = 5.0;           // hyperbolic grid stretching
inline constexpr double dxi   = 1.0 / (N - 1);

// --- absorbing sponge ---
inline constexpr double sponge_start_fraction = 0.75;
inline constexpr double sponge_strength       = 8.0;
inline constexpr int    sponge_power          = 2;

// The hyperbolic map packs points near r = 0 and spreads them toward R_max:
//   r(xi) = R_max sinh(alpha*xi) / sinh(alpha), 0 <= xi <= 1.
inline const double dr_min = R_max * alpha / std::sinh(alpha) * dxi;

inline double radius(int i) {
    return R_max * std::sinh(alpha * i * dxi) / std::sinh(alpha);
}

inline double radius_prime(int i) {
    return R_max * alpha * std::cosh(alpha * i * dxi) / std::sinh(alpha);
}

inline double radius_second(int i) {
    return R_max * alpha * alpha * std::sinh(alpha * i * dxi) / std::sinh(alpha);
}

inline double sponge(double r) {
    const double start = sponge_start_fraction * R_max;
    if (r <= start)
        return 0.0;
    const double x = (r - start) / (R_max - start);
    return sponge_strength * std::pow(x, sponge_power);
}

// --- time integration ---
inline constexpr double cfl       = 0.25;
inline const double dt             = cfl * dr_min;
inline constexpr double t_end     = 60.0;
inline constexpr int    out_every = 100;       // snapshot every this many steps

// --- initial Gaussian profile, placed at the center r = 0 ---
// phi(r, 0) = amp0 * exp(-r^2/(2*width0^2)),  pi(r, 0) = 0
inline constexpr double amp0   = 0.4;
inline constexpr double width0 = 2.0;

// --- potential: V(phi) = 1/2 mass^2 phi^2 - 1/4 lambda phi^4 ---
inline constexpr double mass   = 1.0;
inline constexpr double lambda = 1.0;

inline const double pi = std::acos(-1.0);

} // namespace params
