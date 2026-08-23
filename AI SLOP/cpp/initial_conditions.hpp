#pragma once
#include "fields.hpp"

// Sets phi(r,0) = amp0 * exp(-r^2/width0^2), pi(r,0) = 0 (see params.hpp).
void set_initial_conditions(State& s);
