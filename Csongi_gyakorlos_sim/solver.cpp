#include <vector>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

// list of the tunable parameters
namespace params{
    double m = 1; // mass in the potential function
    double lambda = 1; // parameter of the phi^4 part
    double alpha = 4; // grid stretching parameter 
    double R_max = 10; // maximum distance
    int N = 500; // the number of the spatial points
    double dx = 1.0 / (N - 1); // the computational step
    double T = 10; // full time
    double CFL = 0.2; // CFL constant
    double dt = CFL * (alpha * R_max / (std::sinh(alpha))) * dx;
    double A = 0.4;
    double R0 = 2;

}

// stretching function
// transforms the computational grid to the actual physical grid
double r(double x){
    return params::R_max * std::sinh(params::alpha * x) / std::sinh(params::alpha);
}

// the physical distance differentiated with respect to the computational grid
// it is needed for the chain
double r_x(double x){
    return params::R_max * params::alpha * std::cosh(params::alpha * x) / std::sinh(params::alpha);
}

// second derivative
double r_x_2(double x){
    return params::R_max * params::alpha * params::alpha * std::sinh(params::alpha * x) / std::sinh(params::alpha);
}


// calculate the derivative vector of phi respect to ksi
// todo boundary
std::vector<double> phi_x(const std::vector<double>& phi){
    // todo return N elements
    std::vector<double> res;
    res.push_back(0);
    for(int i = 1; i < params::N - 1; ++i){
        res.push_back( (phi[i+1] - phi[i-1]) / (2*params::dx) );
    }
    return res;
}

// calculate the second derivative vector of phi respect to ksi
// todo boundary
std::vector<double> phi_x_2(const std::vector<double>& phi){
    std::vector<double> res;
    res.push_back(0);
    for(int i = 1; i < params::N - 1; ++i){
        res.push_back( (phi[i+1] - 2 * phi[i] + phi[i-1]) / (params::dx * params::dx) );
    }
    return res;
}




// We create a struct to hold the acutal state
// of the scalar field and it's time derivative
struct State{
    std::vector<double> phi;
    std::vector<double> pi;
};


// potential function (symmetric phi^4 model)
double V(double phi){
    return 0.5*params::m*params::m*phi*phi - 0.25*params::lambda*phi*phi*phi*phi;
}

// the derivative of the potential function with respect to phi
double V_phi(double phi){
    return params::m*params::m*phi - params::lambda*phi*phi*phi;
}

// we claculate the state a little timestep later and put it to dydt
void rhs(const State& y, State& dydt){
    // todo boundary
    std::vector<double> phi_x_ = phi_x(y.phi);
    std::vector<double> phi_x_2_ = phi_x_2(y.phi);

    // clear the previous calculations
    dydt.pi.assign(params::N, 0.0);
    dydt.phi.assign(params::N, 0.0);

    // calcualting the change at the origin
    dydt.phi[0] = y.pi[0];
    dydt.pi[0] = 6 / (r_x(0.0) * r_x(0.0)) * (y.phi[1] -y.phi[0]) / (params::dx * params::dx) -
    V_phi(y.phi[0]);


    for(int i = 1; i < params::N-1; ++i){
        double x = i * params::dx;
        dydt.phi[i] = y.pi[i];
        // todo this is ugly looking
        dydt.pi[i] = 1/(r_x(x) * r_x(x)) * phi_x_2_[i] + (2/(r(x)*r_x(x))
        - r_x_2(x)/(r_x(x) * r_x(x) * r_x(x))) * phi_x_[i] - V_phi(y.phi[i]);
    }

    // dirichlet condition
    dydt.pi[params::N-1] = 0.0;
    dydt.phi[params::N-1] = 0.0;
}


void rk4_step(State& y){

    State k1, k2, k3, k4, stage;
    double dt = params::dt;

    rhs(y, k1);

    stage.phi.resize(params::N);
    stage.pi.resize(params::N);

    for(int i = 0; i < params::N; ++i){
        stage.phi[i] = y.phi[i] + 0.5 * dt * k1.phi[i];
        stage.pi[i] = y.pi[i] + 0.5 * dt * k1.pi[i];
    }

    rhs(stage, k2);

    for(int i = 0; i < params::N; ++i){
        stage.phi[i] = y.phi[i] + 0.5 * dt * k2.phi[i];
        stage.pi[i] = y.pi[i] + 0.5 * dt * k2.pi[i];
    }

    rhs(stage, k3);

    for(int i = 0; i < params::N; ++i){
        stage.phi[i] = y.phi[i] +  dt * k3.phi[i];
        stage.pi[i] = y.pi[i] + dt * k3.pi[i];
    }

    rhs(stage, k4);

    // do one step with the rk4
    for(int i = 0; i < params::N; ++i){
        y.phi[i] += dt / 6 * (k1.phi[i] + 2.0 * k2.phi[i] + 2.0 * k3.phi[i] + k4.phi[i]);
        y.pi[i] += dt / 6 * (k1.pi[i] + 2.0 * k2.pi[i] + 2.0 * k3.pi[i] + k4.pi[i]);
    }

    // force the end to be zero
    y.phi[params::N - 1] = 0.0;
    y.pi[params::N - 1] = 0.0;
}

void write_snapshot(std::ofstream& out, const State& y, double t){
    out << std::setprecision(15);
    for(int i = 0; i < params::N; ++i){
        const double x = i * params::dx;
        out << t << ' ' << r(x) << ' ' << y.phi[i] << '\n';
    }
    out << '\n';
}

int main(){

    constexpr int output_every = 25;

    // initialize the profile
    State y;
    y.phi.assign(params::N, 0.0);
    y.pi.assign(params::N, 0.0);
    for(int i = 0; i < params::N; ++i){
        const double radius = r(i * params::dx);
        y.phi[i] = params::A * std::exp(
            -0.5 * radius * radius / (params::R0 * params::R0));
    }

    y.phi[params::N - 1] = 0.0;

    std::filesystem::create_directories("output");
    std::ofstream snapshots("output/snapshots.dat");
    if(!snapshots){
        std::cerr << "failed to open output/snapshots.dat\n";
        return 1;
    }

    const int n_steps = static_cast<int>(params::T / params::dt);
    for(int i = 0; i <= n_steps; ++i){
        const double t = i * params::dt;

        if(i % output_every == 0 || i == n_steps)
            write_snapshot(snapshots, y, t);

        if(i == n_steps)
            break;

        rk4_step(y);
    }

    std::cout << "wrote output/snapshots.dat\n";


    return 0;
}
