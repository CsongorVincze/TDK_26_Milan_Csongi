#include <vector>
#include <cmath>

// list of the tunable parameters
namespace params{
    double m = 1; // mass in the potential function
    double lambda = 2; // parameter of the phi^4 part
    double alpha = 4; // grid stretching parameter 
    double R_max = 10; // maximum distance
    int N = 500; // the number of the spatial points
    double dx = 1 / (N - 1); // the computational step
    double T = 10; // full time
    double CFL = 0.2; // CFL constant
    double dt = CFL * (alpha * R_max / (std::sinh(alpha))) * dx;

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
std::vector<double> phi_x(std::vector<double>& phi){
    std::vector<double> res;
    res.push_back(0);
    for(int i = 1; i < params::N - 1; ++i){
        res.push_back( (phi[i+1] - phi[i-1]) / (2*params::dx) );
    }
    return res;
}

// calculate the second derivative vector of phi respect to ksi
// todo boundary
std::vector<double> phi_x_2(std::vector<double>& phi){
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
    return 0.5*params::m*params::m*phi*phi + 0.25*params::lambda*phi*phi*phi*phi;
}

// the derivative of the potential function with respect to phi
double V_phi(double phi){
    return params::m*params::m*phi + params::lambda*phi*phi*phi;
}

void rhs(State& y){
    // todo boundary
    std::vector<double> phi_x_ = phi_x(y.phi);
    std::vector<double> phi_x_2_ = phi_x_2(y.phi);
    for(int i = 1; i < params::N; ++i){
        double x = i * params::dx;
        double pi_t = 1/(r_x(x) * r_x(x)) * phi_x_2_[i] + (2/(r(x)*r_x(x)) - r_x_2(x)/(r_x(x) * r_x(x) * r_x(x))) * phi_x_[i] - V_phi(y.phi[i]);
    }
}