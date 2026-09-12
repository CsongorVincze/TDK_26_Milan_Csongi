#include <vector>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

// list of the tunable parameters
namespace params{
    double m = 1; // mass in the potential function
    double lambda = 1; // parameter of the phi^4 part
    double alpha = 4; // grid stretching parameter 
    double R_max = 10; // maximum distance
    double R_sp = 7.5; // start of the sponge layer
    double G_0 = 7.0; // Gamma_0, sponge coeff.
    double S_exp = 2.0; // the exponent in the sponge layer
    int N = 2000; // the number of the spatial points
    double dx = 1.0 / (N - 1); // the computational step
    double T = 10; // full time
    double CFL = 0.2; // CFL constant
    double dt = CFL * (alpha * R_max / (std::sinh(alpha))) * dx;
    double A = 0.4;
    double R0 = 2;
    int output_every = 25;
    int origin_every = 25;

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
    std::vector<double> res(params::N, 0.0);
    for(int i = 1; i < params::N - 1; ++i){
        res[i] = (phi[i+1] - phi[i-1]) / (2*params::dx);
    }
    return res;
}

// calculate the second derivative vector of phi respect to ksi
// todo boundary
std::vector<double> phi_x_2(const std::vector<double>& phi){
    std::vector<double> res(params::N, 0.0);
    for(int i = 1; i < params::N - 1; ++i){
        res[i] = (phi[i+1] - 2 * phi[i] + phi[i-1])
               / (params::dx * params::dx);
    }
    return res;
}




// We create a struct to hold the acutal state
// of the scalar field and it's time derivative
struct State{
    std::vector<double> phi;
    std::vector<double> pi;
};


// potential function (original symmetric phi^4 model)
double V(double phi){
    return 0.5*params::m*params::m*phi*phi - 0.25*params::lambda*phi*phi*phi*phi;
}

// the derivative of the potential function with respect to phi
double V_phi(double phi){
    return params::m*params::m*phi - params::lambda*phi*phi*phi;
}

double integrate_uniform(const std::vector<double>& values){
    if(values.size() < 2)
        return 0.0;

    const int last = static_cast<int>(values.size()) - 1;
    double result = 0.0;

    // Simpson's rule over an even number of intervals.  With the current
    // default N=500, the final interval is integrated with a trapezoid.
    const int simpson_last = (last % 2 == 0) ? last : last - 1;
    result = values[0] + values[simpson_last];
    for(int i = 1; i < simpson_last; ++i)
        result += (i % 2 == 0 ? 2.0 : 4.0) * values[i];
    result *= params::dx / 3.0;

    if(simpson_last != last)
        result += 0.5 * params::dx * (values[last - 1] + values[last]);
    return result;
}

double total_energy(const State& y){
    std::vector<double> energy_density(params::N, 0.0);
    for(int i = 0; i < params::N; ++i){
        const double x = i * params::dx;
        double phi_r = 0.0;
        if(i == params::N - 1){
            phi_r = (y.phi[i] - y.phi[i - 1])
                  / (params::dx * r_x(x));
        } else if(i > 0){
            const double phi_x = (y.phi[i + 1] - y.phi[i - 1])
                               / (2.0 * params::dx);
            phi_r = phi_x / r_x(x);
        }

        const double radius = r(x);
        const double density = 0.5 * y.pi[i] * y.pi[i]
                             + 0.5 * phi_r * phi_r
                             + V(y.phi[i]);
        energy_density[i] = density * radius * radius * r_x(x);
    }

    return 4.0 * std::acos(-1.0) * integrate_uniform(energy_density);
}

//sponge layer
double sponge(double r){
    if(r < params::R_sp){
        return 0.0;
    }
    const double fraction = (r - params::R_sp)
                          / (params::R_max - params::R_sp);
    return params::G_0 * std::pow(fraction, params::S_exp);
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
    dydt.pi[0] = 6.0 / (r_x(0.0) * r_x(0.0))
               * (y.phi[1] - y.phi[0]) / (params::dx * params::dx)
               - V_phi(y.phi[0]);


    for(int i = 1; i < params::N-1; ++i){
        double x = i * params::dx;
        dydt.phi[i] = y.pi[i];
        // todo this is ugly looking
        dydt.pi[i] = 1/(r_x(x) * r_x(x)) * phi_x_2_[i] + (2/(r(x)*r_x(x))
        - r_x_2(x)/(r_x(x) * r_x(x) * r_x(x))) * phi_x_[i] - V_phi(y.phi[i]) - sponge(r(x)) * y.pi[i];
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

std::string run_date() {
    const std::time_t now = std::time(nullptr);
    const std::tm* local = std::localtime(&now);
    if (!local)
        return "unknown";

    std::ostringstream result;
    result << std::put_time(local, "%Y-%m-%d");
    return result.str();
}

std::string run_name(int number, double amplitude, double radius0) {
    std::ostringstream name;
    name << "run_" << std::setw(3) << std::setfill('0') << number
         << "_A" << std::fixed << std::setprecision(2) << amplitude
         << "_R0" << std::fixed << std::setprecision(2) << radius0;

    std::string result = name.str();
    std::replace(result.begin(), result.end(), '.', 'p');
    return result;
}

void write_run_info(std::ofstream& out, const std::string& id,
                    const std::string& date) {
    out << "run_id=" << id << '\n'
        << "date=" << date << '\n'
        << std::setprecision(17)
        << "m=" << params::m << '\n'
        << "lambda=" << params::lambda << '\n'
        << "alpha=" << params::alpha << '\n'
        << "R_max=" << params::R_max << '\n'
        << "R_sp=" << params::R_sp << '\n'
        << "G_0=" << params::G_0 << '\n'
        << "S_exp=" << params::S_exp << '\n'
        << "N=" << params::N << '\n'
        << "dx=" << params::dx << '\n'
        << "T=" << params::T << '\n'
        << "CFL=" << params::CFL << '\n'
        << "dt=" << params::dt << '\n'
        << "A=" << params::A << '\n'
        << "R0=" << params::R0 << '\n'
        << "output_every=" << params::output_every << '\n'
        << "origin_every=" << params::origin_every << '\n'
        << "energy_every=" << params::origin_every << '\n'
        << "initial_pi=0\n";
}

State initial_state() {
    State y;
    y.phi.assign(params::N, 0.0);
    y.pi.assign(params::N, 0.0);

    for(int i = 0; i < params::N; ++i){
        const double radius = r(i * params::dx);
        y.phi[i] = params::A * std::exp(
            -0.5 * radius * radius / (params::R0 * params::R0));
    }

    y.phi[params::N - 1] = 0.0;
    return y;
}

class SolverTimer {
public:
    explicit SolverTimer(const std::string& id)
        : id_(id), start_(std::chrono::steady_clock::now()) {}

    ~SolverTimer() {
        const std::chrono::duration<double> elapsed =
            std::chrono::steady_clock::now() - start_;
        std::cout << id_ << ": solver runtime = "
                  << std::fixed << std::setprecision(3)
                  << elapsed.count() << " s" << std::endl;
    }

private:
    std::string id_;
    std::chrono::steady_clock::time_point start_;
};

bool run_simulation(const std::filesystem::path& directory,
                    const std::string& id, const std::string& date,
                    bool print_progress) {
    const SolverTimer timer(id);
    std::filesystem::create_directories(directory);
    std::ofstream snapshots(directory / "snapshots.dat");
    std::ofstream info(directory / "run_info.txt");
    std::ofstream origin(directory / "origin.dat");
    std::ofstream energies(directory / "energy.dat");
    if(!snapshots || !info || !origin || !energies){
        std::cerr << "failed to open output files in " << directory << '\n';
        return false;
    }
    write_run_info(info, id, date);

    State y = initial_state();
    const int n_steps = static_cast<int>(params::T / params::dt);
    for(int i = 0; i <= n_steps; ++i){
        const double t = i * params::dt;

        if(i % params::origin_every == 0 || i == n_steps)
            origin << std::setprecision(15) << t << ' ' << y.phi[0] << '\n';

        if(i % params::origin_every == 0 || i == n_steps)
            energies << std::setprecision(15) << t << ' '
                     << total_energy(y) << '\n';

        if(i % params::output_every == 0 || i == n_steps)
            write_snapshot(snapshots, y, t);

        if(i == n_steps)
            break;

        rk4_step(y);

        if (!std::all_of(y.phi.begin(), y.phi.end(),
                         [](double value) { return std::isfinite(value); }) ||
            !std::all_of(y.pi.begin(), y.pi.end(),
                         [](double value) { return std::isfinite(value); })) {
            std::cerr << id << ": non-finite value at step " << i << '\n';
            return false;
        }
    }

    if (print_progress)
        std::cout << "wrote " << (directory / "snapshots.dat") << '\n';
    return true;
}

int run_single(const std::filesystem::path& output_directory) {
    return run_simulation(output_directory, "single", run_date(), true) ? 0 : 1;
}

int run_sweep() {
    // Literature-informed exploratory grid.  The numerical parameters remain
    // fixed; only the physical initial amplitude and radius are swept here.
    const std::vector<double> amplitudes = {0.4, 0.8, 1.2, 1.6, 2.0};
    const std::vector<double> radii = {2.0, 2.5, 2.86, 3.5, 5.0};
    const std::string date = run_date();
    const std::filesystem::path sweep_dir = "output/sweep";

    std::filesystem::create_directories(sweep_dir);
    std::ofstream manifest(sweep_dir / "manifest.csv");
    if(!manifest){
        std::cerr << "failed to open " << (sweep_dir / "manifest.csv") << '\n';
        return 1;
    }

    manifest << "run_id,A,R0,m,lambda,alpha,R_max,R_sp,G_0,S_exp,N,dx,T,CFL,dt,"
             << "output_every,status,relative_path\n"
             << std::setprecision(17);

    int run_number = 1;
    bool all_completed = true;
    for (double amplitude : amplitudes) {
        for (double radius0 : radii) {
            params::A = amplitude;
            params::R0 = radius0;
            const std::string id = run_name(run_number, amplitude, radius0);
            const std::filesystem::path directory = sweep_dir / id;
            const bool completed = run_simulation(directory, id, date, false);
            all_completed = all_completed && completed;

            manifest << id << ',' << params::A << ',' << params::R0 << ','
                     << params::m << ',' << params::lambda << ','
                     << params::alpha << ',' << params::R_max << ','
                     << params::R_sp << ',' << params::G_0 << ','
                     << params::S_exp << ','
                     << params::N << ',' << params::dx << ',' << params::T << ','
                     << params::CFL << ',' << params::dt << ','
                     << params::output_every << ','
                     << (completed ? "ok" : "nonfinite") << ',' << id << '\n';
            std::cout << id << ": " << (completed ? "done" : "failed") << '\n';
            ++run_number;
        }
    }

    std::cout << "wrote " << (sweep_dir / "manifest.csv") << '\n';
    if (!all_completed)
        std::cout << "completed with nonfinite cases; see manifest.csv\n";
    return 0;
}

int main(int argc, char** argv){
    if (argc > 1 && std::string(argv[1]) == "sweep")
        return run_sweep();
    if (argc > 1 && std::string(argv[1]) != "single") {
        std::cerr << "usage: " << argv[0]
                  << " [single [A R0 T output_every origin_every output_dir N]]\n";
        return 2;
    }

    try {
        if (argc > 2) params::A = std::stod(argv[2]);
        if (argc > 3) params::R0 = std::stod(argv[3]);
        if (argc > 4) params::T = std::stod(argv[4]);
        if (argc > 5) params::output_every = std::stoi(argv[5]);
        if (argc > 6) params::origin_every = std::stoi(argv[6]);
        if (params::output_every <= 0 || params::origin_every <= 0)
            throw std::invalid_argument("output intervals must be positive");
        if (argc > 7) {
            std::filesystem::path output_directory = argv[7];
            if (argc > 8) params::N = std::stoi(argv[8]);
            if (params::N < 3)
                throw std::invalid_argument("N must be at least 3");
            params::dx = 1.0 / (params::N - 1);
            params::dt = params::CFL
                       * (params::alpha * params::R_max / std::sinh(params::alpha))
                       * params::dx;
            return run_single(output_directory);
        }
    } catch (const std::exception& error) {
        std::cerr << "invalid command-line parameter: " << error.what() << '\n';
        return 2;
    }

    return run_single("output");
}
