#include "Data.h"
#include "MPRVModel.h"
#include "DNest4.h"
#include <cmath>
#include <Eigen/Core>
#include <celerite2/celerite2.h>
#include <celerite2/terms.hpp>

#include <Eigen/Dense>

int MPRVModel::num_planets = 2;

using namespace std;
using namespace DNest4;

MPRVModel::MPRVModel()
{

}

double MPRVModel::kepler_solve(double M, double e) const
{
    // mean anomaly between 0 and 2 pi
    // M = fmod(M, 2.0 * M_PI);
    // if (M < 0.0) M += 2.0 * M_PI;
    wrap(M, 0, 2 * M_PI);

    // guess initially
    // eccentric anomaly is E
    double E = M; 
    if (e > 0.8) E = M_PI; 

    // newton root finding loop
    for (int i = 0; i < 50; i++) 
    {
        // M = E - e sin(E)

        double f = E - e * sin(E) - M;
        double f_prime = 1.0 - e * cos(E);
        double E_next = E - f / f_prime;

        if (std::abs(E_next - E) < 1e-10) {
            return E_next;
        }
        E = E_next;
    }
    return E;
}

void MPRVModel::calculate_mu()
{
    const auto& t = Data::get_instance().get_x();
    if(mu_model.size() != t.size()) 
        mu_model.resize(t.size());
    
    for(size_t i = 0; i < t.size(); i++)    
    {
        mu_model[i] = 0.0;

        for(int j = 0; j < num_planets; j++) {
            // get mean anomaly
            double M_anom = (2.0 * M_PI * t[i] / T[j]) + phi[j];
            
            // solve for ecc anomaly
            double E_anom = kepler_solve(M_anom, ecc[j]);
            
            // convert ecc anomaly to true anomaly, f
            double f_true = 2.0 * atan(sqrt((1.0 + ecc[j]) / (1.0 - ecc[j])) * tan(E_anom / 2.0));
            
            // calcualte keplerian RV
            mu_model[i] += K[j] * (cos(omega[j] + f_true) + ecc[j] * cos(omega[j]));
        }
    }
}

void MPRVModel::from_prior(DNest4::RNG& rng)
{
    int num_instr = Data::get_instance().get_num_instruments();
    double mean_rv = Data::get_instance().get_mean_rv();

    gamma.resize(num_instr);
    sigma.resize(num_instr);

    for (int i = 0; i < num_instr; i++) {
        gamma[i] = mean_rv + 1000.0 * rng.randn();
        sigma[i] = exp(log(0.01) + log(10000.0) * rng.rand());
    }
    K.resize(num_planets);
    T.resize(num_planets);
    phi.resize(num_planets);
    ecc.resize(num_planets);
    omega.resize(num_planets);

    for(int j=0; j<num_planets; j++) {
        K[j] = exp(log(0.1) + log(100.0) * rng.rand());
        T[j] = exp(log(0.1) + log(1000.0) * rng.rand());
        phi[j] = 2.0 * M_PI * rng.rand();
        while(true)
        {
            ecc[j] = 0.1*tan(M_PI*(rng.rand() - 0.5));
            if(ecc[j] >= 0.0 && ecc[j] < 1.0)
                break;
        }
        omega[j] = 2.0 * M_PI * rng.rand();
    }

    log_S0 = -20.0 + rng.rand() * 40.0;
    log_w0 = -5.0 + rng.rand() * 10.0;
    log_Q = -0.34 + rng.rand() * 6.0;
    
    calculate_mu();
}

double MPRVModel::perturb(RNG& rng)
{
    double log_H = 0.0;
    int num_instr = Data::get_instance().get_num_instruments();
    double mean_rv = Data::get_instance().get_mean_rv();
    // cout << to_string(mean_rv);
    // 5 standard parameters: semi-amplitude, period, phase, eccentricity and omega
    // n gamma offsets and n sigmas where n in number of instruments
    int which = rng.rand_int(10); // now 10 to match the 3 new SHO terms
    int j = rng.rand_int(num_planets); // choose which planet's params to perturb

    if(which == 0) // offset
    {
        int idx = rng.rand_int(num_instr);
        log_H -= -0.5 * pow((gamma[idx] - mean_rv) / 1000.0, 2);
        gamma[idx] += 1000.0 * rng.randh();
        log_H += -0.5 * pow((gamma[idx] - mean_rv) / 1000.0, 2);
    }
    else if(which == 1) // semi amplitude
    {
        K[j] = log(K[j]);
        K[j] += log(100.0) * rng.randh();
        wrap(K[j], log(0.1), log(100.0));
        K[j] = exp(K[j]);
    }
    else if(which == 2) // period
    {
        T[j] = log(T[j]); 
        T[j] += log(1000.0) * rng.randh(); 
        wrap(T[j], log(0.1), log(1000.0)); 
        T[j] = exp(T[j]);
        
        // reject if periods cross
        // if(num_planets == 2 && T[0] > T[1]) {
        //     return -1e300;
        // }
    }
    else if(which == 3) // phi
    {
        phi[j] += 2.0 * M_PI * rng.randh();
        wrap(phi[j], 0.0, 2.0 * M_PI);
    }
    else if(which == 4) // eccentricity
    {
        log_H -= -log(1.0 + pow(ecc[j]/0.1, 2));
        ecc[j] += rng.randh();
        wrap(ecc[j], 0.0, 1.0);
        log_H += -log(1.0 + pow(ecc[j]/0.1, 2));
    }
    else if(which == 5) // omega
    {
        omega[j] += 2.0 * M_PI * rng.randh();
        wrap(omega[j], 0.0, 2.0 * M_PI);
    }
    else if(which == 6) // sigma
    {
        int idx = rng.rand_int(num_instr);
        sigma[idx] = log(sigma[idx]);
        sigma[idx] += log(10000.0) * rng.randh(); 
        wrap(sigma[idx], log(0.01), log(100.0));
        sigma[idx] = exp(sigma[idx]);
    }
    else if(which == 7) // log_S0
    {
        log_S0 += 40.0 * rng.randh();
        wrap(log_S0, -20.0, 20.0);
    }
    else if(which == 8) // log_w0
    {
        log_w0 += 10.0 * rng.randh();
        wrap(log_w0, -5.0, 5.0);
    }
    else if(which == 9) // log_Q
    {
        log_Q += 6.0 * rng.randh();
        wrap(log_Q, -0.34, 5.66); 
    }

    if(which < 6)
        calculate_mu();

    return log_H;
}

// double MPRVModel::log_likelihood() const
// {
//     const auto& y = Data::get_instance().get_y();
//     const auto& err = Data::get_instance().get_err();
//     const auto& instr = Data::get_instance().get_instr();

//     double logl = 0.0;
// 	for(size_t i=0; i < y.size(); ++i) {
//         double mu = mu_model[i] + gamma[instr[i]];
//         double var = (err[i] * err[i]) + (sigma[instr[i]] * sigma[instr[i]]);
// 		logl += -0.5*log(2*M_PI*var) - 0.5*pow(y[i] - mu, 2)/var;
//     }
// 	return logl;
// }

double MPRVModel::log_likelihood() const
{
    const auto& t = Data::get_instance().get_x();
    const auto& y = Data::get_instance().get_y();
    const auto& err = Data::get_instance().get_err();
    const auto& instr = Data::get_instance().get_instr();

    int N = t.size();
    
    Eigen::VectorXd t_arr(N);
    Eigen::VectorXd res_arr(N);
    Eigen::VectorXd diag_arr(N);

    for(int i = 0; i < N; i++) {
        t_arr(i) = t[i];
        double mu = mu_model[i] + gamma[instr[i]];
        res_arr(i) = y[i] - mu;
        // White noise variance (instrument error^2 + stellar jitter^2)
        diag_arr(i) = (err[i] * err[i]) + (sigma[instr[i]] * sigma[instr[i]]);
    }

    double S0 = exp(log_S0);
    double w0 = exp(log_w0);
    double Q = exp(log_Q);

    // Initialize the templated SHO term
    celerite2::SHOTerm<double> term(S0, w0, Q);
    
    // This returns a std::tuple containing: {c, a, U, V}
    auto matrices = term.get_celerite_matrices(t_arr, diag_arr);
    
    // Unpack the tuple references cleanly without copying any memory
    const auto& c = std::get<0>(matrices);
    const auto& a = std::get<1>(matrices);
    const auto& U = std::get<2>(matrices);
    const auto& V = std::get<3>(matrices);
    
    // Match the exact Row-Major layout and fixed width (2) provided by celerite2
    Eigen::Matrix<double, Eigen::Dynamic, 2, Eigen::RowMajor> W(N, 2);
    Eigen::VectorXd d(N);
    
    // Define the execution workspaces explicitly as Row-Major
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> O_N_factor_work(N, 4); // J*J = 4
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> O_N_solve_work(N, 2);  // J = 2

    // O(N) Cholesky Factorization
    int flag = celerite2::core::factor(t_arr, c, a, U, V, d, W, O_N_factor_work);
    
    // Catch mathematical anomalies and reject the step safely
    if (flag != 0 || (d.array() <= 0.0).any()) {
        return -1e300; 
    }

    // O(N) Solver: (L * D * L^T) * alpha = res_arr
    Eigen::VectorXd alpha = res_arr;
    celerite2::core::solve_lower(t_arr, c, U, W, alpha, alpha, O_N_solve_work);
    alpha = alpha.array() / d.array();
    celerite2::core::solve_upper(t_arr, c, U, W, alpha, alpha, O_N_solve_work);

    // Final Log-Likelihood combination
    double log_det = d.array().log().sum();
    double chi_squared = res_arr.dot(alpha);
    
    return -0.5 * chi_squared - 0.5 * log_det - 0.5 * N * std::log(2.0 * M_PI);
}

void MPRVModel::print(std::ostream& out) const
{
    // write all 7 parameters to sample.txt
    // offsets
    for (size_t i = 0; i < gamma.size(); i++) {
        out << gamma[i] << " ";
    }

    // params
    for (int j = 0; j < num_planets; j++) {
        out << ' ' << K[j] << ' ' << T[j] << ' ' << phi[j] << ' ' << ecc[j] << ' ' << omega[j] << ' ';
    }

    // sigmas
    for (size_t i = 0; i < (sigma.size()); i++) {
        out << sigma[i] << " ";
    }

    // SHO terms
    out << log_S0 << " " << log_w0 << " " << log_Q << " ";
}

string MPRVModel::description() const
{
    string descr = "";
    int num_instr = Data::get_instance().get_num_instruments();

    // offsets
    for (int i = 0; i < num_instr; i++) {
        descr += "gamma_" + to_string(i) + ", ";
    }

    // params
    for (int j = 0; j < num_planets; j++) {
        descr += "K_" + to_string(j) + ", T_" + to_string(j) + ", phi_" + to_string(j);
        descr += ", ecc_" + to_string(j) + ", omega_" + to_string(j) + ", ";
    }

    // sigmas
    for (int i = 0; i < num_instr; i++) {
        descr += "sigma_" + to_string(i) + ", ";
    }
    
    // SHO terms
    descr += "log_S0, log_w0, log_Q, ";

    return descr;
}