#include "Data.h"
#include "MPRVModel.h"
#include "DNest4.h"
#include <cmath>
#include <string>

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

// OLD caulculate_mu()
// void MPRVModel::calculate_mu()
// {
//     const auto& t = Data::get_instance().get_x();
//     const auto& set_idx = Data::get_instance().get_instrument();
//     if(mu_model.size() != t.size()) 
//         mu_model.resize(t.size());
    
//     for(size_t i = 0; i < t.size(); i++)    
//     {
//         int inst = static_cast<int>(set_idx[i]);
//         mu_model[i] = gamma[inst] + K * cos(2.0 * M_PI * t[i] / T + phi);
//     }
// }

// NEW calculate_mu()
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
    // 2. Initialize Planets
    K.assign(num_planets, 0.0);
    T.assign(num_planets, 0.0);
    phi.assign(num_planets, 0.0);
    ecc.assign(num_planets, 0.0);
    omega.assign(num_planets, 0.0);

    for(int j=0; j<num_planets; j++) {
        K[j] = 100.0 * rng.rand();
        T[j] = exp(log(1.0) + log(10000.0) * rng.rand());
        phi[j] = 2.0 * M_PI * rng.rand();
        ecc[j] = 0.9 * rng.rand();
        omega[j] = 2.0 * M_PI * rng.rand();
    }
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
    int which = rng.rand_int(7);
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
        K[j] += 100.0 * rng.randh();
        wrap(K[j], 0.0, 100.0);
    }
    else if(which == 2) // period
    {
        T[j] = log(T[j]); 
        T[j] += log(10000.0) * rng.randh(); 
        wrap(T[j], log(1.0), log(10000.0)); 
        T[j] = exp(T[j]);
        
        // reject if periods cross
        if(num_planets == 2 && T[0] > T[1]) {
            return -1e300;
        }
    }
    else if(which == 3) // phi
    {
        phi[j] += 2.0 * M_PI * rng.randh();
        wrap(phi[j], 0.0, 2.0 * M_PI);
    }
    else if(which == 4) // eccentricity
    {
        ecc[j] += 0.90 * rng.randh();
        wrap(ecc[j], 0.0, 0.90);
    }
    else if(which == 5) // omega
    {
        omega[j] += 2.0 * M_PI * rng.randh();
        wrap(omega[j], 0.0, 2.0 * M_PI);
    }
    else if(which >= 6) // sigma
    {
        int idx = rng.rand_int(num_instr);
        sigma[idx] = log(sigma[idx]);
        sigma[idx] += log(10000.0) * rng.randh(); 
        wrap(sigma[idx], log(0.01), log(100.0));
        sigma[idx] = exp(sigma[idx]);
    }

    if(which < 6)
        calculate_mu();

    return log_H;
}

double MPRVModel::log_likelihood() const
{
    const auto& y = Data::get_instance().get_y();
    const auto& err = Data::get_instance().get_err();
    const auto& instr = Data::get_instance().get_instr();

    double logl = 0.0;
	for(size_t i=0; i < y.size(); ++i) {
        double mu = mu_model[i] + gamma[instr[i]];
        double var = (err[i] * err[i]) + (sigma[instr[i]] * sigma[instr[i]]);
		logl += -0.5*log(2*M_PI*var) - 0.5*pow(y[i] - mu, 2)/var;
    }
	return logl;
}

void MPRVModel::print(std::ostream& out) const
{
    // write all 7 parameters to sample.txt
    for (size_t i = 0; i < gamma.size(); i++) {
        out << gamma[i] << " ";
    }

    for (int j = 0; j < num_planets; j++) {
        out << ' ' << K[j] << ' ' << T[j] << ' ' << phi[j] << ' ' << ecc[j] << ' ' << omega[j] << ' ';
    }

    for (size_t i = 0; i < (sigma.size()); i++) {
        out << sigma[i] << " ";
    }
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
    return descr;
}