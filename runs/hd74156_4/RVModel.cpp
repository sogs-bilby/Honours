#include "Data.h"
#include "RVModel.h"
#include "DNest4.h"
#include <cmath>
#include <string>

using namespace std;
using namespace DNest4;

RVModel::RVModel()
{

}

double RVModel::kepler_solve(double M, double e) const
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
// void RVModel::calculate_mu()
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
void RVModel::calculate_mu()
{
    const auto& t = Data::get_instance().get_x();
    if(mu_model.size() != t.size()) 
        mu_model.resize(t.size());
    
    for(size_t i = 0; i < t.size(); i++)    
    {
        
        // get mean anomaly
        double M_anom = (2.0 * M_PI * t[i] / T) + phi;
        
        // solve for ecc anomaly
        double E_anom = kepler_solve(M_anom, ecc);
        
        // convert ecc anomaly to true anomaly, f
        double f_true = 2.0 * atan(sqrt((1.0 + ecc) / (1.0 - ecc)) * tan(E_anom / 2.0));
        
        // calcualte keplerian RV
        mu_model[i] = K * (cos(omega + f_true) + ecc * cos(omega));
    }
}

void RVModel::from_prior(DNest4::RNG& rng)
{
    int num_instr = Data::get_instance().get_num_instruments();
    double mean_rv = Data::get_instance().get_mean_rv();

    gamma.resize(num_instr);
    sigma.resize(num_instr);

    for (int i = 0; i < num_instr; i++) {
        gamma[i] = mean_rv + 1000.0 * rng.randn();
        sigma[i] = exp(log(0.01) + log(10000.0) * rng.rand());
    }
    // semi-amplitude
    K = 1e2 * rng.rand();
    // log uniform from 1 to 10 days
    T = exp(log(1.0) + log(1000.0) * rng.rand());
    // phase from 0 to 2 pi
    phi = 2 * M_PI * rng.rand();

    ecc = 0.90 * rng.rand();
    omega = 2.0 * M_PI * rng.rand();

    calculate_mu();
}

double RVModel::perturb(RNG& rng)
{
    double log_H = 0.0;
    int num_instr = Data::get_instance().get_num_instruments();
    double mean_rv = Data::get_instance().get_mean_rv();
    // cout << to_string(mean_rv);
    // 5 standard parameters: semi-amplitude, period, phase, eccentricity and omega
    // n gamma offsets and n sigmas where n in number of instruments
    int which = rng.rand_int(8);

    if(which == 0) // offset
    {
        int idx = rng.rand_int(num_instr);
        log_H -= -0.5 * pow((gamma[idx] - mean_rv) / 1000.0, 2);
        gamma[idx] += 1000.0 * rng.randh();
        log_H += -0.5 * pow((gamma[idx] - mean_rv) / 1000.0, 2);
    }
    else if(which == 1) // semi amplitude
    {
        K += 100.0 * rng.randh();
        wrap(K, 0.0, 100.0);
    }
    else if(which == 2) // period
    {
        T = log(T); 
        T += log(1000.0) * rng.randh(); 
        wrap(T, log(1.0), log(1000.0)); 
        T = exp(T);
    }
    else if(which == 3) // phi
    {
        phi += 2.0 * M_PI * rng.randh();
        wrap(phi, 0.0, 2.0 * M_PI);
    }
    else if(which == 4) // eccentricity
    {
        ecc += 0.90 * rng.randh();
        wrap(ecc, 0.0, 0.90);
    }
    else if(which == 5) // omega
    {
        omega += 2.0 * M_PI * rng.randh();
        wrap(omega, 0.0, 2.0 * M_PI);
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

double RVModel::log_likelihood() const
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

void RVModel::print(std::ostream& out) const
{
    // write all 7 parameters to sample.txt
    for (size_t i = 0; i < gamma.size(); i++) {
        out << gamma[i] << " ";
    }
    out << ' ' << K << ' ' << T << ' ' << phi << ' ' << ecc << ' ' << omega << ' ';

    for (size_t i = 0; i < (sigma.size()); i++) {
        out << sigma[i] << " ";
    }
}

string RVModel::description() const
{
    string descr = "";
    int num_instr = Data::get_instance().get_num_instruments();

    for (int i = 0; i < num_instr; i++) {
        descr += "gamma_" + to_string(i) + ", ";
    }
    descr += "K, T, phi, ecc, omega, ";
    for (int i = 0; i < num_instr; i++) {
        descr += "sigma_" + to_string(i) + ", ";
    }
    return descr;
}