#ifndef DNest4_MPRVModel
#define DNest4_MPRVModel

#include "DNest4.h"
#include <vector>
#include <valarray>
#include <ostream>

class MPRVModel
{
	private:
		int num_planets = 2;
    
		// instrument params
		std::vector<double> gamma;
		std::vector<double> sigma;

		// planetary params
		std::vector<double> K;
		std::vector<double> T;
		std::vector<double> phi;
		std::vector<double> ecc;
		std::vector<double> omega;

		void calculate_mu();

		// kepler solver
		double kepler_solve(double M, double e) const;

	public:
		// Constructor
		MPRVModel();
		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);
		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);
		// Likelihood function
		double log_likelihood() const;
		// Print to stream
		void print(std::ostream& out) const;
		// Return string with column information
		std::string description() const;
};

#endif

