#ifndef DNest4_RVModel
#define DNest4_RVModel

#include "DNest4.h"
#include <vector>
#include <valarray>
#include <ostream>

class RVModel
{
	private:
		// planet params:
		double K;		// semi amplitude
		double T;		// period
		double phi;		// phase
		double ecc;		// eccentricity
		double omega;	// argument of periastron

		std::vector<double> gamma;
		std::vector<double> sigma;

		// Model prediction
		std::valarray<double> mu_model;

		// Compute the RVs
		void calculate_mu();

		// kepler solver
		double kepler_solve(double M, double e) const;

	public:
		// Constructor
		RVModel();
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

