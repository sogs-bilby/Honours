#ifndef DNest4_Data
#define DNest4_Data

#include <valarray>

/*
* An object of this class is a dataset
*/
class Data
{
	private:
		// The data points
		std::valarray<double> x;
		std::valarray<double> y;
		std::valarray<double> err; // measurement error
		std::valarray<int> instr;

		double mean_rv;
		int num_instruments;

	public:
		// Constructor
		Data();

		// Load data from a file
		void load(const char* filename);

		// Access to the data points
		const std::valarray<double>& get_x() const
		{ return x; }
		const std::valarray<double>& get_y() const
		{ return y; }
		const std::valarray<double>& get_err() const
		{ return err; }
		const std::valarray<int>& get_instr() const
		{ return instr; }

		int get_num_instruments() const
		{ return num_instruments; }
		double get_mean_rv() const
		{ return mean_rv; }

	private:
		// Static "global" instance
		static Data instance;

	public:
		// Getter for the global instance
		static Data& get_instance()
		{ return instance; }
};

#endif

