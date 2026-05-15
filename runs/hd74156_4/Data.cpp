#include "Data.h"
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

// The static instance
Data Data::instance;

Data::Data() : num_instruments(0)
{

}

void Data::load(const char* filename)
{
	// Vectors to hold the data
	std::vector<double> _x;
	std::vector<double> _y;
	std::vector<double> _err;
	std::vector<int> _instr;

	// Open the file
	fstream fin(filename, ios::in);

	// Temporary variables
	double temp1, temp2, temp3;
	int temp4;
	int max_index = -1;
	double mean_rv = 0.0;

	// Read until end of file
	while(fin>>temp1 && fin>>temp2 && fin>>temp3 && fin >> temp4)
	{
		_x.push_back(temp1);
		_y.push_back(temp2);
		_err.push_back(temp3);
		_instr.push_back(temp4);

		mean_rv += temp2;

		if(temp4 > max_index) {
			max_index = temp4;
		}
	}

	mean_rv /= _y.size();
	this->mean_rv = mean_rv;
	// Close the file
	fin.close();

	// Copy the data to the valarrays
	x = valarray<double>(&_x[0], _x.size());
	y = valarray<double>(&_y[0], _y.size());
	err = valarray<double>(&_err[0], _err.size());
	instr = valarray<int>(&_instr[0], _instr.size());

	num_instruments = max_index + 1;
}

