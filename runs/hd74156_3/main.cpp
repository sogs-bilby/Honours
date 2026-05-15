#include <iostream>
#include "Data.h"
#include "DNest4.h"
#include "RVModel.h"

using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
	Data::get_instance().load("data.txt");

    cout << "mean rv: " << Data::get_instance().get_mean_rv() << endl;

	start<RVModel>(argc, argv);
	return 0;
}

