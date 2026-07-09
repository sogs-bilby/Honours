#include <iostream>
#include "Data.h"
#include "DNest4.h"
#include "MPRVModel.h"

using namespace std;
using namespace DNest4;

int main(int argc, char** argv)
{
	// Read num_planets from config.json before starting anything else
    ifstream config_file("config.json");
    string line;
    
    if (config_file.is_open()) {
        while (getline(config_file, line)) {
            if (line.find("\"num_planets\"") != string::npos) {
                size_t colon_pos = line.find(":");
                size_t comma_pos = line.find(",");
                
                // If it's the last line in the JSON, there might not be a comma
                if (comma_pos == string::npos) {
                    comma_pos = line.find("}");
                }
                
                string val_str = line.substr(colon_pos + 1, comma_pos - colon_pos - 1);
                MPRVModel::num_planets = stoi(val_str); // std::stoi automatically handles whitespace
                
                cout << "Config: Setting num_planets = " << MPRVModel::num_planets << " from config.json" << endl;
                break;
            }
        }
        config_file.close();
    } else {
        cout << "Warning: Could not open config.json. Defaulting to " << MPRVModel::num_planets << " planets." << endl;
    }

	Data::get_instance().load("data.txt");

    cout << "mean rv: " << Data::get_instance().get_mean_rv() << endl;

	start<MPRVModel>(argc, argv);
	return 0;
}

