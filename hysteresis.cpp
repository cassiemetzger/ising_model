#include "isingModel.h"

int main() {
    // Initialzing constants
    int Nx = 50;
    int Ny = 50;
    double J = 1.0;
    double KbT = 1.0;
    int numSteps = 5'000;
    int startAverage = 3'000;
    double H = 3.0;

    std::ofstream outfile("H-data.csv");
    std::vector<double> output(3);

    isingModel model(Nx, Ny, J, KbT,H);
    while(H > -3.0) {
        model.set_externalMag(H);
        model.run_model(numSteps);

        output = model.output(startAverage);
        // output[0] = average_magnetization
        // output[1] = average_susceptibility
        // output[2] = average_heat_capacity
        outfile << H << ", " << output[0] << std::endl;
        
        H -= 0.01;
        //std::cout << H << std::endl;
    }
    while(H < 3.0) {
        model.set_externalMag(H);
        model.run_model(numSteps);

        output = model.output(startAverage);
        // output[0] = average_magnetization
        // output[1] = average_susceptibility
        // output[2] = average_heat_capacity
        outfile << H << ", " << output[0] << std::endl;
        
        H += 0.01;
        //std::cout << H << std::endl;
    }

    return 0;
}