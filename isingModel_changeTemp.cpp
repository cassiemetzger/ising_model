#include "isingModel.h"

int main() {
    // Initialzing constants
    int Nx = 50;
    int Ny = 50;
    double J = 1;
    double KbT = 0.01;
    int numSteps = 5'000'000;
    int startAverage = 4'000'000;

    std::ofstream outfile("temperature-data.csv");
    std::vector<double> output(3);

    isingModel model(Nx, Ny, J, KbT);
    while(KbT < 4.0) {
        model.set_temperature(KbT);
        model.run_model(numSteps);

        output = model.output(startAverage);
        // output[0] = average_magnetization
        // output[1] = average_susceptibility
        // output[2] = average_heat_capacity
        outfile << KbT << ", " << std::abs(output[0]) << ", " << output[1] << ", " << output[2] << std::endl;
        
        KbT += 0.01;
        model.reset();
        std::cout << KbT << std::endl;
    }

    return 0;
}