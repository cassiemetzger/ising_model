#include "isingModel.h"

int main() {

    // Initialzing constants
    int Nx = 100;
    int Ny = 100;
    double J = 1.0;
    double KbT = 2.0;

    std::vector<double> output(3);

    isingModel model(Nx, Ny, J, KbT);
    model.run_model(5000000, "timestep-data-2.csv");
    output = model.output(1000000, "lattice-2KbT.csv");
    std::cout << "Average Magnetization: " << output[0] << std::endl;
    model.reset();

    KbT = 3.0;

    model.set_temperature(KbT);
    model.run_model(5000000, "timestep-data-3.csv");
    output = model.output(1000000, "lattice-3KbT.csv");
    std::cout << "Average Magnetization: " << output[0] << std::endl;
    model.reset();


    return 0;
}
