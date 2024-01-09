#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "random.h"
#include <tuple>

class isingModel {
public:
    isingModel(int nx, int ny, double J, double KbT,double H=0):nx(nx), ny(ny), J(J), KbT(KbT){
        Nx = nx + 2;
        Ny = ny + 2;
        lattice.resize(Nx * Ny, 1);
        beta = 1.0 / (KbT);
        this->H=H;

        //fill lattice with random spins
        //DONT DO THIS because it is a lot less computationally heavy to just start with the highest energy state
        //since each flip is much more likely to be accepted than rejected
        /*for(int i = 0; i < Nx * Ny; i++) {
            if(rng.boolean()) {
                lattice[i] = 1;
            } else {
                lattice[i] = -1;
            }
        }*/
        periodic_boundary();
        calculate_energy();
        calculate_magnetization();
        std::cout << "Initial Energy: " << energy << std::endl;
        std::cout << "Initial Magnitization: " << magnetization << std::endl;
    }

    ~isingModel() {}

    //If you don't want to save the data, just don't pass a second argument
    void run_model(int timeSteps, std::string filename = "none") {
        energyList.resize(timeSteps);
        magnetizationList.resize(timeSteps); 
        for(int i = 0; i < timeSteps; i++) {
            bool updatedSpin = false;
            int random_i = rng.uniform_int(1, Nx-1);
            int random_j = rng.uniform_int(1, Ny-1);

            //calculate possible new energy of fllipped spin
            int spin = -lattice[random_i + random_j * Nx];
            double newEnergy = energy;
            newEnergy += -2 *J *lattice[(random_i - 1) + Nx * random_j] * spin-2.0*H*spin;
            newEnergy += -2 * J*lattice[(random_i + 1) + Nx * random_j] * spin-2.0*H*spin;
            newEnergy += -2 * J*lattice[random_i + Nx * (random_j - 1)] * spin-2.0*H*spin;
            newEnergy += -2 * J*lattice[random_i + Nx * (random_j + 1)] * spin-2.0*H*spin;
            double deltaE = newEnergy - energy;

            /*std::cout << "random i: " << random_i << std::endl;
            std::cout << "random j: " << random_j << std::endl;
            std::cout << "spin: " << spin << std::endl;
            std::cout << "spin left: " << lattice[(random_i - 1) + Nx * random_j] << std::endl;
            std::cout << "spin right: " << lattice[(random_i + 1) + Nx * random_j] << std::endl;
            std::cout << "spin up: " << lattice[random_i + Nx * (random_j - 1)] << std::endl;
            std::cout << "spin down: " << lattice[random_i + Nx * (random_j + 1)] << std::endl;
            std::cout << "energy: " << energy << std::endl;
            std::cout << "new energy: " << newEnergy << std::endl;
            std::cout << "deltaE: " << deltaE << std::endl;*/
            

            
            //check if we accept new energy
            if(newEnergy < energy) {
                updatedSpin = true;
            } else if (newEnergy > energy) {
                //want to check to accept based on probability
                if(rng.boolean(std::exp(-beta*deltaE))) {
                    updatedSpin = true;
                }
            }

            //std::cout << "updated spin: " << updatedSpin << std::endl;

            //if spin updated, update magnetization
            if(updatedSpin){
                //std::cout << "spin updated" << std::endl;
                magnetization += 2 * (double)spin/(nx * ny);
                lattice[random_i + random_j * Nx] = spin;
                energy = newEnergy;
            }

            //check boundaries to see if they need to be updated
            if(updatedSpin) {
                if(random_i == 1) {
                    lattice[(Nx - 1) + Nx * random_j] = spin;
                }
                if(random_i == Nx - 2) {
                    lattice[0 + Nx * random_j] = spin;
                }
                if(random_j == 1) {
                    lattice[random_i + Nx * (Nx - 1)] = spin;
                }
                if(random_j == Ny - 2) {
                    lattice[random_i + Nx * 0] = spin;
                }
            }

            //std::cout << "---------------------" << std::endl;
            energyList[i] = energy;
            magnetizationList[i] = magnetization;
            
        }

        if(filename != "none" && filename.substr(filename.length() - 4) == ".csv") {
            std::ofstream outfile(filename);
            for(int i = 0; i < timeSteps; i++) {
                outfile << i << ", " << energyList[i] << ", " << magnetizationList[i] << std::endl;
            }
            outfile.close();
        }

    }

    void reset() {
        for(int i = 0; i < Nx * Ny; i++) {
            lattice[i] = 1;
        }
        periodic_boundary();
        calculate_energy();
        calculate_magnetization();
    }

    void set_temperature(double KbT) {
        this->KbT = KbT;
        beta = 1.0 / (KbT);
    }
    void set_externalMag(double H ){
        this-> H=H;
    }

    void calculate_energy() {
        energy = 0.0;
        for(int i = 1; i < Nx - 1; i++) {
            for(int j = 1; j < Ny - 1; j++) {
                energy+=(-0.5)*(J)*(lattice[i+(j+1)*Nx]+lattice[i+(j-1)*Nx]+lattice[(i+1)+j*Nx]+lattice[(i-1)+j*Nx])-H*lattice[i+j*Nx];
            }
        }
    }

    //This is a very inefficient way to calculate magnetization
    //It is much better to just keep track of the magnetization as you go
    //This is only used for initializtion in the constructor
    //Same with calculate_energy and periodic_boundary which isn't even usefull now
    void calculate_magnetization() {
        magnetization = 0;
        for(int i = 1; i < Nx - 1; i++) {
            for(int j = 1; j < Ny - 1; j++) {
                magnetization += lattice[i+j*Nx];
            }
        }
        magnetization = magnetization / (nx * ny);
    }

    void periodic_boundary() {
        for(int j = 0; j < Ny; j++) {
            lattice[0 + j * Nx] = lattice[(Nx - 2) + j * Nx];
            lattice[(Nx - 1) + j * Nx] = lattice[1 + j * Nx];
        }

        for(int i = 0; i < Nx; i++) {
            lattice[i + 0 * Nx] = lattice[i + (Ny - 2) * Nx];
            lattice[i + (Ny - 1) * Nx] = lattice[i + 1 * Nx];
        }
    }

    /*
    All of these functions are used to calculate averages once the model has been run
    */

    std::tuple<double,double> calculate_average_magnetization(int startAverage) {
        double average = 0.0;
        double averagesquare = 0.0;
        for(int i = startAverage; i < magnetizationList.size(); i++) {
            average += magnetizationList[i];
            averagesquare += magnetizationList[i]*magnetizationList[i];
        }
        average = average / (magnetizationList.size() - startAverage);
        averagesquare = averagesquare/ (magnetizationList.size()-startAverage);
        double magnetizationVariance = averagesquare - average*average;
        std::tuple<double,double> result;
        result = std::make_tuple(average,magnetizationVariance);
        return result;
    }

    std::tuple<double,double> calculate_average_energy(int startAverage) {
        double average = 0.0;
        double averagesquare = 0.0;
        for(int i = startAverage; i < energyList.size(); i++) {
            average += energyList[i];
            averagesquare += energyList[i]*energyList[i];
        }
        average = average / (energyList.size() - startAverage);
        averagesquare = averagesquare / (energyList.size()-startAverage);
        //energy variance?
        double energyVariance = averagesquare - average*average;
        std::tuple<double,double> result;
        result = std::make_tuple(average,energyVariance);
        return result;
    }

    double calculate_average_susceptibility(int startAverage) {
        double variance = std::get<1>(calculate_average_magnetization(startAverage));
        // double average = calculate_average_magnetization(startAverage);
        // double variance = 0.0;
        // for(int i = startAverage; i < magnetizationList.size(); i++) {
        //     variance += std::pow(average - magnetizationList[i], 2);
        // }
        // //The -1 is because we are calculating the sample variance
        // variance = variance / (magnetizationList.size() - startAverage - 1);
        double susceptibility = beta * variance;
        return susceptibility;
    }

    double calculate_average_heat_capacity(int startAverage) {
        double variance = std::get<1>(calculate_average_energy(startAverage));
        // double average = calculate_average_energy(startAverage);
        // double variance = 0.0;
        // for(int i = startAverage; i < energyList.size(); i++) {
        //     variance += std::pow(average - energyList[i], 2);
        // }
        // //The -1 is because we are calculating the sample variance
        // variance = variance / (energyList.size() - startAverage - 1);
        double heatCapacity = beta * beta * variance;
        return heatCapacity;
    }
  
  
    std::vector<double> output(int startAverage, std::string filename = "none") {
        if(filename != "none") {
            std::ofstream outfile(filename);
            for (int j = 1; j < Ny - 1; j++) {
                for (int i = 1; i < Nx - 1; i++) {
                int idx = i + j * Nx;
                outfile << lattice[idx];
                if (i != Nx - 2)
                    outfile << ", ";
                else
                    outfile << std::endl;
                }
            }
            outfile.close();
        }
        
        //calculate_magnetization();
        double average_magnetization = std::get<0>(calculate_average_magnetization(startAverage));
        double average_susceptibility = calculate_average_susceptibility(startAverage);
        double average_heat_capacity = calculate_average_heat_capacity(startAverage);

        std::vector<double> output(3);
        output[0] = average_magnetization;
        output[1] = average_susceptibility;
        output[2] = average_heat_capacity;

        return output;
    }



    int Nx, Ny, nx, ny;
    double H;
    double energy, J, KbT, magnetization, beta;
    std::vector<int> lattice;
    std::vector<double> energyList, magnetizationList;
    random_number_generator rng;
};
