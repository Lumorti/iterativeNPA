#include "optim.h"
#include "utils.h"

// Optim
#define OPTIM_ENABLE_EIGEN_WRAPPERS
#include "optim.hpp"

// Import Eigen
#include <Eigen/Dense>

// The structure used to hold the data for the optimisation
struct optimData {
    double penalty;
    Poly objective;
    std::vector<Poly> equalityCons;
};

// Cost/gradient function for optim
static double gradFunction(const Eigen::VectorXd& x, Eigen::VectorXd* gradOut, void* optData) {
    
    // Recast this generic pointer into the correct format
    optimData* optDataRecast = reinterpret_cast<optimData*>(optData);

    std::map<Mon, std::complex<double>> xMap;
    for (size_t i=0; i<x.size(); i++) {
        xMap[Mon("<X" + std::to_string(i) + ">")] = std::complex<double>(x(i), 0.0);
    }

    //std::cout << "X = " << x << std::endl;

    // The cost is a combination of the objective and the equality constraints
    double cost = 0.0;
    cost -= optDataRecast->objective.eval(xMap).real();
    for (size_t i=0; i<optDataRecast->equalityCons.size(); i++) {
        //std::cout << optDataRecast->equalityCons[i] << std::endl;
        //std::cout << optDataRecast->equalityCons[i].eval(xMap) << std::endl;
        cost += optDataRecast->penalty * std::pow(optDataRecast->equalityCons[i].eval(xMap).real(), 2);
    }

    // Calculate the gradient
    if (gradOut) {

        // Reset the gradient
        *gradOut = Eigen::VectorXd::Zero(x.size());

    }

    // Return just the cost
    return cost;

}

double solveOptim(Poly objective, std::vector<Poly> cons) {

    // Get the list of variables
    int highestXInd = 0;
    for (auto const& [mon, coeff] : objective.polynomial) {
        for (auto const& part : mon.monomial) {
            if (part.first == 'X') {
                highestXInd = std::max(highestXInd, part.second);
            }
        }
    }
    for (size_t i=0; i<cons.size(); i++) {
        for (auto const& [mon, coeff] : cons[i].polynomial) {
            for (auto const& part : mon.monomial) {
                if (part.first == 'X') {
                    highestXInd = std::max(highestXInd, part.second);
                }
            }
        }
    }
    int numVars = highestXInd + 1;

    // Set up the optimisation
    Eigen::VectorXd x = Eigen::VectorXd::Zero(numVars);
    optimData optData;
    optData.objective = objective;
    optData.equalityCons = cons;
    optData.penalty = 1e3;

    // Check the initial cost
    Eigen::VectorXd gradData = Eigen::VectorXd::Zero(numVars);
    double startCost = gradFunction(x, &gradData, &optData);
    std::cout << "Starting cost: " << startCost << std::endl;

    // Settings for the optimiser
    optim::algo_settings_t settings;
    settings.print_level = 1;
    settings.iter_max = 10000;
    settings.rel_objfn_change_tol = 1e-10;
    settings.gd_settings.method = 6;
    settings.lbfgs_settings.par_M = 10;
    settings.lbfgs_settings.wolfe_cons_1 = 1e-3;
    settings.lbfgs_settings.wolfe_cons_2 = 0.8;
    settings.vals_bound = true;
    settings.lower_bounds = -1*Eigen::VectorXd::Ones(numVars);
    settings.upper_bounds = 1*Eigen::VectorXd::Ones(numVars);

    //optData.penalty = 1e2;
    //bool success = optim::nm(x, gradFunction, &optData, settings);

    // Optimise
    optData.penalty = 1e3;
    for (int i=0; i<10; i++) {
        std::cout << "Optimising with " << numVars << " variables and penalty " << optData.penalty << std::endl;
        //bool success = optim::lbfgs(x, gradFunction, &optData, settings);
        bool success = optim::nm(x, gradFunction, &optData, settings);
        std::cout << "Result = " << gradFunction(x, &gradData, &optData) << std::endl;
        optData.penalty *= 2;
    }

    optData.penalty = 0;
    double result = gradFunction(x, &gradData, &optData);
    std::cout << "Result = " << result << std::endl;

    return 0.0;

}
