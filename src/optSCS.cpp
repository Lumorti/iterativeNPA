#include "optSCS.h"
#include "utils.h"
#include "printing.h"
#include <iostream>

// Import SCS
#include "scs.h"

// Convert to SCS form and solve
double solveSCS(Poly obj, std::vector<std::vector<std::vector<Poly>>>& psd, std::vector<Poly> constraintsZero, std::vector<Poly> constraintsPositive, int verbosity, std::pair<int,int> varBounds, std::map<Mon, std::complex<double>>* solution) {

    // Get the list of variables
    int oneIndex = 0;
    std::set<Mon> variableSet;
    variableSet.insert(Mon());
    for (size_t i=0; i<psd.size(); i++) {
        addVariables(variableSet, psd[i]);
    }
    for (size_t i=0; i<constraintsZero.size(); i++) {
        addVariables(variableSet, constraintsZero[i]);
    }
    for (size_t i=0; i<constraintsPositive.size(); i++) {
        addVariables(variableSet, constraintsPositive[i]);
    }
    addVariables(variableSet, obj);
    std::vector<Mon> variables = toVector(variableSet);

    // Output the variable list
    if (verbosity >= 2) {
        std::cout << "Num variables: " << variables.size() << std::endl;
    }
    if (verbosity >= 3) {
        std::cout << "Variables:" << std::endl;
        for (size_t i=0; i<variables.size(); i++) {
            std::cout << i << " " << variables[i] << std::endl;
        }
        std::cout << std::endl;
    }

    // Cache the variable locations
    std::map<Mon, int> variableLocs;
    for (size_t i=0; i<variables.size(); i++) {
        variableLocs[variables[i]] = i;
    }

    // The c vector defining the objective
    std::vector<double> c(variables.size());
    for (auto& term : obj.polynomial) {

        // Find the location of this variable
        int varLoc = variableLocs[term.first];

        // Add the real and imaginary terms
        c[varLoc] -= std::real(term.second);

    }

    // Number of variables and constraints
    int n = variables.size();
    int numBoxCons = n;
    int m = 1 + constraintsZero.size() + constraintsPositive.size() + numBoxCons;
    for (size_t i=0; i<psd.size(); i++) {
        m += psd[i].size() * (psd[i].size() + 1) / 2;
    }

    // The list of constants
    std::vector<double> bVals(m, 0.0);

    // The A matrix defining all constraints
    std::vector<int> ARows;
    std::vector<int> ACols;
    std::vector<double> AVals;

    // First the zero constraints
    int numConsSoFar = 0;
    for (size_t i=0; i<constraintsZero.size(); i++) {

        // For each term in the constraint
        for (auto& term : constraintsZero[i]) {

            // Find the location of this variable
            int realLoc = variableLocs[term.first];
            double realVal = std::real(term.second);
            if (std::abs(realVal) > 1e-14) {
                ARows.push_back(i + numConsSoFar);
                ACols.push_back(realLoc);
                AVals.push_back(realVal);
            }

        }

    }
    numConsSoFar += constraintsZero.size();

    // Constraint one element to be one
    ARows.push_back(numConsSoFar);
    ACols.push_back(oneIndex);
    AVals.push_back(1.0);
    bVals[constraintsZero.size()] = 1.0;
    numConsSoFar++;

    // Then the positive constraints
    for (size_t i=0; i<constraintsPositive.size(); i++) {

        // For each term in the constraint
        for (auto& term : constraintsPositive[i]) {

            // Find the location of this variable
            int realLoc = variableLocs[term.first];
            double realVal = std::real(term.second);
            if (std::abs(realVal) > 1e-14) {
                ARows.push_back(numConsSoFar + i);
                ACols.push_back(realLoc);
                AVals.push_back(realVal);
            }

        }

    }
    numConsSoFar += constraintsPositive.size();

    // Then box constraints
    for (int i=0; i<numBoxCons; i++) {
        ARows.push_back(numConsSoFar + i);
        ACols.push_back(i);
        AVals.push_back(1.0);
    }
    numConsSoFar += numBoxCons;

    // Then the PSD constraints
    for (size_t k=0; k<psd.size(); k++) {

        // If the matrix is just one
        if (psd[k].size() == 1) {
            continue;
        }

        // The indices and coefficients for the svec
        int fullMatSize = psd[k].size();
        int sVecSize = fullMatSize * (fullMatSize + 1) / 2;
        for (size_t i=0; i<psd[k].size(); i++) {
            for (size_t j=i; j<psd[k][i].size(); j++) {
                int l = 0;
                for (auto& term : psd[k][i][j]) {
                    
                    // Find this in the variable list
                    int realLoc = variableLocs[term.first];

                    // Locations in the svec
                    int realInd = matLocToVecLoc(i, j, fullMatSize) + l*sVecSize;

                    // Set the value in the matrix
                    ARows.push_back(realInd + numConsSoFar);
                    ACols.push_back(realLoc);
                    if (i == j) {
                        AVals.push_back(std::real(term.second));
                    } else {
                        AVals.push_back(std::real(term.second) * std::sqrt(2.0));
                    }

                    l++;

                }
            }

        }
        numConsSoFar += sVecSize;

    }

    // Verbose output of the A matrix
    if (verbosity >= 3) {
        std::cout << "ARows: " << ARows << std::endl;
        std::cout << "ACols: " << ACols << std::endl;
        std::cout << "AVals: " << AVals << std::endl;
        std::cout << "bVals: " << bVals << std::endl;
        std::vector<std::vector<double>> ATemp2(m, std::vector<double>(n, 0.0));
        for (size_t i=0; i<ACols.size(); i++) {
            ATemp2[ARows[i]][ACols[i]] = AVals[i];
        }
        std::cout << "Reconstructed matrix:" << std::endl;
        std::cout << ATemp2 << std::endl;
    }

    // Convert to compressed sparse column format
    toCC(ARows, ACols, AVals, n);

    // Verbose output of the A matrix
    if (verbosity >= 3) {
        std::cout << "Now in CSC format:" << std::endl;
        std::cout << "ARows: " << ARows << std::endl;
        std::cout << "ACols: " << ACols << std::endl;
        std::cout << "AVals: " << AVals << std::endl;
        std::vector<std::vector<double>> ATemp(m, std::vector<double>(n, 0.0));
        for (size_t i=0; i<ACols.size()-1; i++) {
            int col = i;
            int colStart = ACols[i];
            int colEnd = ACols[i+1];
            for (int j=colStart; j<colEnd; j++) {
                ATemp[ARows[j]][col] = AVals[j];
            }
        }
        std::cout << "Reconstructed matrix:" << std::endl;
        std::cout << ATemp << std::endl;
    }

    // Take the negative of A
    for (size_t i=0; i<AVals.size(); i++) {
        AVals[i] = -AVals[i];
    }
    for (size_t i=0; i<bVals.size(); i++) {
        bVals[i] = -bVals[i];
    }

    // Allocate SCS structs
    ScsCone *k = (ScsCone *)calloc(1, sizeof(ScsCone));
    ScsData *d = (ScsData *)calloc(1, sizeof(ScsData));
    ScsSettings *stgs = (ScsSettings *)calloc(1, sizeof(ScsSettings));
    ScsSolution *sol = (ScsSolution *)calloc(1, sizeof(ScsSolution));
    ScsInfo *info = (ScsInfo *)calloc(1, sizeof(ScsInfo));

    // Fill in data struct 
    ScsMatrix A = {&AVals[0], &ARows[0], &ACols[0], m, n};
    d->m = m;
    d->n = n;
    d->b = &bVals[0];
    d->c = &c[0];
    d->A = &A;
    d->P = SCS_NULL;

    // Box constraint cone
    double* boxConsMin = new double[n-1];
    double* boxConsMax = new double[n-1];
    for (int i=0; i<n-1; i++) {
        boxConsMin[i] = varBounds.first;
        boxConsMax[i] = varBounds.second;
    }
    k->bsize = n;
    k->bl = boxConsMin;
    k->bu = boxConsMax;

    // Zero cone
    k->z = 1 + constraintsZero.size();
    
    // Positivity cone
    k->l = constraintsPositive.size();

    // SDP cones
    k->ssize = psd.size();
    int *s = (int *)malloc(k->ssize * sizeof(int));
    for (int i = 0; i < k->ssize; ++i) {
        s[i] = psd[i].size();
    }
    k->s = s;

    // Utility to set default settings
    scs_set_default_settings(stgs);

    // Modify tolerances
    stgs->eps_abs = 1e-6;
    stgs->eps_rel = 1e-6;
    stgs->verbose = verbosity-1;

    // Initialize SCS workspace
    ScsWork *scs_work = scs_init(d, k, stgs);

    // Solve! 
    int exitflag = scs_solve(scs_work, sol, info, 0);
    double pObj = info->pobj;

    // Free SCS workspace 
    scs_finish(scs_work);

    // Verify that SCS solved the problem 
    if (verbosity >= 3) {
        printf("SCS solved successfully: %i\n", exitflag == SCS_SOLVED);
        printf("SCS took %i iters, using the %s linear solver.\n", info->iter, info->lin_sys_solver);
        printf("Optimal solution vector x*:\n");
        for (int i = 0; i < n; ++i) {
            printf("x[%i] = %4f\n", i, sol->x[i]);
        }
        printf("Optimal dual vector y*:\n");
        for (int i = 0; i < m; ++i) {
            printf("y[%i] = %4f\n", i, sol->y[i]);
        }
    }

    // Get the solution
    std::map<Mon, std::complex<double>> solMap;
    for (int i=0; i<n; i++) {
        solMap[variables[i]] = sol->x[i];
    }

    // Check the PSD matrices
    if (verbosity >= 3) {
        for (size_t k=0; k<psd.size(); k++) {
            Eigen::MatrixXcd mat = replaceVariables(psd[k], solMap);
            std::cout << mat << std::endl;
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(mat);
            Eigen::VectorXd eigenvalues = es.eigenvalues();
            std::cout << "Min eigenvalue: " << eigenvalues.minCoeff() << std::endl;
        }
    }

    // Free allocated memory 
    free(k);
    free(d);
    free(s);
    free(stgs);
    free(info);
    free(sol->x);
    free(sol->y);
    free(sol->s);
    free(sol);
    free(boxConsMin);
    free(boxConsMax);
    free(scs_work);

    return -pObj;

}

