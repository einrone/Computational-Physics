#ifndef SOLVER_H
#define SOLVER_H
#include "quantities.h"
#include "dumpfiles.h"
#include <random>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
using namespace std;

class solver
{
public:
    solver(bool ses, bool arg5);
    ~solver();
    bool save_each_step;
    friend class quantities;

    double *DumpVal;
    int sum_spin, m, n, dE;
    int flipped = 0;
    int Mc_equilibrium;
    bool equilibrium;

    int sum_nabo_spin(quantities &model, int i,int j, int n_spins);
    void metropolis(quantities &model, dumpfiles &sendfiles,string argument,int mc_cycle,
                    int n_spins, string arg, double temp, int rank, int Np);

    //expectation value of quantity X

};

#endif // SOLVER_H
