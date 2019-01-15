#ifndef QUANTITIES_H
#define QUANTITIES_H
#include <iostream>
#include <random>
#include <math.h>
#include <mpi.h>



using namespace std;

class quantities
{
public:
    quantities(int L);
    ~quantities();



    //pointers
    int *ptr_lattice;
    double *ptr_Ediff;


    //variables
    int n_spins;

    void initialize_matrix(string arg);
    void ising(double temp, double E, double M);


    //inline methods
    inline int Idx(int i, int j){
        return i*n_spins + j;
    };

    inline int PeriodicBoundary(int k, int limit, int add) {
      return (k+limit+add) % (limit);
    }

    double initialiaze_energy(double E);
    void initialize_energy_difference(double temp);
    double initialize_magnetization(double M);
};

#endif // QUANTITIES_H
