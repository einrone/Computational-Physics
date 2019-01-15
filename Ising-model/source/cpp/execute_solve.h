#ifndef EXECUTE_SOLVE_H
#define EXECUTE_SOLVE_H
#include <iomanip>
#include <iostream>
#include "mpi.h"
#include "quantities.h"
#include "solver.h"
#include "dumpfiles.h"


class execute_solve
{
public:
    friend class solver;
    execute_solve(int mc_cycle, int NSpins, double Ti, double Tf, double Ts, string arg1,
                  string arg2, bool arg3, string fname, bool arg4);
    ~execute_solve();

    double *final_dumped_value;

    int MCcycle;
    int Ns;
    double InitialTemp;
    double FinalTemp;
    double TempStep;
    string  argument1;
    string  argument2;
    bool  argument3;
    bool  argument4;
    string filename;


    void run_solver(solver &solutions_solver, quantities &info, dumpfiles &files_dumped,int argc, char* argv[]);
};
#endif // EXECUTE_SOLVE_H
