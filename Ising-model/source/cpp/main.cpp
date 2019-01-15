#include <iostream>
#include <mpi.h>
#include "quantities.h"
#include "solver.h"
#include "execute_solve.h"
#include "dumpfiles.h"
#include <iomanip>
using namespace std;

ofstream file;

void dumpfiles_each_step(solver &solve_system, double, double);

int main(int argc, char* argv[])
{
    string filename, arg1, arg2;
    bool arg3,arg4, arg5;
    int NSpins, MCcycles;
    double InitialTemp, FinalTemp, TempStep;
    if (argc <= 5) {
        cout << "Bad Usage: " << argv[0] <<
                " read output file, Number of spins, MC cycles, initial and final temperature and tempurate step" << endl;
        exit(1);
    }

    if (argc > 1) {
        filename = argv[1];
        NSpins = atoi(argv[2]);
        MCcycles = atoi(argv[3]);
        InitialTemp = atof(argv[4]);
        FinalTemp = atof(argv[5]);
        TempStep = atof(argv[6]);
        arg1 = argv[7];
        arg2 = argv[8];
        arg3 = (bool) atoi(argv[9]);
        arg4 = (bool)atoi(argv[10]);
        arg5 = (bool)atoi(argv[11]);

    }
    quantities something(NSpins);
    solver solve_system(arg3, arg5);
    dumpfiles system_dumped_file(MCcycles, NSpins, arg3, filename);
    execute_solve execute(MCcycles, NSpins,InitialTemp, FinalTemp, TempStep, arg1, arg2, arg3,filename, arg4);
    execute.run_solver(solve_system, something, system_dumped_file,argc, argv);


    return 0;
}



