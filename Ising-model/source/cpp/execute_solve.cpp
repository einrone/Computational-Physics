#include "execute_solve.h"

execute_solve::execute_solve(int mc_cycle, int NSpins, double Ti, double Tf, double Ts,
                             string arg1,string arg2, bool arg3, string fname, bool arg4)
{
    MCcycle = mc_cycle;
    Ns = NSpins;
    InitialTemp = Ti;
    FinalTemp = Tf;
    TempStep = Ts;
    argument1 = arg1;
    argument2 = arg2;
    argument3 = arg3;
    argument4 = arg4;
    filename = fname;

    final_dumped_value = new double [5];

}

execute_solve::~execute_solve()
{
    delete[] final_dumped_value;
}
void execute_solve::run_solver(solver &solutions_solver, quantities &info, dumpfiles &files_dumped,int argc, char* argv[])
{
    bool save_each_step = argument3;

    int NProcesses, RankProcess;

    //initializes the mpi
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &NProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &RankProcess);

    files_dumped.setRank(RankProcess);

    MPI_Bcast(&MCcycle, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Ns, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&InitialTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&FinalTemp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&TempStep, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //RankProcess= 0;
    double TimeStart = MPI_Wtime();
    if(RankProcess == 0){
        cout <<"Simulation started for "<< MCcycle << " monte carlo cycles per core. " << NProcesses << " cores are used" << endl;
    }
    //this loops for the initial and final temperature given.
    for(double temperature = InitialTemp; temperature <= FinalTemp; temperature += TempStep){
        solutions_solver.metropolis(info, files_dumped, argument2, MCcycle, Ns,argument1,temperature, RankProcess,
                                    NProcesses);

        if(save_each_step == true){
            //saves all values when running metroplis function
            files_dumped.save_file_each_step(argument2, filename,temperature,RankProcess, NProcesses);
        }

        if(save_each_step == false){
            //saves the last value when running metroplis function
            for(int i = 0; i < 5; i++)final_dumped_value[i] = 0;
            MPI_Reduce(&solutions_solver.DumpVal[0], &final_dumped_value[0],5, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            if(RankProcess == 0) cout << "Dumping thermodynamic quantities for T = " << temperature << " to file"<<endl;
            files_dumped.save_last_step(final_dumped_value, temperature, MCcycle, Ns, RankProcess,NProcesses,
                                        filename, argument4);
        }
    }

    double TimeEnd = MPI_Wtime();


    if(RankProcess == 0 && save_each_step == true ){
        cout << "Simulation done " << filename + "_" + argument2 + ".txt" << " is generated"<< endl;
        cout << "Time spent " << TimeEnd - TimeStart << endl;
        cout << "Program exited" << endl;
    }

    if(RankProcess == 0 && save_each_step == false ){
        cout << "Simulation done " << filename + ".txt" << " is generated with values in last step"<< endl;
        cout << "Time spent " << TimeEnd - TimeStart << endl;
        cout << "Program exited" << endl;
    }

    MPI_Finalize();
}
