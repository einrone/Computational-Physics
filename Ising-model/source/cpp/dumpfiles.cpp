#include "dumpfiles.h"
#include "mpi.h"

dumpfiles::dumpfiles(int mc_cycles,int Nspins, bool ses, string filename):
    MCcycles(mc_cycles)
{   int criterion = 1e5;

    if(MCcycles >= criterion){

        length = MCcycles/100;
    }
    if(MCcycles < criterion){
        length = MCcycles;

    }

    accepted_energy = new double[length];
    expected_value1 = new double[length];
    expected_value2 = new double[length];
    accepted_flipps = new double[length];

    expected_value_sum1 = new double[length];
    expected_value_sum2 = new double[length];
    accepted_flipps_sum = new double[length];
    accepted_energy_sum = new double[length];


    for(int i = 0; i < length; i++){
        accepted_energy[i] = 0;
        expected_value1[i] = 0;
        expected_value2[i] = 0;

        expected_value_sum1[i] = 0;
        expected_value_sum2[i] = 0;

        accepted_flipps[i] = 0;
        accepted_flipps_sum[i] = 0;
        accepted_energy_sum[i] =0;
    }

    total_spin = int(Nspins*Nspins);


    if(ses == false){
        mfile.open(filename + ".txt");
        /*
        mfile << setw(15) <<    "temp"      <<
                 setw(15) <<    "Nspin"     <<
                 setw(15) <<    "MCcycle"   <<
                 setw(15) <<    "mean E"    <<
                 setw(15) <<    "mean E2"   <<
                 setw(15) <<    "Cv"        <<
                 setw(15) <<    "mean M"    <<
                 setw(15) <<    "mean M2"   <<
                 setw(15) <<    "X"         <<
                 setw(15) <<    "mean abs(M)"<< endl;
        */
    }

    /*
     A class which saves values into a .txt file.
     you have the option to either save all calculated result
      or you can save the last calculated value
     a clear description of this is given in the article ising_model.pdf
     */
}

dumpfiles::~dumpfiles(){
    //DELETES ALL POINTERS USED IN THIS CLASS
    mfile.close();
    delete[] expected_value_sum1;
    delete[] expected_value_sum2;
    delete[] expected_value1;
    delete[] expected_value2;
    delete[] accepted_flipps;
    delete[] accepted_flipps_sum;
    delete[] accepted_energy;
    delete[] accepted_energy_sum;
}


void dumpfiles::save_each_step(int index, double val1, double val2, int flipps, int step, int E){
    //SAVES THE VALUES FROM METROPOLIS INTO NEW POINTS. THIS VALUES ARE SENT TO SAVE_FILE_EACH_STEP
    //THE VALUES ARE ALSO NORMALIZED WITH RESPECT TO MONTE CARLO CYCLES AND TOTAL SPIN
    expected_value1[index] = double(val1/double((step)*total_spin));
    expected_value2[index] = double(val2/double((step)*total_spin));
    accepted_flipps[index] = double(flipps/double(total_spin*step));
    accepted_energy[index] = double(double(E)/total_spin);


}

void dumpfiles::save_file_each_step(string argument, string filename, double temp, int rank, int Nprocessors){
    //REDUCES ALL VALUES AND SAVES THEM INTO NEW POINTERS
    MPI_Allreduce(expected_value1, expected_value_sum1,length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(expected_value2, expected_value_sum2,length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(accepted_flipps, accepted_flipps_sum,length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(accepted_energy, accepted_energy_sum,length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if(rank == 0){
        double Np = double(Nprocessors);
        ofstream file;
        file.open(filename + "_" + argument + ".txt");
        //SAVES THE VALUES INTO FILES. THIS SAVES ENERGY, MAGNETIC MOMENT, ACCEPTED FLIPPS AND ACCEPTED ENERGY
        //NOTE THAT IT IS ALSO NORMALIZED WITH RESPECT TO NUMBER OF CORES USED
        for(int k = 0; k < length; k++){
            file << setprecision(8) << setw(15) << expected_value_sum1[k]/Np << setw(15)
                 << expected_value_sum2[k]/Np << setw(15) << k << setw(15) << accepted_flipps_sum[k]/Np
                 <<setw(15) <<accepted_energy_sum[k]/Np<< endl;
        }

        file.close();
    }
}

void dumpfiles::save_last_step(double mypointer[],double temp, int mc, int spin, int rank,
                               int Nprocessors, string filename, bool large_system){
    if(rank == 0){
        int Mc_equilibrium = int(10000/Nprocessors);
        //SAVES THE LAST VALUES IN MONTE CARLO SIMULATION
        double X;
        double e1;
        double e2;
        double e3;
        double e4;
        double e5;
        if(large_system == true){
            double Np = double(Nprocessors);
            //IF YOU WANT TO START AFTER STEADY STATE, THE NUMBER OF CYCLES TO
            //STEADY STATE ARE REMOVED FRPOM THE TOTAL
            //NORMALIZED WITH RESPECT TO MC CYCLE AND NUMBER PROCESSORS

            e1 = mypointer[0]/double((MCcycles-Mc_equilibrium)*Np);
            e2 = mypointer[1]/double((MCcycles-Mc_equilibrium)*Np);
            e3 = mypointer[2]/double((MCcycles-Mc_equilibrium)*Np);
            e4 = mypointer[3]/double((MCcycles-Mc_equilibrium)*Np);
            e5 = mypointer[4]/double((MCcycles-Mc_equilibrium)*Np);
            X = variance(e4, e5,temp*spin*spin);

        }
        if(large_system == false){
            //IF YOU WANT TO INCLUDE ALL CYCLES
            //NORMALIZED WITH RESPECT TO MC CYCLE AND NUMBER PROCESSORS
            double Np = double(Nprocessors);
            e1 = mypointer[0]/(MCcycles*Np);
            e2 = mypointer[1]/(MCcycles*Np);
            e3 = mypointer[2]/(MCcycles*Np);
            e4 = mypointer[3]/(MCcycles*Np);
            e5 = mypointer[4]/(MCcycles*Np);
            X = variance(e4, e3,temp*spin*spin);

        }

        double Cv = variance(e2, e1,temp*temp*spin*spin);
        //NORMALIZES WITH RESPECT TO TOTAL SPIN AND THEN SAVES TO .TXT FILE
        mfile <<    setw(15) << right << setprecision(8)   <<   temp                    <<
                    setw(15) << right <<setprecision(8)    <<   spin                    <<
                    setw(15) << right <<setprecision(8)    <<   mc                      <<
                    setw(15) << right <<setprecision(8)    <<   e1/(double(spin*spin))  <<
                    setw(15) << right <<setprecision(8)    <<   e2/(double(spin*spin))  <<
                    setw(15) << right <<setprecision(8)    <<   Cv                      <<
                    setw(15) << right <<setprecision(8)    <<   e3/(double(spin*spin))  <<
                    setw(15) << right <<setprecision(8)    <<   e4/(double(spin*spin))  <<
                    setw(15) <<right << setprecision(8)    <<   X                       <<
                    setw(15) << right <<setprecision(8)    <<   e5/(double(spin*spin))  << endl;

    }

}
