#include "solver.h"

solver::solver(bool ses, bool arg5):
    save_each_step(ses)
{
    DumpVal = new double [5]; //pointer array of size 5
    for(int i = 0; i < 5; i++){
        DumpVal[i] = 0;
    }
    equilibrium = arg5;


}

solver::~solver(){
    //DELETES ALL POINTERS USED IN THIS CLASS
    delete[] DumpVal;
}

int solver::sum_nabo_spin(quantities &model,int i,int j, int n_spins){
    //THIS FUNCTION CALCULATES THE SUM OF THE NEIGHBOURING SPIN
    //quantities &model;
    //1.mxn left
    //2.mxn right
    //3.mxn left
    //4.mxn right
    int sum = 0;

    sum += (model.ptr_lattice[model.Idx(model.PeriodicBoundary(i,n_spins,-1),j)]
            + model.ptr_lattice[model.Idx(i,model.PeriodicBoundary(j,n_spins,-1))]

            + model.ptr_lattice[model.Idx(model.PeriodicBoundary(i,n_spins,1),j)]
            + model.ptr_lattice[model.Idx(i,model.PeriodicBoundary(j,n_spins,1))]);
    //cout << "sum "<<sum << endl;
    return sum;
}

void solver::metropolis(quantities &model, dumpfiles &sendfiles,string argument ,int mc_cycle,
                        int n_spins, string arg, double temp, int rank, int Np){

    //THIS FUNCTION IS THE METROPOLIS FUNCTION
    int criterion = 1e5;
    int length;
    if(mc_cycle >= criterion){
        length = mc_cycle/100;
        //SAMPLING PER HUNDRED STEP IF YOU CHOOSE CYCLES BIGGER THAN 1E5
    }
    if(mc_cycle < criterion){
        length = mc_cycle;
    }


    Mc_equilibrium = int(10000/Np);
    random_device rd;
    mt19937_64 gen(rd() + rank);
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
    uniform_int_distribution<int> RandomPosition(0, n_spins-1);

    double E = 0;
    double M = 0;
    int Ns = n_spins*n_spins;
    model.initialize_matrix(arg);
    E = model.initialiaze_energy(E);
    M = model.initialize_magnetization(M);
    model.initialize_energy_difference(temp);
    int flipped = 0;

    for(int i = 0; i < 5; i++){
        DumpVal[i] = 0;
    }

    for(int k = 1; k < mc_cycle; k++){
        for(int j = 0; j < Ns;  j++){
            m = (int)(RandomPosition(gen));
            n = (int)(RandomPosition(gen));
            dE = 2*model.ptr_lattice[model.Idx(m,n)]*sum_nabo_spin(model, m,n, n_spins);

            double R = RandomNumberGenerator(gen);
            //CHECK IF THE RANDOM NUMBER IS LOWER THAN THE ACCEPTANCE RATIO
            if(R <= model.ptr_Ediff[dE + 8]){
                model.ptr_lattice[model.Idx(m,n)] *= -1;
                E += dE;
                M += 2*model.ptr_lattice[model.Idx(m,n)];
                flipped += 1;

            }
        }

        if(equilibrium == true && k > Mc_equilibrium){
            //IF YOU WANT TO START SAVING AFTER EQUILIBRIUM THIS VARIABLE IS TRUE
            DumpVal[0] += E;
            DumpVal[1] += E*E;
            DumpVal[2] += M;
            DumpVal[3] += M*M;
            DumpVal[4] += fabs((double)M);
        }

        if(equilibrium == false){
            //IF YOU DONT WANT TO START AFTER EQUILIBRIUM, ALL VALUES ARE SUMMED AND SAVED
            DumpVal[0] += E;
            DumpVal[1] += E*E;
            DumpVal[2] += M;
            DumpVal[3] += M*M;
            DumpVal[4] += fabs((double)M);
        }


        //THESE IF TEST SAMPLES AND SAVES IF YOU WANT TO SAVE FOR EACH STEP
        if(save_each_step == true){
            if(length == mc_cycle && (temp == double(2.4) || temp == double(1.0))){
                sendfiles.save_each_step(int(k), DumpVal[0],DumpVal[4], flipped, int(k), E);
            }

            if(k%100 == 0 && length < mc_cycle && (temp == double(2.4) || temp == double(1.0))){
                sendfiles.save_each_step(int(k/100), DumpVal[0],DumpVal[1], flipped, int(k), E);
            }

        }
    }

}
