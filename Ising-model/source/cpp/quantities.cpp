#include "quantities.h"

quantities::quantities(int L):
    n_spins(L)
{

    ptr_lattice = new int [(int)((n_spins)*(n_spins))];
    ptr_Ediff = new double [17];
}


quantities::~quantities(){
    delete[] ptr_lattice;
    delete[] ptr_Ediff;
}

void quantities::initialize_matrix(string arg){

    if(arg == "random"){
        random_device rd;
        mt19937_64 gen(rd());
        uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0);
        //this gives the lattice a random configurations spins
        for(int i = 0; i < n_spins; i++){
            for(int j = 0; j < n_spins; j++)
                if(RandomNumberGenerator(gen) >= 0.5){
                    ptr_lattice[Idx(i,j)] = 1;
                }
                else{
                    ptr_lattice[Idx(i,j)] = -1;

                }
        }

    }

    if(arg == "up"){
        //this gives the lattice a ordered configuration, where all spins points upward
        for(int i = 0; i < n_spins; i++){
            for(int j = 0; j < n_spins; j++){
                ptr_lattice[Idx(i,j)] = 1;

            }
        }
    }

    if(arg == "down"){
        //this gives the lattice a ordered configuration, where all spins points downwards
        for(int i = 0; i < n_spins; i++){
            for(int j = 0; j < n_spins; j++){
                ptr_lattice[Idx(i,j)] = -1;
            }
        }
    }
}


double quantities::initialiaze_energy(double E){
    //initializes the energy in the system when looking at a spin with four spins around it

    for(int ix = 0; ix < n_spins; ix++){
        for(int jx = 0; jx < n_spins; jx++){
            E -= (double)ptr_lattice[Idx(ix,jx)]*(ptr_lattice[Idx(PeriodicBoundary(ix,n_spins,-1),jx)]
                    + ptr_lattice[Idx(ix,PeriodicBoundary(jx,n_spins,-1))]);
        }
    }

    return E;
}

double quantities::initialize_magnetization(double M){
    //initializes the magnetic moment when looking at a spin with four spins around it
    for(int i = 0; i < n_spins; i++){
        for(int j = 0; j < n_spins; j++){
            M += (double)ptr_lattice[Idx(i,j)];
        }
    }

    return M;
}
void quantities::initialize_energy_difference(double temp){
    //precomputes the acceptance ratio
    for(int dE = -8; dE <= 8; dE+= 4){
        //cout << temp << endl;
        ptr_Ediff[dE + 8] = exp(-(double)dE/temp);

    }
}

