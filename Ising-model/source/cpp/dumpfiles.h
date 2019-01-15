#ifndef DUMPFILE_H
#define DUMPFILE_H
#include <iomanip>
#include <iostream>
#include <fstream>


using namespace std;
class dumpfiles
{
public:
    friend class solver;
    dumpfiles(int mc_cycles, int Nspins, bool ses, string filename);
    ~dumpfiles();

    int rank;
    int MCcycles;

    inline double variance(double X2, double X, double temp_nspin){
        return (double)((X2-X*X)/(temp_nspin));
    }


    double *expected_value1;
    double *expected_value2;
    double *expected_value_sum1;
    double *expected_value_sum2;
    double *accepted_flipps;
    double *accepted_flipps_sum;
    double *accepted_energy;
    double *accepted_energy_sum;

    ofstream mfile;

    void setRank(int new_rank) { rank = new_rank; }


    void save_each_step(int index, double val1, double val2,int flipps, int step, int E);


    int total_spin;
    int length;

    void save_file_each_step(string argument, string filename, double temp, int rank, int Nprocessors);
    void save_last_step(double mypointer[], double temp, int mc, int spin, int rank, int Nprocessors,
                        string filename, bool large_system);
};

#endif // DUMPFILE_H
