//Program to calculate predicted sigma values for loci in a model system, where data is generated by a sampling process
//This code uses a deterministic method to run the underlying population distribution, and samples it at low frequency

#include <iostream>
#include <vector>
#include <list>
#include <deque>
#include <sstream>
using namespace std;

#include "shared.h"

/*#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>*/

vector<trajectory> traj; 

vector<haplo> population;
vector<haplo> pop_temp;
vector<haplo> pop_sample;

int main(int argc, char *argv[]){

	//Set up random number generator
	int seed=atoi(argv[1]);
    int dim=atoi(argv[2]);
    double sigma=atof(argv[3]); /*selection coefficient */
    int N_s=1000;

	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);

	//Initial parameters
	int L=50; /* length of a genome*/
	int N=100000; /* number of individuals in the population */
	int tSim=100000; /* simulation time*/	
	double mu=pow(10,-5);  /* mutation rate */
	//double sigma=0.01;
	int min_sample=100; /*output parameter*/
	int nDrivers=L; /*All loci are drivers - can modify though*/
    int dt=1;
    
    //Get haplotypes
    vector<haplo> haps;
    ImportHapsLinear(dim,haps);

    //Construct haplotyps fitnesses
    GetHapFits(sigma,haps);

    //Construct list of mutations
    vector<mut> muts;
    GetMutList(haps,muts);
    
    //Initialise population
    SetupPopulation(N,haps);

    //Propagation
    ofstream freqs_file;
    freqs_file.open("Freqs.out");
    vector<double> p_select;
    vector<haplo> haps_temp;
    for (int t=0;t<1000;t++) {
        //Mutation step
        CalcMutation(mu,haps,muts,rgen);
    
        //Calculate fitness
        double glob_fit=GetGlobalFit(haps);

        //Sample next generation
        CalcProbs(glob_fit,haps,p_select);
        GetNextGen(N,haps,p_select,rgen);
        
        //Generate output - standard multinomial
        if (t%dt==0) {
            SampleNextGen (N,N_s,haps,p_select,freqs_file,rgen);
        }
        
        //Terminate at 99%
        if ((haps[haps.size()-1].n+0.)/(N+0.)>0.99) {
            break;
        }
    }
    
    return 0;
    
}
	
	
