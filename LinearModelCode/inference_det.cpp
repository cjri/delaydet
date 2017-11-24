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

	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, seed);

    vector<double> fact_store;
    FindLogFact(fact_store,10000);
    
	//Initial parameters
	int L=50; /* length of a genome*/
	int N=100000; /* number of individuals in the population */
    double sigma=gsl_rng_uniform(rgen); /*selection coefficient */
    cout << "Initial sigma " << sigma << "\n";
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
    InitialFreqPopulation(haps,rgen);

    /*for (int i=0;i<haps.size();i++) {
        cout << haps[i].q << " ";
    }
    cout << "\n";*/
    vector<double> initial_freqs;
    for (int i=0;i<haps.size();i++) {
        initial_freqs.push_back(haps[i].q);
    }
    
    //Get Data
    vector< vector<int> > data;
    ImportData(dim,data);

    int N_s=0;
    for (int j=0;j<data[0].size();j++) {
        N_s=N_s+data[0][j];
    }
    
    //cout << data.size() << "\n";
    
    for (int i=0;i<initial_freqs.size();i++) {
        initial_freqs[i]=0;
    }
    initial_freqs[0]=1;
    
    //Iteration
    vector<double> initial_freqs_best=initial_freqs;
    double sigma_best=sigma;
    double lL_opt=-1e10;
    double lL_opt_best=-1e10;
    int first1=1;
    double ds=0.1;
    int max=dt*data.size();

    cout << "Max " << max << "\n";
    
    ofstream inf_file;

    //Implement simple 1-D search for sigma.  Keep high and low values - reject change outside of this.
    double sigma_low=0;
    double sigma_high=10;
    
    
    //First systematic search
    int reject=0;
    double lL_top=-1e10;
    double sigma_top=0;
    for (double sigma=0.01;sigma<=1.91;sigma=sigma+0.05) {
        //Second optimisation
        double lL_best=-1e10;
        double lL=-1e10;
        int first2=1;
        //Reset initial frequencies
        InitialFreqPopulation(haps,rgen);
        initial_freqs.clear();
        cout << "B ";
        for (int i=0;i<haps.size();i++) {
            initial_freqs.push_back(haps[i].q);
            cout << haps[i].q << " ";
        }
        cout << "\n";
        for (int i=0;i<haps.size();i++) {
            cout << haps[i].fitness << " ";
        }
        cout << "\n";
        int fail=0;
        for (int it2=0;it2<5000;it2++) {
            if (first2==0) {
                if (lL>lL_best) {
                    lL_best=lL;
                    initial_freqs_best=initial_freqs;
                    fail=0;
                } else {
                    initial_freqs=initial_freqs_best;
                    fail++;
                }
                
            }
            if (fail>=100) {break;}
            first2=0;
            //Change parameter
            ChangeFreq(dim,ds,initial_freqs,rgen);
            //Propagate the system
            DoPropagation(0,N_s,mu,dt,max,sigma,lL,initial_freqs,muts,fact_store,data,haps);
            lL_opt=lL_best;
        }
        cout << "Likelihood " << sigma << " " << lL_best << "\n";
        if (lL_best>lL_top) {
            lL_top=lL_best;
            sigma_top=sigma;
            reject=0;
        } else {
            reject++;
        }
        if (reject==9) {break;}
    }

    
    
    cout << "Top " << sigma_top << "\n";
    
    lL_top=-1e10;
    if (sigma_top<0.05) {sigma_top=0.05;}
    for (double sigma=sigma_top-0.05;sigma<=0.91;sigma=sigma+0.01) {
        //Second optimisation
        double lL_best=-1e10;
        double lL=-1e10;
        int first2=1;
        //Reset initial frequencies
        InitialFreqPopulation(haps,rgen);
        initial_freqs.clear();
        int fail=0;
        for (int it2=0;it2<5000;it2++) {
            if (first2==0) {
                if (lL>lL_best) {
                    lL_best=lL;
                    initial_freqs_best=initial_freqs;
                    fail=0;
                } else {
                    initial_freqs=initial_freqs_best;
                    fail++;
                }
                
            }
            if (fail>100) {break;}
            first2=0;
            //Change parameter
            ChangeFreq(dim,ds,initial_freqs,rgen);
            //Propagate the system
            DoPropagation(0,N_s,mu,dt,max,sigma,lL,initial_freqs,muts,fact_store,data,haps);
            lL_opt=lL_best;
        }
        cout << "A ";
        for (int i=0;i<initial_freqs_best.size();i++) {
            cout << initial_freqs_best[i] << " ";
        }
        cout << "\n";
        cout << "Likelihood " << sigma << " " << lL_best << "\n";
        if (lL_best>lL_top) {
            lL_top=lL_best;
            sigma_top=sigma;
            reject=0;
        } else {
            reject++;
        }
        if (reject>=3&&sigma>sigma_top+0.05) {break;}
    }
    
    cout << "Top " << sigma_top << "\n";
    
    sigma=sigma_top;
    
    for (int it1=0;it1<1000;it1++) {
        if (first1==0) {
            if (lL_opt>lL_opt_best) {
                lL_opt_best=lL_opt;
                cout << "Better lL_opt " << lL_opt_best << " sigma " << sigma << "\n";
                for (int i=0;i<initial_freqs_best.size();i++) {
                    cout << initial_freqs_best[i] << " ";
                }
                cout << "\n";
                sigma_best=sigma;
                DoPropagation(1,N_s,mu,dt,max,sigma,lL_opt,initial_freqs_best,muts,fact_store,data,haps);

            } else {
                if (sigma<sigma_best) {
                    sigma_low=sigma;
                }
                if (sigma>sigma_best) {
                    sigma_high=sigma;
                }
                sigma=sigma_best;
            }
        }
        
         cout << sigma_low << " " << sigma_best << " " << sigma_high << "\n";
        //Termination step
        if (sigma_high-sigma_low<0.0005) {
            cout << "Final lL_opt " << lL_opt_best << "\n";
            cout << "Initial_freqs ";
            for (int i=0;i<initial_freqs_best.size();i++) {
                cout << initial_freqs_best[i] << " ";
            }
            cout << "\n";
            cout << "Sigma " << sigma_best << "\n";
            break;
        }

        
        first1=0;
        if (it1>0) {
        //Change parameter
        int checked=0;
        do {
            sigma=sigma+(gsl_rng_uniform(rgen)*2*ds)-ds;
            if (sigma>sigma_low&&sigma<sigma_high) {
                checked=1;
            }
        } while (checked==0);
        
        
        if (sigma<0) {
            sigma=-sigma;
        }
        }
        //Second optimisation
        double lL_best=-1e10;
        double lL=-1e10;
        int first2=1;
        //Reset initial frequencies
        InitialFreqPopulation(haps,rgen);
        initial_freqs.clear();
        for (int i=0;i<haps.size();i++) {
            initial_freqs.push_back(haps[i].q);
        }
        for (int it2=0;it2<5000;it2++) {
            if (first2==0) {
                if (lL>lL_best) {
                    lL_best=lL;
                    initial_freqs_best=initial_freqs;
                } else {
                    initial_freqs=initial_freqs_best;
                }

            }
            first2=0;
            //Change parameter
            ChangeFreq(dim,ds,initial_freqs,rgen);
            //Propagate the system
            DoPropagation(0,N_s,mu,dt,max,sigma,lL,initial_freqs,muts,fact_store,data,haps);
            lL_opt=lL_best;
        }
        
    }

    return 0;
    
    
}
	
	