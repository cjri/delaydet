#include "shared.h"
#include <iostream>
#include <sstream>
#include <string>

void ImportHapsLinear (int dim, vector<haplo>& haps) {
    ifstream haps_file;
    stringstream s;
    s << dim;
    string d=s.str();
    string file="../../Haps"+d+".in";
    //cout << d << "\n";
    //cout << file << "\n";
    haps_file.open(file.c_str());
    int x;
    for (int i=0;i<=dim;i++) {
        haplo h;
        for (int j=0;j<dim;j++) {
            if (!(haps_file >> x)) break;
           // cout << x << " ";
            h.seq.push_back(x);
        }
        //cout << "\n";
        haps.push_back(h);
        
    }

}

void GetMutList (vector<haplo>& haps, vector<mut>& muts) {
    for (int i=0;i<haps.size();i++) {
        for (int j=0;j<haps.size();j++) {
            int diff=0;
            for (int k=0;k<haps[0].seq.size();k++) {
                if (haps[i].seq[k]!=haps[j].seq[k]) {
                    diff++;
                }
            }
            if (diff==1) {
                mut m;
                m.a=i;
                m.b=j;
                muts.push_back(m);
            }
        }
    }
}

void GetHapFits (double sigma, vector<haplo>& haps) {
    for (int i=0;i<haps.size();i++) {
        haps[i].fitness=1;
        for (int j=0;j<haps[i].seq.size();j++) {
            if (haps[i].seq[j]==1) {
                haps[i].fitness=haps[i].fitness+sigma;
            }
        }
    }
}

void SetupPopulation (int N, vector<haplo>& haps) {
    haps[0].n=N;
    for (int i=1;i<haps.size();i++) {
        haps[i].n=0;
    }
}

void InitialFreqPopulation (vector<haplo>& haps, gsl_rng *rgen) {
    double tot=0;
    for (int i=0;i<haps.size();i++) {
        double r=gsl_rng_uniform(rgen);
        haps[i].q=r;
        tot=tot+r;
    }
    for (int i=0;i<haps.size();i++) {
        haps[i].q=haps[i].q/tot;
    }
}

void InitialFreqPopulationRapid (vector<haplo>& haps) {
    double tot=0;
    haps[0].q=1;
    for (int i=1;i<haps.size();i++) {
        haps[i].q=0;
    }
}

void InitialFreqPopulationN (int N, vector<haplo>& haps, gsl_rng *rgen) {
    double tot=0;
    for (int i=0;i<haps.size();i++) {
        double r=gsl_rng_uniform(rgen);
        haps[i].q=r;
        tot=tot+r;
    }
    for (int i=0;i<haps.size();i++) {
        haps[i].q=haps[i].q/tot;
    }
    int tot_i=0;
    for (int i=0;i<haps.size()-1;i++) {
        haps[i].n=haps[i].q*N;
        tot_i=tot_i+haps[i].n;
    }
    haps[haps.size()-1].n=N-tot_i;
}


void ImportData (int dim, vector< vector<int> >& data) {
    ifstream freqs_file;
    freqs_file.open("Freqs.out");
    double x;
    int index=0;
    vector<int> d;
    do {
        if (!(freqs_file >> x)) {break;}
        d.push_back(x);
        index++;
        if (index==dim+1) {
            data.push_back(d);
            d.clear();
            index=0;
        }
    } while (1==1);
}

void CalcMutation (double mu, vector<haplo>& haps, vector<mut>& muts, gsl_rng *rgen) {
    //Numbers
    for (int i=0;i<muts.size();i++) {
        muts[i].n=gsl_ran_poisson(rgen,haps[muts[i].a].n*mu);/* Number of mutations along this pathway */
    }
    
    //Implementation
    for (int i=0;i<muts.size();i++) {
        haps[muts[i].b].n=haps[muts[i].b].n+muts[i].n;
        haps[muts[i].a].n=haps[muts[i].a].n-muts[i].n;
    }
}

void DetMutation (double mu, vector<mut>& muts, vector<haplo>& haps) {
    vector<double> mval;
    for (int i=0;i<muts.size();i++) {
        double m=mu*haps[muts[i].a].q;
        mval.push_back(m);
    }
    for (int i=0;i<muts.size();i++) {
        haps[muts[i].b].q=haps[muts[i].b].q+mval[i];
        haps[muts[i].a].q=haps[muts[i].a].q-mval[i];
    }
}

void DelayDetMutation (double mu, double beta, vector<mut>& muts, vector<haplo>& haps) {
    vector<double> mval;
    for (int i=0;i<muts.size();i++) {
        double m=mu*haps[muts[i].a].q;
        mval.push_back(m);
    }
    for (int i=0;i<muts.size();i++) {
       // cout << haps[muts[i].a].q << " " << beta <<" " << mval[i] << "\n";
        if (haps[muts[i].a].q>beta) {
            haps[muts[i].b].q=haps[muts[i].b].q+mval[i];
            haps[muts[i].a].q=haps[muts[i].a].q-mval[i];
        }
    }
}


double GetGlobalFit (vector<haplo> haps) {
    double glob_fit=0;
    for (int i=0;i<haps.size();i++) {
        glob_fit=glob_fit+(haps[i].n*haps[i].fitness);
    }
    return glob_fit;
}

double GetGlobalFitQ (vector<haplo> haps) {
    double glob_fit=0;
    for (int i=0;i<haps.size();i++) {
        glob_fit=glob_fit+(haps[i].q*haps[i].fitness);
    }
    return glob_fit;
}


void CalcProbs (double glob_fit, vector<haplo> haps, vector<double>& p_select) {
    p_select.clear();
    for (int i=0;i<haps.size();i++) {
        double p=(haps[i].n*haps[i].fitness)/glob_fit;
        p_select.push_back(p);
    }
}

void SampleNextGen (int N, int c, vector<haplo> haps, vector<double>& p_select, ofstream& freqs_file, gsl_rng *rgen) {
    p_select.clear();
    for (int i=0;i<haps.size();i++) {
        double prob=(haps[i].n+0.)/(N+0.);
      //  cout << prob << " ";
        p_select.push_back(prob);
    }
    //cout << "\n";
    double prob_select_mirror[haps.size()];
    unsigned int next_gen[haps.size()];
    for (unsigned int i=0;i<haps.size();i++) {
        prob_select_mirror[i]=p_select[i];
    }
    gsl_ran_multinomial(rgen,haps.size(),c,prob_select_mirror,next_gen);
    for (int i=0;i<haps.size();i++) {
        freqs_file << next_gen[i] << " ";
    }
    freqs_file << "\n";
}

void DetNextGen(double glob_fit, vector<haplo>& haps) {
    for (int i=0;i<haps.size();i++) {
        haps[i].q=(haps[i].q*haps[i].fitness)/glob_fit;
    }
}

void GetNextGen (int N, vector<haplo>& haps, vector<double> p_select, gsl_rng *rgen) {
    double prob_select_mirror[haps.size()];
    unsigned int next_gen[haps.size()];
    for (unsigned int i=0;i<haps.size();i++) {
        prob_select_mirror[i]=p_select[i];
    }
    gsl_ran_multinomial(rgen,haps.size(),N,prob_select_mirror,next_gen);
    for (int i=0;i<haps.size();i++) {
        haps[i].n=next_gen[i];
        //cout << next_gen[i] << " ";
    }
    //cout << "\n";
}

void ChangeFreq (int dim, double ds, vector<double>& initial_freqs, gsl_rng *rgen) {
    int r=floor(gsl_rng_uniform(rgen)*dim)+1;
    int s=0;
//    do {
//        s=floor(gsl_rng_uniform(rgen)*dim)+1;
//    } while (s==r);
    double diff=(gsl_rng_uniform(rgen)*2*ds)-ds;
    initial_freqs[r]=initial_freqs[r]+diff;
//    initial_freqs[s]=initial_freqs[s]+diff;
    if (initial_freqs[r]<0) {
        initial_freqs[r]=0;
    }
//    if (initial_freqs[s]<0) {
//        initial_freqs[s]=1e-20;
//    }
    double tot=0;
    for (int i=0;i<initial_freqs.size();i++) {
        tot=tot+initial_freqs[i];
    }
    for (int i=0;i<initial_freqs.size();i++) {
        initial_freqs[i]=initial_freqs[i]/tot;
    }
}

void DoPropagation (int verb, int N_s, double mu, int dt, int max, double sigma, double& lL, vector<double>& initial_freqs, vector<mut>& muts, vector<double>& fact_store, vector< vector<int> >& data, vector<haplo>& haps) {
//    for (int i=0;i<initial_freqs.size();i++) {
  //      cout << initial_freqs[i] << " ";
   // }
   // cout << "\n";
   // cout << sigma << "\n";
    for (int i=0;i<haps.size();i++) {
        haps[i].q=initial_freqs[i];
    }
    GetHapFits(sigma,haps);
    lL=0;
    ofstream inf_file;
    if (verb==1) {
        inf_file.open("Inf.out");
    }
    for (int t=0;t<max;t++) {
        DetMutation(mu,muts,haps);
        double glob_fit=GetGlobalFitQ(haps);
        DetNextGen(glob_fit,haps);
        if (t%dt==0) {
            vector<double> inf;
            for (int i=0;i<haps.size();i++) {
                if (haps[i].q>1e-20) {
                    inf.push_back(haps[i].q);
                } else {
                    inf.push_back(1e-20);
                }
                if (verb==1) {
                    inf_file << haps[i].q << " ";
                }
            }
            if (verb==1) {
                inf_file << "\n";
            }
            int tdt=t/dt;
            double logL=MultiCalc(N_s,data[tdt],inf,fact_store);
            lL=lL+logL;
        }
    }
    //cout << "lL " << lL << "\n";
}

void DoPropagationStoch (int verb, int N_s, int N, double mu, int dt, int max, double sigma, double& lL, vector<double>& initial_freqs, vector<mut>& muts, vector<double>& fact_store, vector< vector<int> >& data, vector<haplo>& haps, gsl_rng *rgen) {
    for (int i=0;i<haps.size();i++) {
        haps[i].q=initial_freqs[i];
    }
    GetHapFits(sigma,haps);
    lL=0;
    ofstream inf_file;
    if (verb==1) {
        inf_file.open("Inf_stoch.out");
    }
    vector<double> p_select;
    for (int t=0;t<max;t++) {
        CalcMutation(mu,haps,muts,rgen);
        double glob_fit=GetGlobalFit(haps);
        //Sample next generation
        CalcProbs(glob_fit,haps,p_select);
        GetNextGen(N,haps,p_select,rgen);
        
        if (t%dt==0) {
            //Get haplotype frequencies
            for (int i=0;i<haps.size();i++) {
                haps[i].q=(haps[i].n+0.)/(N+0.);
            }
            vector<double> inf;
            for (int i=0;i<haps.size();i++) {
                if (haps[i].q>1e-20) {
                    inf.push_back(haps[i].q);
                } else {
                    inf.push_back(1e-20);
                }
                if (verb==1) {
                    inf_file << haps[i].q << " ";
                }
            }
            if (verb==1) {
                inf_file << "\n";
            }
            int tdt=t/dt;
            double logL=MultiCalc(N_s,data[tdt],inf,fact_store);
            //cout << logL << "\n";
            lL=lL+logL;
        }
    }
}


void DoPropagationR (int verb, int N_s, double mu, int dt, int max, double sigma, double& lL, vector<double>& initial_freqs, vector<mut>& muts, vector<double>& fact_store, vector< vector<int> >& data, vector<haplo>& haps) {
    for (int i=0;i<haps.size();i++) {
        haps[i].q=initial_freqs[i];
    }
    GetHapFits(sigma,haps);
    lL=0;
    ofstream inf_file;
    if (verb==1) {
        inf_file.open("Inf_rapid.out");
    }
    for (int t=0;t<max;t++) {
        DetMutation(mu,muts,haps);
        double glob_fit=GetGlobalFitQ(haps);
        DetNextGen(glob_fit,haps);
        if (t%dt==0) {
            vector<double> inf;
            for (int i=0;i<haps.size();i++) {
                if (haps[i].q>1e-20) {
                    inf.push_back(haps[i].q);
                } else {
                    inf.push_back(1e-20);
                }
                if (verb==1) {
                    inf_file << haps[i].q << " ";
                }
            }
            if (verb==1) {
                inf_file << "\n";
            }
            int tdt=t/dt;
            double logL=MultiCalc(N_s,data[tdt],inf,fact_store);
            lL=lL+logL;
        }
    }
    //cout << "lL " << lL << "\n";
}

void DoPropagationDD (int verb, int N_s, double mu, double beta, int dt, int max, double sigma, double& lL, vector<double>& initial_freqs, vector<mut>& muts, vector<double>& fact_store, vector< vector<int> >& data, vector<haplo>& haps) {
    for (int i=0;i<haps.size();i++) {
        haps[i].q=initial_freqs[i];
    }
    GetHapFits(sigma,haps);
    lL=0;
    ofstream inf_file;
    if (verb==1) {
        inf_file.open("Inf_DD.out");
    }
    for (int t=0;t<max;t++) {
        DelayDetMutation(mu,beta,muts,haps);
        double glob_fit=GetGlobalFitQ(haps);
        DetNextGen(glob_fit,haps);
        if (t%dt==0) {
            vector<double> inf;
            for (int i=0;i<haps.size();i++) {
                if (haps[i].q>1e-20) {
                    inf.push_back(haps[i].q);
                } else {
                    inf.push_back(1e-20);
                }
                if (verb==1) {
                    inf_file << haps[i].q << " ";
                }
            }
            if (verb==1) {
                inf_file << "\n";
            }
            int tdt=t/dt;
            double logL=MultiCalc(N_s,data[tdt],inf,fact_store);
            lL=lL+logL;
        }
    }
    //cout << "lL " << lL << "\n";
}


void DoPropagationDDR (int verb, int N_s, double mu, double beta, int dt, int max, double sigma, double& lL, vector<double>& initial_freqs, vector<mut>& muts, vector<double>& fact_store, vector< vector<int> >& data, vector<haplo>& haps) {
    for (int i=0;i<haps.size();i++) {
        haps[i].q=initial_freqs[i];
    }
    GetHapFits(sigma,haps);
    lL=0;
    ofstream inf_file;
    if (verb==1) {
        inf_file.open("Inf_rapid_DD.out");
    }
    for (int t=0;t<max;t++) {
        DelayDetMutation(mu,beta,muts,haps);
        double glob_fit=GetGlobalFitQ(haps);
        DetNextGen(glob_fit,haps);
        if (t%dt==0) {
            vector<double> inf;
            for (int i=0;i<haps.size();i++) {
                if (haps[i].q>1e-20) {
                    inf.push_back(haps[i].q);
                } else {
                    inf.push_back(1e-20);
                }
                if (verb==1) {
                    inf_file << haps[i].q << " ";
                }
            }
            if (verb==1) {
                inf_file << "\n";
            }
            int tdt=t/dt;
            double logL=MultiCalc(N_s,data[tdt],inf,fact_store);
            lL=lL+logL;
        }
    }
    //cout << "lL " << lL << "\n";
}


void DoPropagationN (int verb, int N, int N_s, double mu, int dt, int max, double sigma, double& lL, vector<int>& initial_vals, vector<mut>& muts, vector<double>& fact_store, vector< vector<int> >& data, vector<haplo>& haps,gsl_rng *rgen) {
    
    //Alter this to run 100 times (N.B. may need to store values as likelihood calculated from former values...
    
    for (int i=0;i<haps.size();i++) {
        haps[i].n=initial_vals[i];
    }
    GetHapFits(sigma,haps);
    lL=0;
    ofstream inf_file;
    if (verb==1) {
        inf_file.open("Inf.out");
    }
    for (int t=0;t<max;t++) {
        CalcMutation(mu,haps,muts,rgen);
        double glob_fit=GetGlobalFit(haps);
        vector<double> p_select;
        CalcProbs(glob_fit,haps,p_select);
        GetNextGen(N,haps,p_select,rgen);
        if (t%dt==0) {
            vector<double> inf;
            for (int i=0;i<haps.size();i++) {
                haps[i].q=haps[i].n/N;
                if (haps[i].q>1e-20) {
                    inf.push_back(haps[i].q);
                } else {
                    inf.push_back(1e-20);
                }
                if (verb==1) {
                    inf_file << haps[i].q << " ";
                }
            }
            if (verb==1) {
                inf_file << "\n";
            }
            int tdt=t/dt;
            double logL=MultiCalc(N_s,data[tdt],inf,fact_store);
            lL=lL+logL;
        }
    }
}


void FindLogFact(vector<double>& fact_store,int N){
    double logN=0;
    fact_store.push_back(0);
    for (int i=1;i<=N;i++) {
        logN=logN+log(i);
        fact_store.push_back(logN);
    }
}


double MultiCalc(int N, vector<int> obs, vector<double> inf, vector<double> fact_store) {
    double bin=fact_store[N];
    for (unsigned int i=0;i<obs.size();i++) {
        bin=bin-fact_store[obs[i]];
    }
    for (unsigned int i=0;i<obs.size();i++) {
        bin=bin+(obs[i]*log(inf[i]));
    }
    return(bin);
}

