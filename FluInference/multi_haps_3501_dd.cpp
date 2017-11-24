//Applies linear model to unobserved haplotypes.  Retains loci for which there is no variation.

#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>

using namespace std;

#include "shared_haps_dd.h"

rec data;
int main(int argc, const char **argv) {

	run_params p;
	GetOptionsThree (p,argc,argv);
	
	//Initialise random number generator
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, p.seed);
	
	vector<double> fact_store;
	FindLogFact(fact_store,50000);

	//Observation times
	vector<int> times2;
	Times1(times2);
	
	//Load F3501 data
	vector<const char*> contribs;
	vector<const char*> haps;
	GetF3501Files (contribs,p,haps);
	vector< vector<hapdat> > dat2;
	ImportData (p,contribs,haps,times2,dat2);
	
	//Calculate N for each dataset
	vector< vector<int> > N2;
	CalculateN (times2,dat2,N2);

	//Load model haplotypes and calculate haplotype contributions to observed data
	vector<mod> mod_haps1;
	int dim;
	OrganiseHaps (p,dat2,dim,mod_haps1);

	vector<char> sel_model; //Code describing selection or non-selection on each locus e.g. A C 0 0.  0: Neutral, A: Constant selection for allele A
	vector<int> tds;
	vector<epi> epistat2;
	vector<epi> epistat3;
	vector<epi> epistat4;
	vector<epi> epistat5;
	cout << "Read selection model\n";
	ReadSelectionModel(dim,sel_model,tds,epistat2,epistat3,epistat4,epistat5);

	
	//Construct mutational information.  First check distances
	double mu=p.mu;
	vector< vector<double> > mut;
	if (p.getmatrix==1) {
		ifstream m_file;
		m_file.open("../Data/PartialHaplotype/HAdata/Haps1.mut");
		cout << mod_haps1.size() << "\n";
		for (int i=0;i<=mod_haps1.size();i++) {
			vector<double> m;
			double x;
			for (int j=0;j<=mod_haps1.size();j++) {
				m_file >> x;
				m.push_back(x);
			}
			mut.push_back(m);
		}
	} else {
		MakeMutationMatrix (dim,mu,mut,mod_haps1);
	}
	SquareMutationMatrix(mut);

	cout << "Here\n";

	
	cout << "Squared mutation matrix\n";
	for (int i=0;i<mut.size();i++) {
		for (int j=0;j<mut[i].size();j++) {
			cout << mut[i][j] << " ";
		}
		cout << "\n";
	}

	//Setup initial data for selection model and initial frequency model
	cout << "Initial Frequencies\n";
	vector<double> init_freqs1;
	SetInitialFreqs (init_freqs1,mut,mod_haps1,rgen);
	
	double beta=0.0001;

	if (p.read==1) {
		ifstream ifreqs_file;
		ifreqs_file.open(p.init_freqs1);
		for (int i=0;i<100;i++) {
			double x=0;
			if (!(ifreqs_file >> x)) break;
			init_freqs1[i]=x;
		}
		ifreqs_file.close();
		cout << "Initial freqs read\n";
		for (int i=0;i<init_freqs1.size();i++) {
			cout << init_freqs1[i] << " ";
		}
		cout << "\n";
		ifreqs_file.close();
		ifreqs_file.open("Beta.in");
		for (int i=0;i<100;i++) {
			double x=0;
			if (!(ifreqs_file >> x)) break;
			beta=x;
		}
		
	}

	cout << "Initial selection coefficients\n";
	vector< vector<double> > sigs;
	SetInitialSelection (dim,tds,sel_model,times2,sigs,rgen);
	if (p.read==1) {
		ifstream sc_file;
		sc_file.open(p.sel_coeffs);
		cout << "Selection coefficients read\n";
		
		for (int i=0;i<sigs.size();i++) {
			double x=0;
			for (int j=0;j<sigs[i].size();j++) {
				if (!(sc_file >> x)) break;
				sigs[i][j]=x;
				cout << x << "\n";
			}
		}
		sc_file.close();
	}

	double nullsig=0;
	SetInitialEpistasis (epistat2,epistat3,epistat4,epistat5,rgen);
	if (p.read==1) {
		ifstream ep_file;
		ep_file.open("Epi_stats.in");
		for (int i=0;i<epistat2.size();i++) {
			double x=0;
			if (!(ep_file >> x)) break;
			epistat2[i].x=x;
		}
		ep_file.close();
		
		ep_file.open("Epi_stats3.in");
		for (int i=0;i<epistat3.size();i++) {
			double x=0;
			if (!(ep_file >> x)) break;
			epistat3[i].x=x;
		}
		ep_file.close();
		
		ep_file.open("Epi_stats4.in");
		for (int i=0;i<epistat4.size();i++) {
			double x=0;
			if (!(ep_file >> x)) break;
			epistat4[i].x=x;
		}
		ep_file.close();
		
		ep_file.open("Epi_stats5.in");
		for (int i=0;i<epistat5.size();i++) {
			double x=0;
			if (!(ep_file >> x)) break;
			epistat5[i].x=x;
		}
		ep_file.close();
	}

	//Check epistasis
	int go=1;
	CheckEpistasis(go,dim,sel_model,mod_haps1,epistat2,epistat3,epistat4,epistat5);
	if (go==0) {
		return 0;
	}
	
	vector< vector<double> > hapsigs;
	SigsHapSigs (1,dim,nullsig,p,times2,sel_model,mod_haps1,sigs,epistat2,epistat3,epistat4,epistat5,hapsigs);

	//Do evaluation
	int first=1;
	vector< vector<double> > inf2;
	double L=-1e-10;
	double bestL=-1e10;
	double storeL=-1e10;
	double dL=-1;
	double beta_store=beta;
	double beta_best=beta;
	vector<double> init_freqs1_best;
	vector< vector<double> > sigs_best;
	vector<epi> epistat2_best;
	vector<epi> epistat3_best;
	vector<epi> epistat4_best;
	vector<epi> epistat5_best;
	vector<double> init_freqs1_store;
	vector< vector<double> > sigs_store;
	vector<epi> epistat2_store;
	vector<epi> epistat3_store;
	vector<epi> epistat4_store;
	vector<epi> epistat5_store;
	vector<double> progress;
	double changex=0.05;
	double changes=0.1;
	double c=p.c;
	int acceptx=0;
	int accepts=0;
	int trys=0;
	int tryx=0;
	int movex=-1;
	int maxrep=1;
	if (p.read==1) {
		maxrep=5;
	}
	double a_rate=0;
	
	//Check selection exists
	vector<int> sel_alls;
	for (int i=0;i<sigs.size();i++) {
		if (sigs[i].size()>0) {
			sel_alls.push_back(i);
			cout << "Selection at " << i << "\n";
		}
	}
	
	for (int rep=0;rep<maxrep;rep++) {
		changex=0.05;
		changes=0.01;
		acceptx=0;
		tryx=0;
		for (int it=0;it<p.iterations;it++) {
			
			//if (it%10000==0) {
			//	cout << "Iteration " << it << "\n";
			//}
			if (first==0) {
				//Evaluate likelihood
				dL=L-storeL;
				if (dL>0) {
					progress.push_back(L);
					if ((it>p.iterations-10000)||(it%20==1)) {

						PrintFreqs1(init_freqs1);
						PrintSigEpi (beta,sigs,epistat2,epistat3,epistat4,epistat5);
						cout << "System\n";
						Propagation(1,beta,times2,init_freqs1,mut,hapsigs,inf2);
						cout << "Better L = " << L << "\n";

					}
					if (movex==1) {
						acceptx++;
					}
					if (movex==0) {
						accepts++;
					}
						
					if (progress.size()>200&&progress[progress.size()-51]-progress[progress.size()-1]>-0.0001) {
						PrintFreqs1(init_freqs1);
						PrintSigEpi (beta,sigs,epistat2,epistat3,epistat4,epistat5);
						cout << "System\n";
						Propagation(1,beta,times2,init_freqs1,mut,hapsigs,inf2);
						cout << "Better L = " << L << "\n";

						break;
					}
						
					storeL=L;
					init_freqs1_store=init_freqs1;
					sigs_store=sigs;
					beta_store=beta;
					epistat2_store=epistat2;
					epistat3_store=epistat3;
					epistat4_store=epistat4;
					epistat5_store=epistat5;
					//Check for best so far
					if (L>bestL) {
						bestL=L;
						init_freqs1_best=init_freqs1;
						sigs_best=sigs;
						beta_best=beta;
						epistat2_best=epistat2;
						epistat3_best=epistat3;
						epistat4_best=epistat4;
						epistat5_best=epistat5;
					}
						
				} else {
					init_freqs1=init_freqs1_store;
					sigs=sigs_store;
					beta=beta_store;
					epistat2=epistat2_store;
					epistat3=epistat3_store;
					epistat4=epistat4_store;
					epistat5=epistat5_store;
				}
				
				//cout << "Modify parameters " << "\n";
				//Modify parameter changes
				if (it>10000) {
					if (it%1000==0) {
						a_rate=(acceptx+0.)/(tryx+0.);
						changex=changex*(0.95+a_rate);
						if (trys>0) {
							a_rate=(accepts+0.)/(trys+0.);
							changes=changes*(0.95+a_rate);
						}
						acceptx=0;
						tryx=0;
						accepts=0;
						trys=0;
					}
				}
				
				//Change parameters init_freqs and sigs
				if (it%11==0) {
					beta=log(beta);
					beta=beta+(gsl_rng_uniform(rgen)*changex*20)-(changex*10);
					if (beta>0) {
						beta=-200;
					}
					beta=exp(beta);
				} else {
					if (it%2==0) {
						ChangeFreq (tryx,movex,changex,init_freqs1,rgen);
					}
					if (it%2==1) {
						ChangeSigs (trys,movex,changes,sel_alls,sigs,epistat2,epistat3,epistat4,epistat5,rgen);
					}
				}
			}
			first=0;
			//Convert sigmas to haplotype sigmas
			hapsigs.clear();
			SigsHapSigs (0,dim,nullsig,p,times2,sel_model,mod_haps1,sigs,epistat2,epistat3,epistat4,epistat5,hapsigs);

			//Propagate under mutation/selection
			Propagation(0,beta,times2,init_freqs1,mut,hapsigs,inf2);

		//	cout << "Inference:\n";
		//	for (unsigned int i=0;i<inf.size();i++) {
		//		cout << "Size " << inf[i].size() << "\n";
		//		for (unsigned int j=0;j<inf[i].size();j++) {
		//			cout << inf[i][j] << " ";
		//		}
		//		cout << "\n";
		//	}
			
			//Evaluate likelihood
			
			L=0;
			AllocateHapInf(times2,mod_haps1,inf2,dat2);
			
		//	cout << "Mod values:\n";
		//	for (unsigned int i=0;i<times1.size();i++) {
		//		for (int j=0;j<mod_haps1.size();j++) {
		//			for (int k=0;k<mod_haps1[j].set.size();k++) {
		//				cout << j << " " << k << " " << mod_haps1[j].set[k] << " " << mod_haps1[j].num[k] << "\n";
		//			}
		//			cout << "\n";
		//		}
		//	}
			
	//		cout << "Assigned inference:\n";
	//		for (unsigned int i=0;i<dat.size();i++) {
	//			cout << "i= " << i << "\n";
	//			for (unsigned int j=0;j<dat[i].size();j++) {
	//				cout << "j= " << j << "\n";
	//				for (unsigned int l=0;l<dat[i][j].inf.size();l++) {
	//					cout << dat[i][j].inf[l] << " ";
	//				}
	//				cout << "\n";
	//			}
	//		}

			for (unsigned int i=0;i<dat2.size();i++) {
				for (unsigned int t=0;t<times2.size();t++) {
					double lL=DirichletMultiCalcNew(i,t,c,N2,dat2,fact_store);
					L=L+lL;
				}
			}
			
			//cout << "Log L " << L << "\n";
			
		}
		init_freqs1=init_freqs1_best;
		sigs=sigs_best;
		epistat2=epistat2_best;
		epistat3=epistat3_best;
		epistat4=epistat4_best;
		epistat5=epistat5_best;
	}

	return 0;
	 
}
							
 
						
	
