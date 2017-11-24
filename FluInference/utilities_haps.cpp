#include "shared_haps.h"
#include <iostream>
#include <string>
#include <sstream>

void GetOptions (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mcmc=0;
	p.data_sets=27;
	p.hap_file="../Haps1.dat";
	p.iterations=100000;
	p.mu=0.333333333333e-5;
	p.c=95.9096;
	p.init_freqs="Initial_freqs.in";
	p.sel_coeffs="Sel_coeffs.in";
	p.seed=atoi(argv[1]);
	p.data=argv[2];
	int x=3;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--sel_coeffs")==0) {
			x++;
			p.sel_coeffs=argv[x];
		} else if (p_switch.compare("--init_freqs")==0) {
			x++;
			p.init_freqs=argv[x];
		} else if (p_switch.compare("--hap_file")==0) {
			x++;
			p.hap_file=argv[x];
		} else if (p_switch.compare("--datasets")==0) {
			x++;
			p.data_sets=atoi(argv[x]);
		} else if (p_switch.compare("--mcmc")==0) {
			x++;
			p.mcmc=atoi(argv[x]);
		} else if (p_switch.compare("--it")==0) {
			x++;
			p.iterations=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void GetOptionsThree (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mcmc=0;
	p.getmatrix=0;
	p.data_sets=27;
	p.hap_file="../Haps1.dat";
	p.iterations=100000;
	p.mu=0.333333333333e-5;
	p.c=95.9096;
	p.init_freqs1="Initial_freqs1.in";
	p.init_freqs2="Initial_freqs2.in";
	p.init_freqs3="Initial_freqs3.in";
	p.sel_coeffs="Sel_coeffs.in";
	p.seed=atoi(argv[1]);
	int x=2;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--sel_coeffs")==0) {
			x++;
			p.sel_coeffs=argv[x];
		} else if (p_switch.compare("--init_freqs1")==0) {
			x++;
			p.init_freqs1=argv[x];
		} else if (p_switch.compare("--init_freqs2")==0) {
			x++;
			p.init_freqs2=argv[x];
		} else if (p_switch.compare("--init_freqs3")==0) {
			x++;
			p.init_freqs3=argv[x];
		} else if (p_switch.compare("--hap_file")==0) {
			x++;
			p.hap_file=argv[x];
		} else if (p_switch.compare("--datasets")==0) {
			x++;
			p.data_sets=atoi(argv[x]);
		} else if (p_switch.compare("--mcmc")==0) {
			x++;
			p.mcmc=atoi(argv[x]);
		} else if (p_switch.compare("--it")==0) {
			x++;
			p.iterations=atoi(argv[x]);
		} else if (p_switch.compare("--getmat")==0) {
			x++;
			p.getmatrix=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}


void GetOptionsEight (run_params& p, int argc, const char **argv) {
	string p_switch;
	p.read=0;
	p.unobs=0;
	p.seed=0;
	p.timefile=0;
	p.mcmc=0;
	p.data_sets=61;
	p.hap_file="Haps1.dat";
	p.iterations=100000;
	p.mu=0.333333333333e-5;
	p.c=95.9096;
	p.init_freqs1="Initial_freqs1.in";
	p.init_freqs2="Initial_freqs2.in";
	p.init_freqs3="Initial_freqs3.in";
	p.init_freqs4="Initial_freqs4.in";
	p.init_freqs5="Initial_freqs5.in";
	p.init_freqs6="Initial_freqs6.in";
	p.init_freqs7="Initial_freqs7.in";
	p.init_freqs8="Initial_freqs8.in";
	p.sel_coeffs="Sel_coeffs.in";
	p.seed=atoi(argv[1]);
	int x=2;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--mu")==0) {
			x++;
			p.mu=atof(argv[x]);
		} else if (p_switch.compare("--c")==0) {
			x++;
			p.c=atof(argv[x]);
		} else if (p_switch.compare("--read")==0) {
			x++;
			p.read=atoi(argv[x]);
		} else if (p_switch.compare("--sel_coeffs")==0) {
			x++;
			p.sel_coeffs=argv[x];
		} else if (p_switch.compare("--init_freqs1")==0) {
			x++;
			p.init_freqs1=argv[x];
		} else if (p_switch.compare("--init_freqs2")==0) {
			x++;
			p.init_freqs2=argv[x];
		} else if (p_switch.compare("--init_freqs3")==0) {
			x++;
			p.init_freqs3=argv[x];
		} else if (p_switch.compare("--hap_file")==0) {
			x++;
			p.hap_file=argv[x];
		} else if (p_switch.compare("--datasets")==0) {
			x++;
			p.data_sets=atoi(argv[x]);
		} else if (p_switch.compare("--mcmc")==0) {
			x++;
			p.mcmc=atoi(argv[x]);
		} else if (p_switch.compare("--it")==0) {
			x++;
			p.iterations=atoi(argv[x]);
		} else if (p_switch.compare("--getmat")==0) {
			x++;
			p.getmatrix=atoi(argv[x]);
		} else if (p_switch.compare("--unobs")==0) {
			x++;
			p.unobs=atoi(argv[x]);
		} else if (p_switch.compare("--readtimes")==0) {
			x++;
			p.timefile=atoi(argv[x]);
		}
		else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
}

void FindLogFact(vector<double>& fact_store,int N){
	double logN=0;
	fact_store.push_back(0);
	for (int i=1;i<=N;i++) {
		logN=logN+log(i);
		fact_store.push_back(logN);
		//cout << "fact_store "<<i<<" "<<gsl_vector_get(fact_store,i)<<"\n";
	}
}

void Times0 (vector<int>& times) {
	times.push_back(1);
}

void Times1 (vector<int>& times) {
	times.push_back(0);
	times.push_back(1);
	times.push_back(3);
	times.push_back(5);
}

void Times2 (vector<int>& times) {
	times.push_back(1);
	times.push_back(3);
}

void Times3 (vector<int>& times) {
	times.push_back(1);
	times.push_back(3);
	times.push_back(5);
}

void ImportData (run_params p, vector<const char*> contribs, vector<const char*> haps, vector<int> times, vector< vector<hapdat> >& dat) {
	ifstream in_file;
	ifstream con_file;
	string s;
	string line;
	for (int i=0;i<p.data_sets;i++) {
		vector<hapdat> dset;
		in_file.open(haps[i]);
		con_file.open(contribs[i]);
		do {
			hapdat h;
			if (!(in_file >> s)) break;
			//		cout << s << "\n";
			h.st=s;
			int k;
			for (int i=0;i<times.size();i++) {
				in_file >> k;
				h.obs.push_back(k);
				h.inf.push_back(0);
			}
			for (int i=0;i<s.size();i++) {
				h.seq.push_back(s[i]);
			}
			getline(con_file,line);
			istringstream is(line);
			//		cout << line << "\n";
			for (int n; is >> n; ) {
				h.contribs.push_back(n);
			}
			dset.push_back(h);
		} while (1==1);
		dat.push_back(dset);
		in_file.close();
		con_file.close();
	}

}

void CalculateN (vector<int> times, vector< vector<hapdat> > dat, vector< vector<int> >& N) {
	for (int i=0;i<dat.size();i++) {
		vector<int> n;
		for (unsigned int t=0;t<times.size();t++) {
			n.push_back(0);
		}
		for (int j=0;j<dat[i].size();j++) {
			for (int k=0;k<dat[i][j].obs.size();k++) {
				n[k]=n[k]+dat[i][j].obs[k];
			}
		}
		N.push_back(n);
	}
}

void OrganiseHaps (run_params p, vector< vector<hapdat> > dat, int& dim, vector<mod>& mod_haps) {
	ifstream haps_file;
	cout << p.hap_file << "\n";
	string s;
	haps_file.open(p.hap_file);
	do {
		mod m;
		if (!(haps_file >> s)) break;
		cout << s << "\n";
		m.st=s;
		for (int i=0;i<s.size();i++) {
			m.seq.push_back(s[i]);
		}
		mod_haps.push_back(m);
	} while (1==1);
	dim=mod_haps[0].seq.size();

	for (unsigned int i=0;i<mod_haps.size();i++) {
		int count=0;
		for (unsigned int j=0;j<dat.size();j++) {
			for (unsigned int k=0;k<dat[j].size();k++) {
				for (unsigned int l=0;l<dat[j][k].contribs.size();l++) {
					if (dat[j][k].contribs[l]==i+1) {
						mod_haps[i].set.push_back(j);
						mod_haps[i].num.push_back(k);
						count++;
						//cout << i << " " << j << " " << k << "\n";
					}
					
					//cout << j << " " << k << " " << dat[j][k].contribs[l] << "\n";
				}
			}
		}
	}
}

void ReadSelectionModel(int dim, vector<char>& sel_model, vector<int>& tds, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5) {
	//Identify which alleles are under selection / epistasis
	ifstream sel_file;
	sel_file.open("Selection.in");
	for (int i=0;i<dim;i++) {
		char k='0'; //Default of no selection
		sel_file >> k;
		sel_model.push_back(k);
		cout << k << " ";
	}
	cout << "\n";
	sel_file.close();
	ifstream tds_file;
	tds_file.open("TDS.in");
	for (int i=0;i<dim;i++) {
		int k=0; //Default of no selection
		tds_file >> k;
		tds.push_back(k);
		cout << k << " ";
	}
	cout << "\n";
	tds_file.close();
	
	ifstream epi_file;
	epi_file.open("TwoWay.in");
	for (int i=0;i<1000;i++) {
		epi e;
		int l1=0;
		int l2=0;
		if (!(epi_file >> l1)) break;
		if (!(epi_file >> l2)) break;
		l1--;
		l2--;
		e.loc.push_back(l1);
		e.loc.push_back(l2);
		e.x=0;
		epistat2.push_back(e);
	}
	epi_file.close();
	
	epi_file.open("ThreeWay.in");
	for (int i=0;i<1000;i++) {
		epi e;
		int l1=0;
		int l2=0;
		int l3=0;
		if (!(epi_file >> l1)) break;
		if (!(epi_file >> l2)) break;
		if (!(epi_file >> l3)) break;
		l1--;
		l2--;
		l3--;
		e.loc.push_back(l1);
		e.loc.push_back(l2);
		e.loc.push_back(l3);
		e.x=0;
		epistat3.push_back(e);
	}
	epi_file.close();
	
	epi_file.open("FourWay.in");
	for (int i=0;i<1000;i++) {
		epi e;
		int l1=0;
		int l2=0;
		int l3=0;
		int l4=0;
		if (!(epi_file >> l1)) break;
		if (!(epi_file >> l2)) break;
		if (!(epi_file >> l3)) break;
		if (!(epi_file >> l4)) break;
		l1--;
		l2--;
		l3--;
		l4--;
		e.loc.push_back(l1);
		e.loc.push_back(l2);
		e.loc.push_back(l3);
		e.loc.push_back(l4);
		e.x=0;
		epistat4.push_back(e);
	}
	epi_file.close();
	
	epi_file.open("FiveWay.in");
	for (int i=0;i<1000;i++) {
		epi e;
		int l1=0;
		int l2=0;
		int l3=0;
		int l4=0;
		int l5=0;
		if (!(epi_file >> l1)) break;
		if (!(epi_file >> l2)) break;
		if (!(epi_file >> l3)) break;
		if (!(epi_file >> l4)) break;
		if (!(epi_file >> l5)) break;
		l1--;
		l2--;
		l3--;
		l4--;
		l5--;
		e.loc.push_back(l1);
		e.loc.push_back(l2);
		e.loc.push_back(l3);
		e.loc.push_back(l4);
		e.loc.push_back(l5);
		e.x=0;
		epistat5.push_back(e);
	}
	epi_file.close();
	
}

void MakeMutationMatrix (int dim, double mu, vector< vector<double> >& mut, vector<mod> mod_haps) {
	//Check distances between haplotypes
	vector< vector<int> > m;
	for (int i=0;i<mod_haps.size();i++) {
		vector<int> row;
		for (int j=0;j<mod_haps.size();j++) {
			int dist=0;
			for (int k=0;k<mod_haps[0].seq.size();k++) {
				if (mod_haps[i].seq[k]!=mod_haps[j].seq[k]) {
					dist++;
				}
			}
			if (dist==1) {
				if (i<j) {
					cout << i << " " << j << " " << dist << "\n";
				}
				row.push_back(1);
			} else {
				row.push_back(0);
			}
		}
		m.push_back(row);
	}
	
	vector<int> tot;
	for (int i=0;i<m.size();i++) {
		int t=0;
		for (int j=0;j<m[i].size();j++) {
			//                      cout << m[i][j] << " ";
			if (m[i][j]==1) {
				t++;
			}
		}
		cout << t << " " << (3*dim)-t << " ";
		tot.push_back(t);
		cout << "\n";
	}
	
	for (int i=0;i<m.size();i++) {
		vector<double> mrow;
		for (int j=0;j<m[i].size();j++) {
			if (m[i][j]==1) {
				mrow.push_back(mu);
			} else if (i==j) {
				mrow.push_back(pow(1-mu,3*dim));  //Incorporate mutations to other haplotypes
			} else {
				mrow.push_back(0);
			}
		}
		mrow.push_back(0);
		mut.push_back(mrow);
	}
	vector<double> mrow;
	for (int i=0;i<m[0].size();i++) {
		mrow.push_back(mu*((3*dim)-tot[i]));
	}
	mrow.push_back(1);
	mut.push_back(mrow);
	
	cout << "Mutation matrix\n";
	for (int i=0;i<mut.size();i++) {
		for (int j=0;j<mut[i].size();j++) {
			cout << mut[i][j] << " ";
		}
		cout << "\n";
	}
	
}

void SquareMutationMatrix (vector<vector<double> > &m ) {
	//	cout << "M size " << m.size() << "\n";
	vector<vector<double> > m2;
	for (unsigned int i=0;i<m.size();i++) {
		vector<double> m2i;
		m2i.assign(m.size(),0);
		m2.push_back(m2i);
	}
	//	cout << "Now multiplying...\n";
	for (unsigned int i=0;i<m.size();i++) {
		for (unsigned int j=0;j<m[i].size();j++) {
			for (unsigned int k=0;k<m[i].size();k++) {
				m2[i][j]=m2[i][j]+(m[i][k]*m[k][j]);
			}
		}
	}
	m=m2;
}

void SetInitialFreqs (vector<double>& init_freqs, vector< vector<double> > mut, vector<mod> mod_haps, gsl_rng *rgen) {
	for (int i=0;i<mut.size();i++) {
		double x=1e-20;
		if (i<=mod_haps.size()) {  //Equality includes one additional point for non-included haplotypes
			x=gsl_rng_uniform(rgen);
		}
		init_freqs.push_back(x);
		//	cout << x << " ";
	}
	NormaliseFreqs(init_freqs);
	//cout << "\n";
}

void NormaliseFreqs (vector<double>& init_freqs) {
	double tot=0;
	for (int i=0;i<init_freqs.size();i++) {
		tot=tot+init_freqs[i];
	}
	for (int i=0;i<init_freqs.size();i++) {
		init_freqs[i]=init_freqs[i]/tot;
		//	cout << init_freqs[i] << " ";
	}
	//cout << "\n";
}

void SetInitialSelection (int dim, vector<int> tds, vector<char> sel_model, vector<int> times, vector< vector<double> >& sigs, gsl_rng *rgen) {
	for (int i=0;i<dim;i++) {
		vector<double> sig;
		if (sel_model[i]!='0') {
			if (tds[i]==0) {
				double s=(2*gsl_rng_uniform(rgen))-0.2;
				sig.push_back(s);
			} else {
				for (int t=0;t<times.size()-1;t++) {
					double s=(2*gsl_rng_uniform(rgen))-0.2;
					sig.push_back(s);
				}
			}
		}
		sigs.push_back(sig);
	}
	for (int i=0;i<sigs.size();i++) {
		cout << i << " ";
		for (int j=0;j<sigs[i].size();j++) {
			cout << sigs[i][j] << " ";
		}
		cout << "\n";
	}
}

void SetInitialEpistasis (vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5, gsl_rng *rgen) {
	for (int i=0;i<epistat2.size();i++) {
		double s=gsl_rng_uniform(rgen)-0.5;
		epistat2[i].x=s;
		cout << "Epi2 " << epistat2[i].loc[0] << " " << epistat2[i].loc[1] << " " << epistat2[i].x << "\n";
	}
	for (int i=0;i<epistat3.size();i++) {
		double s=gsl_rng_uniform(rgen)-0.5;
		epistat3[i].x=s;
		cout << "Epi3 " << epistat3[i].loc[0] << " " << epistat3[i].loc[1] << " " << epistat3[i].loc[2] << " " << epistat3[i].x << "\n";
	}
	for (int i=0;i<epistat4.size();i++) {
		double s=gsl_rng_uniform(rgen)-0.5;
		epistat4[i].x=s;
		cout << "Epi4 " << epistat4[i].loc[0] << " " << epistat4[i].loc[1] << " " << epistat4[i].loc[2] << " " << epistat4[i].loc[3] << " " << epistat4[i].x << "\n";
	}
	for (int i=0;i<epistat5.size();i++) {
		double s=gsl_rng_uniform(rgen)-0.5;
		epistat5[i].x=s;
		cout << "Epi5 " << epistat5[i].loc[0] << " " << epistat5[i].loc[1] << " " << epistat5[i].loc[2] << " " << epistat5[i].loc[3] << " " << epistat5[i].loc[4] << " " << epistat5[i].x << "\n";
	}
}

void CheckEpistasis (int& go, int dim, vector<char> sel_model, vector<mod> mod_haps, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5) {
	//	mod_haps[i].seq
	
	cout << "Check epistasis\n";
	int sugg=0;
	int seen=0;
	for (int i=0;i<epistat2.size();i++) {
		sugg++;
		int f2=0;
		for (int j=0;j<dim;j++) {
			for (int k=0;k<dim;k++) {
				if (epistat2[i].loc[0]==j&&epistat2[i].loc[1]==k) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
					//cout << "j= " << j << " k= " << k << " " << sel_model[j] << " " << sel_model[k] << "\n";
					int d=0;
					for (unsigned int s=0;s<mod_haps.size();s++) {
						if (mod_haps[s].seq[j]==sel_model[j]&&mod_haps[s].seq[k]==sel_model[k]) {
							d=1;
						}
						if (f2==0&&d==1) {
							cout << "Have two-way interaction " << i << " " << epistat2[i].loc[0] << " " << epistat2[i].loc[1] << "\n";
							f2=1;
							seen++;
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat3.size();i++) {
		sugg++;
		int f3=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					if (epistat3[i].loc[0]==k1&&epistat3[i].loc[1]==k2&&epistat3[i].loc[2]==k3) { //Have an epistatic effect of epistat3[i].x on triple mutants k1,k2,k3
						int d=0;
						for (unsigned int s=0;s<mod_haps.size();s++) {
							if (mod_haps[s].seq[k1]==sel_model[k1]&&mod_haps[s].seq[k2]==sel_model[k2]&&mod_haps[s].seq[k3]==sel_model[k3]) {
								d=1;
							}
							if (f3==0&&d==1) {
								cout << "Have three-way interaction " << i << " " << epistat3[i].loc[0] << " " << epistat3[i].loc[1] << " " << epistat3[i].loc[2] << "\n";
								f3=1;
								seen++;
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat4.size();i++) {
		sugg++;
		int f4=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						if (epistat4[i].loc[0]==k1&&epistat4[i].loc[1]==k2&&epistat4[i].loc[2]==k3&&epistat4[i].loc[3]==k4) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
							int d=0;
							for (unsigned int s=0;s<mod_haps.size();s++) {
								if (mod_haps[s].seq[k1]==sel_model[k1]&&mod_haps[s].seq[k2]==sel_model[k2]&&mod_haps[s].seq[k3]==sel_model[k3]&&mod_haps[s].seq[k4]==sel_model[k4]) {
									d=1;
								}
								if (f4==0&&d==1) {
									cout << "Have four-way interaction " << i << " " << epistat4[i].loc[0] << " " << epistat4[i].loc[1] << " " << epistat4[i].loc[2] << " " << epistat4[i].loc[3] << "\n";
									f4=1;
									seen++;
								}
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat5.size();i++) {
		sugg++;
		int f5=0;
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						for (int k5=0;k5<dim;k5++) {
							if (epistat5[i].loc[0]==k1&&epistat5[i].loc[1]==k2&&epistat5[i].loc[2]==k3&&epistat5[i].loc[3]==k4&&epistat5[i].loc[4]==k5) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
								int d=0;
								for (unsigned int s=0;s<mod_haps.size();s++) {
									if (mod_haps[s].seq[k1]==sel_model[k1]&&mod_haps[s].seq[k2]==sel_model[k2]&&mod_haps[s].seq[k3]==sel_model[k3]&&mod_haps[s].seq[k4]==sel_model[k4]&&mod_haps[s].seq[k5]==sel_model[k5]) {
										d=1;
									}
									if (f5==0&&d==1) {
										cout << "Have five-way interaction " << i << " " << epistat5[i].loc[0] << " " << epistat5[i].loc[1] << " " << epistat5[i].loc[2] << " " << epistat5[i].loc[3] << " " << epistat5[i].loc[4] << "\n";
										f5=1;
										seen++;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if (seen<sugg) {
		cout << "Not all epistatic interactions were observed\n";
		go=0;
	}
}

void SigsHapSigs (int verb, int dim, double nullsig, run_params p, vector<int> times, vector<char> sel_model, vector<mod> mod_haps, vector< vector<double> > sigs, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5, vector< vector<double> >& hapsigs) {
	for (int i=0;i<mod_haps.size();i++) {
		vector<double> hs;
		for (int j=0;j<times.size()-1;j++) {
			hs.push_back(0);
		}
		for (int j=0;j<dim;j++) {
			if (sigs[j].size()==1) {
				if (mod_haps[i].seq[j]==sel_model[j]) {
					for (int k=0;k<times.size()-1;k++) {
						hs[k]=hs[k]+sigs[j][0];
					}
				}
			}
			if (sigs[j].size()>1) {
				if (mod_haps[i].seq[j]==sel_model[j]) {
					for (int k=0;k<times.size()-1;k++) {
						hs[k]=hs[k]+sigs[j][k];
					}
				}
			}
		}
		hapsigs.push_back(hs);
	}
	
	//Epistasis 2-way
	for (int i=0;i<epistat2.size();i++) {
		for (int j=0;j<dim;j++) {
			for (int k=0;k<dim;k++) {
				if (epistat2[i].loc[0]==j&&epistat2[i].loc[1]==k) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
					//cout << "j= " << j << " k= " << k << " " << sel_model[j] << " " << sel_model[k] << "\n";
					for (int s=0;s<mod_haps.size();s++) {
						if (mod_haps[s].seq[j]==sel_model[j]&&mod_haps[s].seq[k]==sel_model[k]) {
							//cout << s << " " << mod_haps[s].seq[j] << " " << mod_haps[s].seq[k] << "\n";
							for (int t=0;t<times.size()-1;t++) {
								hapsigs[s][t]=hapsigs[s][t]+epistat2[i].x;
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat3.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					if (epistat3[i].loc[0]==k1&&epistat3[i].loc[1]==k2&&epistat3[i].loc[2]==k3) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
						for (int s=0;s<mod_haps.size();s++) {
							if (mod_haps[s].seq[k1]==sel_model[k1]&&mod_haps[s].seq[k2]==sel_model[k2]&&mod_haps[s].seq[k3]==sel_model[k3]) {
								for (int t=0;t<times.size()-1;t++) {
									hapsigs[s][t]=hapsigs[s][t]+epistat3[i].x;
								}
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat4.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						if (epistat4[i].loc[0]==k1&&epistat4[i].loc[1]==k2&&epistat4[i].loc[2]==k3&&epistat4[i].loc[3]==k4) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
							for (int s=0;s<mod_haps.size();s++) {
								if (mod_haps[s].seq[k1]==sel_model[k1]&&mod_haps[s].seq[k2]==sel_model[k2]&&mod_haps[s].seq[k3]==sel_model[k3]&&mod_haps[s].seq[k4]==sel_model[k4]) {
									for (int t=0;t<times.size()-1;t++) {
										hapsigs[s][t]=hapsigs[s][t]+epistat4[i].x;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	for (int i=0;i<epistat5.size();i++) {
		for (int k1=0;k1<dim;k1++) {
			for (int k2=0;k2<dim;k2++) {
				for (int k3=0;k3<dim;k3++) {
					for (int k4=0;k4<dim;k4++) {
						for (int k5=0;k5<dim;k5++) {
							if (epistat4[i].loc[0]==k1&&epistat5[i].loc[1]==k2&&epistat5[i].loc[2]==k3&&epistat5[i].loc[3]==k4&&epistat5[i].loc[4]==k5) { //Have an epistatic effect of epistat2[i].x on double mutants j,k
								for (int s=0;s<mod_haps.size();s++) {
									if (mod_haps[s].seq[k1]==sel_model[k1]&&mod_haps[s].seq[k2]==sel_model[k2]&&mod_haps[s].seq[k3]==sel_model[k3]&&mod_haps[s].seq[k4]==sel_model[k4]&&mod_haps[s].seq[k5]==sel_model[k5]) {
										for (int t=0;t<times.size()-1;t++) {
											hapsigs[s][t]=hapsigs[s][t]+epistat5[i].x;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	
	//Non-observed state
	vector<double> hs;
	for (int j=0;j<times.size()-1;j++) {
		//if (p.unobs==1) {
		hs.push_back(nullsig);
		//} else if (p.unobs==2) {
		//	hs.push_back(-10);
		//}
	}
	hapsigs.push_back(hs);
	
	for (int i=0;i<hapsigs.size();i++) {
		if (hapsigs[i].size()==0) {
			hapsigs[i].push_back(0);
		}
	}
	
	if (verb==1) {
		cout << "Hapsigs\n";
		for (int i=0;i<hapsigs.size();i++) {
			cout << "i " << i << " ";
			for (int j=0;j<hapsigs[i].size();j++) {
				cout << hapsigs[i][j] << " ";
			}
			cout << "\n";
		}
	}
	
}

void PrintFreqs1(vector<double> init_freqs) {
	cout << "New freqs \n";
	for (int i=0;i<init_freqs.size();i++) {
		cout << init_freqs[i] << " ";
	}
	cout << "\n";
}


void PrintFreqs3(vector<double> init_freqs1,vector<double> init_freqs2, vector<double> init_freqs3) {
	cout << "New freqs 1\n";
	for (int i=0;i<init_freqs1.size();i++) {
		cout << init_freqs1[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 2\n";
	for (int i=0;i<init_freqs2.size();i++) {
		cout << init_freqs2[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 3\n";
	for (int i=0;i<init_freqs3.size();i++) {
		cout << init_freqs3[i] << " ";
	}
	cout << "\n";
}

void PrintFreqs8(vector<double> init_freqs1,vector<double> init_freqs2, vector<double> init_freqs3, vector<double> init_freqs4, vector<double> init_freqs5, vector<double> init_freqs6, vector<double> init_freqs7, vector<double> init_freqs8) {
	cout << "New freqs 1\n";
	for (int i=0;i<init_freqs1.size();i++) {
		cout << init_freqs1[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 2\n";
	for (int i=0;i<init_freqs2.size();i++) {
		cout << init_freqs2[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 3\n";
	for (int i=0;i<init_freqs3.size();i++) {
		cout << init_freqs3[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 4\n";
	for (int i=0;i<init_freqs4.size();i++) {
		cout << init_freqs4[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 5\n";
	for (int i=0;i<init_freqs5.size();i++) {
		cout << init_freqs5[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 6\n";
	for (int i=0;i<init_freqs6.size();i++) {
		cout << init_freqs6[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 7\n";
	for (int i=0;i<init_freqs7.size();i++) {
		cout << init_freqs7[i] << " ";
	}
	cout << "\n";
	cout << "New freqs 8\n";
	for (int i=0;i<init_freqs8.size();i++) {
		cout << init_freqs8[i] << " ";
	}
	cout << "\n";
}

void PrintSigEpi (vector< vector<double> >& sigs, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5) {
	cout << "New sigma\n";
	for	(int i=0;i<sigs.size();i++) {
		cout << i << " ";
		for (int j=0;j<sigs[i].size();j++) {
			cout << sigs[i][j] << " ";
		}
		cout << "\n";
	}
	
	for (int i=0;i<epistat2.size();i++) {
		cout << "Epistasis2 " << epistat2[i].loc[0] << " " << epistat2[i].loc[1] << " " << epistat2[i].x << "\n";
	}
	for (int i=0;i<epistat3.size();i++) {
		cout << "Epistasis3 " << epistat3[i].loc[0] << " " << epistat3[i].loc[1] << " "  << epistat3[i].loc[2] << " " << epistat3[i].x << "\n";
	}
	for (int i=0;i<epistat4.size();i++) {
		cout << "Epistasis4 " << epistat4[i].loc[0] << " " << epistat4[i].loc[1] << " "  << epistat4[i].loc[2] << " " << epistat4[i].loc[3] << " " << epistat4[i].x << "\n";
	}
	for (int i=0;i<epistat5.size();i++) {
		cout << "Epistasis5 " << epistat5[i].loc[0] << " " << epistat5[i].loc[1] << " "  << epistat5[i].loc[2] << " " << epistat5[i].loc[3] << " " << epistat5[i].loc[4] << " " << epistat5[i].x << "\n";
	}

}

void Propagation(int verb, vector<int> times, vector<double> init_freqs, vector< vector<double> > mut, vector< vector<double> > hapsigs, vector< vector<double> >& inf) {
	//cout << "Size " << init_freqs.size() << "\n";
//	if (verb==1) {
//		cout << "Times\n";
//		for (int i=0;i<times.size();i++) {
//			cout << times[i] << " ";
//		}
//		cout << "\n";
//	}
	inf.clear();
	//	int start=times[0];
	int start=0;  //Infection is assumed to begin at time point zero.  Code outputs inferred frequencies at time-points in the times vector.
	int fin=times[times.size()-1];
	int index=0;
	vector<double> si;
	GetSI(index,hapsigs,si);
	vector<double> q=init_freqs;
	
	
	
	for (int j=0;j<=fin;j++) {

	//	if (verb==1) {
	//		cout << "j= " << j << "\n";
	//	}
		
		if (j==0&&j==times[index]) {
			//	cout << "Time " << times[index] << "\n";
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf.push_back(q);
			index++;
		}
		if (j>start) {
			//if (verb==1) {
			//	cout << "Grow Generation\n";
			//}
			GrowGeneration(mut,q,si);
		}
		if (j==times[index]) {
			//	cout << "Time " << times[index] << "\n";
			GetSI(index,hapsigs,si);
			index++;
			if (verb==1) {
				for (int i=0;i<q.size();i++) {
					cout << q[i] << " ";
				}
				cout << "\n";
			}
			inf.push_back(q);
			
		}
	}
}

void GetSI (int index, vector< vector<double> > hapsigs, vector<double>& si) {
	si.clear();
	for (int i=0;i<hapsigs.size();i++) {
		si.push_back(hapsigs[i][index]);
	}
}


void GrowGeneration(vector<vector<double> >& m, vector<double>& q, vector<double>& sigma) {
	vector<double> temp;
	MatMult(m,q,temp);  //Now using m^2
	q=temp;
//	cout << "Print q post mut \n";
//	for (unsigned int i=0;i<q.size();i++) {
//		cout << q[i] << " ";
//	}
//	cout << "\n\n";
	GrowSig(q,sigma);
	MatMult(m,q,temp);
	q=temp;
//	cout << "Print q post mut \n";
//	for (unsigned int i=0;i<q.size();i++) {
//		cout << q[i] << " ";
//	}
//	cout << "\n\n";

	
	GrowSig(q,sigma);
}

void MatMult (vector<vector<double> >& m, vector<double>& v, vector<double>& mv) {
	mv.clear();
	mv.assign(v.size(),0);
	for (unsigned int i=0;i<v.size();i++) {
		for (unsigned int j=0;j<v.size();j++) {
			mv[i]=mv[i]+m[i][j]*v[j];
		}
	}
}

void GrowSig (vector<double>& q, vector<double>& sigma) {
	double tot=0;
	//cout << "Print q post sig\n";
	//for (unsigned int i=0;i<q.size();i++) {
	//	cout << q[i] << " ";
	//}
	//cout << "\n";
	//for (unsigned int i=0;i<q.size();i++) {
	//	cout << exp(sigma[i]) << " ";
	//}
	//cout << "\n";

	
	for (unsigned int i=0;i<q.size();i++) {
		q[i]=q[i]*exp(sigma[i]);
		tot=tot+q[i];
	}
	for (unsigned int i=0;i<q.size();i++) {
		q[i]=q[i]/tot;
	}
	//for (unsigned int i=0;i<q.size();i++) {
	//	cout << q[i] << " ";
	//}
	//cout << "\n\n";

}

void ChangeFreq (int& tryx, int& movex, double changex, vector<double>& init_freqs, gsl_rng *rgen) {
	int j=floor(gsl_rng_uniform(rgen)*(init_freqs.size()));
	init_freqs[j]=init_freqs[j]+(gsl_rng_uniform(rgen)*changex)-(changex/2);
	if(init_freqs[j]<0) {
		init_freqs[j]=0;
	}
	NormaliseFreqs(init_freqs);
	tryx++;
	movex=1;
}

void ChangeSigs (int& trys, int& movex, double changes, vector<int> sel_alls, vector< vector<double> >& sigs, vector<epi>& epistat2, vector<epi>& epistat3, vector<epi>& epistat4, vector<epi>& epistat5, gsl_rng *rgen) {
	if (sel_alls.size()>0) {
		int s=sel_alls.size()+epistat2.size()+epistat3.size()+epistat4.size()+epistat5.size();
		int r=floor(gsl_rng_uniform(rgen)*s);
		//	cout << r << " " << sel_alls.size() << " " << epistat2.size() << " "  << epistat3.size() << "\n";
		if (r<sel_alls.size()) {
			int i=sel_alls[r];
			for (int j=0;j<sigs[i].size();j++) {
				sigs[i][j]=sigs[i][j]+(gsl_rng_uniform(rgen)*changes)-(changes/2);
			}
		} else if (r<sel_alls.size()+epistat2.size()) {
			int i=r-sel_alls.size();
			epistat2[i].x=epistat2[i].x+(gsl_rng_uniform(rgen)*changes)-(changes/2);
			
		} else if (r<sel_alls.size()+epistat2.size()+epistat3.size()) {
			int i=r-sel_alls.size()-epistat2.size();
			epistat3[i].x=epistat3[i].x+(gsl_rng_uniform(rgen)*changes)-(changes/2);
			
		} else if (r<sel_alls.size()+epistat2.size()+epistat3.size()+epistat4.size()) {
			int i=r-sel_alls.size()-epistat2.size()-epistat3.size();
			epistat4[i].x=epistat4[i].x+(gsl_rng_uniform(rgen)*changes)-(changes/2);
			
		} else {
			int i=r-sel_alls.size()-epistat2.size()-epistat3.size()-epistat4.size();
			epistat5[i].x=epistat5[i].x+(gsl_rng_uniform(rgen)*changes)-(changes/2);
		}
		trys++;
		movex=0;
	}
}

void AllocateHapInf (vector<int> times, vector<mod> mod_haps, vector< vector<double> > inf, vector< vector<hapdat> >& dat) {
	//Create inference vector
	for (int j=0;j<dat.size();j++) {
		for (int k=0;k<dat[j].size();k++) {
			fill(dat[j][k].inf.begin(),dat[j][k].inf.end(),1e-40);
		}
	}
	
	//Sum inferred full haplotypes into partial haplotype observations
	for (int i=0;i<times.size();i++) {
		for (int j=0;j<mod_haps.size();j++) {
			for (int k=0;k<mod_haps[j].set.size();k++) {
				dat[mod_haps[j].set[k]][mod_haps[j].num[k]].inf[i]=dat[mod_haps[j].set[k]][mod_haps[j].num[k]].inf[i]+inf[i][j];
			}
		}
	}
	
}

double DirichletMultiCalcNew (int i, int t, double c, vector< vector<int> >& N, vector< vector<hapdat> >& dat, vector<double>& fact_store) {
	vector<int> obs;
	vector<double> inf;
	double bin=0;
	double icorr=1;
	for (int j=0;j<dat[i].size()-1;j++) {
		obs.push_back(dat[i][j].obs[t]);
		inf.push_back(dat[i][j].inf[t]);
		icorr=icorr-dat[i][j].inf[t];
	}
	obs.push_back(dat[i][dat[i].size()-1].obs[t]);
	if (icorr==0) {
		icorr=1e-20;
	}
	inf.push_back(icorr);
	
	/*cout << "Obs\n";
	 for (unsigned int i=0;i<obs.size();i++) {
	 cout << obs[i] << " ";
	 }
	 cout << "\n";
	 cout << "Inf\n";
	 for (unsigned int i=0;i<inf.size();i++) {
		cout << inf[i] << " ";
	 }
	 cout << "\n";
	 if (obs.size()!=inf.size()) {
		cout << "Size match error " << obs.size() << " " << inf.size() << "\n";
	 }*/
	
	if (N[i][t]>0) {
		bin=fact_store[N[i][t]];
		for (unsigned int i=0;i<obs.size();i++) {
			bin=bin-fact_store[obs[i]];
		}
		vector<double> alpha;
		for (unsigned int i=0;i<inf.size();i++) {
			alpha.push_back(c*inf[i]);
		}
		double a=0;
		for (unsigned int i=0;i<alpha.size();i++) {
			a=a+alpha[i];
			bin=bin-gsl_sf_lngamma(alpha[i]);
		}
		bin=bin+gsl_sf_lngamma(a);
		a=0;
		for (unsigned int i=0;i<alpha.size();i++) {
			double b=alpha[i]+obs[i];
			a=a+b;
			bin=bin+gsl_sf_lngamma(b);
		}
		bin=bin-gsl_sf_lngamma(a);
	} else {
		bin=0;
	}
	//	cout << "L " << bin << "\n";
	//	cout << "\n";
	
	return(bin);
}



