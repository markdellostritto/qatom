#include "nlo_atom.hpp"

//***********************************************************************************************************************************
//Profile
//***********************************************************************************************************************************

std::ostream& operator<<(std::ostream& out, const Profile& p){
	return out<<"profile("<<p.sigma()<<","<<p.b1()<<","<<p.b2()<<")";
}

//***********************************************************************************************************************************
//NLO class
//***********************************************************************************************************************************

//constants
const double NLO::mevPerThz=0.24180;
const double NLO::cmiPerThz=0.02998;

std::ostream& operator<<(std::ostream& out, const NLO& nlo){
	out<<"*************************************************\n";
	out<<"****************** NLO2D_ATOM ******************\n";
	out<<"FREQ_UNIT          = "<<nlo.freqUnit_<<"\n";
	out<<"\tMIN_FREQ         = "<<nlo.minFreq_<<"\n";
	out<<"\tMAX_FREQ         = "<<nlo.maxFreq_<<"\n";
	out<<"\tFREQ_CUT         = "<<nlo.freqCut_<<"\n";
	out<<"\tFREQ_RES         = "<<nlo.freqRes_<<"\n";
	out<<"\tSTRIDE_CHG       = "<<nlo.strideChg_<<"\n";
	out<<"\tSTRIDE_ALPHA     = "<<nlo.strideAlpha_<<"\n";
	out<<"\tWINDOW           = "<<nlo.windowType_<<"\n";
	out<<"\tFILE_SPECTRUM    = "<<nlo.fileSpectrum_<<"\n";
	out<<"\tCALC_CHG         = "<<nlo.calcChg_<<"\n";
	out<<"\tCALC_ALPHA       = "<<nlo.calcAlpha_<<"\n";
	out<<"\tCALC_SPECTRUM    = "<<nlo.calcSpectrum_<<"\n";
	out<<"\tPRINT_ALPHA_TOT  = "<<nlo.printAlphaT_<<"\n";
	out<<"\tPRINT_DIPOLE_TOT = "<<nlo.printDipoleT_<<"\n";
	out<<"\tPRINT_CHG_TOT    = "<<nlo.printChgT_<<"\n";
	out<<"\tPRINT_CHG        = "<<nlo.printChg_<<"\n";
	out<<"\tNORMALIZE        = "<<nlo.normalize_<<"\n";
	out<<"\tPROFILE_CALC     = "<<nlo.profileCalc_<<"\n";
	out<<"\tPROFILE_LOAD     = "<<nlo.profileLoad_<<"\n";
	out<<"****************** NLO2D_ATOM ******************\n";
	out<<"*************************************************";
	return out;
}

void NLO::defaults(){
	if(DEBUG_NLO_ATOM>0) log<<"defaults()\n";
	//fft parameters
		freqUnit_=fourier::FreqUnit::THZ;
		freqCut_=0;
		freqN_=0;
		minFreq_=0;
		maxFreq_=0;
		freqRes_=0;
		fourierN_=0;
	//window
		sigma_=0;
		window_=window::Identity();
		windowType_=window::WINDOW_FUNC::IDENTITY;
	//stride
		strideAlpha_=1;//include every timestep
		nStepsAlpha_=0;
	//i/o
		fileSpectrum_=std::string("nlo.dat");
	//calculation flags
		calcChg_=true;
		calcAlpha_=true;
		calcSpectrum_=true;
		normalize_=true;
	//i/o flags
		printAlphaT_=false;
		printDipoleT_=false;
		printChgT_=false;
		printChg_=false;
	//chi2 spectrum
		chi2_.clear();
}

void NLO::load(const char* file){
	if(DEBUG_NLO_ATOM>0) log<<"load(const char*):\n";
	//local function variables
	FILE* reader=NULL;
	char* input=(char*)malloc(sizeof(char)*string::M);
	char* temp=(char*)malloc(sizeof(char)*string::M);
	bool error=false;
	
	try{
		//open the parameter file
		if(DEBUG_NLO_ATOM>1) log<<"Opening parameter file...\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open parameter file.");
		
		//set the defaults
		if(DEBUG_NLO_ATOM>1) log<<"Setting defaults...\n";
		defaults();
		
		//read in the parameters
		if(DEBUG_NLO_ATOM>1) log<<"Reading in parameters...\n";
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_all(string::trim_right(input,string::COMMENT));
			string::to_upper(string::copy_left(temp,input,"="));
			if(std::strcmp(temp,"FILE_SPECTRUM")==0){
				fileSpectrum_=std::string(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"CALC_SPECTRUM")==0){
				calcSpectrum_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"CALC_ALPHA")==0){
				calcAlpha_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"CALC_CHG")==0){
				calcChg_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"PRINT_ALPHA_TOT")==0){
				printAlphaT_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"PRINT_DIPOLE_TOT")==0){
				printDipoleT_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"PRINT_CHG_TOT")==0){
				printChgT_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"PRINT_CHG")==0){
				printChg_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"NORMALIZE")==0){
				normalize_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"FREQ_CUT")==0){
				freqCut_=std::atof(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"FREQ_UNIT")==0){
				string::to_upper(std::strcpy(temp,std::strpbrk(input,"=")+1));
				freqUnit_=fourier::FreqUnit::load(temp);
			} else if(std::strcmp(temp,"RESOLUTION")==0){
				freqRes_=std::atof(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"STRIDE_CHARGE")==0){
				strideChg_=std::atoi(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"STRIDE_ALPHA")==0){
				strideAlpha_=std::atoi(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"WINDOW")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				windowType_=window::WINDOW_FUNC::load(temp);
			} else if(std::strcmp(temp,"PROFILE_CALC")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				if(string::substrN(temp,",")!=3) throw std::invalid_argument("Invalid profile format: must be comma-delimited.");
				profileCalc_.sigma()=std::atof(std::strtok(temp,","));
				profileCalc_.b1()=std::atof(std::strtok(NULL,","));
				profileCalc_.b2()=std::atof(std::strtok(NULL,","));
			} else if(std::strcmp(temp,"PROFILE_LOAD")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				if(string::substrN(temp,",")!=3) throw std::invalid_argument("Invalid profile format: must be comma-delimited.");
				profileLoad_.sigma()=std::atof(std::strtok(temp,","));
				profileLoad_.b1()=std::atof(std::strtok(NULL,","));
				profileLoad_.b2()=std::atof(std::strtok(NULL,","));
			}
		}
		
		//close the parameter file
		if(DEBUG_NLO_ATOM>1) log<<"Closing parameter file...\n";
		fclose(reader);
		reader=NULL;
		
		//check the calculation parameters
		if(DEBUG_NLO_ATOM>1) log<<"Checking parameters...\n";
		if(strideAlpha_==0 || strideChg_==0) throw std::invalid_argument("Invalid stride.");
		if(freqUnit_==fourier::FreqUnit::UNKNOWN) throw std::invalid_argument("Invalid frequency unit.");
	}catch(std::exception& e){
		log<<"ERROR in load(const char*):\n";
		log<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("ERROR in AlphaEff::load(const char*): Failed to load.");
}

void NLO::init(SimAtomic<AtomT>& sim){
	if(DEBUG_NLO_ATOM>0) log<<"init(const SimAtomic<AtomT>&):\n";
	bool error=false;
	
	try{
		//set the fft parameters
		if(DEBUG_NLO_ATOM>0) log<<"Assigning the FFT parameters...\n";
		nStepsChg_=sim.timesteps()/strideChg_;
		nStepsAlpha_=sim.timesteps()/strideAlpha_;
		minFreq_=0.5*(1000.0/sim.timestep())/sim.timesteps();//min freq in THz (1000 for ps, 1/2 due to padding)
		maxFreq_=0.5*(1000.0/sim.timestep());//max freq in THz (1000 for ps, 1/2 due to padding)
		if(freqRes_==0) freqRes_=5.0*minFreq_;
		if(freqCut_==0) freqCut_=maxFreq_;
		
		//set the fourier parameters in the appropriate units
		if(freqUnit_==fourier::FreqUnit::MEV){
			freqCut_*=mevPerThz;
			minFreq_*=mevPerThz;
			maxFreq_*=mevPerThz;
			freqRes_*=mevPerThz;
		} else if(freqUnit_==fourier::FreqUnit::CMI){
			freqCut_*=cmiPerThz;
			minFreq_*=cmiPerThz;
			maxFreq_*=cmiPerThz;
			freqRes_*=cmiPerThz;
		}
		freqN_=(unsigned int)std::ceil((freqCut_+1.0)/minFreq_);
		
		//set the window_
		if(windowType_==window::WINDOW_FUNC::IDENTITY) window_=window::Identity();
		else if(windowType_==window::WINDOW_FUNC::BLACKMANHARRIS) window_=window::BlackmanHarris(sim.timesteps());
		else if(windowType_==window::WINDOW_FUNC::GAUSSIAN){
			sigma_=1/(2*num_const::PI*freqRes_);//from fourier transform
			window_=window::Gaussian(sim.timesteps(),sigma_);
		}
		else if(windowType_==window::WINDOW_FUNC::KAISERBESSEL){
			sigma_=freqRes_*freqRes_;//scale to get sensible values
			window_=window::KaiserBessel(sim.timesteps(),sigma_);
		}
		
		//check the fft parameters
		if(freqCut_-minFreq_<num_const::ZERO) throw std::invalid_argument("Invalid freqency cutoff.");
		
		//allocate space for chi2
		if(DEBUG_NLO_ATOM>0) log<<"Allocating space for chi2...\n";
		chi2_.resize(3);
		for(unsigned int i=0; i<3; ++i){
			chi2_[i].resize(3);
			for(unsigned int j=0; j<3; ++j){
				chi2_[i][j].resize(3);
				for(unsigned int k=0; k<3; ++k){
					chi2_[i][j][k].resize(freqN_);
				}
			}
		}
		
	}catch(std::exception& e){
		log<<"ERROR in load(const char*):\n";
		log<<e.what()<<"\n";
		error=true;
	}
}

void NLO::calcAlpha(SimAtomic<AtomT>& sim, const Thole& thole, const Ewald3D::Dipole& ewald){
	if(DEBUG_NLO_ATOM>0) log<<"calcAlpha(SimAtomic<AtomT>&,const Thole&,const Ewald3D::Dipole&):\n";
	//local function variables
	//parallel
		unsigned int nThreads=1;
		#ifdef _OPENMP
			nThreads=omp_get_max_threads();
		#endif
	//timing
		clock_t start,stop;
		double time;
	//atomic charge
		std::vector<Thole> thole_(nThreads);
	
	//initialize the effective alphas
	if(DEBUG_NLO_ATOM>1) log<<"Initializing Thole objects...\n";
	for(unsigned int i=0; i<nThreads; ++i) thole_[i]=thole;
	
	//calculate the atomic alphas
	if(DEBUG_NLO_ATOM>1) log<<"Calculating Atomic Alphas...\n";
	start=std::clock();
	#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
	for(unsigned int t=0; t<sim.timesteps(); t+=strideAlpha_){
		unsigned int TN=0;
		#ifdef _OPENMP
			TN=omp_get_thread_num();
		#endif
		if(DEBUG_NLO_ATOM>-1) log<<"Timestep: "<<sim.beg()+1+t<<"\n";
		else if(DEBUG_NLO_ATOM>0 && t%1000==0) log<<"Timestep: "<<sim.beg()+1+t<<"\n";
		thole_[TN].alpha_cont(sim,t,ewald);
	}
	stop=std::clock();
	time=((double)(stop-start))/CLOCKS_PER_SEC;
	log<<"time_alpha(s) = "<<time<<"\n";
	
	//interpolate the polarizabilities
	if(strideAlpha_>1){
		if(DEBUG_NLO_ATOM>0) log<<"Interpolating Effective Alphas...\n";
		std::vector<std::vector<Interp::Data> > alphaData;
		alphaData.resize(3);
		for(unsigned int i=0; i<3; ++i){
			alphaData[i].resize(3);
			for(unsigned int j=0; j<3; ++j){
				alphaData[i][j].resize(nStepsAlpha_+1);
			}
		}
		//loop over all atoms included in the calculation
		for(unsigned int n=0; n<sim.nAtoms(); ++n){
			//loop over all timesteps for which we calculated effective polarizabilities
			for(unsigned int t=0; t*strideAlpha_<sim.timesteps(); ++t){
				for(unsigned int i=0; i<3; ++i){
					for(unsigned int j=0; j<3; ++j){
						alphaData[i][j].x(t)=t*strideAlpha_;
						alphaData[i][j].y(t)=sim.atom(t*strideAlpha_,n).alpha()(i,j);
					}
				}
			}
			//interpolate the rest of the timesteps
			for(unsigned int t=0; t*strideAlpha_<sim.timesteps()-strideAlpha_; ++t){
				for(unsigned int tp=1; tp<strideAlpha_; ++tp){
					for(unsigned int i=0; i<3; ++i){
						for(unsigned int j=0; j<3; ++j){
							sim.atom(t*strideAlpha_+tp,n).alpha()(i,j)=Interp::interpAkima(t*strideAlpha_+tp,alphaData[i][j]);
						}
					}
				}
			}
		}
		//fix the final timesteps
		double m;
		//int end=this->fourierN*stride;
		unsigned int end=nStepsAlpha_*strideAlpha_;
		//loop over all molecules
		for(unsigned int n=0; n<sim.nAtoms(); ++n){
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					//calculate the slope at the end
					m=1.0/sim.timestep()*(sim.atom(end,n).alpha()(i,j)-sim.atom(end-1,n).alpha()(i,j));
					//perform linear extrapolation for the last few timesteps
					for(unsigned int t=end; t<sim.timesteps(); ++t){
						sim.atom(t,n).alpha()(i,j)=sim.atom(end,n).alpha()(i,j)+m*(t-end)*sim.timestep();
					}
				}
			}
		}
	}
}

void NLO::calcChg(SimAtomic<AtomT>& sim, const QEQ3& qeq, const Ewald3D::Coulomb& ewald){
	if(DEBUG_NLO_ATOM>0) log<<"calcChg(SimAtomic<AtomT>&,const QEQ3&,const Ewald3D::Coulomb&):\n";
	//local function variables
	//parallel
		unsigned int nThreads=1;
		#ifdef _OPENMP
			nThreads=omp_get_max_threads();
		#endif
	//timing
		clock_t start,stop;
		double time;
	//atomic charge
		std::vector<QEQ3> qeq_(nThreads);
	
	//initialize the effective alphas
	if(DEBUG_NLO_ATOM>1) log<<"Initializing QEQ Objects...\n";
	for(unsigned int i=0; i<nThreads; ++i) qeq_[i]=qeq;
	
	//calculate the atomic charges
	if(DEBUG_NLO_ATOM>1) log<<"Calculating Atomic Charges...\n";
	start=std::clock();
	#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
	for(unsigned int t=0; t<sim.timesteps(); t+=strideChg_){
		unsigned int TN=0;
		#ifdef _OPENMP
			TN=omp_get_thread_num();
		#endif
		if(DEBUG_NLO_ATOM>-1) log<<"Timestep: "<<sim.beg()+1+t<<"\n";
		else if(DEBUG_NLO_ATOM>0 && t%1000==0) log<<"Timestep: "<<sim.beg()+1+t<<"\n";
		qeq_[TN].qt_jZero_cont(t,sim,ewald);
	}
	stop=std::clock();
	time=((double)(stop-start))/CLOCKS_PER_SEC;
	log<<"time_chg(s) = "<<time<<"\n";
	
	if(strideChg_>1){
		if(DEBUG_NLO_ATOM>0) log<<"Interpolating atomic charge...\n";
		Interp::Data chgData(nStepsChg_+1);
		//loop over all atoms
		for(unsigned int n=0; n<sim.nAtoms(); n++){
			//loop over all timesteps for which we calculated atomic charges
			for(int t=0; t*strideChg_<sim.timesteps(); t++){
				chgData.x(t)=t*strideChg_;
				chgData.y(t)=sim.atom(t*strideChg_,n).charge();
			}
			//interpolate the rest of the timesteps
			for(unsigned int t=0; t*strideChg_<sim.timesteps()-strideChg_; ++t){
				for(unsigned int tp=1; tp<strideChg_; ++tp){
					sim.atom(t*strideChg_+tp,n).charge()=Interp::interpAkima(t*strideChg_+tp,chgData);
				}
			}
		}
		//fix the final timesteps
		double m;
		unsigned int end=nStepsChg_*strideChg_;
		//loop over all molecules
		for(unsigned int n=0; n<sim.nAtoms(); ++n){
			//calculate the slope at the end
			m=1.0/sim.timestep()*(sim.atom(end,n).charge()-sim.atom(end-1,n).charge());
			//perform linear extrapolation for the last few timesteps
			for(unsigned int t=end; t<sim.timesteps(); ++t){
				sim.atom(t,n).charge()=sim.atom(end,n).charge()+m*(t-end)*sim.timestep();
			}
		}
	}
}

void NLO::calcSpectrum(SimAtomic<AtomT>& sim, Bonding& bonding){
	if(DEBUG_NLO_ATOM>0) log<<"calcSpectrum(SimAtomic<AtomT>&):\n";
	//local function variables
	fourier::FFT_R2C fftMu[3];//the FFT of the dipole data
	fourier::FFT_R2C fftAlpha[3][3];//the FFT of the polarizability data
	fourier::FFT_C2C fft[3][3][3];//the FFT of chi2
	for(unsigned int i=0; i<3; ++i){
		fftMu[i]=fourier::FFT_R2C(2*sim.timesteps());
		for(unsigned int j=0; j<3; ++j){
			fftAlpha[i][j]=fourier::FFT_R2C(2*sim.timesteps());
			for(unsigned int k=0; k<3; ++k){
				fft[i][j][k]=fourier::FFT_C2C(2*sim.timesteps());
			}
		}
	}
	//transpose of alpha delta
	std::vector<std::vector<std::vector<double> > > alphaTrans(3);
	for(unsigned int i=0; i<3; ++i) alphaTrans[i].resize(3,std::vector<double>(sim.timesteps(),0));
	std::vector<std::vector<double> > dipoleTrans(3,std::vector<double>(sim.timesteps(),0));
	//molecular simulation for calculating bond dipoles
	SimMol<AtomT,MolT> simMol;
	static_cast<SimAtomic<AtomT>&>(simMol).resize(1,sim.nAtomsVec(),sim.atomNames());
	std::vector<std::vector<unsigned int> > bOrder(sim.nSpecies());
	for(unsigned int n=0; n<sim.nSpecies(); ++n) bOrder[n].resize(sim.nAtoms(n));
	if(DEBUG_NLO_ATOM>-1) log<<"SimMol = "<<simMol<<"\n";
	
		//No truncation of correlation function
		if(DEBUG_NLO_ATOM>0) log<<"No truncation of correlation function.\n";
		
		//calculate the transpose of the total dipole moment
		if(DEBUG_NLO_ATOM>0) log<<"Calculating the transpose of the total dipole...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			//local variables
			Eigen::Vector3d muTot=Eigen::Vector3d::Zero();
			//copy over the atoms
			if(DEBUG_NLO_ATOM>1) log<<"Copying over atoms...\n";
			for(unsigned int n=0; n<sim.nSpecies(); ++n){
				for(unsigned int m=0; m<sim.nAtoms(n); ++m){
					simMol.atom(0,n,m).posn().noalias()=sim.atom(t,n,m).posn();
					simMol.atom(0,n,m).charge()=sim.atom(t,n,m).charge();
				}
			}
			simMol.cell(0)=sim.cell(t);
			//calculate bonds
			if(DEBUG_NLO_ATOM>1) log<<"Calculating bonds...\n";
			sim_util::assignBonds(simMol,static_cast<SimAtomic<AtomT>&>(simMol),bonding);
			if(DEBUG_NLO_ATOM>1) log<<"NMol = "<<simMol.nMol()<<"\n";
			//find the bond order
			if(DEBUG_NLO_ATOM>1) log<<"Finding the bond order...\n";
			for(unsigned int n=0; n<sim.nSpecies(); ++n){
				for(unsigned int m=0; m<sim.nAtoms(n); ++m){
					bOrder[n][m]=Bonding::bondOrder(t,sim.atom(t,n,m),sim,bonding);
				}
			}
			if(DEBUG_NLO_ATOM>1) log<<"Divide charge by bond order...\n";
			for(unsigned int n=0; n<simMol.nSpecies(); ++n){
				for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
					simMol.atom(0,n,m).charge()/=((double)bOrder[n][m]);
				}
			}
			if(DEBUG_NLO_ATOM>1) log<<"Find total molecular charge...\n";
			for(unsigned int n=0; n<simMol.nMol(); ++n){
				simMol.molecule(0,n).charge()=0;
				for(unsigned int m=0; m<simMol.molecule(0,n).nAtoms(); ++m){
					simMol.molecule(0,n).charge()+=simMol.molecule(0,n).atom(m).charge();
				}
			}
			if(DEBUG_NLO_ATOM>1) log<<"Find center of charge...\n";
			for(unsigned int n=0; n<simMol.nMol(); ++n){
				coc(simMol.molecule(0,n),simMol.cell(0),simMol.molecule(0,n).posn());
			}
			if(DEBUG_NLO_ATOM>1) log<<"Find molecular dipole moment...\n";
			for(unsigned int n=0; n<simMol.nMol(); ++n){
				dipole(simMol.molecule(0,n),simMol.cell(0),simMol.molecule(0,n).dipole());
			}
			//find the total dipole moment
			if(profileCalc_.sigma()>0){
				for(unsigned int n=0; n<simMol.nMol(); ++n){
					muTot.noalias()+=simMol.molecule(0,n).dipole()*profileCalc_(simMol.molecule(0,n).posn()[2]);
				}
			} else {
				for(unsigned int n=0; n<simMol.nMol(); ++n){
					muTot.noalias()+=simMol.molecule(0,n).dipole();
				}
			}
			/*muTot.setZero();
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				muTot.noalias()+=sim.atom(t,n).charge()*sim.atom(t,n).posn();
			}*/
			for(unsigned int i=0; i<3; ++i) dipoleTrans[i][t]=muTot[i];
		}
		
		FILE* writer=fopen("dipole_tot.dat","w");
		if(writer!=NULL){
			fprintf(writer,"T X Y Z\n");
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				fprintf(writer,"%i %f %f %f\n",t,
					dipoleTrans[0][t],
					dipoleTrans[1][t],
					dipoleTrans[2][t]
				);
			}
		}
		
		//calculate the velocities of the dipole 
		if(DEBUG_NLO_ATOM>0) log<<"Calculating dipole velocity...\n";
		for(unsigned int i=0; i<3; ++i){
			fftMu[i].in(0)=gradient::df1o1(dipoleTrans[i],sim.timestep(),0);
			fftMu[i].in(1)=gradient::dc1o2(dipoleTrans[i],sim.timestep(),1);
			fftMu[i].in(2)=gradient::dc1o4(dipoleTrans[i],sim.timestep(),2);
			fftMu[i].in(3)=gradient::dc1o6(dipoleTrans[i],sim.timestep(),3);
			for(int t=4; t<sim.timesteps()-4; ++t){
				fftMu[i].in(t)=gradient::dc1o8(dipoleTrans[i],sim.timestep(),t);
			}
			fftMu[i].in(sim.timesteps()-4)=gradient::dc1o6(dipoleTrans[i],sim.timestep(),sim.timesteps()-4);
			fftMu[i].in(sim.timesteps()-3)=gradient::dc1o4(dipoleTrans[i],sim.timestep(),sim.timesteps()-3);
			fftMu[i].in(sim.timesteps()-2)=gradient::dc1o2(dipoleTrans[i],sim.timestep(),sim.timesteps()-2);
			fftMu[i].in(sim.timesteps()-1)=gradient::db1o1(dipoleTrans[i],sim.timestep(),sim.timesteps()-1);
		}
		
		//calculate the Fourier transform
		if(DEBUG_NLO_ATOM>0) log<<"Calculating the forward dipole transform...\n";
		for(unsigned int i=0; i<3; ++i) fftMu[i].transformf();
		
		//calculate the transpose of the polarizability
		if(DEBUG_NLO_ATOM>0) log<<"Calculating the transpose of the total polarizability...\n";
		if(profileCalc_.sigma()>0){
			if(DEBUG_NLO_ATOM>0) log<<"Applying profile...\n";
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					double mean=0.0;
					for(unsigned int i=0; i<3; ++i){
						for(unsigned int j=0; j<3; ++j){
							alphaTrans[i][j][t]+=sim.atom(t,n).alpha()(i,j)*profileCalc_(sim.atom(t,n).posn()[2]);
						}
					}
				}
			}
		} else {
			if(DEBUG_NLO_ATOM>0) log<<"No profile...\n";
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					double mean=0.0;
					for(unsigned int i=0; i<3; ++i){
						for(unsigned int j=0; j<3; ++j){
							alphaTrans[i][j][t]+=sim.atom(t,n).alpha()(i,j);
						}
					}
				}
			}
		}
		
		//calculate the polarizability velocity
		if(DEBUG_NLO_ATOM>0) log<<"Calculating the polarizability velocity...\n";
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				fftAlpha[i][j].in(0)=gradient::df1o1(alphaTrans[i][j],sim.timestep(),0);
				fftAlpha[i][j].in(1)=gradient::dc1o2(alphaTrans[i][j],sim.timestep(),1);
				fftAlpha[i][j].in(2)=gradient::dc1o4(alphaTrans[i][j],sim.timestep(),2);
				fftAlpha[i][j].in(3)=gradient::dc1o6(alphaTrans[i][j],sim.timestep(),3);
				for(unsigned int t=4; t<sim.timesteps()-4; t++){
					fftAlpha[i][j].in(t)=gradient::dc1o8(alphaTrans[i][j],sim.timestep(),t);
				}
				fftAlpha[i][j].in(sim.timesteps()-4)=gradient::dc1o6(alphaTrans[i][j],sim.timestep(),sim.timesteps()-4);
				fftAlpha[i][j].in(sim.timesteps()-3)=gradient::dc1o4(alphaTrans[i][j],sim.timestep(),sim.timesteps()-3);
				fftAlpha[i][j].in(sim.timesteps()-2)=gradient::dc1o2(alphaTrans[i][j],sim.timestep(),sim.timesteps()-2);
				fftAlpha[i][j].in(sim.timesteps()-1)=gradient::db1o1(alphaTrans[i][j],sim.timestep(),sim.timesteps()-1);
			}
		}
		
		//tranform into the frequency domain
		if(DEBUG_NLO_ATOM>0) log<<"Calculating the forward polarizability transform...\n";
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				fftAlpha[i][j].transformf();
			}
		}
		
		//calculate the correlation function in the frequency domain
		if(DEBUG_NLO_ATOM>0) log<<"Calculating the correlation function in frequency space...\n";
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				for(unsigned int k=0; k<3; ++k){
					for(unsigned int t=0; t<2*sim.timesteps(); ++t){
						fft[i][j][k].in(t)[0]=(fftAlpha[i][j].out(t)[0]*fftMu[k].out(t)[0]+fftAlpha[i][j].out(t)[1]*fftMu[k].out(t)[1]);//alpha*.mu
						fft[i][j][k].in(t)[1]=(fftAlpha[i][j].out(t)[0]*fftMu[k].out(t)[1]-fftAlpha[i][j].out(t)[1]*fftMu[k].out(t)[0]);//alpha*.mu
					}
				}
			}
		}
	
	//calculate the corerelation function in time-space, then transform to find chi2
	if(DEBUG_NLO_ATOM>0) log<<"Calculating chi2...\n";
	double norm=1.0;
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int k=0; k<3; ++k){
				//calculate the time-space correlation function
				fft[i][j][k].transformr();
				//record the correlation function, normalized and windowed, as input
				for(unsigned int t=0; t<sim.timesteps(); ++t){
					fft[i][j][k].in(t)[0]=fft[i][j][k].out(t)[0]/(sim.timesteps()-t)*window_(t+1);
					fft[i][j][k].in(t)[1]=fft[i][j][k].out(t)[1]/(sim.timesteps()-t)*window_(t+1);
				}
				for(unsigned int t=sim.timesteps(); t<2*sim.timesteps(); ++t){
					fft[i][j][k].in(t)[0]=fft[i][j][k].out(t)[0]/(t-sim.timesteps()+1)*window_(2*sim.timesteps()-t);
					fft[i][j][k].in(t)[1]=fft[i][j][k].out(t)[1]/(t-sim.timesteps()+1)*window_(2*sim.timesteps()-t);
				}
				//calculate the forward transform
				fft[i][j][k].transformf();
				//record chi2
				for(unsigned int t=0; t<freqN_; ++t){
					chi2_[i][j][k][t][0]=1.0/(minFreq_*(t+1))*norm*fft[i][j][k].out(t)[1]*(-1.0);
					chi2_[i][j][k][t][1]=1.0/(minFreq_*(t+1))*norm*fft[i][j][k].out(t)[0];
				}
			}
		}
	}
}

void NLO::printSpectrum(const char* file) const{
	if(DEBUG_NLO_ATOM>0) log<<"printSpectrum():\n";
	//local function variables
	FILE* writer=NULL;
	
	writer=fopen(file,"w");
	if(writer==NULL) throw std::runtime_error("I/O Exception Occured.");
	
	//print the header
	//print the header
	fprintf(writer, "Freq");
	if(freqUnit_==fourier::FreqUnit::MEV) fprintf(writer, "(meV) ");
	else if(freqUnit_==fourier::FreqUnit::CMI) fprintf(writer, "(cm^-1) ");
	else fprintf(writer, "(THz) ");
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int k=0; k<3; ++k){
				fprintf(writer, "%i:%i:%i ", i,j,k);
			}
		}
	}
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int k=0; k<3; ++k){
				fprintf(writer, "%i:%i:%i ", i,j,k);
			}
		}
	}
	fprintf(writer, "\n");
	
	//print the power spectrum
	for(unsigned int t=0; t<freqN_; ++t){
		fprintf(writer, "%f ", minFreq_*t);
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				for(unsigned int k=0; k<3; ++k){
					fprintf(writer, "%f ", chi2_[i][j][k][t][0]);
				}
			}
		}
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				for(unsigned int k=0; k<3; ++k){
					fprintf(writer, "%f ", chi2_[i][j][k][t][1]);
				}
			}
		}
		fprintf(writer, "\n");
	}
	
	//close the file
	fclose(writer);
	writer=NULL;
}

//***********************************************************************************************************************************
//MAIN
//***********************************************************************************************************************************

int main(int argc, char* argv[]){
	/*Parameters:
		argc==2
		argv[0] = Program Execution
		argv[1] = Parameter file
	*/
	
	//local variables
	//simulation
		int beg,end;
		double zFactor=1;
		Eigen::Vector3d offset=Eigen::Vector3d::Zero();
		SimAtomic<AtomT> sim;
		Bonding bonding;
	//ewald
		Ewald3D::Dipole dipole;
		Ewald3D::Coulomb coul;
	//calculation
		NLO nlo;
		Thole thole;
		QEQ3 qeq;
		std::vector<double> jzero;
		std::vector<double> alpha;
		std::vector<std::string> atoms;
	//input/output
		char* paramFile=(char*)malloc(sizeof(char)*string::M);//parameter file for all required parameters
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
		char* simstr=(char*)malloc(sizeof(char)*string::M);
		FILE* reader=NULL;
		FILE_FORMAT::type fileFormat;
	//flags
		bool error=false;
		bool periodic=true;
	//units
		units::System::type unitsys=units::System::METAL;
		units::consts::init(unitsys);
	//logging
		logging::Log& logger=logging::Log::get();//log object
		logger.addSink(logging::Sink(&std::cout));//add console sink
		logger.addSink(logging::Sink(new std::ofstream("nlo_atom.log")));
		logging::DebugLogger log("nlo_atom");
	
	//begin a do loop for easy breaking
	try{
		/* check the number of arguments */
		if(argc!=2){
			log<<"ERROR in main(int,char**):\n";
			log<<"Incorrect number of arguments.\n";
			log<<"\t1. Program Execution\n";
			log<<"\t2. Parameter File\n";
			throw std::invalid_argument("Invalid command-line arguments.");
		}
		
		/* open the parameter file */
		log<<"Opening the parameter file...\n";
		std::strcpy(paramFile,argv[1]);
		reader=fopen(paramFile,"r");
		if(reader==NULL) throw std::invalid_argument("I/O Error: Could not open parameter file.");
		
		/* read in the parameters */
		log<<"Loading parameters...\n";
		while(std::fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);
			string::trim_all(string::to_upper(string::copy_left(temp,input,"=")));
			if(std::strcmp(temp,"SIMULATION")==0){
				std::strcpy(simstr,string::trim_all(std::strpbrk(input,"=")+1));
			} else if(std::strcmp(temp,"PERIODIC")==0){
				periodic=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"OFFSET")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				if(string::substrN(temp,",")!=3) throw std::invalid_argument("Invalid offset format: must be comma-delimited.");
				offset[0]=std::atof(std::strtok(temp,","));
				offset[1]=std::atof(std::strtok(NULL,","));
				offset[2]=std::atof(std::strtok(NULL,","));
			} else if(std::strcmp(temp,"INTERVAL")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				string::to_upper(string::trim_all(temp));
				if(std::strpbrk(temp,",")==NULL) throw std::invalid_argument("Invalid interval format: must be comma-delimited.");
				std::strcpy(input,std::strtok(temp,","));
				if(std::strcmp(input,"BEG")==0) beg=1;
				else beg=std::atoi(input);
				std::strcpy(input,std::strtok(NULL,","));
				if(std::strcmp(input,"END")==0) end=-1;
				else end=std::atoi(input);
			} else if(std::strcmp(temp,"FILE_FORMAT")==0){
				string::to_upper(string::trim(std::strcpy(temp,std::strpbrk(input,"=")+1)));
				fileFormat=FILE_FORMAT::load(temp);
			} else if(std::strcmp(temp,"UNITS")==0){
				string::to_upper(string::trim(std::strcpy(temp,std::strpbrk(input,"=")+1)));
				unitsys=units::System::load(temp);
			} else if(std::strcmp(temp,"ATOMS")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				if(string::substrN(temp,string::WS)==0) throw std::invalid_argument("Invalid number of atoms.");
				atoms.resize(string::substrN(temp,string::WS));
				atoms[0]=std::string(std::strtok(temp,string::WS));
				for(unsigned int n=1; n<atoms.size(); ++n) atoms[n]=std::string(std::strtok(NULL,string::WS));
			} else if(std::strcmp(temp,"JZERO")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				if(string::substrN(temp,string::WS)==0) throw std::invalid_argument("Invalid number of jzeros.");
				jzero.resize(string::substrN(temp,string::WS),0);
				jzero[0]=std::atof(std::strtok(temp,string::WS));
				for(unsigned int n=1; n<jzero.size(); ++n) jzero[n]=std::atof(std::strtok(NULL,string::WS));
			} else if(std::strcmp(temp,"ALPHA")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				if(string::substrN(temp,string::WS)==0) throw std::invalid_argument("Invalid number of alphas.");
				alpha.resize(string::substrN(temp,string::WS),0);
				alpha[0]=std::atof(std::strtok(temp,string::WS));
				for(unsigned int n=1; n<alpha.size(); ++n) alpha[n]=std::atof(std::strtok(NULL,string::WS));
			} 
		}
		fclose(reader);
		reader=NULL;
		
		//check the parameters
		if(unitsys==units::System::UNKNOWN) throw std::invalid_argument("Invalid unit system");
		
		//initialize the unit system
		std::cout<<"Initializing the unit system...\n";
		units::consts::init(unitsys);
		
		//print the parameters to screen
		log<<"SIMULATION PARAMETERS:\n";
		log<<"\tUNITS    = "<<unitsys<<"\n";
		log<<"\tSIM      = "<<simstr<<"\n";
		log<<"\tFORMAT   = "<<fileFormat<<"\n";
		log<<"\tOFFSET   = ("<<offset[0]<<","<<offset[1]<<","<<offset[2]<<")\n";
		log<<"\tINTERVAL = "<<beg<<", "<<end<<"\n";
		log<<"\tATOMS    = "; for(unsigned int i=0; i<atoms.size(); ++i) log.log()<<atoms[i]<<" "; log.log()<<"\n";
		log<<"\tJZERO    = "; for(unsigned int i=0; i<jzero.size(); ++i) log.log()<<jzero[i]<<" "; log.log()<<"\n";
		log<<"\tALPHA    = "; for(unsigned int i=0; i<alpha.size(); ++i) log.log()<<alpha[i]<<" "; log.log()<<"\n";
		
		//load the simulation
		if(DEBUG_NLO_ATOM>0) std::cout<<"Loading simulation...\n";
		if(fileFormat==FILE_FORMAT::POSCAR){
			VASP::POSCAR::load(simstr,sim);
		} else if(fileFormat==FILE_FORMAT::XDATCAR){
			VASP::XDATCAR::load(simstr,sim,beg,end);
		} else if(fileFormat==FILE_FORMAT::LAMMPS){
			if(string::substrN(simstr,",")!=3) throw std::invalid_argument("Invalid LAMMPS file format.");
			std::vector<std::string> substr(3);
			substr[0]=std::strtok(simstr,",");
			substr[1]=std::strtok(NULL,",");
			substr[2]=std::strtok(NULL,",");
			LAMMPS::load(substr[0].c_str(),substr[1].c_str(),substr[2].c_str(),sim,beg,end);
		} else if(fileFormat==FILE_FORMAT::QE){
			if(string::substrN(simstr,",")!=3) throw std::invalid_argument("Invalid QE file format.");
			std::vector<std::string> substr(3);
			substr[0]=std::strtok(simstr,",");
			substr[1]=std::strtok(NULL,",");
			substr[2]=std::strtok(NULL,",");
			QE::load(substr[0].c_str(),substr[1].c_str(),substr[2].c_str(),sim,beg,end);
		} else throw std::invalid_argument("Invalid file format.");
		//set the periodicity
		sim.periodic()=periodic;
		
		/* Set the offset */
		if(offset.norm()>num_const::ZERO){
			//apply the offset 
			log<<"Applying offset...\n";
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int n=0; n<sim.nSpecies(); ++n){
					for(unsigned int m=0; m<sim.nAtoms(n); ++m){
						sim.atom(t,n,m).posn().noalias()+=offset;
						Cell::returnToCell(sim.atom(t,n,m).posn(),sim.atom(t,n,m).posn(),sim.cell(t).R(),sim.cell(t).RInv());
					}
				}
			}
		}
		
		/* assign jzero */
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				for(unsigned int i=0; i<atoms.size(); ++i){
					if(atoms[i]==sim.atom(t,n).name()){
						sim.atom(t,n).jzero()=jzero[i];
						break;
					}
				}
			}
		}
		
		/* assign alpha */
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				for(unsigned int i=0; i<atoms.size(); ++i){
					if(atoms[i]==sim.atom(t,n).name()){
						sim.atom(t,n).alpha().noalias()=Eigen::Matrix3d::Identity()*alpha[i];
						break;
					}
				}
			}
		}
		
		//print the simulation
		log<<"SIM = "<<sim<<"\n";
		
		/* Create the auxiliary objects */
		log<<"Creating the auxiliary objects...\n";
		//bonding object
			log<<"\tBonding object...\n";
			bonding.loadBondParams(paramFile,sim.atomNames());
			log<<"\tBonding = \n"<<bonding<<"\n";
		//ewald object
			log<<"\tEwald objects...\n";
			dipole.init(sim.cell(0),1e-5);
			log<<"Ewald::Dipole = \n"<<dipole<<"\n";
			coul.init(sim,1e-5);
			log<<"Ewald::Coulomb = \n"<<coul<<"\n";
		//qeq object
			log<<"\tQEQ object...\n";
			qeq.load(paramFile);
			qeq.init(sim);
			log<<"QEQ = \n"<<qeq<<"\n";
		//thole object
			log<<"\tThole object...\n";
			Thole::load(paramFile,thole);
			thole.init(sim,dipole);
			log<<"Thole = \n"<<thole<<"\n";
		//nlo object
			log<<"\tNLO object...\n";
			nlo.load(paramFile);
			nlo.init(sim);
			log<<"NLO = \n"<<nlo<<"\n";
		
		/* Perform the calculation */
		log<<"Performing the calculation...\n";
		if(nlo.calcAlpha()){
			log<<"Calculating atomic alphas...\n";
			nlo.calcAlpha(sim,thole,dipole);
		}
		if(nlo.calcChg()){
			log<<"Calculating atomic charges...\n";
			nlo.calcChg(sim,qeq,coul);
		}
		if(nlo.calcSpectrum()){
			log<<"Calculating the spectrum...\n";
			nlo.calcSpectrum(sim,bonding);
			nlo.printSpectrum(nlo.fileSpectrum().c_str());
		}
		
		/* Print the Relevant Data */
		std::cout<<"printing relevant data...\n";
		if(nlo.printAlphaT()){
			std::cout<<"alpha-tot\n";
			Electrostatics::Print::alpha_tot("alpha_tot.dat",sim);
		}
		if(nlo.printDipoleT()){
			std::cout<<"dipole-tot\n";
			Electrostatics::Print::dipole_tot("dipole_tot.dat",sim,bonding);
		}
		if(nlo.printChgT()){
			std::cout<<"chg-tot\n";
			Electrostatics::Print::chg_tot("chg_tot.dat",sim);
		}
		if(nlo.printChg()){
			std::cout<<"chg\n";
			Electrostatics::Print::chg("chg.dat",sim);
		}
	}catch(std::exception& e){
		log<<"ERROR in main(int,char**):\n";
		log<<e.what()<<"\n";
		error=true;
	}
	
	//free all local variables
	log<<"Freeing local variables...\n";
	free(input);
	free(temp);
	free(paramFile);
	free(simstr);
	if(reader!=NULL) fclose(reader);
	
	log<<"Exiting the program...\n";
	if(!error) return 0;
	else return 1;
}
