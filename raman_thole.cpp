#include "raman_thole.hpp"

//***********************************************************************************************************************************
//Profile
//***********************************************************************************************************************************

std::ostream& operator<<(std::ostream& out, const Profile& p){
	return out<<"profile("<<p.sigma()<<","<<p.b1()<<","<<p.b2()<<")";
}

//***********************************************************************************************************************************
//Raman3D class
//***********************************************************************************************************************************

//constants
const double Raman3D::mevPerThz=0.24180;
const double Raman3D::cmiPerThz=0.02998;

std::ostream& operator<<(std::ostream& out, const Raman3D& raman3d){
	out<<"*************************************************\n";
	out<<"****************** RAMAN3D_THOLEV ******************\n";
	out<<"FREQ_UNIT       = "<<raman3d.freqUnit_<<"\n";
	out<<"MIN_FREQ        = "<<raman3d.minFreq_<<"\n";
	out<<"MAX_FREQ        = "<<raman3d.maxFreq_<<"\n";
	out<<"FREQ_CUT        = "<<raman3d.freqCut_<<"\n";
	out<<"FREQ_RES        = "<<raman3d.freqRes_<<"\n";
	out<<"FREQ_VIS        = "<<raman3d.freqVis_<<"\n";
	out<<"STRIDE_ALPHA    = "<<raman3d.strideAlpha_<<"\n";
	out<<"WINDOW          = "<<raman3d.windowType_<<"\n";
	out<<"TEMP            = "<<raman3d.T_<<"\n";
	out<<"FILE_SPECTRUM   = "<<raman3d.fileSpectrum_<<"\n";
	out<<"CALC_ALPHA      = "<<raman3d.calcAlpha_<<"\n";
	out<<"CALC_SPECTRUM   = "<<raman3d.calcSpectrum_<<"\n";
	out<<"PRINT_ALPHA_TOT = "<<raman3d.printAlphaT_<<"\n";
	out<<"NORMALIZE       = "<<raman3d.normalize_<<"\n";
	out<<"PROFILE_CALC    = "<<raman3d.profileCalc_<<"\n";
	out<<"PROFILE_LOAD    = "<<raman3d.profileLoad_<<"\n";
	out<<"****************** RAMAN3D_THOLEV ******************\n";
	out<<"*************************************************";
	return out;
}

void Raman3D::defaults(){
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"defaults()\n";
	//fft parameters
		freqUnit_=fourier::FreqUnit::THZ;
		freqCut_=0;
		freqN_=0;
		minFreq_=0;
		maxFreq_=0;
		freqRes_=0;
		fourierN_=0;
		freqVis_=600;//THz
	//window
		sigma_=0;
		window_=window::Identity();
		windowType_=window::WINDOW_FUNC::IDENTITY;
	//stride
		strideAlpha_=1;//include every timestep
		nStepsAlpha_=0;
	//i/o
		fileSpectrum_=std::string("raman3d.dat");
	//calculation flags
		calcAlpha_=true;
		calcSpectrum_=true;
		normalize_=true;
	//i/o flags
		printAlphaT_=false;
	//temperature
		T_=300.0;
	//raman spectrum
		ramanp_.clear();
		ramans_.clear();
		ramant_.clear();
}

void Raman3D::load(const char* file){
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"load(const char*):\n";
	//local function variables
	FILE* reader=NULL;
	char* input=(char*)malloc(sizeof(char)*string::M);
	char* temp=(char*)malloc(sizeof(char)*string::M);
	bool error=false;
	
	try{
		//open the parameter file
		if(DEBUG_RAMAN3D_THOLEV>1) log<<"Opening parameter file...\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open parameter file.");
		
		//set the defaults
		if(DEBUG_RAMAN3D_THOLEV>1) log<<"Setting defaults...\n";
		defaults();
		
		//read in the parameters
		if(DEBUG_RAMAN3D_THOLEV>1) log<<"Reading in parameters...\n";
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_all(string::trim_right(input,string::COMMENT));
			string::to_upper(string::copy_left(temp,input,"="));
			if(std::strcmp(temp,"FILE_SPECTRUM")==0){
				fileSpectrum_=std::string(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"CALC_SPECTRUM")==0){
				calcSpectrum_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"CALC_ALPHA")==0){
				calcAlpha_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"PRINT_ALPHA_TOT")==0){
				printAlphaT_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"NORMALIZE")==0){
				normalize_=string::boolean(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"FREQ_CUT")==0){
				freqCut_=std::atof(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"FREQ_UNIT")==0){
				string::to_upper(std::strcpy(temp,std::strpbrk(input,"=")+1));
				freqUnit_=fourier::FreqUnit::load(temp);
			} else if(std::strcmp(temp,"RESOLUTION")==0){
				freqRes_=std::atof(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"FREQ_VIS")==0){
				freqVis_=std::atof(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"T")==0){
				T_=std::atof(std::strpbrk(input,"=")+1);
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
		if(DEBUG_RAMAN3D_THOLEV>1) log<<"Closing parameter file...\n";
		fclose(reader);
		reader=NULL;
		
		//check the calculation parameters
		if(DEBUG_RAMAN3D_THOLEV>1) log<<"Checking parameters...\n";
		if(strideAlpha_==0) throw std::invalid_argument("Invalid stride.");
		if(freqUnit_==fourier::FreqUnit::UNKNOWN) throw std::invalid_argument("Invalid frequency unit.");
		if(freqVis_<=0) throw std::invalid_argument("Invalid visible frequency.");
		if(T_<=0) throw std::invalid_argument("Invalid temperature.");
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

void Raman3D::init(SimAtomic<AtomT>& sim){
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"init(const SimAtomic<AtomT>&):\n";
	bool error=false;
	
	try{
		//set the fft parameters
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Assigning the FFT parameters...\n";
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
			freqVis_*=mevPerThz;
		} else if(freqUnit_==fourier::FreqUnit::CMI){
			freqCut_*=cmiPerThz;
			minFreq_*=cmiPerThz;
			maxFreq_*=cmiPerThz;
			freqRes_*=cmiPerThz;
			freqVis_*=cmiPerThz;
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
		
		//allocate space for raman3d
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Allocating space for ir spectrum...\n";
		ramanp_.resize(freqN_);
		ramans_.resize(freqN_);
		ramant_.resize(freqN_);
		
	}catch(std::exception& e){
		log<<"ERROR in load(const char*):\n";
		log<<e.what()<<"\n";
		error=true;
	}
}

void Raman3D::calcAlpha(SimAtomic<AtomT>& sim, const Thole& thole, const Ewald3D::Dipole& ewald){
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"calcAlpha(SimAtomic<AtomT>&,const Thole&,const Ewald3D::Dipole&):\n";
	//local function variables
	//parallel
		unsigned int nThreads=1;
		#ifdef _OPENMP
			nThreads=omp_get_max_threads();
		#endif
	//timing
		std::chrono::high_resolution_clock::time_point start;
		std::chrono::high_resolution_clock::time_point stop;
		std::chrono::duration<double> time;
	//atomic charge
		std::vector<Thole> thole_(nThreads);
	
	//initialize the effective alphas
	if(DEBUG_RAMAN3D_THOLEV>1) log<<"Initializing Thole objects...\n";
	for(unsigned int i=0; i<nThreads; ++i) thole_[i]=thole;
	
	//calculate the atomic alphas
	if(DEBUG_RAMAN3D_THOLEV>1) log<<"Calculating Atomic Alphas...\n";
	if(!sim.cellFixed()){
		std::vector<Ewald3D::Dipole> ewald_(nThreads);
		for(unsigned int i=0; i<nThreads; ++i) ewald_[i]=ewald;
		start=std::chrono::high_resolution_clock::now();
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideAlpha_){
			unsigned int TN=0;
			#ifdef _OPENMP
				TN=omp_get_thread_num();
			#endif
			if(DEBUG_RAMAN3D_THOLEV>-1) log<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_RAMAN3D_THOLEV>0 && t%1000==0) log<<"Timestep: "<<sim.beg()+1+t<<"\n";
			ewald_[TN].init(sim.cell(t),1e-5);
			thole_[TN].alpha(sim,t,ewald);
		}
		stop=std::chrono::high_resolution_clock::now();
		time=std::chrono::duration_cast<std::chrono::duration<double> >(stop-start);
		log<<"alpha time = "<<time.count()<<"\n";
	} else {
		start=std::chrono::high_resolution_clock::now();
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideAlpha_){
			unsigned int TN=0;
			#ifdef _OPENMP
				TN=omp_get_thread_num();
			#endif
			if(DEBUG_RAMAN3D_THOLEV>-1) log<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_RAMAN3D_THOLEV>0 && t%1000==0) log<<"Timestep: "<<sim.beg()+1+t<<"\n";
			thole_[TN].alpha_cont(sim,t,ewald);
		}
		stop=std::chrono::high_resolution_clock::now();
		time=std::chrono::duration_cast<std::chrono::duration<double> >(stop-start);
		log<<"alpha time = "<<time.count()<<"\n";
	}
	
	//interpolate the polarizabilities
	if(strideAlpha_>1){
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Interpolating Effective Alphas...\n";
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

void Raman3D::calcSpectrum(SimAtomic<AtomT>& sim){
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"calcSpectrum(SimAtomic<AtomT>&):\n";
	//local function variables
	fourier::FFT_R2C fftAlpha[3][3];//the FFT of the polarizability data
	fourier::FFT_C2C fft=fourier::FFT_C2C(2*sim.timesteps());//the FFT of chi2
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			fftAlpha[i][j]=fourier::FFT_R2C(2*sim.timesteps());
		}
	}
	//transpose of alpha delta
	std::vector<std::vector<std::vector<double> > > deltaTrans(3);
	for(unsigned int i=0; i<3; ++i){
		deltaTrans[i].resize(3,std::vector<double>(sim.timesteps(),0));
	}
	
		//No truncation of correlation function
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"No truncation of correlation function.\n";
		
		//calculate the transpose of the anisotropy
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating the transpose of the total anisotropy...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				double mean=0;
				for(unsigned int i=0; i<3; ++i){
					mean+=sim.atom(t,n).alpha()(i,i);
					for(unsigned int j=0; j<3; ++j){
						deltaTrans[i][j][t]+=sim.atom(t,n).alpha()(i,j);
					}
				}
				for(unsigned int i=0; i<3; ++i) deltaTrans[i][i][t]-=1.0/3.0*mean;
			}
		}
		
		//calculate the anisotropy velocity
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating the anisotropy velocity...\n";
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				fftAlpha[i][j].in(0)=gradient::df1o1(deltaTrans[i][j],sim.timestep(),0);
				fftAlpha[i][j].in(1)=gradient::dc1o2(deltaTrans[i][j],sim.timestep(),1);
				fftAlpha[i][j].in(2)=gradient::dc1o4(deltaTrans[i][j],sim.timestep(),2);
				fftAlpha[i][j].in(3)=gradient::dc1o6(deltaTrans[i][j],sim.timestep(),3);
				for(int t=4; t<sim.timesteps()-4; t++){
					fftAlpha[i][j].in(t)=gradient::dc1o8(deltaTrans[i][j],sim.timestep(),t);
				}
				fftAlpha[i][j].in(sim.timesteps()-4)=gradient::dc1o6(deltaTrans[i][j],sim.timestep(),sim.timesteps()-4);
				fftAlpha[i][j].in(sim.timesteps()-3)=gradient::dc1o4(deltaTrans[i][j],sim.timestep(),sim.timesteps()-3);
				fftAlpha[i][j].in(sim.timesteps()-2)=gradient::dc1o2(deltaTrans[i][j],sim.timestep(),sim.timesteps()-2);
				fftAlpha[i][j].in(sim.timesteps()-1)=gradient::db1o1(deltaTrans[i][j],sim.timestep(),sim.timesteps()-1);
			}
		}
		
		//tranform into the frequency domain
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating the forward Fourier transform...\n";
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				fftAlpha[i][j].transformf();
			}
		}
		
		//record the correlation function in the frequency domain
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating the correlation function in the frequency domain...\n";
		for(unsigned int t=0; t<2*sim.timesteps(); ++t){
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					fft.in(t)[0]+=fftAlpha[i][j].out(t)[0]*fftAlpha[i][j].out(t)[0]+fftAlpha[i][j].out(t)[1]*fftAlpha[i][j].out(t)[1];
				}
			}
		}
		
		//tranform into the frequency domain
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating the forward Fourier transform...\n";
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				fftAlpha[i][j].transformf();
			}
		}
		
		//record the correlation function in the frequency domain
		if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating the correlation function in the frequency domain...\n";
		for(unsigned int t=0; t<2*sim.timesteps(); ++t){
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					fft.in(t)[0]+=fftAlpha[i][j].out(t)[0]*fftAlpha[i][j].out(t)[0]+fftAlpha[i][j].out(t)[1]*fftAlpha[i][j].out(t)[1];
				}
			}
		}
	
	//transform the correlation function back into the time domain
	fft.transformr();
	
	//normalize and window the correlation function
	for(unsigned int t=0; t<sim.timesteps(); ++t){
		fft.in(t)[0]=fft.out(t)[0]*window_(t+1)/(sim.timesteps()-t);
		fft.in(t)[1]=fft.out(t)[1]*window_(t+1)/(sim.timesteps()-t);
	}
	for(unsigned int t=sim.timesteps(); t<2*sim.timesteps(); ++t){
		fft.in(t)[0]=fft.out(t)[0]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
		fft.in(t)[1]=fft.out(t)[1]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
	}
	
	//trasform to find the windowed, normalized correlation function in the frequency domain
	fft.transformf();
	
	//record the raman spectrum
	double norm=1.0/((2*sim.timesteps())*(2*sim.timesteps()));
	for(unsigned int t=0; t<freqN_; ++t){
		ramanp_[t]=std::sqrt(fft.out(t)[0]*fft.out(t)[0]+fft.out(t)[1]*fft.out(t)[1]);
	}
}

void Raman3D::calcSpectrum2(SimAtomic<AtomT>& sim){
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"calcSpectrum2(SimAtomic<AtomT>&):\n";
	//local function variables
	fourier::FFT_R2C fftAlpha1=fourier::FFT_R2C(2*sim.timesteps());//the FFT of the polarizability data
	fourier::FFT_R2C fftAlpha2=fourier::FFT_R2C(2*sim.timesteps());//the FFT of the polarizability data
	fourier::FFT_C2C fft=fourier::FFT_C2C(2*sim.timesteps());//the FFT of chi2
	//diagonal correlation functions
	std::vector<std::vector<std::vector<double> > > corrDiag(3);
	for(unsigned int i=0; i<3; ++i){
		corrDiag[i].resize(3,std::vector<double>(2*sim.timesteps(),0));
	}
	//off-diagonal correlation functions
	std::vector<std::vector<std::vector<double> > > corrODiag(3);
	for(unsigned int i=0; i<3; ++i){
		corrODiag[i].resize(3,std::vector<double>(2*sim.timesteps(),0));
	}
	//transpose of polarizability
	std::vector<std::vector<std::vector<double> > > alphaTrans(3);
	for(unsigned int i=0; i<3; ++i){
		alphaTrans[i].resize(3,std::vector<double>(sim.timesteps(),0));
	}
	
	//No truncation of correlation function
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"No truncation of correlation function.\n";
	
	//calculate the transpose of the total polarizability
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating the transpose of the total polarizability...\n";
	for(unsigned int t=0; t<sim.timesteps(); ++t){
		Eigen::Matrix3d alpha=Eigen::Matrix3d::Zero();
		for(unsigned int n=0; n<sim.nAtoms(); ++n) alpha.noalias()+=sim.atom(t,n).alpha();
		alpha/=sim.nAtoms();
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				alphaTrans[i][j][t]=alpha(i,j);
			}
		}
	}
	
	//calculate the diagonal correlation functions - diagonal
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating correlation functions - diagonal...\n";
	for(unsigned int i=0; i<3; ++i){
		//calculate velocity
		fftAlpha1.in(0)=gradient::df1o1(alphaTrans[i][i],sim.timestep(),0);
		fftAlpha1.in(1)=gradient::dc1o2(alphaTrans[i][i],sim.timestep(),1);
		fftAlpha1.in(2)=gradient::dc1o4(alphaTrans[i][i],sim.timestep(),2);
		fftAlpha1.in(3)=gradient::dc1o6(alphaTrans[i][i],sim.timestep(),3);
		for(int t=4; t<sim.timesteps()-4; t++){
			fftAlpha1.in(t)=gradient::dc1o8(alphaTrans[i][i],sim.timestep(),t);
		}
		fftAlpha1.in(sim.timesteps()-4)=gradient::dc1o6(alphaTrans[i][i],sim.timestep(),sim.timesteps()-4);
		fftAlpha1.in(sim.timesteps()-3)=gradient::dc1o4(alphaTrans[i][i],sim.timestep(),sim.timesteps()-3);
		fftAlpha1.in(sim.timesteps()-2)=gradient::dc1o2(alphaTrans[i][i],sim.timestep(),sim.timesteps()-2);
		fftAlpha1.in(sim.timesteps()-1)=gradient::db1o1(alphaTrans[i][i],sim.timestep(),sim.timesteps()-1);
		//tranform into the frequency domain
		fftAlpha1.transformf();
		for(unsigned int j=0; j<3; ++j){
			//calculate velocity
			fftAlpha2.in(0)=gradient::df1o1(alphaTrans[j][j],sim.timestep(),0);
			fftAlpha2.in(1)=gradient::dc1o2(alphaTrans[j][j],sim.timestep(),1);
			fftAlpha2.in(2)=gradient::dc1o4(alphaTrans[j][j],sim.timestep(),2);
			fftAlpha2.in(3)=gradient::dc1o6(alphaTrans[j][j],sim.timestep(),3);
			for(int t=4; t<sim.timesteps()-4; t++){
				fftAlpha2.in(t)=gradient::dc1o8(alphaTrans[j][j],sim.timestep(),t);
			}
			fftAlpha2.in(sim.timesteps()-4)=gradient::dc1o6(alphaTrans[j][j],sim.timestep(),sim.timesteps()-4);
			fftAlpha2.in(sim.timesteps()-3)=gradient::dc1o4(alphaTrans[j][j],sim.timestep(),sim.timesteps()-3);
			fftAlpha2.in(sim.timesteps()-2)=gradient::dc1o2(alphaTrans[j][j],sim.timestep(),sim.timesteps()-2);
			fftAlpha2.in(sim.timesteps()-1)=gradient::db1o1(alphaTrans[j][j],sim.timestep(),sim.timesteps()-1);
			//tranform into the frequency domain
			fftAlpha2.transformf();
			//record the correlation function in the frequency domain
			for(unsigned int t=0; t<2*sim.timesteps(); ++t){
				fft.in(t)[0]=fftAlpha1.out(t)[0]*fftAlpha2.out(t)[0]+fftAlpha1.out(t)[1]*fftAlpha2.out(t)[1];//real 
				fft.in(t)[1]=fftAlpha1.out(t)[0]*fftAlpha2.out(t)[1]-fftAlpha1.out(t)[1]*fftAlpha2.out(t)[0];//imag
			}
			//transform the correlation function back into the time domain
			fft.transformr();
			//normalize and window the correlation function
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(t+1)/(sim.timesteps()-t);
				fft.in(t)[1]=fft.out(t)[1]*window_(t+1)/(sim.timesteps()-t);
			}
			for(unsigned int t=sim.timesteps(); t<2*sim.timesteps(); ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
				fft.in(t)[1]=fft.out(t)[1]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
			}
			//transform to find the windowed, normalized correlation function in the frequency domain
			fft.transformf();
			//record the correlation function in the frequency domain
			double norm=1;
			if(normalize_) norm=2.0*sim.timesteps();
			for(unsigned int t=0; t<2*sim.timesteps(); ++t){
				corrDiag[i][j][t]=(fft.out(t)[0]*fft.out(t)[0]+fft.out(t)[1]*fft.out(t)[1])/norm;
			}
		}
	}
	
	//calculate the diagonal correlation functions - off-diagonal
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating correlation functions - off-diagonal...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			//calculate velocity
			fftAlpha1.in(0)=gradient::df1o1(alphaTrans[i][j],sim.timestep(),0);
			fftAlpha1.in(1)=gradient::dc1o2(alphaTrans[i][j],sim.timestep(),1);
			fftAlpha1.in(2)=gradient::dc1o4(alphaTrans[i][j],sim.timestep(),2);
			fftAlpha1.in(3)=gradient::dc1o6(alphaTrans[i][j],sim.timestep(),3);
			for(int t=4; t<sim.timesteps()-4; t++){
				fftAlpha1.in(t)=gradient::dc1o8(alphaTrans[i][j],sim.timestep(),t);
			}
			fftAlpha1.in(sim.timesteps()-4)=gradient::dc1o6(alphaTrans[i][j],sim.timestep(),sim.timesteps()-4);
			fftAlpha1.in(sim.timesteps()-3)=gradient::dc1o4(alphaTrans[i][j],sim.timestep(),sim.timesteps()-3);
			fftAlpha1.in(sim.timesteps()-2)=gradient::dc1o2(alphaTrans[i][j],sim.timestep(),sim.timesteps()-2);
			fftAlpha1.in(sim.timesteps()-1)=gradient::db1o1(alphaTrans[i][j],sim.timestep(),sim.timesteps()-1);
			//tranform into the frequency domain
			fftAlpha1.transformf();
			//calculate velocity
			fftAlpha2.in(0)=gradient::df1o1(alphaTrans[i][j],sim.timestep(),0);
			fftAlpha2.in(1)=gradient::dc1o2(alphaTrans[i][j],sim.timestep(),1);
			fftAlpha2.in(2)=gradient::dc1o4(alphaTrans[i][j],sim.timestep(),2);
			fftAlpha2.in(3)=gradient::dc1o6(alphaTrans[i][j],sim.timestep(),3);
			for(int t=4; t<sim.timesteps()-4; t++){
				fftAlpha2.in(t)=gradient::dc1o8(alphaTrans[i][j],sim.timestep(),t);
			}
			fftAlpha2.in(sim.timesteps()-4)=gradient::dc1o6(alphaTrans[i][j],sim.timestep(),sim.timesteps()-4);
			fftAlpha2.in(sim.timesteps()-3)=gradient::dc1o4(alphaTrans[i][j],sim.timestep(),sim.timesteps()-3);
			fftAlpha2.in(sim.timesteps()-2)=gradient::dc1o2(alphaTrans[i][j],sim.timestep(),sim.timesteps()-2);
			fftAlpha2.in(sim.timesteps()-1)=gradient::db1o1(alphaTrans[i][j],sim.timestep(),sim.timesteps()-1);
			//tranform into the frequency domain
			fftAlpha2.transformf();
			//record the correlation function in the frequency domain 
			for(unsigned int t=0; t<2*sim.timesteps(); ++t){
				fft.in(t)[0]=fftAlpha1.out(t)[0]*fftAlpha2.out(t)[0]+fftAlpha1.out(t)[1]*fftAlpha2.out(t)[1];//real 
				fft.in(t)[1]=fftAlpha1.out(t)[0]*fftAlpha2.out(t)[1]-fftAlpha1.out(t)[1]*fftAlpha2.out(t)[0];//imag
			}
			//transform the correlation function back into the time domain
			fft.transformr();
			//normalize and window the correlation function
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(t+1)/(sim.timesteps()-t);
				fft.in(t)[1]=fft.out(t)[1]*window_(t+1)/(sim.timesteps()-t);
			}
			for(unsigned int t=sim.timesteps(); t<2*sim.timesteps(); ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
				fft.in(t)[1]=fft.out(t)[1]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
			}
			//transform to find the windowed, normalized correlation function in the frequency domain
			fft.transformf();
			//record the correlation function in the frequency domain
			double norm=1;
			if(normalize_) norm=2.0*sim.timesteps();
			for(unsigned int t=0; t<2*sim.timesteps(); ++t){
				corrODiag[i][j][t]=(fft.out(t)[0]*fft.out(t)[0]+fft.out(t)[1]*fft.out(t)[1])/norm;
			}
		}
	}
	
	//reset raman spectra
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Resetting Raman spectra...\n";
	for(unsigned int t=0; t<freqN_; ++t) ramanp_[t]=0.0;
	for(unsigned int t=0; t<freqN_; ++t) ramans_[t]=0.0;
	for(unsigned int t=0; t<freqN_; ++t) ramant_[t]=0.0;
	
	//record the raman spectrum - parallel
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating raman spectrum - parallel...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramanp_[t]+=corrDiag[i][j][t];
			}
		}
	}
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramanp_[t]+=2*corrODiag[i][j][t];
			}
		}
	}
	for(unsigned int t=0; t<freqN_; ++t) ramanp_[t]/=15.0;
	for(unsigned int t=0; t<freqN_; ++t){
		double w=minFreq_*(t+0.5),hbar=units::metal::hbar,kb=units::metal::kb;
		if(freqUnit_==fourier::FreqUnit::MEV) hbar/=mevPerThz;
		else if(freqUnit_==fourier::FreqUnit::CMI) hbar/=cmiPerThz;
		//ramanp_[t]*=std::pow(w/freqVis_-1.0,4)*hbar*w/(w*w)/(kb*T_)/(1.0-std::exp(-hbar*w/(kb*T_)));
		ramanp_[t]*=std::tanh(hbar*w/(kb*T_))/(w*std::pow(w/freqVis_-1.0,4));
	}
	
	//record the raman spectrum - perpendicular
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating raman spectrum - perpendicular...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramans_[t]+=3*corrODiag[i][j][t];
			}
		}
	}
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramans_[t]-=corrDiag[i][j][t];
			}
		}
	}
	for(unsigned int t=0; t<freqN_; ++t) ramans_[t]/=30.0;
	for(unsigned int t=0; t<freqN_; ++t){
		double w=minFreq_*(t+0.5),hbar=units::metal::hbar,kb=units::metal::kb;
		if(freqUnit_==fourier::FreqUnit::MEV) hbar/=mevPerThz;
		else if(freqUnit_==fourier::FreqUnit::CMI) hbar/=cmiPerThz;
		//ramans_[t]*=std::pow(w/freqVis_-1.0,4)*hbar*w/(w*w)/(kb*T_)/(1.0-std::exp(-hbar*w/(kb*T_)));
		ramans_[t]*=std::tanh(hbar*w/(kb*T_))/(w*std::pow(w/freqVis_-1.0,4));
	}
	
	//record the raman spectrum - total
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating raman spectrum - total...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramant_[t]+=corrODiag[i][j][t];
			}
		}
	}
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=i+1; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramant_[t]+=corrDiag[i][j][t];
			}
		}
	}
	for(unsigned int t=0; t<freqN_; ++t){
		double w=minFreq_*(t+0.5),hbar=units::metal::hbar,kb=units::metal::kb;
		if(freqUnit_==fourier::FreqUnit::MEV) hbar/=mevPerThz;
		else if(freqUnit_==fourier::FreqUnit::CMI) hbar/=cmiPerThz;
		//ramant_[t]*=std::pow(w/freqVis_-1.0,4)*hbar*w/(w*w)/(kb*T_)/(1.0-std::exp(-hbar*w/(kb*T_)));
		ramant_[t]*=std::tanh(hbar*w/(kb*T_))/(w*std::pow(w/freqVis_-1.0,4));
	}
}

void Raman3D::calcSpectrum3(SimAtomic<AtomT>& sim){
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"calcSpectrum3(SimAtomic<AtomT>&):\n";
	//local function variables
	fourier::FFT_R2C fftAlpha1=fourier::FFT_R2C(2*sim.timesteps());//the FFT of the polarizability data
	fourier::FFT_R2C fftAlpha2=fourier::FFT_R2C(2*sim.timesteps());//the FFT of the polarizability data
	fourier::FFT_C2C fft=fourier::FFT_C2C(2*sim.timesteps());//the FFT of chi2
	//diagonal correlation functions
	std::vector<std::vector<std::vector<double> > > corrDiag(3);
	for(unsigned int i=0; i<3; ++i){
		corrDiag[i].resize(3,std::vector<double>(2*sim.timesteps(),0));
	}
	//off-diagonal correlation functions
	std::vector<std::vector<std::vector<double> > > corrODiag(3);
	for(unsigned int i=0; i<3; ++i){
		corrODiag[i].resize(3,std::vector<double>(2*sim.timesteps(),0));
	}
	//transpose of polarizability
	std::vector<std::vector<std::vector<double> > > alphaTrans(3);
	for(unsigned int i=0; i<3; ++i){
		alphaTrans[i].resize(3,std::vector<double>(sim.timesteps(),0));
	}
	
	//No truncation of correlation function
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"No truncation of correlation function.\n";
	
	//calculate the transpose of the total polarizability
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating the transpose of the total polarizability...\n";
	for(unsigned int t=0; t<sim.timesteps(); ++t){
		Eigen::Matrix3d alpha=Eigen::Matrix3d::Zero();
		for(unsigned int n=0; n<sim.nAtoms(); ++n) alpha.noalias()+=sim.atom(t,n).alpha();
		alpha/=sim.nAtoms();
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				alphaTrans[i][j][t]=alpha(i,j);
			}
		}
	}
	
	//calculate the diagonal correlation functions - diagonal
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating correlation functions - diagonal...\n";
	for(unsigned int i=0; i<3; ++i){
		//calculate velocity
		for(unsigned int t=0; t<sim.timesteps(); ++t) fftAlpha1.in(t)=alphaTrans[i][i][t];
		//tranform into the frequency domain
		fftAlpha1.transformf();
		for(unsigned int j=0; j<3; ++j){
			double norm=1;
			if(normalize_) norm=2.0*sim.timesteps();
			//calculate velocity
			for(unsigned int t=0; t<sim.timesteps(); ++t) fftAlpha2.in(t)=alphaTrans[j][j][t];
			//tranform into the frequency domain
			fftAlpha2.transformf();
			//record the correlation function in the frequency domain
			for(unsigned int t=0; t<2*sim.timesteps(); ++t){
				fft.in(t)[0]=(fftAlpha1.out(t)[0]*fftAlpha2.out(t)[0]+fftAlpha1.out(t)[1]*fftAlpha2.out(t)[1]);//real 
				fft.in(t)[1]=(fftAlpha1.out(t)[0]*fftAlpha2.out(t)[1]-fftAlpha1.out(t)[1]*fftAlpha2.out(t)[0]);//imag
			}
			//transform the correlation function back into the time domain
			fft.transformr();
			//normalize and window the correlation function
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(t+1)/(sim.timesteps()-t);
				fft.in(t)[1]=fft.out(t)[1]*window_(t+1)/(sim.timesteps()-t);
			}
			for(unsigned int t=sim.timesteps(); t<2*sim.timesteps(); ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
				fft.in(t)[1]=fft.out(t)[1]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
			}
			//transform to find the windowed, normalized correlation function in the frequency domain
			fft.transformf();
			//record the correlation function in the frequency domain
			for(unsigned int t=0; t<2*sim.timesteps(); ++t){
				corrDiag[i][j][t]=(fft.out(t)[0]*fft.out(t)[0]+fft.out(t)[1]*fft.out(t)[1])/(norm*norm);
			}
		}
	}
	
	//calculate the diagonal correlation functions - off-diagonal
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating correlation functions - off-diagonal...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			double norm=1;
			if(normalize_) norm=2.0*sim.timesteps();
			//calculate velocity
			for(unsigned int t=0; t<sim.timesteps(); ++t) fftAlpha1.in(t)=alphaTrans[i][j][t];
			//tranform into the frequency domain
			fftAlpha1.transformf();
			//calculate velocity
			for(unsigned int t=0; t<sim.timesteps(); ++t) fftAlpha2.in(t)=alphaTrans[i][j][t];
			//tranform into the frequency domain
			fftAlpha2.transformf();
			//record the correlation function in the frequency domain 
			for(unsigned int t=0; t<2*sim.timesteps(); ++t){
				fft.in(t)[0]=(fftAlpha1.out(t)[0]*fftAlpha2.out(t)[0]+fftAlpha1.out(t)[1]*fftAlpha2.out(t)[1]);//real 
				fft.in(t)[1]=(fftAlpha1.out(t)[0]*fftAlpha2.out(t)[1]-fftAlpha1.out(t)[1]*fftAlpha2.out(t)[0]);//imag
			}
			//transform the correlation function back into the time domain
			fft.transformr();
			//normalize and window the correlation function
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(t+1)/(sim.timesteps()-t);
				fft.in(t)[1]=fft.out(t)[1]*window_(t+1)/(sim.timesteps()-t);
			}
			for(unsigned int t=sim.timesteps(); t<2*sim.timesteps(); ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
				fft.in(t)[1]=fft.out(t)[1]*window_(2*sim.timesteps()-t)/(t-sim.timesteps()+1);
			}
			//transform to find the windowed, normalized correlation function in the frequency domain
			fft.transformf();
			//record the correlation function in the frequency domain
			for(unsigned int t=0; t<2*sim.timesteps(); ++t){
				corrODiag[i][j][t]=(fft.out(t)[0]*fft.out(t)[0]+fft.out(t)[1]*fft.out(t)[1])/(norm*norm);
			}
		}
	}
	
	//reset raman spectra
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Resetting Raman spectra...\n";
	for(unsigned int t=0; t<freqN_; ++t) ramanp_[t]=0.0;
	for(unsigned int t=0; t<freqN_; ++t) ramans_[t]=0.0;
	for(unsigned int t=0; t<freqN_; ++t) ramant_[t]=0.0;
	
	//record the raman spectrum - parallel
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating raman spectrum - parallel...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramanp_[t]+=corrDiag[i][j][t];
			}
		}
	}
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramanp_[t]+=2*corrODiag[i][j][t];
			}
		}
	}
	for(unsigned int t=0; t<freqN_; ++t) ramanp_[t]/=15.0;
	for(unsigned int t=0; t<freqN_; ++t){
		double w=minFreq_*(t+0.5),hbar=units::metal::hbar,kb=units::metal::kb;
		if(freqUnit_==fourier::FreqUnit::MEV) hbar/=mevPerThz;
		else if(freqUnit_==fourier::FreqUnit::CMI) hbar/=cmiPerThz;
		//ramanp_[t]*=std::pow(w/freqVis_-1.0,4)*hbar*w/(kb*T_)/(1.0-std::exp(-hbar*w/(kb*T_)));
		ramanp_[t]*=w*std::tanh(hbar*w/(kb*T_))/(std::pow(w/freqVis_-1.0,4));
	}
	
	//record the raman spectrum - perpendicular
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating raman spectrum - perpendicular...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramans_[t]+=3*corrODiag[i][j][t];
			}
		}
	}
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramans_[t]-=corrDiag[i][j][t];
			}
		}
	}
	for(unsigned int t=0; t<freqN_; ++t) ramans_[t]/=30.0;
	for(unsigned int t=0; t<freqN_; ++t){
		double w=minFreq_*(t+0.5),hbar=units::metal::hbar,kb=units::metal::kb;
		if(freqUnit_==fourier::FreqUnit::MEV) hbar/=mevPerThz;
		else if(freqUnit_==fourier::FreqUnit::CMI) hbar/=cmiPerThz;
		//ramans_[t]*=std::pow(w/freqVis_-1.0,4)*hbar*w/(kb*T_)/(1.0-std::exp(-hbar*w/(kb*T_)));
		ramans_[t]*=w*std::tanh(hbar*w/(kb*T_))/(std::pow(w/freqVis_-1.0,4));
	}
	
	//record the raman spectrum - total
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"Calculating raman spectrum - total...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramant_[t]+=corrODiag[i][j][t];
			}
		}
	}
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=i+1; j<3; ++j){
			for(unsigned int t=0; t<freqN_; ++t){
				ramant_[t]+=corrDiag[i][j][t];
			}
		}
	}
	for(unsigned int t=0; t<freqN_; ++t){
		double w=minFreq_*(t+0.5),hbar=units::metal::hbar,kb=units::metal::kb;
		if(freqUnit_==fourier::FreqUnit::MEV) hbar/=mevPerThz;
		else if(freqUnit_==fourier::FreqUnit::CMI) hbar/=cmiPerThz;
		//ramant_[t]*=std::pow(w/freqVis_-1.0,4)*hbar*w/(kb*T_)/(1.0-std::exp(-hbar*w/(kb*T_)));
		ramant_[t]*=w*std::tanh(hbar*w/(kb*T_))/(std::pow(w/freqVis_-1.0,4));
	}
}

void Raman3D::printSpectrum(const char* file) const{
	if(DEBUG_RAMAN3D_THOLEV>0) log<<"printSpectrum():\n";
	//local function variables
	FILE* writer=NULL;
	
	writer=fopen(file,"w");
	if(writer==NULL) throw std::runtime_error("I/O Exception Occured.");
	
	//print the header
	fprintf(writer, "Freq");
	if(freqUnit_==fourier::FreqUnit::MEV) fprintf(writer, "(meV) ");
	else if(freqUnit_==fourier::FreqUnit::CMI) fprintf(writer, "(cm^-1) ");
	else fprintf(writer, "(THz) ");
	fprintf(writer, "RamanP RamanS RamanT\n");
	
	//print the power spectrum
	for(unsigned int t=0; t<freqN_; ++t) fprintf(writer, "%f %f %f %f\n", minFreq_*t, ramanp_[t], ramans_[t], ramant_[t]);
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
		unsigned int version=0;
		int beg,end;
		double zFactor=1;
		Eigen::Vector3d offset=Eigen::Vector3d::Zero();
		SimAtomic<AtomT> sim;
	//ewald
		Ewald3D::Dipole dipole;
	//calculation
		Raman3D raman3d;
		Thole thole;
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
		units::System::type unitsys;
	//logging
		logging::Log& logger=logging::Log::get();//log object
		logger.addSink(logging::Sink(&std::cout));//add console sink
		logger.addSink(logging::Sink(new std::ofstream("raman3d_tholev.log")));
		logging::DebugLogger log("raman3d_tholev");
	
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
			} else if(std::strcmp(temp,"VERSION")==0){
				version=std::atoi(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"ATOMS")==0){
				std::strcpy(temp,std::strpbrk(input,"=")+1);
				if(string::substrN(temp,string::WS)==0) throw std::invalid_argument("Invalid number of atoms.");
				atoms.resize(string::substrN(temp,string::WS));
				atoms[0]=std::string(std::strtok(temp,string::WS));
				for(unsigned int n=1; n<atoms.size(); ++n) atoms[n]=std::string(std::strtok(NULL,string::WS));
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
		if(unitsys==units::System::UNKNOWN) throw std::invalid_argument("Invalid unit system.");
		if(!(version==2 || version==3)) throw std::invalid_argument("Invalid version.");
		
		//print the parameters to screen
		log<<"SIMULATION PARAMETERS:\n";
		log<<"\tVERSION  = "<<version<<"\n";
		log<<"\tUNITS    = "<<unitsys<<"\n";
		log<<"\tSIM      = "<<simstr<<"\n";
		log<<"\tFORMAT   = "<<fileFormat<<"\n";
		log<<"\tOFFSET   = ("<<offset[0]<<","<<offset[1]<<","<<offset[2]<<")\n";
		log<<"\tINTERVAL = "<<beg<<", "<<end<<"\n";
		log<<"\tATOMS    = "; for(unsigned int i=0; i<atoms.size(); ++i) log.log()<<atoms[i]<<" "; log.log()<<"\n";
		log<<"\tALPHA    = "; for(unsigned int i=0; i<alpha.size(); ++i) log.log()<<alpha[i]<<" "; log.log()<<"\n";
		
		//initialize the unit system
		log<<"Initializing the unit system...\n";
		units::consts::init(unitsys);
		
		//load the simulation
		if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Loading simulation...\n";
		if(fileFormat==FILE_FORMAT::POSCAR) VASP::POSCAR::load(simstr,sim);
		else if(fileFormat==FILE_FORMAT::XDATCAR) VASP::XDATCAR::load(simstr,sim,beg,end);
		else if(fileFormat==FILE_FORMAT::LAMMPS){
			if(string::substrN(simstr,",")!=3) throw std::invalid_argument("Invalid LAMMPS file format.");
			std::vector<std::string> substr(3);
			substr[0]=std::strtok(simstr,",");
			substr[1]=std::strtok(NULL,",");
			substr[2]=std::strtok(NULL,",");
			LAMMPS::load(substr[0].c_str(),substr[1].c_str(),substr[2].c_str(),sim,beg,end);
		} else if(fileFormat==FILE_FORMAT::GAUSSIAN){
			GAUSSIAN::Format formatGauss;
			formatGauss.formatCalc=GAUSSIAN::FormatCalc::ADMP;
			formatGauss.formatVersion=GAUSSIAN::FormatVersion::g16;
			formatGauss.log=simstr;
			GAUSSIAN::load(formatGauss,sim);
		} else throw std::invalid_argument("Invalid file format.");
		//set the periodicity
		sim.periodic()=periodic;
		
		/* Set the offset */
		if(offset.norm()>num_const::ZERO){
			log<<"Applying offset...\n";
			//apply the offset 
			if(sim.cell(0).R().determinant()>0){
				for(unsigned int t=0; t<sim.timesteps(); ++t){
					for(unsigned int n=0; n<sim.nSpecies(); ++n){
						for(unsigned int m=0; m<sim.nAtoms(n); ++m){
							sim.atom(t,n,m).posn().noalias()+=offset;
							Cell::returnToCell(sim.atom(t,n,m).posn(),sim.atom(t,n,m).posn(),sim.cell(t).R(),sim.cell(t).RInv());
						}
					}
				}
			} else {
				for(unsigned int t=0; t<sim.timesteps(); ++t){
					for(unsigned int n=0; n<sim.nSpecies(); ++n){
						for(unsigned int m=0; m<sim.nAtoms(n); ++m){
							sim.atom(t,n,m).posn().noalias()+=offset;
						}
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
		log<<"CELL = \n"<<sim.cell(0).R()<<"\n";
		
		/* Create the auxiliary objects */
		log<<"Creating the auxiliary objects...\n";
		//ewald object
			log<<"\tEwald objects...\n";
			dipole.init(sim.cell(0),1e-5);
			log<<"Ewald::Dipole = \n"<<dipole<<"\n";
		//thole object
			log<<"\tThole object...\n";
			Thole::load(paramFile,thole);
			thole.init(sim,dipole);
			log<<"Thole = \n"<<thole<<"\n";
		//raman3d object
			log<<"\tRaman3D object...\n";
			raman3d.load(paramFile);
			raman3d.init(sim);
			log<<"Raman3D = \n"<<raman3d<<"\n";
		
		/* Perform the calculation */
		log<<"Performing the calculation...\n";
		if(raman3d.calcAlpha()){
			log<<"Calculating atomic alphas...\n";
			raman3d.calcAlpha(sim,thole,dipole);
		}
		if(raman3d.calcSpectrum()){
			log<<"Calculating the spectrum...\n";
			if(version==2) raman3d.calcSpectrum2(sim);
			else if(version==3) raman3d.calcSpectrum3(sim);
			raman3d.printSpectrum(raman3d.fileSpectrum().c_str());
		}
		
		/* Print the Relevant Data */
		std::cout<<"printing relevant data...\n";
		if(raman3d.printAlphaT()) Electrostatics::Print::alpha_tot("alpha_tot.dat",sim);
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
