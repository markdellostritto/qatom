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
	out<<"FREQ_MIN         = "<<raman3d.minFreq_<<"\n";
	out<<"FREQ_MAX         = "<<raman3d.maxFreq_<<"\n";
	out<<"FREQ_CUT         = "<<raman3d.freqCut_<<"\n";
	out<<"FREQ_VIS         = "<<raman3d.freqVis_<<"\n";
	out<<"FREQ_RES         = "<<raman3d.freqRes_<<"\n";
	out<<"WINDOW           = "<<raman3d.windowType_<<"\n";
	out<<"STRIDE_ALPHA     = "<<raman3d.strideAlpha_<<"\n";
	out<<"CALC_ALPHA       = "<<raman3d.calcAlpha_<<"\n";
	out<<"CALC_SPECTRUM    = "<<raman3d.calcSpectrum_<<"\n";
	out<<"WRITE_ALPHA_TOT  = "<<raman3d.writeAlphaT_<<"\n";
	out<<"WRITE_ALPHA_ATOM = "<<raman3d.writeAlphaA_<<"\n";
	out<<"READ_ALPHA_ATOM  = "<<raman3d.readAlphaA_<<"\n";
	out<<"PROFILE_CALC     = "<<raman3d.profileCalc_<<"\n";
	out<<"PROFILE_LOAD     = "<<raman3d.profileLoad_<<"\n";
	out<<"SUBSET           = "<<raman3d.subsetStr_<<"\n";
	out<<"TEMP             = "<<raman3d.T_<<"\n";
	out<<"NORMALIZE        = "<<raman3d.normalize_<<"\n";
	out<<"FILE_SPECTRUM    = "<<raman3d.fileSpectrum_<<"\n";
	out<<"FILE_ALPHA_TOT   = "<<raman3d.fileAlphaT_<<"\n";
	out<<"FILE_ALPHA_ATOM  = "<<raman3d.fileAlphaA_<<"\n";
	out<<"****************** RAMAN3D_THOLEV ******************\n";
	out<<"*************************************************";
	return out;
}

void Raman3D::defaults(){
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"defaults()\n";
	//fft parameters
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
		strideAlpha_=1;
		nStepsAlpha_=0;
	//subset
		subsetStr_.clear();
		subset_.clear();
	//i/o
		fileSpectrum_=std::string("raman3d.dat");
		fileAlphaT_=std::string("alpha_tot.dat");
		fileAlphaA_=std::string("alpha_atom.dat");
	//calculation flags
		calcAlpha_=true;
		calcSpectrum_=true;
		normalize_=true;
	//i/o flags
		writeAlphaT_=false;
		writeAlphaA_=false;
	//temperature
		T_=300.0;
	//raman spectrum
		ramanp_.clear();
		ramans_.clear();
		ramant_.clear();
}

void Raman3D::read(const char* file){
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"read(const char*):\n";
	//local function variables
	FILE* reader=NULL;
	char* input=new char[string::M];
	char* temp=new char[string::M];
	std::vector<std::string> strlist;
	bool error=false;
	
	try{
		//open the parameter file
		if(DEBUG_RAMAN3D_THOLEV>1) std::cout<<"Opening parameter file...\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open parameter file.");
		
		//set the defaults
		if(DEBUG_RAMAN3D_THOLEV>1) std::cout<<"Setting defaults...\n";
		defaults();
		
		//read in the parameters
		if(DEBUG_RAMAN3D_THOLEV>1) std::cout<<"Reading in parameters...\n";
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);
			if(string::split(input,string::WS,strlist)==0) continue;
			string::to_upper(strlist.at(0));
			if(strlist.at(0)=="FILE_SPECTRUM"){
				fileSpectrum_=strlist.at(1);
			} else if(strlist.at(0)=="CALC_SPECTRUM"){
				calcSpectrum_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="CALC_ALPHA"){
				calcAlpha_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="WRITE_ALPHA_TOT"){
				writeAlphaT_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="WRITE_ALPHA_ATOM"){
				writeAlphaA_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="READ_ALPHA_ATOM"){
				readAlphaA_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="NORMALIZE"){
				normalize_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="FREQ_CUT"){
				freqCut_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="RESOLUTION"){
				freqRes_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="FREQ_VIS"){
				freqVis_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="T"){
				T_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="STRIDE_ALPHA"){
				strideAlpha_=std::atoi(strlist.at(1).c_str());
			} else if(strlist.at(0)=="WINDOW"){
				windowType_=window::WINDOW_FUNC::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="PROFILE_CALC"){
				profileCalc_.sigma()=std::atof(strlist.at(1).c_str());
				profileCalc_.b1()=std::atof(strlist.at(2).c_str());
				profileCalc_.b2()=std::atof(strlist.at(3).c_str());
			} else if(strlist.at(0)=="PROFILE_LOAD"){
				profileLoad_.sigma()=std::atof(strlist.at(1).c_str());
				profileLoad_.b1()=std::atof(strlist.at(2).c_str());
				profileLoad_.b2()=std::atof(strlist.at(3).c_str());
			} else if(strlist.at(0)=="SUBSET"){
				subsetStr_=strlist.at(1);
			}
		}
		
		//close the parameter file
		if(DEBUG_RAMAN3D_THOLEV>1) std::cout<<"Closing parameter file...\n";
		fclose(reader);
		reader=NULL;
		
		//check the calculation parameters
		if(DEBUG_RAMAN3D_THOLEV>1) std::cout<<"Checking parameters...\n";
		if(strideAlpha_==0) throw std::invalid_argument("Invalid stride.");
		if(freqVis_<=0) throw std::invalid_argument("Invalid visible frequency.");
		if(T_<=0) throw std::invalid_argument("Invalid temperature.");
	}catch(std::exception& e){
		std::cout<<"ERROR in read(const char*):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("ERROR in Raman3D::read(const char*): Failed to read.");
}

void Raman3D::init(Simulation& sim){
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"init(const Simulation&):\n";
	bool error=false;
	
	try{
		//set the fft parameters
		if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Assigning the FFT parameters...\n";
		nStepsAlpha_=sim.timesteps()/strideAlpha_;
		minFreq_=0.5*(1000.0/sim.timestep())/sim.timesteps();//min freq in THz (1000 for ps, 1/2 due to padding)
		maxFreq_=0.5*(1000.0/sim.timestep());//max freq in THz (1000 for ps, 1/2 due to padding)
		if(freqRes_==0) freqRes_=5.0*minFreq_;
		if(freqCut_==0) freqCut_=maxFreq_;
		
		//set the fourier parameters in the appropriate units
		freqN_=(unsigned int)std::ceil((freqCut_+1.0)/minFreq_);
		
		//set the window_
		if(windowType_==window::WINDOW_FUNC::IDENTITY) window_=window::Identity();
		else if(windowType_==window::WINDOW_FUNC::BLACKMANHARRIS) window_=window::BlackmanHarris(sim.timesteps());
		else if(windowType_==window::WINDOW_FUNC::GAUSSIAN){
			sigma_=1/(2*num_const::PI*freqRes_);//from fourier transform
			window_=window::Gaussian(sim.timesteps(),sigma_);
		}
		/*else if(windowType_==window::WINDOW_FUNC::KAISERBESSEL){
			sigma_=freqRes_*freqRes_;//scale to get sensible values
			window_=window::KaiserBessel(sim.timesteps(),sigma_);
		}*/
		
		//check the fft parameters
		if(freqCut_-minFreq_<num_const::ZERO) throw std::invalid_argument("Invalid freqency cutoff.");
		
		//allocate space for raman3d
		if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Allocating space for ir spectrum...\n";
		ramanp_.resize(freqN_);
		ramans_.resize(freqN_);
		ramant_.resize(freqN_);
		
		//read the subset
		if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Loading the subset...\n";
		if(subsetStr_.length()>0) StructureI::read_atoms(subsetStr_.c_str(),subset_,sim.frame(0));
		
	}catch(std::exception& e){
		std::cout<<"ERROR in read(const char*):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
}

void Raman3D::calcAlpha(Simulation& sim, const Thole& thole, const Ewald3D::Dipole& ewald){
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"calcAlpha(Simulation&,const Thole&,const Ewald3D::Dipole&):\n";
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
		std::vector<Ewald3D::Dipole> ewald_(nThreads);
		
	//initialize the effective alphas
	if(DEBUG_RAMAN3D_THOLEV>1) std::cout<<"Initializing Thole objects...\n";
	for(unsigned int i=0; i<nThreads; ++i) thole_[i]=thole;
	for(unsigned int i=0; i<nThreads; ++i) ewald_[i]=ewald;
	
	//calculate the atomic alphas
	if(DEBUG_RAMAN3D_THOLEV>1) std::cout<<"Calculating Atomic Alphas...\n";
	start=std::chrono::high_resolution_clock::now();
	if(!sim.cell_fixed()){
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideAlpha_){
			unsigned int TN=0;
			#ifdef _OPENMP
				TN=omp_get_thread_num();
			#endif
			if(DEBUG_RAMAN3D_THOLEV>-1) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_RAMAN3D_THOLEV>0 && t%1000==0) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			ewald_[TN].init(sim.frame(t).cell(),1e-5);
			thole_[TN].alpha(sim.frame(t),ewald_[TN]);
		}
	} else {
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideAlpha_){
			unsigned int TN=0;
			#ifdef _OPENMP
				TN=omp_get_thread_num();
			#endif
			if(DEBUG_RAMAN3D_THOLEV>-1) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_RAMAN3D_THOLEV>0 && t%1000==0) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			thole_[TN].alpha_cont(sim.frame(t),ewald);
		}
	}
	stop=std::chrono::high_resolution_clock::now();
	time=std::chrono::duration_cast<std::chrono::duration<double> >(stop-start);
	std::cout<<"alpha time = "<<time.count()<<"\n";
	
	//interpolate the polarizabilities
	if(strideAlpha_>1){
		if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Interpolating Effective Alphas...\n";
		std::vector<std::vector<Interp::Data> > alphaData;
		alphaData.resize(3);
		for(unsigned int i=0; i<3; ++i){
			alphaData[i].resize(3);
			for(unsigned int j=0; j<3; ++j){
				alphaData[i][j].resize(nStepsAlpha_+1);
			}
		}
		//loop over all atoms included in the calculation
		for(unsigned int n=0; n<sim.frame(0).nAtoms(); ++n){
			//loop over all timesteps for which we calculated effective polarizabilities
			for(unsigned int t=0; t*strideAlpha_<sim.timesteps(); ++t){
				for(unsigned int i=0; i<3; ++i){
					for(unsigned int j=0; j<3; ++j){
						alphaData[i][j].x(t)=t*strideAlpha_;
						alphaData[i][j].y(t)=sim.frame(t*strideAlpha_).alpha(n)(i,j);
					}
				}
			}
			//interpolate the rest of the timesteps
			for(unsigned int t=0; t*strideAlpha_<sim.timesteps()-strideAlpha_; ++t){
				for(unsigned int tp=1; tp<strideAlpha_; ++tp){
					for(unsigned int i=0; i<3; ++i){
						for(unsigned int j=0; j<3; ++j){
							sim.frame(t*strideAlpha_+tp).alpha(n)(i,j)=Interp::interpAkima(t*strideAlpha_+tp,alphaData[i][j]);
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
		for(unsigned int n=0; n<sim.frame(0).nAtoms(); ++n){
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					//calculate the slope at the end
					m=1.0/sim.timestep()*(sim.frame(end).alpha(n)(i,j)-sim.frame(end-1).alpha(n)(i,j));
					//perform linear extrapolation for the last few timesteps
					for(unsigned int t=end; t<sim.timesteps(); ++t){
						sim.frame(t).alpha(n)(i,j)=sim.frame(end).alpha(n)(i,j)+m*(t-end)*sim.timestep();
					}
				}
			}
		}
	}
}

void Raman3D::calcSpectrum(Simulation& sim){
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"calcSpectrum2(Simulation&):\n";
	//constants
	const double TS=sim.timesteps();
	const double ts=sim.timestep();
	double hbar=0,kb=0;
	if(units::consts::system()==units::System::AU){
		hbar=units::au::hbar;
		kb=units::au::kb;
	} else if(units::consts::system()==units::System::METAL){
		hbar=units::metal::hbar;
		kb=units::metal::kb;
	}
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
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"No truncation of correlation function.\n";
	
	//calculate the transpose of the total polarizability
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Calculating the transpose of the total polarizability...\n";
	if(subset_.size()==0){
		for(unsigned int t=0; t<TS; ++t){
			Eigen::Matrix3d alpha=Eigen::Matrix3d::Zero();
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n) alpha.noalias()+=sim.frame(t).alpha(n);
			alpha/=sim.frame(t).nAtoms();
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					alphaTrans[i][j][t]=alpha(i,j);
				}
			}
		}
	} else {
		for(unsigned int t=0; t<TS; ++t){
			Eigen::Matrix3d alpha=Eigen::Matrix3d::Zero();
			for(unsigned int n=0; n<subset_.size(); ++n) alpha.noalias()+=sim.frame(t).alpha(subset_[n]);
			alpha/=sim.frame(t).nAtoms();
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					alphaTrans[i][j][t]=alpha(i,j);
				}
			}
		}
	}
	
	//calculate the diagonal correlation functions - diagonal
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Calculating correlation functions - diagonal...\n";
	for(unsigned int i=0; i<3; ++i){
		//calculate velocity
		fftAlpha1.in(0)=gradient::df1o1(alphaTrans[i][i],ts,0);
		fftAlpha1.in(1)=gradient::dc1o2(alphaTrans[i][i],ts,1);
		fftAlpha1.in(2)=gradient::dc1o4(alphaTrans[i][i],ts,2);
		fftAlpha1.in(3)=gradient::dc1o6(alphaTrans[i][i],ts,3);
		for(int t=4; t<TS-4; t++){
			fftAlpha1.in(t)=gradient::dc1o8(alphaTrans[i][i],ts,t);
		}
		fftAlpha1.in(TS-4)=gradient::dc1o6(alphaTrans[i][i],ts,TS-4);
		fftAlpha1.in(TS-3)=gradient::dc1o4(alphaTrans[i][i],ts,TS-3);
		fftAlpha1.in(TS-2)=gradient::dc1o2(alphaTrans[i][i],ts,TS-2);
		fftAlpha1.in(TS-1)=gradient::db1o1(alphaTrans[i][i],ts,TS-1);
		//tranform into the frequency domain
		fftAlpha1.transformf();
		for(unsigned int j=0; j<3; ++j){
			//calculate velocity
			fftAlpha2.in(0)=gradient::df1o1(alphaTrans[j][j],ts,0);
			fftAlpha2.in(1)=gradient::dc1o2(alphaTrans[j][j],ts,1);
			fftAlpha2.in(2)=gradient::dc1o4(alphaTrans[j][j],ts,2);
			fftAlpha2.in(3)=gradient::dc1o6(alphaTrans[j][j],ts,3);
			for(int t=4; t<TS-4; t++){
				fftAlpha2.in(t)=gradient::dc1o8(alphaTrans[j][j],ts,t);
			}
			fftAlpha2.in(TS-4)=gradient::dc1o6(alphaTrans[j][j],ts,TS-4);
			fftAlpha2.in(TS-3)=gradient::dc1o4(alphaTrans[j][j],ts,TS-3);
			fftAlpha2.in(TS-2)=gradient::dc1o2(alphaTrans[j][j],ts,TS-2);
			fftAlpha2.in(TS-1)=gradient::db1o1(alphaTrans[j][j],ts,TS-1);
			//tranform into the frequency domain
			fftAlpha2.transformf();
			//record the correlation function in the frequency domain
			for(unsigned int t=0; t<2*TS; ++t){
				fft.in(t)[0]=fftAlpha1.out(t)[0]*fftAlpha2.out(t)[0]+fftAlpha1.out(t)[1]*fftAlpha2.out(t)[1];//real 
				fft.in(t)[1]=fftAlpha1.out(t)[0]*fftAlpha2.out(t)[1]-fftAlpha1.out(t)[1]*fftAlpha2.out(t)[0];//imag
			}
			//transform the correlation function back into the time domain
			fft.transformr();
			//normalize and window the correlation function
			for(unsigned int t=0; t<TS; ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(t+1)/(TS-t);
				fft.in(t)[1]=fft.out(t)[1]*window_(t+1)/(TS-t);
			}
			for(unsigned int t=TS; t<2*TS; ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(2*TS-t)/(t-TS+1);
				fft.in(t)[1]=fft.out(t)[1]*window_(2*TS-t)/(t-TS+1);
			}
			//transform to find the windowed, normalized correlation function in the frequency domain
			fft.transformf();
			//record the correlation function in the frequency domain
			double norm=1;
			if(normalize_) norm=2.0*TS;
			for(unsigned int t=0; t<2*TS; ++t){
				corrDiag[i][j][t]=(fft.out(t)[0]*fft.out(t)[0]+fft.out(t)[1]*fft.out(t)[1])/norm;
			}
		}
	}
	
	//calculate the diagonal correlation functions - off-diagonal
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Calculating correlation functions - off-diagonal...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			//calculate velocity
			fftAlpha1.in(0)=gradient::df1o1(alphaTrans[i][j],ts,0);
			fftAlpha1.in(1)=gradient::dc1o2(alphaTrans[i][j],ts,1);
			fftAlpha1.in(2)=gradient::dc1o4(alphaTrans[i][j],ts,2);
			fftAlpha1.in(3)=gradient::dc1o6(alphaTrans[i][j],ts,3);
			for(int t=4; t<TS-4; t++){
				fftAlpha1.in(t)=gradient::dc1o8(alphaTrans[i][j],ts,t);
			}
			fftAlpha1.in(TS-4)=gradient::dc1o6(alphaTrans[i][j],ts,TS-4);
			fftAlpha1.in(TS-3)=gradient::dc1o4(alphaTrans[i][j],ts,TS-3);
			fftAlpha1.in(TS-2)=gradient::dc1o2(alphaTrans[i][j],ts,TS-2);
			fftAlpha1.in(TS-1)=gradient::db1o1(alphaTrans[i][j],ts,TS-1);
			//tranform into the frequency domain
			fftAlpha1.transformf();
			//calculate velocity
			fftAlpha2.in(0)=gradient::df1o1(alphaTrans[i][j],ts,0);
			fftAlpha2.in(1)=gradient::dc1o2(alphaTrans[i][j],ts,1);
			fftAlpha2.in(2)=gradient::dc1o4(alphaTrans[i][j],ts,2);
			fftAlpha2.in(3)=gradient::dc1o6(alphaTrans[i][j],ts,3);
			for(int t=4; t<TS-4; t++){
				fftAlpha2.in(t)=gradient::dc1o8(alphaTrans[i][j],ts,t);
			}
			fftAlpha2.in(TS-4)=gradient::dc1o6(alphaTrans[i][j],ts,TS-4);
			fftAlpha2.in(TS-3)=gradient::dc1o4(alphaTrans[i][j],ts,TS-3);
			fftAlpha2.in(TS-2)=gradient::dc1o2(alphaTrans[i][j],ts,TS-2);
			fftAlpha2.in(TS-1)=gradient::db1o1(alphaTrans[i][j],ts,TS-1);
			//tranform into the frequency domain
			fftAlpha2.transformf();
			//record the correlation function in the frequency domain 
			for(unsigned int t=0; t<2*TS; ++t){
				fft.in(t)[0]=fftAlpha1.out(t)[0]*fftAlpha2.out(t)[0]+fftAlpha1.out(t)[1]*fftAlpha2.out(t)[1];//real 
				fft.in(t)[1]=fftAlpha1.out(t)[0]*fftAlpha2.out(t)[1]-fftAlpha1.out(t)[1]*fftAlpha2.out(t)[0];//imag
			}
			//transform the correlation function back into the time domain
			fft.transformr();
			//normalize and window the correlation function
			for(unsigned int t=0; t<TS; ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(t+1)/(TS-t);
				fft.in(t)[1]=fft.out(t)[1]*window_(t+1)/(TS-t);
			}
			for(unsigned int t=TS; t<2*TS; ++t){
				fft.in(t)[0]=fft.out(t)[0]*window_(2*TS-t)/(t-TS+1);
				fft.in(t)[1]=fft.out(t)[1]*window_(2*TS-t)/(t-TS+1);
			}
			//transform to find the windowed, normalized correlation function in the frequency domain
			fft.transformf();
			//record the correlation function in the frequency domain
			double norm=1;
			if(normalize_) norm=2.0*TS;
			for(unsigned int t=0; t<2*TS; ++t){
				corrODiag[i][j][t]=(fft.out(t)[0]*fft.out(t)[0]+fft.out(t)[1]*fft.out(t)[1])/norm;
			}
		}
	}
	
	//reset raman spectra
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Resetting Raman spectra...\n";
	for(unsigned int t=0; t<freqN_; ++t) ramanp_[t]=0.0;
	for(unsigned int t=0; t<freqN_; ++t) ramans_[t]=0.0;
	for(unsigned int t=0; t<freqN_; ++t) ramant_[t]=0.0;
	
	//record the raman spectrum - parallel
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Calculating raman spectrum - parallel...\n";
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
		const double w=minFreq_*(t+0.5);
		//ramanp_[t]*=std::tanh(hbar*w/(kb*T_))/(1.0-std::exp(-hbar*w/(kb*T_)))*std::pow(w/freqVis_-1.0,-4);
		//ramanp_[t]*=std::tanh(hbar*w/(kb*T_))/(hbar*w/(kb*T_)*std::pow(w/freqVis_-1.0,-4));
		ramanp_[t]*=1.0/std::pow(w/freqVis_-1.0,-4);
	}
	
	//record the raman spectrum - perpendicular
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Calculating raman spectrum - perpendicular...\n";
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
		const double w=minFreq_*(t+0.5);
		//ramans_[t]*=std::tanh(hbar*w/(kb*T_))/(1.0-std::exp(-hbar*w/(kb*T_)))*std::pow(w/freqVis_-1.0,-4);
		//ramans_[t]*=std::tanh(hbar*w/(kb*T_))/(hbar*w/(kb*T_)*std::pow(w/freqVis_-1.0,-4));
		ramans_[t]*=1.0/std::pow(w/freqVis_-1.0,-4);
	}
	
	//record the raman spectrum - total
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"Calculating raman spectrum - total...\n";
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
		const double w=minFreq_*(t+0.5);
		//ramant_[t]*=std::tanh(hbar*w/(kb*T_))/(1.0-std::exp(-hbar*w/(kb*T_)))*std::pow(w/freqVis_-1.0,-4);
		//ramant_[t]*=std::tanh(hbar*w/(kb*T_))/(hbar*w/(kb*T_)*std::pow(w/freqVis_-1.0,-4));
		ramant_[t]*=1.0/std::pow(w/freqVis_-1.0,-4);
	}
}

void Raman3D::printSpectrum(const char* file) const{
	if(DEBUG_RAMAN3D_THOLEV>0) std::cout<<"printSpectrum():\n";
	//local function variables
	FILE* writer=NULL;
	
	writer=fopen(file,"w");
	if(writer==NULL) throw std::runtime_error("I/O Exception Occured.");
	
	//print the header
	fprintf(writer, "#Freq(THz) RamanP RamanS RamanT\n");
	
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
	
	//======== local variables ========
	//simulation
		Eigen::Vector3d offset=Eigen::Vector3d::Zero();
		Simulation sim;
		Interval interval;
		AtomType atomT;
		atomT.name=true; atomT.an=true; atomT.specie=true; atomT.index=true;
		atomT.posn=true; atomT.alpha=true;
	//ewald
		Ewald3D::Dipole dipole;
	//calculation
		Raman3D raman3d;
		Thole thole;
		std::vector<double> alpha;
		std::vector<std::string> names;
	//input/output
		char* input=new char[string::M];
		std::string paramfile;
		std::string simstr;
		FILE* reader=NULL;
		FILE_FORMAT::type fileFormat;
	//flags
		bool error=false;
	//units
		units::System::type unitsys;
		std::vector<std::string> strlist;
	
	try{
		//======== check the number of arguments ========
		if(argc!=2){
			std::cout<<"ERROR in main(int,char**):\n";
			std::cout<<"Incorrect number of arguments.\n";
			std::cout<<"\t1. Program Execution\n";
			std::cout<<"\t2. Parameter File\n";
			throw std::invalid_argument("Invalid command-line arguments.");
		}
		
		//======== open parameter file ========
		std::cout<<"opening parameter file...\n";
		paramfile=argv[1];
		reader=fopen(paramfile.c_str(),"r");
		if(reader==NULL) throw std::invalid_argument("I/O Error: Could not open parameter file.");
		
		//======== read parameters ========
		std::cout<<"reading parameters...\n";
		while(std::fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);
			if(string::split(input,string::WS,strlist)==0) continue;
			string::to_upper(strlist.at(0));
			if(strlist.at(0)=="SIM"){
				simstr=strlist.at(1);
			} else if(strlist.at(0)=="TIMESTEP"){
				sim.timestep()=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="INTERVAL"){
				interval=Interval::read(strlist.at(1).c_str());
			} else if(strlist.at(0)=="OFFSET"){
				offset[0]=std::atof(strlist.at(1).c_str());
				offset[1]=std::atof(strlist.at(2).c_str());
				offset[2]=std::atof(strlist.at(3).c_str());
			} else if(strlist.at(0)=="FORMAT"){
				fileFormat=FILE_FORMAT::read(strlist.at(1));
			} else if(strlist.at(0)=="UNITS"){
				unitsys=units::System::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="ALPHA"){
				names.push_back(strlist.at(1));
				alpha.push_back(std::atof(strlist.at(2).c_str()));
			}
		}
		fclose(reader);
		reader=NULL;
		
		//======== check the parameters ========
		if(unitsys==units::System::UNKNOWN) throw std::invalid_argument("Invalid unit system.");
		
		//======== print the parameters to screen ========
		std::cout<<"SIMULATION PARAMETERS:\n";
		std::cout<<"\tUNITS    = "<<unitsys<<"\n";
		std::cout<<"\tSIM      = "<<simstr<<"\n";
		std::cout<<"\tFORMAT   = "<<fileFormat<<"\n";
		std::cout<<"\tOFFSET   = "<<offset[0]<<" "<<offset[1]<<" "<<offset[2]<<"\n";
		std::cout<<"\tINERVAL  = "<<interval<<"\n";
		std::cout<<"\tATOMS    = "; for(unsigned int i=0; i<names.size(); ++i) std::cout<<names[i]<<" "; std::cout<<"\n";
		std::cout<<"\tALPHA    = "; for(unsigned int i=0; i<alpha.size(); ++i) std::cout<<alpha[i]<<" "; std::cout<<"\n";
		
		//======== initialize the unit system ========
		std::cout<<"initializing the unit system...\n";
		units::consts::init(unitsys);
		
		//======== read the simulation ========
		std::cout<<"reading simulation...\n";
		if(fileFormat==FILE_FORMAT::XDATCAR){
			VASP::XDATCAR::read(simstr.c_str(),interval,atomT,sim);
		} else if(fileFormat==FILE_FORMAT::LAMMPS){
			LAMMPS::Format format;
			std::strcpy(input,simstr.c_str());
			if(string::substrN(input,",")!=3) throw std::invalid_argument("Invalid LAMMPS file format.");
			format.in=std::strtok(input,",");
			format.data=std::strtok(NULL,",");
			format.dump=std::strtok(NULL,",");
			LAMMPS::read(format,interval,atomT,sim);
		} else if(fileFormat==FILE_FORMAT::QE){
			QE::Format format;
			std::strcpy(input,simstr.c_str());
			if(string::substrN(input,",")!=3) throw std::invalid_argument("Invalid QE file format.");
			format.fileIn=std::strtok(input,",");
			format.fileCel=std::strtok(NULL,",");
			format.filePos=std::strtok(NULL,",");
			QE::read(format,interval,atomT,sim);
		} else throw std::invalid_argument("Invalid file format.");
		
		//======== set offset ========
		if(offset.norm()>num_const::ZERO){
			std::cout<<"setting offset...\n";
			if(sim.frame(0).cell().R().determinant()>0){
				for(unsigned int t=0; t<sim.timesteps(); ++t){
					for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
						sim.frame(t).posn(n).noalias()+=offset;
						Cell::returnToCell(
							sim.frame(t).posn(n),sim.frame(t).posn(n),
							sim.frame(t).cell().R(),sim.frame(t).cell().RInv()
						);
					}
				}
			} else {
				for(unsigned int t=0; t<sim.timesteps(); ++t){
					for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
						sim.frame(t).posn(n).noalias()+=offset;
					}
				}
			}
		}
		
		//======== assign alpha ========
		std::cout<<"assigning alpha...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				for(unsigned int i=0; i<names.size(); ++i){
					if(names[i]==sim.frame(t).name(n)){
						sim.frame(t).alpha(n)=Eigen::Matrix3d::Identity()*alpha[i];
						break;
					}
				}
			}
		}
		
		//======== print the simulation ========
		std::cout<<sim<<"\n";
		std::cout<<sim.frame(0).cell().R()<<"\n";
		
		//read auxiliary objects
		std::cout<<"reading auxiliary objects...\n";
		//ewald object
			std::cout<<"Ewald objects...\n";
			dipole.init(sim.frame(0).cell(),1e-5);
			std::cout<<dipole<<"\n";
		//thole object
			std::cout<<"Thole object...\n";
			thole.read(paramfile.c_str());
			thole.init(sim.frame(0),dipole);
			std::cout<<thole<<"\n";
		//raman3d object
			std::cout<<"Raman3D object...\n";
			raman3d.read(paramfile.c_str());
			raman3d.init(sim);
			std::cout<<raman3d<<"\n";
		
		//======== read alpha ========
		if(raman3d.readAlphaA()){
			std::cout<<"reading alpha data...\n";
			Utility::Read::alpha_atom(raman3d.fileAlphaA().c_str(),sim);
		}
		
		//======== execute calculation ========
		std::cout<<"executing calculation...\n";
		if(raman3d.calcAlpha()){
			std::cout<<"computing atomic alphas...\n";
			raman3d.calcAlpha(sim,thole,dipole);
		}
		if(raman3d.calcSpectrum()){
			std::cout<<"computing raman spectrum...\n";
			raman3d.calcSpectrum(sim);
			raman3d.printSpectrum(raman3d.fileSpectrum().c_str());
		}
		
		//======== write relevant data ========
		std::cout<<"writing relevant data...\n";
		if(raman3d.writeAlphaT()) Utility::Write::alpha_tot(raman3d.fileAlphaT().c_str(),sim);
		if(raman3d.writeAlphaA()) Utility::Write::alpha_atom(raman3d.fileAlphaA().c_str(),sim);
	}catch(std::exception& e){
		std::cout<<"ERROR in main(int,char**):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//======== free all local variables ========
	std::cout<<"freeing local variables...\n";
	delete[] input;
	if(reader!=NULL) fclose(reader);
	
	std::cout<<"exiting the program...\n";
	if(!error) return 0;
	else return 1;
}
