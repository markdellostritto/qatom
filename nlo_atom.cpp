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
	out<<"******************* NLO_ATOM *******************\n";
	out<<"T                = "<<nlo.T_<<"\n";
	out<<"FREQ_MIN         = "<<nlo.minFreq_<<"\n";
	out<<"FREQ_MAX         = "<<nlo.maxFreq_<<"\n";
	out<<"FREQ_CUT         = "<<nlo.freqCut_<<"\n";
	out<<"FREQ_RES         = "<<nlo.freqRes_<<"\n";
	out<<"WINDOW           = "<<nlo.windowType_<<"\n";
	out<<"STRIDE_CHARGE    = "<<nlo.strideChg_<<"\n";
	out<<"STRIDE_ALPHA     = "<<nlo.strideAlpha_<<"\n";
	out<<"CALC_CHG         = "<<nlo.calcChg_<<"\n";
	out<<"CALC_ALPHA       = "<<nlo.calcAlpha_<<"\n";
	out<<"CALC_SPECTRUM    = "<<nlo.calcSpectrum_<<"\n";
	out<<"WRITE_ALPHA_TOT  = "<<nlo.writeAlphaT_<<"\n";
	out<<"WRITE_DIPOLE_TOT = "<<nlo.writeDipoleT_<<"\n";
	out<<"WRITE_CHG_TOT    = "<<nlo.writeChgT_<<"\n";
	out<<"WRITE_CHG_ATOM   = "<<nlo.writeChgA_<<"\n";
	out<<"WRITE_ALPHA_ATOM = "<<nlo.writeAlphaA_<<"\n";
	out<<"READ_CHG_ATOM    = "<<nlo.readChgA_<<"\n";
	out<<"READ_ALPHA_ATOM  = "<<nlo.readAlphaA_<<"\n";
	out<<"NORMALIZE        = "<<nlo.normalize_<<"\n";
	out<<"PROFILE_CALC     = "<<nlo.profileCalc_<<"\n";
	out<<"FILE_SPECTRUM    = "<<nlo.fileSpectrum_<<"\n";
	out<<"FILE_DIPOLE_TOT  = "<<nlo.fileDipoleT_<<"\n";
	out<<"FILE_ALPHA_TOT   = "<<nlo.fileAlphaT_<<"\n";
	out<<"FILE_CHG_TOT     = "<<nlo.fileChgT_<<"\n";
	out<<"FILE_ALPHA_ATOM  = "<<nlo.fileAlphaA_<<"\n";
	out<<"FILE_CHG_ATOM    = "<<nlo.fileChgA_<<"\n";
	out<<"******************* NLO_ATOM *******************\n";
	out<<"*************************************************";
	return out;
}

void NLO::defaults(){
	if(DEBUG_NLO_ATOM>0) std::cout<<"defaults()\n";
	//fft parameters
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
		strideChg_=1;
		nStepsChg_=0;
		strideAlpha_=1;
		nStepsAlpha_=0;
	//i/o
		fileSpectrum_=std::string("nlo.dat");
		fileDipoleT_=std::string("dipole_tot.dat");
		fileAlphaT_=std::string("alpha_tot.dat");
		fileChgT_=std::string("chg_tot.dat");
		fileAlphaA_=std::string("alpha_atom.dat");
		fileChgA_=std::string("chg_atom.dat");
	//calculation flags
		calcChg_=true;
		calcAlpha_=true;
		calcSpectrum_=true;
		normalize_=true;
	//i/o flags
		writeAlphaT_=false;
		writeDipoleT_=false;
		writeChgT_=false;
		writeChgA_=false;
		readAlphaA_=false;
		readChgA_=false;
	//chi2 spectrum
		chi2_.clear();
}

void NLO::read(const char* file){
	if(DEBUG_NLO_ATOM>0) std::cout<<"read(const char*):\n";
	//local function variables
	FILE* reader=NULL;
	char* input=new char[string::M];
	std::vector<std::string> strlist;
	bool error=false;
	
	try{
		//open the parameter file
		if(DEBUG_NLO_ATOM>1) std::cout<<"Opening parameter file...\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open parameter file.");
		
		//set the defaults
		if(DEBUG_NLO_ATOM>1) std::cout<<"Setting defaults...\n";
		defaults();
		
		//read in the parameters
		if(DEBUG_NLO_ATOM>1) std::cout<<"Reading in parameters...\n";
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);
			if(string::split(input,string::WS,strlist)==0) continue;
			string::to_upper(strlist.at(0));
			if(strlist.at(0)=="T"){
				T_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="FILE_SPECTRUM"){
				fileSpectrum_=strlist.at(1);
			} else if(strlist.at(0)=="FILE_DIPOLE_TOT"){
				fileDipoleT_=strlist.at(1);
			} else if(strlist.at(0)=="FILE_ALPHA_TOT"){
				fileAlphaT_=strlist.at(1);
			} else if(strlist.at(0)=="FILE_CHG_TOT"){
				fileChgT_=strlist.at(1);
			} else if(strlist.at(0)=="FILE_ALPHA_ATOM"){
				fileAlphaA_=strlist.at(1);
			} else if(strlist.at(0)=="FILE_CHG_ATOM"){
				fileChgA_=strlist.at(1);
			} else if(strlist.at(0)=="CALC_SPECTRUM"){
				calcSpectrum_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="CALC_ALPHA"){
				calcAlpha_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="CALC_CHG"){
				calcChg_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="WRITE_DIPOLE_TOT"){
				writeDipoleT_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="WRITE_ALPHA_TOT"){
				writeAlphaT_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="WRITE_ALPHA_ATOM"){
				writeAlphaT_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="WRITE_CHG_TOT"){
				writeChgT_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="WRITE_CHG_ATOM"){
				writeChgA_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="READ_CHG_ATOM"){
				readChgA_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="READ_ALPHA_ATOM"){
				readAlphaA_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="NORMALIZE"){
				normalize_=string::boolean(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="FREQ_CUT"){
				freqCut_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="RESOLUTION"){
				freqRes_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="STRIDE_CHARGE"){
				strideChg_=std::atoi(strlist.at(1).c_str());
			} else if(strlist.at(0)=="STRIDE_ALPHA"){
				strideAlpha_=std::atoi(strlist.at(1).c_str());
			} else if(strlist.at(0)=="WINDOW"){
				windowType_=window::WINDOW_FUNC::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="PROFILE_CALC"){
				profileCalc_.sigma()=std::atof(strlist.at(1).c_str());
				profileCalc_.b1()=std::atof(strlist.at(2).c_str());
				profileCalc_.b2()=std::atof(strlist.at(3).c_str());
			} 
		}
		
		//close the parameter file
		if(DEBUG_NLO_ATOM>1) std::cout<<"Closing parameter file...\n";
		fclose(reader);
		reader=NULL;
		
		//check the calculation parameters
		if(DEBUG_NLO_ATOM>1) std::cout<<"Checking parameters...\n";
		if(strideAlpha_==0 || strideChg_==0) throw std::invalid_argument("Invalid stride.");
	}catch(std::exception& e){
		std::cout<<"ERROR in read(const char*):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	delete[] input;
	
	if(error) throw std::runtime_error("ERROR in AlphaEff::read(const char*): Failed to read.");
}

void NLO::init(Simulation& sim){
	if(DEBUG_NLO_ATOM>0) std::cout<<"init(const Simulation&):\n";
	bool error=false;
	
	try{
		//set the fft parameters
		if(DEBUG_NLO_ATOM>0) std::cout<<"Assigning the FFT parameters...\n";
		nStepsChg_=sim.timesteps()/strideChg_;
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
		
		//allocate space for chi2
		if(DEBUG_NLO_ATOM>0) std::cout<<"Allocating space for chi2...\n";
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
		std::cout<<"ERROR in read(const char*):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
}

void NLO::calcAlpha(Simulation& sim, const Thole& thole, const Ewald3D::Dipole& ewald){
	if(DEBUG_NLO_ATOM>0) std::cout<<"calcAlpha(Simulation&,const Thole&,const Ewald3D::Dipole&):\n";
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
		std::vector<Ewald3D::Dipole> ewald_(nThreads);
		
	//initialize the effective alphas
	if(DEBUG_NLO_ATOM>1) std::cout<<"Initializing Thole objects...\n";
	for(unsigned int i=0; i<nThreads; ++i) thole_[i]=thole;
	for(unsigned int i=0; i<nThreads; ++i) ewald_[i]=ewald;
		
	//calculate the atomic alphas
	if(DEBUG_NLO_ATOM>1) std::cout<<"Calculating Atomic Alphas...\n";
	start=std::clock();
	if(!sim.cell_fixed()){
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideAlpha_){
			unsigned int TN=0;
			#ifdef _OPENMP
				TN=omp_get_thread_num();
			#endif
			if(DEBUG_NLO_ATOM>-1) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_NLO_ATOM>0 && t%1000==0) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
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
			if(DEBUG_NLO_ATOM>-1) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_NLO_ATOM>0 && t%1000==0) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			thole_[TN].alpha_cont(sim.frame(t),ewald_[TN]);
		}
	}
	stop=std::clock();
	time=((double)(stop-start))/CLOCKS_PER_SEC;
	std::cout<<"alpha-time = "<<time<<"\n";
	
	//interpolate the polarizabilities
	if(strideAlpha_>1){
		if(DEBUG_NLO_ATOM>0) std::cout<<"Interpolating Effective Alphas...\n";
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

void NLO::calcChg(Simulation& sim, const QEQ& qeq, const Ewald3D::Coulomb& ewald){
	if(DEBUG_NLO_ATOM>0) std::cout<<"calcChg(Simulation&,const QEQ&,const Ewald3D::Coulomb&):\n";
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
		std::vector<QEQ> qeq_(nThreads);
		std::vector<Ewald3D::Coulomb> ewald_(nThreads);
		
	//initialize the effective alphas
	if(DEBUG_NLO_ATOM>1) std::cout<<"initializing QEQ objects...\n";
	for(unsigned int i=0; i<nThreads; ++i) qeq_[i]=qeq;
	for(unsigned int i=0; i<nThreads; ++i) ewald_[i]=ewald;
	
	//calculate the atomic charges
	if(DEBUG_NLO_ATOM>1) std::cout<<"computing atomic charges...\n";
	start=std::clock();
	if(sim.cell_fixed()){
		if(DEBUG_NLO_ATOM>1) std::cout<<"cell_fixed true\n";
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideChg_){
			unsigned int TN=0;
			#ifdef _OPENMP
				TN=omp_get_thread_num();
			#endif
			if(DEBUG_NLO_ATOM>-1) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_NLO_ATOM>0 && t%1000==0) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			qeq_[TN].qt_jZero_cont(sim.frame(t),ewald_[TN]);
		}
	} else {
		if(DEBUG_NLO_ATOM>1) std::cout<<"cell_fixed false\n";
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideChg_){
			unsigned int TN=0;
			#ifdef _OPENMP
				TN=omp_get_thread_num();
			#endif
			if(DEBUG_NLO_ATOM>-1) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_NLO_ATOM>0 && t%1000==0) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			ewald_[TN].init(sim.frame(t),1e-5);
			qeq_[TN].qt_jZero(sim.frame(t),ewald_[TN]);
		}
	}
	stop=std::clock();
	time=((double)(stop-start))/CLOCKS_PER_SEC;
	std::cout<<"chg-time = "<<time<<"\n";
	
	if(strideChg_>1){
		if(DEBUG_NLO_ATOM>0) std::cout<<"Interpolating atomic charge...\n";
		Interp::Data chgData(nStepsChg_+1);
		//loop over all atoms
		for(int n=0; n<sim.frame(0).nAtoms(); n++){
			//loop over all timesteps for which we calculated atomic charges
			for(int t=0; t*strideChg_<sim.timesteps(); t++){
				chgData.x(t)=t*strideChg_;
				chgData.y(t)=sim.frame(t*strideChg_).charge(n);
			}
			//interpolate the rest of the timesteps
			for(unsigned int t=0; t*strideChg_<sim.timesteps()-strideChg_; ++t){
				for(unsigned int tp=1; tp<strideChg_; ++tp){
					sim.frame(t*strideChg_+tp).charge(n)=Interp::interpAkima(t*strideChg_+tp,chgData);
				}
			}
		}
		//fix the final timesteps
		double m;
		unsigned int end=nStepsChg_*strideChg_;
		//loop over all molecules
		for(unsigned int n=0; n<sim.frame(0).nAtoms(); ++n){
			//calculate the slope at the end
			m=1.0/sim.timestep()*(sim.frame(end).charge(n)-sim.frame(end-1).charge(n));
			//perform linear extrapolation for the last few timesteps
			for(unsigned int t=end; t<sim.timesteps(); ++t){
				sim.frame(t).charge(n)=sim.frame(end).charge(n)+m*(t-end)*sim.timestep();
			}
		}
	}
}

void NLO::calcSpectrum(Simulation& sim){
	if(DEBUG_NLO_ATOM>0) std::cout<<"calcSpectrum(Simulation&):\n";
	//==== constants ====
	double hbar=0,kb=0;
	if(units::consts::system()==units::System::AU){
		hbar=units::au::hbar;
		kb=units::au::kb;
	} else if(units::consts::system()==units::System::METAL){
		hbar=units::metal::hbar;
		kb=units::metal::kb;
	}
	//==== local function variables ====
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
	//==== transpose of alpha delta ====
	std::vector<std::vector<std::vector<double> > > alphaTrans(3);
	for(unsigned int i=0; i<3; ++i) alphaTrans[i].resize(3,std::vector<double>(sim.timesteps(),0));
	//==== dipole utilities ====
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > qdotr(sim.timesteps(),Eigen::Vector3d::Zero());
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > rdotq(sim.timesteps(),Eigen::Vector3d::Zero());
	std::vector<double> q(sim.timesteps());
	std::vector<double> qv(sim.timesteps());
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > p(sim.timesteps(),Eigen::Vector3d::Zero());
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > v(sim.timesteps(),Eigen::Vector3d::Zero());
	
	//==== no truncation of correlation function ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"No truncation of correlation function.\n";
		
	//==== compute qdotr ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"computing qdotr...\n";
	for(unsigned int n=0; n<sim.frame(0).nAtoms(); ++n){
		unsigned int tt;
		double ts=sim.timestep();
		//read the charge
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			q[t]=sim.frame(t).charge(n);
		}
		//compute the charge velocity
		tt=0; qv[tt]=(q[tt+1]-q[tt])/ts;
		tt=1; qv[tt]=0.5*(q[tt+1]-q[tt-1])/ts;
		tt=2; qv[tt]=(-1.0/12.0*(q[tt+2]-q[tt-2])+2.0/3.0*(q[tt+1]-q[tt-1]))/ts;
		tt=3; qv[tt]=(1.0/60.0*(q[tt+3]-q[tt-3])-3.0/20.0*(q[tt+2]-q[tt-2])+3.0/4.0*(q[tt+1]-q[tt-1]))/ts;
		for(unsigned int t=4; t<sim.timesteps()-4; ++t){
			qv[t]=(-1.0/280.0*(q[t+4]-q[t-4])+4.0/105.0*(q[t+3]-q[t-3])-1.0/5.0*(q[t+2]-q[t-2])+4.0/5.0*(q[t+1]-q[t-1]))/ts;
		}
		tt=sim.timesteps()-4; qv[tt]=(1.0/60.0*(q[tt+3]-q[tt-3])-3.0/20.0*(q[tt+2]-q[tt-2])+3.0/4.0*(q[tt+1]-q[tt-1]))/ts;
		tt=sim.timesteps()-3; qv[tt]=(-1.0/12.0*(q[tt+2]-q[tt-2])+2.0/3.0*(q[tt+1]-q[tt-1]))/ts;
		tt=sim.timesteps()-2; qv[tt]=0.5*(q[tt+1]-q[tt-1])/ts;
		tt=sim.timesteps()-1; qv[tt]=(q[tt]-q[tt-1])/ts;
		//compute qdotr
		if(profileCalc_.sigma()>0){
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				qdotr[t].noalias()+=qv[t]*sim.frame(t).posn(n)*profileCalc_(sim.frame(t).posn(n)[2]);
			}
		} else {
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				qdotr[t].noalias()+=qv[t]*sim.frame(t).posn(n);
			}
		}
	}
	
	//==== compute qdotq ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"computing rdotq...\n";
	for(unsigned int n=0; n<sim.frame(0).nAtoms(); ++n){
		unsigned int tt;
		double ts=sim.timestep();
		//read the position - fractional coordinates
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			p[t].noalias()=sim.frame(t).cell().RInv()*sim.frame(t).posn(n);
		}
		//compute the velocities - fractional coordinates
		tt=0;
		v[tt][0]=special::mod(p[tt+1][0]-p[tt][0],-0.5,0.5)/ts;
		v[tt][1]=special::mod(p[tt+1][1]-p[tt][1],-0.5,0.5)/ts;
		v[tt][2]=special::mod(p[tt+1][2]-p[tt][2],-0.5,0.5)/ts;
		tt=1;
		v[tt][0]=0.5*special::mod(p[tt+1][0]-p[tt-1][0],-0.5,0.5)/ts;
		v[tt][1]=0.5*special::mod(p[tt+1][1]-p[tt-1][1],-0.5,0.5)/ts;
		v[tt][2]=0.5*special::mod(p[tt+1][2]-p[tt-1][2],-0.5,0.5)/ts;
		tt=2;
		v[tt][0]=(-1.0/12.0*special::mod(p[tt+2][0]-p[tt-2][0],-0.5,0.5)+2.0/3.0*special::mod(p[tt+1][0]-p[tt-1][0],-0.5,0.5))/ts;
		v[tt][1]=(-1.0/12.0*special::mod(p[tt+2][1]-p[tt-2][1],-0.5,0.5)+2.0/3.0*special::mod(p[tt+1][1]-p[tt-1][1],-0.5,0.5))/ts;
		v[tt][2]=(-1.0/12.0*special::mod(p[tt+2][2]-p[tt-2][2],-0.5,0.5)+2.0/3.0*special::mod(p[tt+1][2]-p[tt-1][2],-0.5,0.5))/ts;
		tt=3;
		v[tt][0]=(1.0/60.0*special::mod(p[tt+3][0]-p[tt-3][0],-0.5,0.5)-3.0/20.0*special::mod(p[tt+2][0]-p[tt-2][0],-0.5,0.5)+3.0/4.0*special::mod(p[tt+1][0]-p[tt-1][0],-0.5,0.5))/ts;
		v[tt][1]=(1.0/60.0*special::mod(p[tt+3][1]-p[tt-3][1],-0.5,0.5)-3.0/20.0*special::mod(p[tt+2][1]-p[tt-2][1],-0.5,0.5)+3.0/4.0*special::mod(p[tt+1][1]-p[tt-1][1],-0.5,0.5))/ts;
		v[tt][2]=(1.0/60.0*special::mod(p[tt+3][2]-p[tt-3][2],-0.5,0.5)-3.0/20.0*special::mod(p[tt+2][2]-p[tt-2][2],-0.5,0.5)+3.0/4.0*special::mod(p[tt+1][2]-p[tt-1][2],-0.5,0.5))/ts;
		for(unsigned int t=4; t<sim.timesteps()-4; ++t){
			v[t][0]=(-1.0/280.0*special::mod(p[t+4][0]-p[t-4][0],-0.5,0.5)+4.0/105.0*special::mod(p[t+3][0]-p[t-3][0],-0.5,0.5)-1.0/5.0*special::mod(p[t+2][0]-p[t-2][0],-0.5,0.5)+4.0/5.0*special::mod(p[t+1][0]-p[t-1][0],-0.5,0.5))/ts;
			v[t][1]=(-1.0/280.0*special::mod(p[t+4][1]-p[t-4][1],-0.5,0.5)+4.0/105.0*special::mod(p[t+3][1]-p[t-3][1],-0.5,0.5)-1.0/5.0*special::mod(p[t+2][1]-p[t-2][1],-0.5,0.5)+4.0/5.0*special::mod(p[t+1][1]-p[t-1][1],-0.5,0.5))/ts;
			v[t][2]=(-1.0/280.0*special::mod(p[t+4][2]-p[t-4][2],-0.5,0.5)+4.0/105.0*special::mod(p[t+3][2]-p[t-3][2],-0.5,0.5)-1.0/5.0*special::mod(p[t+2][2]-p[t-2][2],-0.5,0.5)+4.0/5.0*special::mod(p[t+1][2]-p[t-1][2],-0.5,0.5))/ts;
		}
		tt=sim.timesteps()-4;
		v[tt][0]=(1.0/60.0*special::mod(p[tt+3][0]-p[tt-3][0],-0.5,0.5)-3.0/20.0*special::mod(p[tt+2][0]-p[tt-2][0],-0.5,0.5)+3.0/4.0*special::mod(p[tt+1][0]-p[tt-1][0],-0.5,0.5))/ts;
		v[tt][1]=(1.0/60.0*special::mod(p[tt+3][1]-p[tt-3][1],-0.5,0.5)-3.0/20.0*special::mod(p[tt+2][1]-p[tt-2][1],-0.5,0.5)+3.0/4.0*special::mod(p[tt+1][1]-p[tt-1][1],-0.5,0.5))/ts;
		v[tt][2]=(1.0/60.0*special::mod(p[tt+3][2]-p[tt-3][2],-0.5,0.5)-3.0/20.0*special::mod(p[tt+2][2]-p[tt-2][2],-0.5,0.5)+3.0/4.0*special::mod(p[tt+1][2]-p[tt-1][2],-0.5,0.5))/ts;
		tt=sim.timesteps()-3;
		v[tt][0]=(-1.0/12.0*special::mod(p[tt+2][0]-p[tt-2][0],-0.5,0.5)+2.0/3.0*special::mod(p[tt+1][0]-p[tt-1][0],-0.5,0.5))/ts;
		v[tt][1]=(-1.0/12.0*special::mod(p[tt+2][1]-p[tt-2][1],-0.5,0.5)+2.0/3.0*special::mod(p[tt+1][1]-p[tt-1][1],-0.5,0.5))/ts;
		v[tt][2]=(-1.0/12.0*special::mod(p[tt+2][2]-p[tt-2][2],-0.5,0.5)+2.0/3.0*special::mod(p[tt+1][2]-p[tt-1][2],-0.5,0.5))/ts;
		tt=sim.timesteps()-2;
		v[tt][0]=0.5*special::mod(p[tt+1][0]-p[tt-1][0],-0.5,0.5)/ts;
		v[tt][1]=0.5*special::mod(p[tt+1][1]-p[tt-1][1],-0.5,0.5)/ts;
		v[tt][2]=0.5*special::mod(p[tt+1][2]-p[tt-1][2],-0.5,0.5)/ts;
		tt=sim.timesteps()-1;
		v[tt][0]=special::mod(p[tt][0]-p[tt-1][0],-0.5,0.5)/ts;
		v[tt][1]=special::mod(p[tt][1]-p[tt-1][1],-0.5,0.5)/ts;
		v[tt][2]=special::mod(p[tt][2]-p[tt-1][2],-0.5,0.5)/ts;
		//transform the velocities from fractional to Cartesian
		for(unsigned int t=0; t<sim.timesteps(); ++t) v[t]=sim.frame(t).cell().R()*v[t];
		//compute rdotq
		if(profileCalc_.sigma()>0){
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				rdotq[t].noalias()+=v[t]*sim.frame(t).charge(n)*profileCalc_(sim.frame(t).posn(n)[2]);
			}
		} else {
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				rdotq[t].noalias()+=v[t]*sim.frame(t).charge(n);
			}
		}
	}
	
	//==== sum qdotr and rdotq to get the total dipole moment time derivative ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"computing dipole time derivative (sum of qdotr and rdotq)...\n";
	for(unsigned int t=0; t<sim.timesteps(); ++t){
		for(unsigned int i=0; i<3; ++i){
			fftMu[i].in(t)=qdotr[t][i]+rdotq[t][i];
		}
	}
	/*FILE* reader=fopen("dipole-v.dat","w");
	if(reader!=NULL){
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			fprintf(reader,"%i %f %f %f\n",t,fftMu[0].in(t),fftMu[1].in(t),fftMu[2].in(t));
		}
		fclose(reader);
		reader=NULL;
	}*/
		
	//==== compute the Fourier transform ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"Calculating the forward dipole transform...\n";
	for(unsigned int i=0; i<3; ++i) fftMu[i].transformf();
	
	//==== compute the transpose of the polarizability ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"Calculating the transpose of the total polarizability...\n";
	if(profileCalc_.sigma()>0){
		if(DEBUG_NLO_ATOM>0) std::cout<<"Applying profile...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				double mean=0.0;
				for(unsigned int i=0; i<3; ++i){
					for(unsigned int j=0; j<3; ++j){
						alphaTrans[i][j][t]+=sim.frame(t).alpha(n)(i,j)*profileCalc_(sim.frame(t).posn(n)[2]);
					}
				}
			}
		}
	} else {
		if(DEBUG_NLO_ATOM>0) std::cout<<"No profile...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				double mean=0.0;
				for(unsigned int i=0; i<3; ++i){
					for(unsigned int j=0; j<3; ++j){
						alphaTrans[i][j][t]+=sim.frame(t).alpha(n)(i,j);
					}
				}
			}
		}
	}
	
	//==== compute the polarizability velocity ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"Calculating the polarizability velocity...\n";
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
	
	//==== tranform into the frequency domain ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"Calculating the forward polarizability transform...\n";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			fftAlpha[i][j].transformf();
		}
	}
	
	//==== compute the correlation function in the frequency domain ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"Calculating the correlation function in frequency space...\n";
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
	
	//==== compute the corerelation function in time-space, then transform to find chi2 ====
	if(DEBUG_NLO_ATOM>0) std::cout<<"Calculating chi2...\n";
	double norm=1.0/(2.0*sim.timesteps());
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			for(unsigned int k=0; k<3; ++k){
				//compute the time-space correlation function
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
				//compute the forward transform
				fft[i][j][k].transformf();
				//record chi2
				for(unsigned int t=0; t<freqN_; ++t){
					const double w=minFreq_*(t+0.5);
					//chi2_[i][j][k][t][0]=1.0/w*norm*fft[i][j][k].out(t)[1]*(-1.0);
					//chi2_[i][j][k][t][1]=1.0/w*norm*fft[i][j][k].out(t)[0];
					chi2_[i][j][k][t][0]=norm*fft[i][j][k].out(t)[1]*(-1.0);
					chi2_[i][j][k][t][1]=norm*fft[i][j][k].out(t)[0];
				}
			}
		}
	}
}

void NLO::printSpectrum(const char* file) const{
	if(DEBUG_NLO_ATOM>0) std::cout<<"printSpectrum():\n";
	//local function variables
	FILE* writer=NULL;
	
	writer=fopen(file,"w");
	if(writer==NULL) throw std::runtime_error("I/O Exception Occured.");
	
	//print the header
	fprintf(writer, "#Freq(THz) ");
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
	
	//======== local variables ========
	//simulation
		double timestep=0.5;
		Interval interval;
		Eigen::Vector3d offset=Eigen::Vector3d::Zero();
		AtomType atomT;
		atomT.name=true; atomT.an=true; atomT.specie=true; atomT.index=true;
		atomT.posn=true; atomT.charge=true; atomT.alpha=true; atomT.jzero=true;
		Simulation sim;
	//ewald
		Ewald3D::Dipole dipole;
		Ewald3D::Coulomb coul;
	//calculation
		NLO nlo;
		Thole thole;
		QEQ qeq;
		std::vector<double> jzero;
		std::vector<double> alpha;
		std::vector<std::string> atoms_jzero;
		std::vector<std::string> atoms_alpha;
	//input/output
		std::vector<std::string> strlist;
		char* input=new char[string::M];
		std::string paramfile,simstr;
		FILE* reader=NULL;
		FILE_FORMAT::type fileFormat;
	//flags
		bool error=false;
	//units
		units::System::type unitsys=units::System::METAL;
		units::consts::init(unitsys);
	//charges
		std::vector<std::string> atoms_charge_str;
		std::vector<std::vector<unsigned int> > atoms_charge;
		std::vector<double> charge;
	
	try{
		//======== check number of arguments ========
		if(argc!=2){
			std::cout<<"ERROR in main(int,char**):\n";
			std::cout<<"Incorrect number of arguments.\n";
			std::cout<<"\t1. Program Execution\n";
			std::cout<<"\t2. Parameter File\n";
			throw std::invalid_argument("Invalid command-line arguments.");
		}
		
		//======== open parameter file ========
		std::cout<<"opening the parameter file...\n";
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
			} else if(strlist.at(0)=="INTERVAL"){
				interval=Interval::read(strlist.at(1).c_str());
			} else if(strlist.at(0)=="OFFSET"){
				offset[0]=std::atof(strlist.at(1).c_str());
				offset[1]=std::atof(strlist.at(2).c_str());
				offset[2]=std::atof(strlist.at(3).c_str());
			} else if(strlist.at(0)=="OFFSET"){
				interval=Interval::read(strlist.at(1).c_str());
			} else if(strlist.at(0)=="FORMAT"){
				fileFormat=FILE_FORMAT::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="UNITS"){
				unitsys=units::System::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="JZERO"){
				atoms_jzero.push_back(strlist.at(1));
				jzero.push_back(std::atof(strlist.at(2).c_str()));
			} else if(strlist.at(0)=="ALPHA"){
				atoms_alpha.push_back(strlist.at(1));
				alpha.push_back(std::atof(strlist.at(2).c_str()));
			} else if(strlist.at(0)=="CHARGE"){
				atoms_charge_str.push_back(strlist.at(1));
				charge.push_back(std::atof(strlist.at(2).c_str()));
			} else if(strlist.at(0)=="TIMESTEP"){
				timestep=std::atof(strlist.at(1).c_str());
			}
		}
		fclose(reader);
		reader=NULL;
		
		//======== check the parameters ========
		if(unitsys==units::System::UNKNOWN) throw std::invalid_argument("Invalid unit system");
		
		//======== initialize the unit system ========
		std::cout<<"initializing the unit system...\n";
		units::consts::init(unitsys);
		
		//======== print the parameters to screen ========
		std::cout<<"SIMULATION PARAMETERS:\n";
		std::cout<<"\tUNITS    = "<<unitsys<<"\n";
		std::cout<<"\tSIM      = "<<simstr<<"\n";
		std::cout<<"\tTIMESTEP = "<<timestep<<"\n";
		std::cout<<"\tFORMAT   = "<<fileFormat<<"\n";
		std::cout<<"\tOFFSET   = "<<offset[0]<<","<<offset[1]<<","<<offset[2]<<"\n";
		std::cout<<"\tINTERVAL = "<<interval<<"\n";
		std::cout<<"\tATOMS    = "; for(unsigned int i=0; i<atoms_alpha.size(); ++i) std::cout<<atoms_alpha[i]<<" "; std::cout<<"\n";
		std::cout<<"\tJZERO    = "; for(unsigned int i=0; i<jzero.size(); ++i) std::cout<<jzero[i]<<" "; std::cout<<"\n";
		std::cout<<"\tALPHA    = "; for(unsigned int i=0; i<alpha.size(); ++i) std::cout<<alpha[i]<<" "; std::cout<<"\n";
		
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
		
		//======== set charges ========
		if(atoms_charge_str.size()>0){
			std::cout<<"setting charges...\n";
			for(unsigned int i=0; i<atoms_charge_str.size(); ++i){
				std::cout<<atoms_charge_str[i]<<" "<<charge[i]<<"\n";
			}
			atoms_charge.resize(atoms_charge_str.size());
			for(unsigned int i=0; i<atoms_charge_str.size(); ++i){
				StructureI::read_atoms(atoms_charge_str[i].c_str(),atoms_charge[i],sim.frame(0));
			}
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int i=0; i<atoms_charge.size(); ++i){
					for(unsigned int j=0; j<atoms_charge[i].size(); ++j){
						sim.frame(t).charge(atoms_charge[i][j])=charge[i];
					}
				}
			}
		}
		
		//======== set offset ========
		if(offset.norm()>num_const::ZERO){
			std::cout<<"setting offset...\n";
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
					sim.frame(t).posn(n).noalias()+=offset;
					Cell::returnToCell(sim.frame(t).posn(n),sim.frame(t).posn(n),sim.frame(t).cell().R(),sim.frame(t).cell().RInv());
				}
			}
		}
		
		//======== assign jzero ========
		std::cout<<"assigning jzero...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				for(unsigned int i=0; i<atoms_jzero.size(); ++i){
					if(atoms_jzero[i]==sim.frame(t).name(n)){
						sim.frame(t).jzero(n)=jzero[i];
						break;
					}
				}
			}
		}
		
		//======== assign alpha ========
		std::cout<<"assigning alpha...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				for(unsigned int i=0; i<atoms_alpha.size(); ++i){
					if(atoms_alpha[i]==sim.frame(t).name(n)){
						sim.frame(t).alpha(n)=Eigen::Matrix3d::Identity()*alpha[i];
						break;
					}
				}
			}
		}
		
		//======== print the simulation ========
		sim.timestep()=timestep;
		std::cout<<sim<<"\n";
		std::cout<<sim.frame(0).cell().R()<<"\n";
		
		//======== read auxiliary objects ========
		std::cout<<"reading the auxiliary objects...\n";
		//ewald object
			std::cout<<"Ewald objects...\n";
			dipole.init(sim.frame(0).cell(),1e-5);
			std::cout<<dipole<<"\n";
			coul.init(sim.frame(0),1e-5);
			std::cout<<coul<<"\n";
		//qeq object
			std::cout<<"QEQ object...\n";
			qeq.read(paramfile.c_str());
			qeq.init(sim.frame(0));
			std::cout<<qeq<<"\n";
		//thole object
			std::cout<<"Thole object...\n";
			thole.read(paramfile.c_str());
			thole.init(sim.frame(0),dipole);
			std::cout<<thole<<"\n";
		//nlo object
			std::cout<<"NLO object...\n";
			nlo.read(paramfile.c_str());
			nlo.init(sim);
			std::cout<<nlo<<"\n";
		
		//======== reading chg ========
		if(nlo.readChgA()){
			std::cout<<"reading chg data...\n";
			Utility::Read::chg_atom(nlo.fileChgA().c_str(),sim);
		}
		
		//======== read alpha ========
		if(nlo.readAlphaA()){
			std::cout<<"reading alpha data...\n";
			Utility::Read::alpha_atom(nlo.fileAlphaA().c_str(),sim);
		}
		
		//======== execute calculation ========
		std::cout<<"executing calculation...\n";
		if(nlo.calcAlpha()){
			std::cout<<"computing atomic alphas...\n";
			nlo.calcAlpha(sim,thole,dipole);
		}
		if(nlo.calcChg()){
			std::cout<<"computing atomic charges...\n";
			nlo.calcChg(sim,qeq,coul);
		}
		if(nlo.calcSpectrum()){
			std::cout<<"computing the spectrum...\n";
			nlo.calcSpectrum(sim);
			nlo.printSpectrum(nlo.fileSpectrum().c_str());
		}
		
		//======== write relevant data ========
		std::cout<<"writing relevant data...\n";
		if(nlo.writeDipoleT()) Utility::Write::dipole_tot(nlo.fileDipoleT().c_str(),sim);
		if(nlo.writeAlphaT()) Utility::Write::alpha_tot(nlo.fileAlphaT().c_str(),sim);
		if(nlo.writeAlphaA()) Utility::Write::alpha_atom(nlo.fileAlphaA().c_str(),sim);
		if(nlo.writeChgT()) Utility::Write::chg_tot(nlo.fileChgT().c_str(),sim);
		if(nlo.writeChgA()) Utility::Write::chg_atom(nlo.fileChgA().c_str(),sim);
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
