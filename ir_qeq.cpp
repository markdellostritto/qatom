#include "ir_qeq.hpp"

//***********************************************************************************************************************************
//Profile
//***********************************************************************************************************************************

std::ostream& operator<<(std::ostream& out, const Profile& p){
	return out<<"profile("<<p.sigma()<<","<<p.b1()<<","<<p.b2()<<")";
}

//***********************************************************************************************************************************
//IRQEQ class
//***********************************************************************************************************************************

//constants
const double IRQEQ::mevPerThz=0.24180;
const double IRQEQ::cmiPerThz=0.02998;

std::ostream& operator<<(std::ostream& out, const IRQEQ& ir3d){
	out<<"*************************************************\n";
	out<<"******************** IR_QEQ ********************\n";
	out<<"FREQ_UNIT        = "<<ir3d.freqUnit_<<"\n";
	out<<"FREQ_MIN         = "<<ir3d.minFreq_<<"\n";
	out<<"FREQ_MAX         = "<<ir3d.maxFreq_<<"\n";
	out<<"FREQ_CUT         = "<<ir3d.freqCut_<<"\n";
	out<<"FREQ_RES         = "<<ir3d.freqRes_<<"\n";
	out<<"WINDOW           = "<<ir3d.windowType_<<"\n";
	out<<"STRIDE_CHG       = "<<ir3d.strideChg_<<"\n";
	out<<"CALC_CHG         = "<<ir3d.calcChg_<<"\n";
	out<<"CALC_SPECTRUM    = "<<ir3d.calcSpectrum_<<"\n";
	out<<"WRITE_DIPOLE_TOT = "<<ir3d.writeDipoleT_<<"\n";
	out<<"WRITE_CHG_TOT    = "<<ir3d.writeChgT_<<"\n";
	out<<"WRITE_CHG_ATOM   = "<<ir3d.writeChgA_<<"\n";
	out<<"READ_CHG_ATOM    = "<<ir3d.readChgA_<<"\n";
	out<<"NORMALIZE        = "<<ir3d.normalize_<<"\n";
	out<<"PROFILE_CALC     = "<<ir3d.profile_<<"\n";
	out<<"FILE_SPECTRUM    = "<<ir3d.fileSpectrum_<<"\n";
	out<<"FILE_DIPOLE_TOT  = "<<ir3d.fileDipoleT_<<"\n";
	out<<"FILE_CHG_TOT     = "<<ir3d.fileChgT_<<"\n";
	out<<"FILE_CHG_ATOM    = "<<ir3d.fileChgA_<<"\n";
	out<<"T                = "<<ir3d.T_<<"\n";
	if(ir3d.subset_.size()>0){
		out<<"SUBSET           = ";
		for(unsigned int i=0; i<ir3d.subset_.size(); ++i) std::cout<<ir3d.subset_[i]<<" ";
		std::cout<<"\n";
	}
	out<<"******************** IR_QEQ ********************\n";
	out<<"*************************************************";
	return out;
}

void IRQEQ::defaults(){
	if(DEBUG_IR_QEQ>0) std::cout<<"defaults()\n";
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
		strideChg_=1;//include every timestep
		nStepsChg_=0;
	//i/o
		fileSpectrum_=std::string("ir3d_qeq3.dat");
		fileDipoleT_=std::string("dipole_tot.dat");
		fileChgT_=std::string("chg_tot.dat");
		fileChgA_=std::string("chg_atom.dat");
	//calculation flags
		calcChg_=true;
		calcSpectrum_=true;
		normalize_=true;
	//i/o flags
		writeDipoleT_=false;
		writeChgT_=false;
		writeChgA_=false;
		readChgA_=false;
	//temp
		T_=300;
	//ir spectrum
		irSpectrum_.clear();
	//subset
		subsetstr_.clear();
		subset_.clear();
}

void IRQEQ::read(const char* file){
	if(DEBUG_IR_QEQ>0) std::cout<<"read(const char*):\n";
	//local function variables
	FILE* reader=NULL;
	char* input=new char[string::M];
	std::vector<std::string> strlist;
	bool error=false;
	
	try{
		//open the parameter file
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open parameter file.");
		
		//set the defaults
		defaults();
		
		//read in the parameters
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);
			if(string::split(input,string::WS,strlist)==0) continue;
			string::to_upper(strlist.at(0));
			if(strlist.at(0)=="FILE_SPECTRUM"){
				fileSpectrum_=strlist.at(1);
			} else if(strlist.at(0)=="FILE_DIPOLE_TOT"){
				fileDipoleT_=strlist.at(1);
			} else if(strlist.at(0)=="FILE_CHG_TOT"){
				fileChgT_=strlist.at(1);
			} else if(strlist.at(0)=="FILE_CHG_ATOM"){
				fileChgA_=strlist.at(1);
			} else if(strlist.at(0)=="CALC_SPECTRUM"){
				calcSpectrum_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="CALC_CHARGE"){
				calcChg_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="WRITE_DIPOLE_TOT"){
				writeDipoleT_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="WRITE_CHG_TOT"){
				writeChgT_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="WRITE_CHG_ATOM"){
				writeChgA_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="READ_CHG_ATOM"){
				readChgA_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="NORMALIZE"){
				normalize_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="FREQ_CUT"){
				freqCut_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="FREQ_UNIT"){
				freqUnit_=fourier::FreqUnit::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="RESOLUTION"){
				freqRes_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="STRIDE_CHARGE"){
				strideChg_=std::atoi(strlist.at(1).c_str());
			} else if(strlist.at(0)=="WINDOW"){
				windowType_=window::WINDOW_FUNC::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="T"){
				T_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="T"){
				subsetstr_=strlist.at(1);
			} else if(strlist.at(0)=="PROFILE_CALC"){
				profile_.sigma()=std::atof(strlist.at(1).c_str());
				profile_.b1()=std::atof(strlist.at(2).c_str());
				profile_.b2()=std::atof(strlist.at(3).c_str());
			}
		}
		//close the parameter file
		fclose(reader);
		reader=NULL;
		
		//check the calculation parameters
		if(strideChg_==0) throw std::invalid_argument("Invalid stride.");
		if(freqUnit_==fourier::FreqUnit::UNKNOWN) throw std::invalid_argument("Invalid frequency unit.");
	}catch(std::exception& e){
		std::cout<<"ERROR in read(const char*):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	delete[] input;
	
	if(error) throw std::runtime_error("ERROR in AlphaEff::read(const char*): Failed to read.");
}

void IRQEQ::init(Simulation& sim){
	if(DEBUG_IR_QEQ>0) std::cout<<"init(const Simulation&):\n";
	bool error=false;
	
	try{
		//set the fft parameters
		if(DEBUG_IR_QEQ>0) std::cout<<"Assigning the FFT parameters...\n";
		nStepsChg_=sim.timesteps()/strideChg_;
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
		/*else if(windowType_==window::WINDOW_FUNC::KAISERBESSEL){
			sigma_=freqRes_*freqRes_;//scale to get sensible values
			window_=window::KaiserBessel(sim.timesteps(),sigma_);
		}*/
		
		//check the fft parameters
		if(freqCut_-minFreq_<num_const::ZERO) throw std::invalid_argument("Invalid freqency cutoff.");
		
		//allocate space for ir3d
		if(DEBUG_IR_QEQ>0) std::cout<<"Allocating space for ir spectrum...\n";
		irSpectrum_.resize(freqN_);
		
		//read the subset
		if(DEBUG_NLO_ATOM>0) std::cout<<"Reading the subset...\n";
		if(subsetstr_.size()>0) Structure::read_atoms(subsetstr_.c_str(),subset_,sim.frame(0));
		
	}catch(std::exception& e){
		std::cout<<"ERROR in read(const char*):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
}

void IRQEQ::calcChg(Simulation& sim, const QEQ& qeq, const Ewald3D::Coulomb& ewald){
	if(DEBUG_IR_QEQ>0) std::cout<<"calcChg(Simulation&,const QEQ&,const Ewald3D::Coulomb&):\n";
	//local function variables
	//parallel
		#ifdef _OPENMP
			const unsigned int nThreads=omp_get_max_threads();
		#else
			const unsigned int nThreads=1;
		#endif
	//atomic charge
		std::vector<QEQ> qeq_(nThreads);
		std::vector<Ewald3D::Coulomb> ewald_(nThreads);
		
	//initialize the effective alphas
	if(DEBUG_IR_QEQ>1) std::cout<<"initializing QEQ objects...\n";
	for(unsigned int i=0; i<nThreads; ++i) qeq_[i]=qeq;
	for(unsigned int i=0; i<nThreads; ++i) ewald_[i]=ewald;
	
	//calculate the atomic charges
	if(DEBUG_IR_QEQ>1) std::cout<<"computing atomic charges...\n";
	const clock_t start=std::clock();
	if(sim.cell_fixed()){
		if(DEBUG_IR_QEQ>1) std::cout<<"cell_fixed true\n";
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideChg_){
			#ifdef _OPENMP
				const unsigned int TN=omp_get_thread_num();
			#else
				const unsigned int TN=0;
			#endif
			if(DEBUG_IR_QEQ>-1) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_IR_QEQ>0 && t%1000==0) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			qeq_[TN].qt_jZero_cont(sim.frame(t),ewald_[TN]);
		}
	} else {
		if(DEBUG_IR_QEQ>1) std::cout<<"cell_fixed false\n";
		#pragma omp parallel for schedule(static) num_threads(omp_get_max_threads())
		for(unsigned int t=0; t<sim.timesteps(); t+=strideChg_){
			#ifdef _OPENMP
				const unsigned int TN=omp_get_thread_num();
			#else
				const unsigned int TN=0;
			#endif
			if(DEBUG_IR_QEQ>-1) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			else if(DEBUG_IR_QEQ>0 && t%1000==0) std::cout<<"Timestep: "<<sim.beg()+1+t<<"\n";
			ewald_[TN].init(sim.frame(t),1e-5);
			qeq_[TN].qt_jZero(sim.frame(t),ewald_[TN]);
		}
	}
	const clock_t stop=std::clock();
	const double time=((double)(stop-start))/CLOCKS_PER_SEC;
	std::cout<<"chg-time = "<<time<<"\n";
	
	if(strideChg_>1){
		if(DEBUG_IR_QEQ>0) std::cout<<"Interpolating atomic charge...\n";
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

void IRQEQ::calcSpectrum(Simulation& sim){
	if(DEBUG_IR_QEQ>0) std::cout<<"calcSpectrum(Simulation&):\n";
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
	fourier::FFT_C2C fft(2*sim.timesteps());//the FFT of the dipole spectrum
	for(unsigned int i=0; i<3; ++i) fftMu[i]=fourier::FFT_R2C(2*sim.timesteps());
	std::vector<std::vector<double> > dipoleTrans(3,std::vector<double>(sim.timesteps(),0));
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > qdotr(sim.timesteps(),Eigen::Vector3d::Zero());
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > rdotq(sim.timesteps(),Eigen::Vector3d::Zero());
	std::vector<double> q(sim.timesteps());
	std::vector<double> qv(sim.timesteps());
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > p(sim.timesteps(),Eigen::Vector3d::Zero());
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > v(sim.timesteps(),Eigen::Vector3d::Zero());
	//==== atoms ====
	std::vector<unsigned int> atoms;
	if(subset_.size()>0) atoms=subset_;
	else {
		atoms.resize(sim.frame(0).nAtoms());
		for(unsigned int i=0; i<atoms.size(); ++i) atoms[i]=i;
	}
	
	//compute qdotr
	if(DEBUG_IR_QEQ>0) std::cout<<"computing qdotr...\n";
	for(unsigned int n=0; n<atoms.size(); ++n){
		unsigned int tt;
		const double ts=sim.timestep();
		//read the charge
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			q[t]=sim.frame(t).charge(atoms[n]);
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
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			qdotr[t].noalias()+=qv[t]*sim.frame(t).posn(atoms[n]);
		}
	}
	
	//compute qdotq
	if(DEBUG_IR_QEQ>0) std::cout<<"computing rdotq...\n";
	for(unsigned int n=0; n<atoms.size(); ++n){
		unsigned int tt;
		const double ts=sim.timestep();
		//read the position - fractional coordinates
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			p[t].noalias()=sim.frame(t).cell().RInv()*sim.frame(t).posn(atoms[n]);
		}
		//calculate the velocities - fractional coordinates
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
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			rdotq[t].noalias()+=v[t]*sim.frame(t).charge(atoms[n]);
		}
	}
	
	//sum qdotr and rdotq to get the total dipole moment time derivative
	if(DEBUG_IR_QEQ>0) std::cout<<"computing dipole time derivative (sum of qdotr and rdotq)...\n";
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
	
	//compute the Fourier transform
	if(DEBUG_IR_QEQ>0) std::cout<<"computing the forward Fourier transforms...\n";
	for(unsigned int i=0; i<3; ++i) fftMu[i].transformf();
	
	//record the correlation function in frequency domain
	if(DEBUG_IR_QEQ>0) std::cout<<"recording the correlation function in the frequency domain...\n";
	for(unsigned int t=0; t<2*sim.timesteps(); ++t){
		fft.in(t)[0]=1.0/(2*sim.timesteps())*(
			fftMu[0].out(t)[0]*fftMu[0].out(t)[0]+fftMu[0].out(t)[1]*fftMu[0].out(t)[1]
			+fftMu[1].out(t)[0]*fftMu[1].out(t)[0]+fftMu[1].out(t)[1]*fftMu[1].out(t)[1]
			+fftMu[2].out(t)[0]*fftMu[2].out(t)[0]+fftMu[2].out(t)[1]*fftMu[2].out(t)[1]
		);
		fft.in(t)[1]=0.0;
	}
	
	/*
		//No truncation of correlation function
		if(DEBUG_IR_QEQ>0) std::cout<<"No truncation of correlation function.\n";
		
		//calculate the transpose of the total dipole moment
		if(DEBUG_IR_QEQ>0) std::cout<<"Calculating the transpose of the total dipole...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			Eigen::Vector3d muTot=Eigen::Vector3d::Zero();
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				muTot.noalias()+=(sim.frame(t).posn(n)-origin)*sim.frame(t).charge(n);
			}
			for(unsigned int i=0; i<3; ++i) dipoleTrans[i][t]=muTot[i];
		}
		
		//calculate the velocities of the dipole 
		if(DEBUG_IR_QEQ>0) std::cout<<"Calculating dipole velocity...\n";
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
		if(DEBUG_IR_QEQ>0) std::cout<<"Calculating the forward Fourier transforms...\n";
		for(unsigned int i=0; i<3; ++i) fftMu[i].transformf();
		
		//record the correlation function in frequency domain
		if(DEBUG_IR_QEQ>0) std::cout<<"Record the correlation function in the frequency domain...\n";
		for(unsigned int t=0; t<2*sim.timesteps(); ++t){
			fft.in(t)[0]=1.0/(2*sim.timesteps())*(
				fftMu[0].out(t)[0]*fftMu[0].out(t)[0]+fftMu[0].out(t)[1]*fftMu[0].out(t)[1]
				+fftMu[1].out(t)[0]*fftMu[1].out(t)[0]+fftMu[1].out(t)[1]*fftMu[1].out(t)[1]
				+fftMu[2].out(t)[0]*fftMu[2].out(t)[0]+fftMu[2].out(t)[1]*fftMu[2].out(t)[1]
			);
			fft.in(t)[1]=0.0;
		}
	*/
	
	//transform back into the time domain
	fft.transformr();
	
	//normalize and window the time-domain data
	for(unsigned int t=0; t<sim.timesteps(); ++t){
		fft.in(t)[0]=fft.out(t)[0]*window_(t)/(sim.timesteps()-t);
		fft.in(t)[1]=fft.out(t)[1]*window_(t)/(sim.timesteps()-t);
	}
	for(unsigned int t=sim.timesteps(); t<2*sim.timesteps(); ++t){
		fft.in(t)[0]=fft.out(t)[0]*window_(t-sim.timesteps())/(t-sim.timesteps()+1);
		fft.in(t)[1]=fft.out(t)[1]*window_(t-sim.timesteps())/(t-sim.timesteps()+1);
	}
	
	//transform into the frequency domain
	fft.transformf();
	
	//calculate the ir spectrum
	double norm=1.0/((2*sim.timesteps()));
	for(unsigned int t=0; t<freqN_; ++t){
		const double w=minFreq_*(t+0.5);
		irSpectrum_[t]=std::tanh(hbar*w/(kb*T_))/(1.0-std::exp(-hbar*w/(kb*T_)))*(fft.out(t)[0]*fft.out(t)[0]+fft.out(t)[1]*fft.out(t)[1]);
	}
}

void IRQEQ::printSpectrum(const char* file) const{
	if(DEBUG_IR_QEQ>0) std::cout<<"printSpectrum():\n";
	//local function variables
	FILE* writer=NULL;
	
	writer=fopen(file,"w");
	if(writer==NULL){
		std::cout<<"ERROR in IR_I::printSpectrum(const char*):\n";
		std::cout<<"Could not open parameter file.\n";
		throw std::runtime_error("I/O Exception Occured.");
	}
	
	//print the header
	fprintf(writer, "#Freq");
	if(freqUnit_==fourier::FreqUnit::MEV) fprintf(writer, "(meV) ");
	else if(freqUnit_==fourier::FreqUnit::CMI) fprintf(writer, "(cm^-1) ");
	else fprintf(writer, "(THz) ");
	fprintf(writer, "IR\n");
	
	//print the power spectrum
	for(unsigned int t=0; t<freqN_; ++t){
		fprintf(writer, "%f %f\n", minFreq_*t, irSpectrum_[t]);
	}
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
		Interval interval;
		Eigen::Vector3d offset=Eigen::Vector3d::Zero();
		Simulation sim;
		std::vector<double> jzero;
		std::vector<std::string> names;
		AtomType atomT;
		atomT.name=true; atomT.an=true; atomT.specie=true; atomT.index=true;
		atomT.posn=true; atomT.charge=true; atomT.jzero=true;
	//ewald
		Ewald3D::Coulomb coul;
		double precEwald=1e-5;
	//calculation
		IRQEQ ir3d;
		QEQ qeq;
	//input/output
		std::vector<std::string> strlist;
		std::string paramfile;
		std::string simstr;
		char* input=new char[string::M];
		FILE* reader=NULL;
		FILE_FORMAT::type fileFormat;
	//flags
		bool error=false;
		bool periodic=true;
	//units
		units::System::type unitsys;
	
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
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open parameter file.");
		
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
			} else if(strlist.at(0)=="OFFSET"){
				offset[0]=std::atof(strlist.at(1).c_str());
				offset[1]=std::atof(strlist.at(2).c_str());
				offset[2]=std::atof(strlist.at(3).c_str());
			} else if(strlist.at(0)=="INTERVAL"){
				interval=Interval::read(strlist.at(1).c_str());
			} else if(strlist.at(0)=="FORMAT"){
				fileFormat=FILE_FORMAT::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="UNITS"){
				unitsys=units::System::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="JZERO"){
				names.push_back(strlist.at(1));
				jzero.push_back(std::atof(strlist.at(2).c_str()));
			} else if(strlist.at(0)=="PREC_EWALD"){
				precEwald=std::atof(strlist.at(1).c_str());
			}
		}
		fclose(reader);
		reader=NULL;
		
		//======== check the parameters ========
		if(unitsys==units::System::UNKNOWN) throw std::invalid_argument("Invalid unit system");
		
		//======== print the parameters to screen ========
		std::cout<<"**********************************************\n";
		std::cout<<"************* GENERAL PARAMETERS *************\n";
		std::cout<<"SIMULATION PARAMETERS:\n";
		std::cout<<"\tUNITS    = "<<unitsys<<"\n";
		std::cout<<"\tSIM      = "<<simstr<<"\n";
		std::cout<<"\tFORMAT   = "<<fileFormat<<"\n";
		std::cout<<"\tOFFSET   = ("<<offset[0]<<","<<offset[1]<<","<<offset[2]<<")\n";
		std::cout<<"\tINTERVAL = "<<interval<<"\n";
		std::cout<<"\tATOMS    = "; for(unsigned int n=0; n<names.size(); ++n) std::cout<<names[n]<<" "; std::cout<<"\n";
		std::cout<<"\tJZERO    = "; for(unsigned int n=0; n<jzero.size(); ++n) std::cout<<jzero[n]<<" "; std::cout<<"\n";
		std::cout<<"************* GENERAL PARAMETERS *************\n";
		std::cout<<"**********************************************\n";
		
		//======== initialize the unit system ========
		std::cout<<"Initializing the unit system...\n";
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
			std::cout<<"structure = "<<sim.frame(0)<<"\n";
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
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
					sim.frame(t).posn(n).noalias()+=offset;
					Cell::returnToCell(
						sim.frame(t).posn(n),sim.frame(t).posn(n),
						sim.frame(t).cell().R(),sim.frame(t).cell().RInv()
					);
				}
			}
		}
		
		//======== print the simulation ========
		std::cout<<sim<<"\n";
		std::cout<<sim.frame(0)<<"\n";
		
		//======== read auxiliary objects ========
		std::cout<<"reading auxiliary objects...\n";
		//ewald object
			std::cout<<"Ewald objects...\n";
			coul.init(sim.frame(0),precEwald);
			std::cout<<coul<<"\n";
		//qeq object
			std::cout<<"QEQ object...\n";
			qeq.read(paramfile.c_str());
			qeq.init(sim.frame(0));
			std::cout<<qeq<<"\n";
		//ir3d object
			std::cout<<"IR3D object...\n";
			ir3d.read(paramfile.c_str());
			ir3d.init(sim);
			std::cout<<ir3d<<"\n";
		
		//======== assign jzero ========
		std::cout<<"assigning jzero...\n";
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				for(unsigned int i=0; i<names.size(); ++i){
					if(names[i]==sim.frame(t).name(n)){
						sim.frame(t).jzero(n)=jzero[i];
						break;
					}
				}
			}
		}
		
		//======== read chg data ========
		if(ir3d.readChgA()){
			std::cout<<"reading chg data...\n";
			Utility::Read::chg_atom(ir3d.fileChgA().c_str(),sim);
		}
		
		//======== exectute calculation ========
		std::cout<<"executing calculation...\n";
		if(ir3d.calcChg()){
			std::cout<<"computing atomic charges...\n";
			ir3d.calcChg(sim,qeq,coul);
		}
		if(ir3d.calcSpectrum()){
			std::cout<<"computing the spectrum...\n";
			ir3d.calcSpectrum(sim);
			ir3d.printSpectrum(ir3d.fileSpectrum().c_str());
		}
		
		//======== write simulation data ========
		std::cout<<"writing simulation data...\n";
		if(ir3d.writeDipoleT()) Utility::Write::dipole_tot(ir3d.fileDipoleT().c_str(),sim);
		if(ir3d.writeChgT()) Utility::Write::chg_tot(ir3d.fileChgT().c_str(),sim);
		if(ir3d.writeChgA()) Utility::Write::chg_atom(ir3d.fileChgA().c_str(),sim);
	}catch(std::exception& e){
		std::cout<<"ERROR in main(int,char**):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//======== free all local variables ========
	std::cout<<"freeing local variables...\n";
	delete[] input;
	if(reader!=NULL) fclose(reader);
	
	std::cout<<"exiting program...\n";
	if(!error) return 0;
	else return 1;
}
