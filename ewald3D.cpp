#include "ewald3D.hpp"

namespace Ewald3D{

//**********************************************************************************************************
//Utility Class
//**********************************************************************************************************

//operators

std::ostream& operator<<(std::ostream& out, const Utility& u){
	out<<"**********************************************\n";
	out<<"************** EWALD3D::UTILITY **************\n";
	out<<"PREC    = "<<u.prec_<<"\n";
	out<<"ALPHA   = "<<u.alpha_<<"\n";
	out<<"EPSILON = "<<u.eps_<<"\n";
	out<<"WEIGHT  = "<<u.weight_<<"\n";
	out<<"R_MAX   = "<<u.rMax_<<"\n";
	out<<"K_MAX   = "<<u.kMax_<<"\n";
	out<<"************** EWALD3D::UTILITY **************\n";
	out<<"**********************************************";
}

//member functions

void Utility::defaults(){
	if(DEBUG_EWALD_3D>0) std::cout<<"defaults():\n";
	prec_=1E-4;
	rMax_=0;
	kMax_=0;
	weight_=1;
	alpha_=0;
	eps_=1;//vacuum boundary conditions
}

void Utility::init(double prec){
	if(DEBUG_EWALD_3D>0) std::cout<<"init(const Cell&, int, double):\n";
	
	//set the precision
	if(DEBUG_EWALD_3D>1) std::cout<<"Setting the precision...\n";
	if(prec-1>=num_const::ZERO || prec<=0) throw std::invalid_argument("Precision must be in (0,1)");
	else prec_=prec;
}

//**********************************************************************************************************
//Coulomb Class
//**********************************************************************************************************

//operators

std::ostream& operator<<(std::ostream& out, const Coulomb& c){
	return out<<static_cast<const Utility&>(c);
}

//member functions

void Coulomb::defaults(){
	if(DEBUG_EWALD_3D>0) std::cout<<"defaults():\n";
	Utility::defaults();
	vSelfR_=0;
	vSelfK_=0;
	vSelfC_=0;
}

void Coulomb::init(const Structure& struc, double prec){
	if(DEBUG_EWALD_3D>0) std::cout<<"init(const Cell&, int, double):\n";
	Utility::init(prec);
	//local function variables
	const double pi=num_const::PI;
	
	//calculate the limits
	if(DEBUG_EWALD_3D>0) std::cout<<"Calculating the limits...\n";
	alpha_=std::pow(struc.nAtoms()*weight_*pi*pi*pi/(struc.cell().vol()*struc.cell().vol()),1.0/6.0);
	rMax_=std::sqrt(-1.0*std::log(prec))/alpha_;
	kMax_=2.0*alpha_*std::sqrt(-1.0*std::log(prec));
	if(DEBUG_EWALD_3D>0) std::cout<<"ALPHA = "<<alpha_<<" R_MAX = "<<rMax_<<" K_MAX = "<<kMax_<<"\n";
	
	//find the lattice vectors
	if(DEBUG_EWALD_3D>0) std::cout<<"Finding the lattice vectors...\n";
	int rNMaxX=std::ceil(rMax_/struc.cell().R().col(0).norm());
	int rNMaxY=std::ceil(rMax_/struc.cell().R().col(1).norm());
	int rNMaxZ=std::ceil(rMax_/struc.cell().R().col(2).norm());
	int kNMaxX=std::ceil(kMax_/struc.cell().K().col(0).norm());
	int kNMaxY=std::ceil(kMax_/struc.cell().K().col(1).norm());
	int kNMaxZ=std::ceil(kMax_/struc.cell().K().col(2).norm());
	if(DEBUG_EWALD_3D>0) std::cout<<"RN = ("<<rNMaxX<<","<<rNMaxY<<","<<rNMaxZ<<") - "<<(2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1)<<"\n";
	if(DEBUG_EWALD_3D>0) std::cout<<"KN = ("<<kNMaxX<<","<<kNMaxY<<","<<kNMaxZ<<") - "<<(2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1)<<"\n";
	R.resize((2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1));
	K.resize((2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1));
	int rCount=0;
	int kCount=0;
	const double f=1.05;
	for(int i=-rNMaxX; i<=rNMaxX; ++i){
		for(int j=-rNMaxY; j<=rNMaxY; ++j){
			for(int k=-rNMaxZ; k<=rNMaxZ; ++k){
				//note: don't skip R=0
				double norm=(i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2)).norm();
				if(norm<f*rMax_) R[rCount++].noalias()=i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2);
			}
		}
	}
	R.resize(rCount);
	for(int i=-kNMaxX; i<=kNMaxX; ++i){
		for(int j=-kNMaxY; j<=kNMaxY; ++j){
			for(int k=-kNMaxZ; k<=kNMaxZ; ++k){
				//note: skip K=0
				double norm=(i*struc.cell().K().col(0)+j*struc.cell().K().col(1)+k*struc.cell().K().col(2)).norm();
				if(norm>0 && norm<f*kMax_) K[kCount++].noalias()=i*struc.cell().K().col(0)+j*struc.cell().K().col(1)+k*struc.cell().K().col(2);
			}
		}
	}
	K.resize(kCount);
	if(DEBUG_EWALD_3D>0) std::cout<<"rCount = "<<rCount<<", kCount = "<<kCount<<"\n";
	
	//calculate the k sum amplitudes
	if(DEBUG_EWALD_3D>0) std::cout<<"computing the K-amplitudes...\n";
	kAmp.resize(kCount);
	for(int i=kCount-1; i>=0; --i){
		kAmp[i]=4.0*num_const::PI/struc.cell().vol()*std::exp(-K[i].dot(K[i])/(4.0*alpha_*alpha_))/(K[i].dot(K[i]));
	}
	
	//calculate the self-interaction strength
	if(DEBUG_EWALD_3D>0) std::cout<<"computing self-interaction strength...\n";
	vSelfR_=0; vSelfK_=0;
	for(int i=R.size()-1; i>=0; --i){
		if(R[i].norm()!=0) vSelfR_+=std::erfc(alpha_*R[i].norm())/R[i].norm();
	}
	for(int i=kAmp.size()-1; i>=0; --i){
		vSelfK_+=kAmp[i];
	}
	vSelfC_=-2.0*alpha_/num_const::RadPI;
	if(DEBUG_EWALD_3D>0) std::cout<<"vSelf = "<<vSelfR_+vSelfK_+vSelfC_<<"\n";
}

void Coulomb::init_alpha(const Structure& struc, double prec){
	if(DEBUG_EWALD_3D>0) std::cout<<"Coulomb::init_alpha(const Structure&,double):\n";
	if(prec>0) prec_=prec;
	double alpha=std::pow(struc.nAtoms()*weight_*num_const::PI*num_const::PI*num_const::PI/(struc.cell().vol()*struc.cell().vol()),1.0/6.0);
	rMax_=std::sqrt(-1.0*std::log(prec_))/alpha_;
	kMax_=2.0*alpha_*std::sqrt(-1.0*std::log(prec_));
	if(DEBUG_EWALD_3D>0) std::cout<<"ALPHA = "<<alpha_<<" R_MAX = "<<rMax_<<" K_MAX = "<<kMax_<<"\n";
}

//calculation - energy

double Coulomb::energy(const Structure& struc)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"energy(const Structure&,int,int*):\n";
	//local function variables
	double energyT=0,q2s=0;
	const std::complex<double> I(0,1);
	Eigen::Vector3d vec;
	if(DEBUG_EWALD_3D==0){
		//r-space
		for(int i=struc.nAtoms()-1; i>=0; --i){
			for(int j=i-1; j>=0; --j){
				Cell::diff(struc.posn(i),struc.posn(j),vec,struc.cell().R(),struc.cell().RInv());
				double energyS=0;
				for(int n=R.size()-1; n>=0; --n){
					double dist=(vec+R[n]).norm();
					energyS+=std::erfc(alpha_*dist)/dist;
				}
				energyT+=struc.charge(i)*struc.charge(j)*energyS;
			}
		}
		//k-space
		vec.setZero();
		for(int i=struc.nAtoms()-1; i>=0; --i){
			q2s+=struc.charge(i)*struc.charge(i);
			vec.noalias()+=struc.charge(i)*struc.posn(i);
			std::complex<double> sf(0,0);
			for(int n=K.size()-1; n>=0; --n){
				sf+=kAmp[n]*std::exp(-I*K[n].dot(struc.posn(i)));
			}
			energyT+=std::norm(struc.charge(i)*sf);
		}
		//self-energy
		energyT+=0.5*q2s*(vSelfR_+vSelfC_);
		//polarization energy
		if(DEBUG_EWALD_3D>0) std::cout<<"dipole = "<<vec.transpose()<<"\n";
		energyT+=2.0*num_const::PI/(2.0*eps_+1.0)*vec.dot(vec)/struc.cell().vol();
	} else {
		double energyR=0,energyK=0,energyS,energyP,q2s=0;
		//r-space
		for(int i=struc.nAtoms()-1; i>=0; --i){
			for(int j=i-1; j>=0; --j){
				Cell::diff(struc.posn(i),struc.posn(j),vec,struc.cell().R(),struc.cell().RInv());
				energyS=0;
				for(int n=R.size()-1; n>=0; --n){
					double dist=(vec+R[n]).norm();
					energyS+=std::erfc(alpha_*dist)/dist;
				}
				energyR+=struc.charge(i)*struc.charge(j)*energyS;
			}
		}
		//k-space
		vec.setZero();
		for(int i=struc.nAtoms()-1; i>=0; --i){
			q2s+=struc.charge(i)*struc.charge(i);
			vec.noalias()+=struc.charge(i)*struc.posn(i);
			std::complex<double> sf(0,0);
			for(int n=K.size()-1; n>=0; --n){
				sf+=kAmp[n]*std::exp(-I*K[n].dot(struc.posn(i)));
			}
			energyK+=struc.charge(i)*std::norm(sf);
		}
		//self-energy
		energyS=0.5*q2s*(vSelfR_+vSelfC_);
		//polarization energy
		if(DEBUG_EWALD_3D>0) std::cout<<"dipole = "<<vec.transpose()<<"\n";
		energyP=2.0*num_const::PI/(2.0*eps_+1.0)*vec.dot(vec)/struc.cell().vol();
	
		std::cout<<"energy-r = "<<energyR<<"\n";
		std::cout<<"energy-k = "<<energyK<<"\n";
		std::cout<<"energy-p = "<<energyP<<"\n";
		std::cout<<"energy-s = "<<energyS<<"\n";
		std::cout<<"energy-t = "<<energyR+energyK+energyS+energyP<<"\n";
		energyT=energyR+energyK+energyS+energyP;
	}
	
	return energyT;
}

double Coulomb::energy_single(const Structure& struc){
	if(DEBUG_EWALD_3D>0) std::cout<<"energy(const Structure&,int,int*):\n";
	//local function variables
	int shellx,shelly,shellz,count=0;
	double energyR=0,energyK=0,energyS,energyP,q2s=0;
	const std::complex<double> I(0,1);
	Eigen::Vector3d vec;
	
	//init alpha
	init_alpha(struc);
	
	//r-space
	shellx=std::ceil(rMax_/struc.cell().R().row(0).lpNorm<Eigen::Infinity>());//max of row - max abs x the lattice vectors
	shelly=std::ceil(rMax_/struc.cell().R().row(1).lpNorm<Eigen::Infinity>());//max of row - max abs y the lattice vectors
	shellz=std::ceil(rMax_/struc.cell().R().row(2).lpNorm<Eigen::Infinity>());//max of row - max abs z the lattice vectors
	if(DEBUG_EWALD_3D>0) std::cout<<"R_SHELL = ("<<shellx<<","<<shelly<<","<<shellz<<")\n";
	count=0;
	for(int i=-shellx; i<=shellx; ++i){
		for(int j=-shelly; j<=shelly; ++j){
			for(int k=-shellz; k<=shellz; ++k){
				//note: don't skip R=0
				double norm=(i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2)).norm();
				if(norm<1.05*rMax_) R[count++].noalias()=i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2);
			}
		}
	}
	R.resize(count);
	//r-space
	for(int i=struc.nAtoms()-1; i>=0; --i){
		for(int j=i-1; j>=0; --j){
			Cell::diff(struc.posn(i),struc.posn(j),vec,struc.cell().R(),struc.cell().RInv());
			energyS=0;
			for(int n=R.size()-1; n>=0; --n){
				double dist=(vec+R[n]).norm();
				energyS+=std::erfc(alpha_*dist)/dist;
			}
			energyR+=struc.charge(i)*struc.charge(j)*energyS;
		}
	}
	
	//k-space
	shellx=std::ceil(kMax_/struc.cell().K().row(0).lpNorm<Eigen::Infinity>());//max of row - max abs x the lattice vectors
	shelly=std::ceil(kMax_/struc.cell().K().row(1).lpNorm<Eigen::Infinity>());//max of row - max abs y the lattice vectors
	shellz=std::ceil(kMax_/struc.cell().K().row(2).lpNorm<Eigen::Infinity>());//max of row - max abs z the lattice vectors
	if(DEBUG_EWALD_3D>0) std::cout<<"K_SHELL = ("<<shellx<<","<<shelly<<","<<shellz<<")\n";
	count=0;
	for(int i=-shellx; i<=shellx; ++i){
		for(int j=-shelly; j<=shelly; ++j){
			for(int k=-shellz; k<=shellz; ++k){
				double norm=(i*struc.cell().K().col(0)+j*struc.cell().K().col(1)+k*struc.cell().K().col(2)).norm();//note: skip K=0
				if(norm>0 && norm<1.05*kMax_) K[count++].noalias()=i*struc.cell().K().col(0)+j*struc.cell().K().col(1)+k*struc.cell().K().col(2);
			}
		}
	}
	K.resize(count);
	kAmp.resize(count);
	for(unsigned int i=0; i<kAmp.size(); ++i){
		kAmp[i]=4.0*num_const::PI/struc.cell().vol()*std::exp(-K[i].dot(K[i])/(4.0*alpha_*alpha_))/(K[i].dot(K[i]));
	}
	for(int i=struc.nAtoms()-1; i>=0; --i){
		std::complex<double> sf(0,0);
		for(int n=K.size()-1; n>=0; --n){
			sf+=kAmp[n]*std::exp(-I*K[n].dot(struc.posn(i)));
		}
		energyK+=std::norm(struc.charge(i)*sf);
	}
	
	//total charge, dipole
	vec.setZero();
	for(int i=struc.nAtoms()-1; i>=0; --i){
		q2s+=struc.charge(i)*struc.charge(i);
		vec.noalias()+=struc.charge(i)*struc.posn(i);
	}
	
	//self-energy
	double vSelfR=0;
	for(int i=-shellx; i<=shellx; ++i){
		for(int j=-shelly; j<=shelly; ++j){
			for(int k=-shellz; k<=shellz; ++k){
				double norm=(i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2)).norm();
				if(norm>0) vSelfR+=std::erfc(alpha_*norm)/norm;
			}
		}
	}
	energyS=0.5*q2s*(vSelfR_-2.0*alpha_/num_const::RadPI);
	
	//polarization energy
	if(DEBUG_EWALD_3D>0) std::cout<<"dipole = "<<vec.transpose()<<"\n";
	energyP=2.0*num_const::PI/(2.0*eps_+1.0)*vec.dot(vec)/struc.cell().vol();
	
	if(DEBUG_EWALD_3D>0){
		std::cout<<"energy-r = "<<energyR<<"\n";
		std::cout<<"energy-k = "<<energyK<<"\n";
		std::cout<<"energy-s = "<<energyS<<"\n";
		std::cout<<"energy-p = "<<energyP<<"\n";
		std::cout<<"energy-t = "<<energyR+energyK+energyS+energyP<<"\n";
	}
	
	return energyR+energyK+energyS+energyP;
}

double Coulomb::energyBrute(const Structure& struc, int N)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"energyBrute(const Structure&,int,int*):\n";
	//local function variables
	double interEnergy=0,energySelf=0,chargeSum=0;
	double chgProd;
	Eigen::Vector3d dr;
	
	if(DEBUG_EWALD_3D>0) std::cout<<"Interaction Energy...\n";
	for(unsigned int n=0; n<struc.nAtoms(); ++n){
		for(unsigned int m=n+1; m<struc.nAtoms(); ++m){
			if(DEBUG_EWALD_3D>1) std::cout<<"Pair ("<<n<<","<<m<<")\n";
			//find the zeroth cell distance and energy
			Cell::diff(struc.posn(n),struc.posn(m),dr,struc.cell().R(),struc.cell().RInv());
			chgProd=struc.charge(n)*struc.charge(m);
			for(int i=-N; i<=N; ++i){
				for(int j=-N; j<=N; ++j){
					for(int k=-N; k<=N; ++k){
						interEnergy+=chgProd/(dr+i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2)).norm();
					}
				}
			}
		}
	}
	
	if(DEBUG_EWALD_3D>0) std::cout<<"Self-Energy...\n";
	for(unsigned int n=0; n<struc.nAtoms(); ++n){
		chargeSum+=struc.charge(n)*struc.charge(n);
	}
	energySelf=0;
	for(int i=-N; i<=N; ++i){
		for(int j=-N; j<=N; ++j){
			for(int k=-N; k<=N; ++k){
				double norm=(i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2)).norm();
				energySelf+=(norm>0)?1.0/norm:0;
			}
		}
	}
	energySelf*=0.5*chargeSum;
	
	if(DEBUG_EWALD_3D>0){
		std::cout<<"energy-s = "<<energySelf<<"\n";
		std::cout<<"energy-i = "<<interEnergy<<"\n";
		std::cout<<"energy-t = "<<interEnergy+energySelf<<"\n";
	}
	
	return interEnergy+energySelf;
}

//calculation - potential

double Coulomb::potential(const Structure& struc, unsigned int nn)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"potential(const Structure&,int)const:\n";
	double vTot=0,r2=0,chgsum=struc.charge(nn);
	Eigen::Vector3d dr;
	if(DEBUG_EWALD_3D==0){
		for(int i=struc.nAtoms()-1; i>=0; --i){
			if(i==nn) continue;
			double vLoc=0;
			Cell::diff(struc.posn(i),struc.posn(nn),dr,struc.cell().R(),struc.cell().RInv());
			for(int n=R.size()-1; n>=0; --n){
				double dist=(dr+R[n]).norm();
				vLoc+=std::erfc(alpha_*dist)/dist;
			}
			for(int n=K.size()-1; n>=0; --n){
				vLoc+=kAmp[n]*std::cos(K[n].dot(dr));
			}
			vTot+=vLoc*struc.charge(i);
			r2+=dr.squaredNorm()*struc.charge(i);
			chgsum+=struc.charge(i);
		}
		
		vTot+=-2.0*num_const::PI/(2.0*eps_+1.0)*r2/struc.cell().vol();
		vTot+=-num_const::PI/(alpha_*alpha_*struc.cell().vol())*chgsum;
		vTot+=(vSelfR_+vSelfK_+vSelfC_)*struc.charge(nn);
	} else {
		double vR=0,vK=0,vP,vC,vS;
		Eigen::Vector3d dr;
		
		for(int i=struc.nAtoms()-1; i>=0; --i){
			if(i==nn) continue;
			Cell::diff(struc.posn(i),struc.posn(nn),dr,struc.cell().R(),struc.cell().RInv());
			for(int n=R.size()-1; n>=0; --n){
				double dist=(dr+R[n]).norm();
				vR+=std::erfc(alpha_*dist)/dist*struc.charge(i);
			}
			for(int n=K.size()-1; n>=0; --n){
				vK+=kAmp[n]*std::cos(K[n].dot(dr))*struc.charge(i);
			}
			r2+=dr.squaredNorm()*struc.charge(i);
			chgsum+=struc.charge(i);
		}
		
		vP=-2.0*num_const::PI/(2.0*eps_+1.0)*r2/struc.cell().vol();
		vC=-num_const::PI/(alpha_*alpha_*struc.cell().vol())*chgsum;
		vS=(vSelfR_+vSelfK_+vSelfC_)*struc.charge(nn);
		
		std::cout<<"chg-t = "<<chgsum<<"\n";
		std::cout<<"pot-r = "<<vR<<"\n";
		std::cout<<"pot-k = "<<vK<<"\n";
		std::cout<<"pot-p = "<<vP<<"\n";
		std::cout<<"pot-c = "<<vC<<"\n";
		std::cout<<"pot-s = "<<vS<<"\n";
		std::cout<<"pot-t = "<<vR+vK+vP+vC+vS<<"\n";
		vTot=vR+vK+vP+vC+vS;
	}
	
	return vTot;
}

double Coulomb::potentialBrute(const Structure& struc, unsigned int n, int N)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"potentialBrute(const Structure&,int,const Eigen::Vector3d&)const:\n";
	double v=0;
	Eigen::Vector3d dr;
	
	for(unsigned int ii=0; ii<struc.nAtoms(); ++ii){
		if(ii==n) continue;
		Cell::diff(struc.posn(ii),struc.posn(n),dr,struc.cell().R(),struc.cell().RInv());
		for(int i=-N; i<=N; ++i){
			for(int j=-N; j<=N; ++j){
				for(int k=-N; k<=N; ++k){
					v+=1.0/(dr+i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2)).norm()*struc.charge(ii);
				}
			}
		}
	}
	
	return v;
}

double Coulomb::phi(const Structure& struc, const Eigen::Vector3d& dr)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"phi(const Structure&,int)const:\n";
	double vTot=0;
	if(DEBUG_EWALD_3D==0){
		for(int n=R.size()-1; n>=0; --n){
			double dist=(dr+R[n]).norm();
			vTot+=std::erfc(alpha_*dist)/dist;
		}
		for(int n=K.size()-1; n>=0; --n){
			vTot+=kAmp[n]*std::cos(K[n].dot(dr));
		}
		vTot+=-2.0*num_const::PI/(2.0*eps_+1.0)*dr.squaredNorm()/struc.cell().vol();
		vTot+=-num_const::PI/(alpha_*alpha_*struc.cell().vol());
		vTot+=(vSelfR_+vSelfK_+vSelfC_);
	} else {
		double vR=0,vK=0,vP,vC,vS;
		for(int n=R.size()-1; n>=0; --n){
			double dist=(dr+R[n]).norm();
			vR+=std::erfc(alpha_*dist)/dist;
		}
		for(int n=K.size()-1; n>=0; --n){
			vK+=kAmp[n]*std::cos(K[n].dot(dr));
		}
		
		vP=-2.0*num_const::PI/(2.0*eps_+1.0)*dr.squaredNorm()/struc.cell().vol();
		vC=-num_const::PI/(alpha_*alpha_*struc.cell().vol());
		vS=(vSelfR_+vSelfK_+vSelfC_);
		
		std::cout<<"pot-r = "<<vR<<"\n";
		std::cout<<"pot-k = "<<vK<<"\n";
		std::cout<<"pot-p = "<<vP<<"\n";
		std::cout<<"pot-c = "<<vC<<"\n";
		std::cout<<"pot-s = "<<vS<<"\n";
		std::cout<<"pot-t = "<<vR+vK+vP+vC+vS<<"\n";
		vTot=vR+vK+vP+vC+vS;
	}
	
	return vTot;
}

//calculation - electric field

Eigen::Vector3d& Coulomb::efield(const Structure& struc, unsigned int n, Eigen::Vector3d& field)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"field(const Structure&,int,Eigen::Vector3d&)const:\n";
	const double aa=2.0*alpha_/std::sqrt(num_const::PI);
	const double bb=4.0*num_const::PI/(2.0*eps_+1.0);
	Eigen::Vector3d dr;
	
	//zero field 
	field.setZero();
	
	//compute real/reciprocal space contributions
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		if(i==n) continue;
		Cell::diff(struc.posn(i),struc.posn(n),dr,struc.cell().R(),struc.cell().RInv());
		field.noalias()+=bb*dr*struc.charge(i);
		for(unsigned int n=0; n<R.size(); ++n){
			double dist=(dr+R[n]).norm();
			field.noalias()+=(dr+R[n])*(aa*std::exp(-alpha_*alpha_*dist*dist)
				+std::erfc(alpha_*dist)/dist)*struc.charge(i)/(dist*dist);
		}
		for(unsigned int n=0; n<K.size(); ++n){
			field.noalias()+=K[n]*kAmp[n]*std::sin(K[n].dot(dr))*struc.charge(i);
		}
	}
	
	return field;
}

Eigen::Vector3d& Coulomb::efieldBrute(const Structure& struc, unsigned int n, Eigen::Vector3d& field, int N)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"potentialBrute(const Structure&,int,Eigen::Vector3d&,unsigned int)const:\n";
	Eigen::Vector3d dr;
	
	//zero field
	field.setZero();
	
	for(unsigned int ii=0; ii<struc.nAtoms(); ++ii){
		Cell::diff(struc.posn(ii),struc.posn(n),dr,struc.cell().R(),struc.cell().RInv());
		for(int i=-N; i<=N; ++i){
			for(int j=-N; j<=N; ++j){
				for(int k=-N; k<=N; ++k){
					double drn=(dr+i*struc.cell().R().col(0)+j*struc.cell().R().col(1)+k*struc.cell().R().col(2)).norm();
					if(drn>num_const::ZERO) field.noalias()+=dr/(drn*drn*drn)*struc.charge(ii);
				}
			}
		}
	}
	
	return field;
}

//potential

double Coulomb::phiSelf()const{
	return vSelfR_+vSelfK_+vSelfC_;
}

//**********************************************************************************************************
//Dipole Class
//**********************************************************************************************************

//operators

std::ostream& operator<<(std::ostream& out, const Dipole& d){
	return out<<static_cast<const Utility&>(d);
}

//member functions

void Dipole::init(const Cell& cell, double prec){
	if(DEBUG_EWALD_3D>0) std::cout<<"init(const Cell&,double):\n";
	Utility::init(prec);
	double pi=num_const::PI;
	
	Eigen::Vector3d Rx=cell.R().col(0);
	Eigen::Vector3d Ry=cell.R().col(1);
	Eigen::Vector3d Rz=cell.R().col(2);
	Eigen::Vector3d Kx=cell.K().col(0);
	Eigen::Vector3d Ky=cell.K().col(1);
	Eigen::Vector3d Kz=cell.K().col(2);
	
	//calculate the limits
	if(DEBUG_EWALD_3D>0) std::cout<<"Calculating the limits...\n";
	alpha_=std::pow(weight_*pi*pi*pi/(cell.vol()*cell.vol()),1.0/6.0);
	rMax_=std::sqrt(-1.0*std::log(prec))/alpha_;
	kMax_=2*alpha_*std::sqrt(-1.0*std::log(prec));
	if(DEBUG_EWALD_3D>0) std::cout<<"Alpha = "<<alpha_<<", R_MAX = "<<rMax_<<", K_MAX = "<<kMax_<<"\n";
	
	//find the lattice vectors
	if(DEBUG_EWALD_3D>0) std::cout<<"Finding the lattice vectors...\n";
	int rNMaxX=std::ceil(rMax_/Rx.norm());
	int rNMaxY=std::ceil(rMax_/Ry.norm());
	int rNMaxZ=std::ceil(rMax_/Rz.norm());
	int kNMaxX=std::ceil(kMax_/Kx.norm());
	int kNMaxY=std::ceil(kMax_/Ky.norm());
	int kNMaxZ=std::ceil(kMax_/Kz.norm());
	if(DEBUG_EWALD_3D>0) std::cout<<"RN = ("<<rNMaxX<<","<<rNMaxY<<","<<rNMaxZ<<") - "<<(2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1)<<"\n";
	if(DEBUG_EWALD_3D>0) std::cout<<"KN = ("<<kNMaxX<<","<<kNMaxY<<","<<kNMaxZ<<") - "<<(2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1)<<"\n";
	R.resize((2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1));
	K.resize((2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1));
	int rCount=0;
	int kCount=0;
	double f=1.05;
	for(int i=-rNMaxX; i<=rNMaxX; ++i){
		for(int j=-rNMaxY; j<=rNMaxY; ++j){
			for(int k=-rNMaxZ; k<=rNMaxZ; ++k){
				//note: don't skip R=0
				double norm=(i*Rx+j*Ry+k*Rz).norm();
				if(norm<f*rMax_) R[rCount++].noalias()=i*Rx+j*Ry+k*Rz;
			}
		}
	}
	R.resize(rCount);
	for(int i=-kNMaxX; i<=kNMaxX; ++i){
		for(int j=-kNMaxY; j<=kNMaxY; ++j){
			for(int k=-kNMaxZ; k<=kNMaxZ; ++k){
				//note: skip K=0
				double norm=(i*Kx+j*Ky+k*Kz).norm();
				if(norm>0 && norm<f*kMax_) K[kCount++].noalias()=i*Kx+j*Ky+k*Kz;
			}
		}
	}
	K.resize(kCount);
	if(DEBUG_EWALD_3D>0) std::cout<<"rCount = "<<rCount<<", kCount = "<<kCount<<"\n";
	
	//calculate the K-matrices
	kMats.resize(K.size());
	for(unsigned int i=0; i<K.size(); ++i){
		kMats[i].noalias()=4.0*num_const::PI/cell.vol()*K[i]*K[i].transpose()*std::exp(-K[i].dot(K[i])/(4*alpha_*alpha_))/(K[i].dot(K[i]));
	}
	
	matS.noalias()=4.0*num_const::PI/(3.0*cell.vol())*Eigen::Matrix3d::Identity();
}

void Dipole::init(const Cell& cell, int nAtoms, double prec){
	if(DEBUG_EWALD_3D>0) std::cout<<"init(const Cell&,double):\n";
	Utility::init(prec);
	double pi=num_const::PI;
	
	Eigen::Vector3d Rx=cell.R().col(0);
	Eigen::Vector3d Ry=cell.R().col(1);
	Eigen::Vector3d Rz=cell.R().col(2);
	Eigen::Vector3d Kx=cell.K().col(0);
	Eigen::Vector3d Ky=cell.K().col(1);
	Eigen::Vector3d Kz=cell.K().col(2);
	
	//calculate the limits
	if(DEBUG_EWALD_3D>0) std::cout<<"Calculating the limits...\n";
	alpha_=std::pow(nAtoms*weight_*pi*pi*pi/(cell.vol()*cell.vol()),1.0/6.0);
	rMax_=std::sqrt(-1.0*std::log(prec))/alpha_;
	kMax_=2*alpha_*std::sqrt(-1.0*std::log(prec));
	if(DEBUG_EWALD_3D>0) std::cout<<"Alpha = "<<alpha_<<", R_MAX = "<<rMax_<<", K_MAX = "<<kMax_<<"\n";
	
	//find the lattice vectors
	if(DEBUG_EWALD_3D>0) std::cout<<"Finding the lattice vectors...\n";
	int rNMaxX=std::ceil(rMax_/Rx.norm());
	int rNMaxY=std::ceil(rMax_/Ry.norm());
	int rNMaxZ=std::ceil(rMax_/Rz.norm());
	int kNMaxX=std::ceil(kMax_/Kx.norm());
	int kNMaxY=std::ceil(kMax_/Ky.norm());
	int kNMaxZ=std::ceil(kMax_/Kz.norm());
	if(DEBUG_EWALD_3D>0) std::cout<<"RN = ("<<rNMaxX<<","<<rNMaxY<<","<<rNMaxZ<<") - "<<(2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1)<<"\n";
	if(DEBUG_EWALD_3D>0) std::cout<<"KN = ("<<kNMaxX<<","<<kNMaxY<<","<<kNMaxZ<<") - "<<(2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1)<<"\n";
	R.resize((2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1));
	K.resize((2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1));
	int rCount=0;
	int kCount=0;
	double f=1.05;
	for(int i=-rNMaxX; i<=rNMaxX; ++i){
		for(int j=-rNMaxY; j<=rNMaxY; ++j){
			for(int k=-rNMaxZ; k<=rNMaxZ; ++k){
				//note: don't skip R=0
				double norm=(i*Rx+j*Ry+k*Rz).norm();
				if(norm<f*rMax_) R[rCount++].noalias()=i*Rx+j*Ry+k*Rz;
			}
		}
	}
	R.resize(rCount);
	for(int i=-kNMaxX; i<=kNMaxX; ++i){
		for(int j=-kNMaxY; j<=kNMaxY; ++j){
			for(int k=-kNMaxZ; k<=kNMaxZ; ++k){
				//note: skip K=0
				double norm=(i*Kx+j*Ky+k*Kz).norm();
				if(norm>0 && norm<f*kMax_) K[kCount++].noalias()=i*Kx+j*Ky+k*Kz;
			}
		}
	}
	K.resize(kCount);
	if(DEBUG_EWALD_3D>0) std::cout<<"rCount = "<<rCount<<", kCount = "<<kCount<<"\n";
	
	//calculate the K-matrices
	kMats.resize(K.size());
	for(int i=0; i<K.size(); ++i){
		kMats[i].noalias()=4*num_const::PI/cell.vol()*K[i]*K[i].transpose()*std::exp(-K[i].dot(K[i])/(4*alpha_*alpha_))/(K[i].dot(K[i]));
	}
	
	matS.noalias()=4.0*num_const::PI/(3.0*cell.vol())*Eigen::Matrix3d::Identity();
	if(DEBUG_EWALD_3D>0) std::cout<<"matS = \n"<<matS<<"\n";
}

Eigen::Matrix3d& Dipole::interMat(const Eigen::Vector3d& r, Eigen::Matrix3d& mat)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"interMat(const Eigen::Vector3d&,Eigen::Matrix3d&):\n";
	Eigen::Vector3d dr;
	mat.setZero();
	
	//evaluate real-space part
	for(int i=R.size()-1; i>=0; --i){
		dr.noalias()=r+R[i];
		double norm=dr.norm();
		double exp=2*alpha_/num_const::RadPI*std::exp(-alpha_*alpha_*norm*norm);
		double erfc=(exp+std::erfc(alpha_*norm)/norm)/(norm*norm);
		mat.noalias()+=dr*dr.transpose()*(2*alpha_*alpha_*exp+3*erfc)/(norm*norm);
		mat.noalias()-=Eigen::Matrix3d::Identity()*erfc;
	}
	
	//evaluate k-space part
	for(int i=K.size()-1; i>=0; --i){
		mat.noalias()+=-std::cos(K[i].dot(r))*kMats[i];
	}
	
	//add surface term(?)
	//mat.noalias()+=matS*qtot;//only if qtot!=0
	
	return mat;
}

Eigen::Matrix3d& Dipole::interMatBrute(const Eigen::Vector3d& r, Eigen::Matrix3d& mat, int N)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"interMatBrute(const Eigen::Vector3d&,Eigen::Matrix3d&,int):\n";
	Eigen::Vector3d dr;
	mat.setZero();
	
	for(int i=-N; i<=N; ++i){
		for(int j=-N; j<=N; ++j){
			for(int k=-N; k<=N; ++k){
				//dr.noalias()=r+i*cell_.R().col(0)+j*cell_.R().col(1)+k*cell_.R().col(2);
				double n=dr.norm();
				//mat.noalias()+=3*dr*dr.transpose()/(n*n*n*n*n)-Eigen::Matrix3d::Identity()*1.0/(n*n*n);
				mat.noalias()+=3*dr*dr.transpose()/(n*n*n*n*n)-Eigen::Matrix3d::Identity()*1.0/(n*n*n);
				mat.noalias()-=Eigen::Matrix3d::Identity()*1.0/(n*n*n);
			}
		}
	}
	
	return mat;
}

Eigen::Matrix3d& Dipole::interMatSelf(Eigen::Matrix3d& mat)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"interMatSelf(Eigen::Matrix3d&):\n";
	mat.setZero();
	return mat;
}

Eigen::Matrix3d& Dipole::interMatSelfR(Eigen::Matrix3d& mat, int N)const{
	if(DEBUG_EWALD_3D>0) std::cout<<"interMatSelfR(Eigen::Matrix3d&):\n";
	Eigen::Vector3d r;
	mat.setZero();
	
	for(int i=-N; i<=N; ++i){
		for(int j=-N; j<=N; ++j){
			for(int k=-N; k<=N; ++k){
				//r=i*cell_.R().col(0)+j*cell_.R().col(1)+k*cell_.R().col(2);
				double n=r.norm();
				//if(n!=0) mat.noalias()+=3*r*r.transpose()/(n*n*n*n*n)-Eigen::Matrix3d::Identity()*1.0/(n*n*n);
				if(n!=0){
					mat.noalias()+=3*r*r.transpose()/(n*n*n*n*n);
					mat.noalias()-=Eigen::Matrix3d::Identity()*1.0/(n*n*n);
				}
			}
		}
	}
	
	return mat;
}

}