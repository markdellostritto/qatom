#include "ewald3D.hpp"

namespace Ewald3D{

//**********************************************************************************************************
//Utility Class
//**********************************************************************************************************

//operators

Utility& Utility::operator=(const Utility& u){
	init(u.cell(),u.nAtoms(),u.prec());
	return *this;
}

std::ostream& operator<<(std::ostream& out, const Utility& u){
	out<<"**********************************************\n";
	out<<"************** EWALD3D::UTILITY **************\n";
	out<<"PREC    = "<<u.prec_<<"\n";
	out<<"ALPHA   = "<<u.alpha_<<"\n";
	out<<"EPSILON = "<<u.epsilon_<<"\n";
	out<<"R_MAX   = "<<u.rMax_<<"\n";
	out<<"K_MAX   = "<<u.kMax_<<"\n";
	out<<"N_ATOMS = "<<u.nAtoms_<<"\n";
	out<<"W       = "<<u.weight_<<"\n";
	out<<"************** EWALD3D::UTILITY **************\n";
	out<<"**********************************************";
}

//member functions

void Utility::defaults(){
	if(DEBUG_EWALD_3D>0) log<<"defaults():\n";
	prec_=1E-4;
	rMax_=0;
	kMax_=0;
	weight_=1;
	alpha_=0;
	nAtoms_=0;
	epsilon_=1;//vacuum boundary conditions
}

void Utility::init(const Cell& cell, int nAtoms, double prec){
	if(DEBUG_EWALD_3D>0) log<<"init(const Cell&, int, double):\n";
	
	//set the precision
	if(DEBUG_EWALD_3D>1) log<<"Setting the precision...\n";
	if(prec-1>=num_const::ZERO || prec<=0) throw std::invalid_argument("Precision must be in (0,1)");
	else prec_=prec;
	if(DEBUG_EWALD_3D>1) log<<"prec = "<<prec_<<"\n";
	
	//set the unit cell
	if(DEBUG_EWALD_3D>1) log<<"Setting the unit cell...\n";
	if(DEBUG_EWALD_3D>1) log<<"Cell = \n"<<cell<<"\n";
	cell_=cell;
	nAtoms_=nAtoms;
	if(DEBUG_EWALD_3D>1) log<<"Cell = \n"<<cell_<<"\n";
}

//**********************************************************************************************************
//Coulomb Class
//**********************************************************************************************************

//operators

Coulomb& Coulomb::operator=(const Coulomb& c){
	init(c.cell(),c.nAtoms(),c.prec());
	return *this;
}

std::ostream& operator<<(std::ostream& out, const Coulomb& c){
	return out<<static_cast<const Utility&>(c);
}

//member functions

void Coulomb::defaults(){
	if(DEBUG_EWALD_3D>0) log<<"defaults():\n";
	Utility::defaults();
	vSelf_=0;
}

void Coulomb::init(const Cell& cell, int nAtoms, double prec){
	if(DEBUG_EWALD_3D>0) log<<"init(const Cell&, int, double):\n";
	Utility::init(cell,nAtoms,prec);
	//local function variables
	double pi=num_const::PI;
	
	if(DEBUG_EWALD_3D>0) log<<"Setting the lattice info...\n";
	Eigen::Vector3d Rx=cell_.R().col(0);
	Eigen::Vector3d Ry=cell_.R().col(1);
	Eigen::Vector3d Rz=cell_.R().col(2);
	Eigen::Vector3d Kx=cell_.K().col(0);
	Eigen::Vector3d Ky=cell_.K().col(1);
	Eigen::Vector3d Kz=cell_.K().col(2);
	double rx=Rx[0]; double ry=Ry[1]; double rz=Rz[2];
	Iyz=8*std::atan(ry*rz/(rx*std::sqrt(rx*rx+ry*ry+rz*rz)));
	Izx=8*std::atan(rz*rx/(ry*std::sqrt(rx*rx+ry*ry+rz*rz)));
	Ixy=8*std::atan(rx*ry/(rz*std::sqrt(rx*rx+ry*ry+rz*rz)));
	
	//calculate the limits
	if(DEBUG_EWALD_3D>0) log<<"Calculating the limits...\n";
	double vol=cell_.vol();
	alpha_=std::pow(nAtoms_*weight_*pi*pi*pi/(vol*vol),1.0/6.0);
	rMax_=std::sqrt(-1.0*std::log(prec))/alpha_;
	kMax_=2*alpha_*std::sqrt(-1.0*std::log(prec));
	if(DEBUG_EWALD_3D>0) log<<"Alpha = "<<alpha_<<", R_MAX = "<<rMax_<<", K_MAX = "<<kMax_<<"\n";
	
	//find the lattice vectors
	if(DEBUG_EWALD_3D>0) log<<"Finding the lattice vectors...\n";
	int rNMaxX=std::ceil(rMax_/Rx.norm());
	int rNMaxY=std::ceil(rMax_/Ry.norm());
	int rNMaxZ=std::ceil(rMax_/Rz.norm());
	int kNMaxX=std::ceil(kMax_/Kx.norm());
	int kNMaxY=std::ceil(kMax_/Ky.norm());
	int kNMaxZ=std::ceil(kMax_/Kz.norm());
	if(DEBUG_EWALD_3D>0) log<<"RN = ("<<rNMaxX<<","<<rNMaxY<<","<<rNMaxZ<<") - "<<(2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1)<<"\n";
	if(DEBUG_EWALD_3D>0) log<<"KN = ("<<kNMaxX<<","<<kNMaxY<<","<<kNMaxZ<<") - "<<(2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1)<<"\n";
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
	if(DEBUG_EWALD_3D>0) log<<"rCount = "<<rCount<<", kCount = "<<kCount<<"\n";
	
	//calculate the k sum amplitudes
	if(DEBUG_EWALD_3D>0) log<<"Calculating the K-amplitudes...\n";
	kAmp.resize(kCount);
	for(unsigned int i=0; i<kAmp.size(); ++i){
		kAmp[i]=4*num_const::PI/vol*std::exp(-K[i].dot(K[i])/(4*alpha_*alpha_))/(K[i].dot(K[i]));
	}
	
	//calculate the self-interaction strength
	if(DEBUG_EWALD_3D>0) log<<"Calculating self-interaction strength...\n";
	double vSelfR=0;
	double vSelfK=0;
	for(unsigned int i=0; i<R.size(); ++i){
		//if(R[i].norm()!=0) vSelfR+=boost::math::erfc(alpha_*R[i].norm())/R[i].norm();
		if(R[i].norm()!=0) vSelfR+=std::erfc(alpha_*R[i].norm())/R[i].norm();
	}
	for(unsigned int i=0; i<kAmp.size(); ++i){
		vSelfK+=kAmp[i];
	}
	vSelf_=vSelfR+vSelfK-2*alpha_/num_const::RadPI;
	if(DEBUG_EWALD_3D>0) log<<"vSelf = "<<vSelf_<<"\n";
}

//potential

double Coulomb::phi(const Eigen::Vector3d& dr)const{
	double phiR=0,phiK=0,phiP,phiC,dist;
	
	for(unsigned int n=0; n<R.size(); ++n){
		dist=(dr+R[n]).norm();
		phiR+=std::erfc(alpha_*dist)/dist;
	}
	for(unsigned int n=0; n<K.size(); ++n){
		phiK+=kAmp[n]*std::cos(K[n].dot(dr));
	}
	phiK*=4*num_const::PI/cell_.vol();
	//phiP=-2*num_const::PI*dr.dot(dr)/(3*cell_.vol());
	phiP=-0.5/cell_.vol()*(dr[0]*dr[0]*Iyz+dr[1]*dr[1]*Izx+dr[2]*dr[2]*Ixy);
	phiC=-num_const::PI/(alpha_*alpha_*cell_.vol());
	
	if(DEBUG_EWALD_3D>0){
		log<<"R-space phi      = "<<phiR<<"\n";
		log<<"K-space phi      = "<<phiK<<"\n";
		log<<"Polarization phi = "<<phiP<<"\n";
		log<<"Constant phi     = "<<phiC<<"\n";
		log<<"Total potential  = "<<phiR+phiK+phiP+phiC<<"\n";
	}
	
	return phiR+phiK+phiP+phiC;
}

double Coulomb::phiSelf()const{
	return vSelf_-num_const::PI/(alpha_*alpha_*cell_.vol());
}

//**********************************************************************************************************
//Dipole Class
//**********************************************************************************************************

//operators

Dipole& Dipole::operator=(const Dipole& d){
	init(d.cell(),d.prec());
	return *this;
}

std::ostream& operator<<(std::ostream& out, const Dipole& d){
	return out<<static_cast<const Utility&>(d);
}

//member functions

void Dipole::init(const Cell& cell, double prec){
	if(DEBUG_EWALD_3D>0) log<<"init(const Cell&,double):\n";
	Utility::init(cell,0,prec);
	double pi=num_const::PI;
	
	Eigen::Vector3d Rx=cell_.R().col(0);
	Eigen::Vector3d Ry=cell_.R().col(1);
	Eigen::Vector3d Rz=cell_.R().col(2);
	Eigen::Vector3d Kx=cell_.K().col(0);
	Eigen::Vector3d Ky=cell_.K().col(1);
	Eigen::Vector3d Kz=cell_.K().col(2);
	
	//calculate the limits
	if(DEBUG_EWALD_3D>0) log<<"Calculating the limits...\n";
	double vol=cell_.vol();
	alpha_=std::pow(weight_*pi*pi*pi/(vol*vol),1.0/6.0);
	rMax_=std::sqrt(-1.0*std::log(prec))/alpha_;
	kMax_=2*alpha_*std::sqrt(-1.0*std::log(prec));
	if(DEBUG_EWALD_3D>0) log<<"Alpha = "<<alpha_<<", R_MAX = "<<rMax_<<", K_MAX = "<<kMax_<<"\n";
	
	//find the lattice vectors
	if(DEBUG_EWALD_3D>0) log<<"Finding the lattice vectors...\n";
	int rNMaxX=std::ceil(rMax_/Rx.norm());
	int rNMaxY=std::ceil(rMax_/Ry.norm());
	int rNMaxZ=std::ceil(rMax_/Rz.norm());
	int kNMaxX=std::ceil(kMax_/Kx.norm());
	int kNMaxY=std::ceil(kMax_/Ky.norm());
	int kNMaxZ=std::ceil(kMax_/Kz.norm());
	if(DEBUG_EWALD_3D>0) log<<"RN = ("<<rNMaxX<<","<<rNMaxY<<","<<rNMaxZ<<") - "<<(2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1)<<"\n";
	if(DEBUG_EWALD_3D>0) log<<"KN = ("<<kNMaxX<<","<<kNMaxY<<","<<kNMaxZ<<") - "<<(2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1)<<"\n";
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
	if(DEBUG_EWALD_3D>0) log<<"rCount = "<<rCount<<", kCount = "<<kCount<<"\n";
	
	//calculate the K-matrices
	kMats.resize(K.size());
	for(int i=0; i<K.size(); ++i){
		kMats[i].noalias()=4*num_const::PI/vol*K[i]*K[i].transpose()*std::exp(-K[i].dot(K[i])/(4*alpha_*alpha_))/(K[i].dot(K[i]));
	}
	
	matS.noalias()=4*num_const::PI/(3*cell_.vol())*Eigen::Matrix3d::Identity();
}

void Dipole::init(const Cell& cell, int nAtoms, double prec){
	if(DEBUG_EWALD_3D>0) log<<"init(const Cell&,double):\n";
	Utility::init(cell,0,prec);
	double pi=num_const::PI;
	nAtoms_=nAtoms;
	
	Eigen::Vector3d Rx=cell_.R().col(0);
	Eigen::Vector3d Ry=cell_.R().col(1);
	Eigen::Vector3d Rz=cell_.R().col(2);
	Eigen::Vector3d Kx=cell_.K().col(0);
	Eigen::Vector3d Ky=cell_.K().col(1);
	Eigen::Vector3d Kz=cell_.K().col(2);
	
	//calculate the limits
	if(DEBUG_EWALD_3D>0) log<<"Calculating the limits...\n";
	double vol=cell_.vol();
	alpha_=std::pow(nAtoms*weight_*pi*pi*pi/(vol*vol),1.0/6.0);
	rMax_=std::sqrt(-1.0*std::log(prec))/alpha_;
	kMax_=2*alpha_*std::sqrt(-1.0*std::log(prec));
	if(DEBUG_EWALD_3D>0) log<<"Alpha = "<<alpha_<<", R_MAX = "<<rMax_<<", K_MAX = "<<kMax_<<"\n";
	
	//find the lattice vectors
	if(DEBUG_EWALD_3D>0) log<<"Finding the lattice vectors...\n";
	int rNMaxX=std::ceil(rMax_/Rx.norm());
	int rNMaxY=std::ceil(rMax_/Ry.norm());
	int rNMaxZ=std::ceil(rMax_/Rz.norm());
	int kNMaxX=std::ceil(kMax_/Kx.norm());
	int kNMaxY=std::ceil(kMax_/Ky.norm());
	int kNMaxZ=std::ceil(kMax_/Kz.norm());
	if(DEBUG_EWALD_3D>0) log<<"RN = ("<<rNMaxX<<","<<rNMaxY<<","<<rNMaxZ<<") - "<<(2*rNMaxX+1)*(2*rNMaxY+1)*(2*rNMaxZ+1)<<"\n";
	if(DEBUG_EWALD_3D>0) log<<"KN = ("<<kNMaxX<<","<<kNMaxY<<","<<kNMaxZ<<") - "<<(2*kNMaxX+1)*(2*kNMaxY+1)*(2*kNMaxZ+1)<<"\n";
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
	if(DEBUG_EWALD_3D>0) log<<"rCount = "<<rCount<<", kCount = "<<kCount<<"\n";
	
	//calculate the K-matrices
	kMats.resize(K.size());
	for(int i=0; i<K.size(); ++i){
		kMats[i].noalias()=4*num_const::PI/vol*K[i]*K[i].transpose()*std::exp(-K[i].dot(K[i])/(4*alpha_*alpha_))/(K[i].dot(K[i]));
	}
	
	matS.noalias()=4.0*num_const::PI/(3.0*cell_.vol())*Eigen::Matrix3d::Identity();
	if(DEBUG_EWALD_3D>0) log<<"matS = \n"<<matS<<"\n";
}

Eigen::Matrix3d& Dipole::interMat(const Eigen::Vector3d& r, Eigen::Matrix3d& mat)const{
	if(DEBUG_EWALD_3D>0) log<<"interMat(const Eigen::Vector3d&,Eigen::Matrix3d&):\n";
	Eigen::Vector3d dr;
	//Eigen::Matrix3d rMat=Eigen::Matrix3d::Zero();
	//Eigen::Matrix3d kMat=Eigen::Matrix3d::Zero();
	mat.setZero();
	
	//evaluate real-space part
	for(unsigned int i=0; i<R.size(); ++i){
		dr.noalias()=r+R[i];
		double norm=dr.norm();
		double exp=2*alpha_/num_const::RadPI*std::exp(-alpha_*alpha_*norm*norm);
		double erfc=(exp+std::erfc(alpha_*norm)/norm)/(norm*norm);
		//mat.noalias()+=dr*dr.transpose()*(2*alpha_*alpha_*exp+3*erfc)/(norm*norm)-Eigen::Matrix3d::Identity()*erfc;
		mat.noalias()+=dr*dr.transpose()*(2*alpha_*alpha_*exp+3*erfc)/(norm*norm);
		mat.noalias()-=Eigen::Matrix3d::Identity()*erfc;
		//rMat.noalias()+=dr*dr.transpose()*(2*alpha_*alpha_*exp+3*erfc)/(norm*norm)-Eigen::Matrix3d::Identity()*erfc;
	}
	
	//evaluate k-space part
	for(unsigned int i=0; i<K.size(); ++i){
		mat.noalias()+=-std::cos(K[i].dot(r))*kMats[i];
		//kMat.noalias()+=-std::cos(K[i].dot(r))*kMats[i];
	}
	
	//add surface term(?)
	//mat.noalias()+=matS*qtot;//only if qtot!=0
	
	/*log<<"R-mat = ";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			log.log()<<rMat(i,j)<<" ";
		}
	}
	log.log()<<"\n";
	log<<"K-mat = ";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			log.log()<<kMat(i,j)<<" ";
		}
	}
	log.log()<<"\n";
	log<<"S-mat = ";
	for(unsigned int i=0; i<3; ++i){
		for(unsigned int j=0; j<3; ++j){
			log.log()<<matS(i,j)<<" ";
		}
	}
	log.log()<<"\n";*/
	
	//mat.noalias()=rMat+kMat+matS;
	//mat.noalias()=rMat+kMat;
	return mat;
}

Eigen::Matrix3d& Dipole::interMatBrute(const Eigen::Vector3d& r, Eigen::Matrix3d& mat, int N)const{
	if(DEBUG_EWALD_3D>0) log<<"interMatBrute(const Eigen::Vector3d&,Eigen::Matrix3d&,int):\n";
	Eigen::Vector3d dr;
	mat.setZero();
	
	for(int i=-N; i<=N; ++i){
		for(int j=-N; j<=N; ++j){
			for(int k=-N; k<=N; ++k){
				dr.noalias()=r+i*cell_.R().col(0)+j*cell_.R().col(1)+k*cell_.R().col(2);
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
	if(DEBUG_EWALD_3D>0) log<<"interMatSelf(Eigen::Matrix3d&):\n";
	mat.setZero();
	return mat;
}

Eigen::Matrix3d& Dipole::interMatSelfR(Eigen::Matrix3d& mat, int N)const{
	if(DEBUG_EWALD_3D>0) log<<"interMatSelfR(Eigen::Matrix3d&):\n";
	Eigen::Vector3d r;
	mat.setZero();
	
	for(int i=-N; i<=N; ++i){
		for(int j=-N; j<=N; ++j){
			for(int k=-N; k<=N; ++k){
				r=i*cell_.R().col(0)+j*cell_.R().col(1)+k*cell_.R().col(2);
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