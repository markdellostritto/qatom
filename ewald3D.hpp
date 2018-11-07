#ifndef EWALD_3D_HPP
#define EWALD_3D_HPP

//c libraries
#include <cmath>
#include <cstdlib>
//Eigen
#include <Eigen/Dense>
#include <Eigen/StdVector>
//ame
#include "math_const.hpp"
//local
#include "log.hpp"
#include "cell.hpp"
#include "sim.hpp"

#ifndef DEBUG_EWALD_3D
#define DEBUG_EWALD_3D 0
#endif

namespace Ewald3D{

//**********************************************************************************************************
//Utility Class
//**********************************************************************************************************

class Utility{
protected:
	//calculation parameters
	double prec_;//precision for the calculation
	double rMax_;//the maximum length of lattice vector included in the real-space sum
	double kMax_;//the maximum length of the lattice vector included in the reciprocal-space sum
	double alpha_;//the integral cutoff separating the real- and reciprocal-space sums
	double weight_;//weighting of real space calculations
	int nAtoms_;
	
	//unit cell
	Cell cell_;//the unit cell
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > R;//the lattice vectors to sum over
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > K;//the reciprocal lattice vectors to sum over
	
	//electrostatics
	double epsilon_;//the dielectric constant of the outer material
	
	//log
	logging::DebugLogger log;
public:
	Utility():log("Ewald3D::Utility"){defaults();};
	Utility(const Cell& cell, int nAtoms, double prec):log("Ewald3D::Utility"){init(cell,nAtoms,prec);};
	Utility(const Utility& u);
	~Utility(){};
	
	//operators
	Utility& operator=(const Utility& u);
	friend std::ostream& operator<<(std::ostream& out, const Utility& u);
	
	//access
	const Cell& cell()const{return cell_;};
	double prec()const{return prec_;};
	double rMax()const{return rMax_;};
	double kMax()const{return kMax_;};
	double alpha()const{return alpha_;};
	double nAtoms()const{return nAtoms_;};
	const double& epsilon()const{return epsilon_;};
	double& epsilon(){return epsilon_;};
	const double& weight()const{return weight_;};
	double& weight(){return weight_;};
	
	//member functions
	void defaults();
	void clear(){defaults();};
	void init(const Cell& cell, int nAtoms, double prec);
};

//**********************************************************************************************************
//Coulomb Class
//**********************************************************************************************************

class Coulomb: public Utility{
private:
	//calculation parameters
	double vSelf_;//the self-interaction strength for an ion in a periodic lattice
	
	//unit cell
	std::vector<double> kAmp;//the reciprocal space sum amplitudes
	double Ixy,Iyz,Izx;
	
	//log
	logging::DebugLogger log;
public:
	Coulomb():log("Ewald3D::Coulomb"){defaults();};
	Coulomb(const Cell& cell, int nAtoms, double prec):log("Ewald3D::Coulomb"){init(cell,nAtoms,prec);};
	Coulomb(const Coulomb& c):log("Ewald3D::Coulomb"){init(c.cell(),c.nAtoms(),c.prec());};
	~Coulomb(){};
	
	//operators
	Coulomb& operator=(const Coulomb& c);
	friend std::ostream& operator<<(std::ostream& out, const Coulomb& c);
	
	//member functions
	void defaults();
	void clear(){defaults();};
	void init(const Cell& cell, int nAtoms, double prec);
	template <class AtomT> void init(const SimAtomic<AtomT>& sim, double prec);
	
	//calculation - energy
	template <class AtomT> double energy(const SimAtomic<AtomT>& sim, int t)const;
	template <class AtomT> double energyBrute(const SimAtomic<AtomT>& sim, int t, int N)const;
	
	//calculation - potential
	double phi(const Eigen::Vector3d& dr)const;
	double phiSelf()const;
	template <class AtomT> double potential(const SimAtomic<AtomT>& sim, int t, const AtomT& atom)const;
	template <class AtomT> double potentialBrute(const SimAtomic<AtomT>& sim, int t, const AtomT& atom, int N)const;
};

//member functions

template <class AtomT> 
void Coulomb::init(const SimAtomic<AtomT>& sim, double prec){
	if(DEBUG_EWALD_3D>0) log<<"init(const SimAtomic<AtomT>&, double):\n";
	init(sim.cell(0),sim.nAtoms(),prec);
}

//calculation - energy

template <class AtomT>
double Coulomb::energy(const SimAtomic<AtomT>& sim, int t)const{
	if(DEBUG_EWALD_3D>0) log<<"energy(const SimAtomic<AtomT>&,int,int*):\n";
	//local function variables
	double energyR=0,energyK=0,energySelf=0,energyP=0;
	double dist,chgProd,chargeSum=0;
	Eigen::Vector3d dr,dipole=Eigen::Vector3d::Zero();
	
	for(unsigned int i=0; i<sim.nAtoms(); ++i){
		chargeSum+=sim.atom(t,i).charge()*sim.atom(t,i).charge();
		dipole.noalias()+=sim.atom(t,i).charge()*sim.atom(t,i).posn();
		for(unsigned int j=i+1; j<sim.nAtoms(); ++j){
			Cell::diff(sim.atom(t,i).posn(),sim.atom(t,j).posn(),dr,sim.R(),sim.RInv());
			chgProd=sim.atom(t,i).charge()*sim.atom(t,j).charge();
			for(unsigned int n=0; n<R.size(); ++n){
				dist=(dr+R[n]).norm();
				//energyR+=chgProd*boost::math::erfc(alpha_*dist)/dist;
				energyR+=chgProd*std::erfc(alpha_*dist)/dist;
			}
			for(unsigned int n=0; n<K.size(); ++n){
				energyK+=chgProd*kAmp[n]*std::cos(K[n].dot(dr));
			}
		}
	}
	energySelf=0.5*chargeSum*vSelf_;
	if(DEBUG_EWALD_3D>0) log<<"Dipole = "<<dipole.transpose()<<"\n";
	//energyP=2*num_const::PI/(2*EPSILON+1)*dipole.dot(dipole)/vol;
	energyP=0.5/sim.cell(t).vol()*(dipole[0]*dipole[0]*Iyz+dipole[1]*dipole[1]*Izx+dipole[2]*dipole[2]*Ixy);
	
	if(DEBUG_EWALD_3D>0){
		log<<"R-space energy = "<<energyR<<"\n";
		log<<"K-space energy = "<<energyK<<"\n";
		log<<"Self-energy = "<<energySelf<<"\n";
		log<<"Polarization energy = "<<energyP<<"\n";
		log<<"Total energy = "<<energyR+energyK+energySelf+energyP<<"\n";
	}
	
	return energyR+energyK+energySelf+energyP;
}

template <class AtomT>
double Coulomb::energyBrute(const SimAtomic<AtomT>& sim, int t, int N)const{
	if(DEBUG_EWALD_3D>0) log<<"energyBrute(const SimAtomic<AtomT>&,int,int*):\n";
	//local function variables
	double interEnergy=0,energySelf=0,chargeSum=0;
	double chgProd;
	Eigen::Vector3d dr;
	
	if(DEBUG_EWALD_3D>0) log<<"Interaction Energy...\n";
	for(unsigned int n=0; n<sim.nAtoms(); ++n){
		for(unsigned int m=n+1; m<sim.nAtoms(); ++m){
			if(DEBUG_EWALD_3D>1) log<<"Pair ("<<n<<","<<m<<")\n";
			//find the zeroth cell distance and energy
			Cell::diff(sim.atom(t,n).posn(),sim.atom(t,m).posn(),dr,sim.R(),sim.RInv());
			chgProd=sim.atom(t,n).charge()*sim.atom(t,m).charge();
			for(int i=-N; i<=N; ++i){
				for(int j=-N; j<=N; ++j){
					for(int k=-N; k<=N; ++k){
						interEnergy+=chgProd/(dr+i*sim.cell(t).R().col(0)+j*sim.cell(t).R().col(1)+k*sim.cell(t).R().col(2)).norm();
					}
				}
			}
		}
	}
	
	log<<"Self-Energy...\n";
	for(unsigned int n=0; n<sim.nAtoms(); ++n){
		chargeSum+=sim.atom(t,n).charge()*sim.atom(t,n).charge();
	}
	energySelf=0;
	for(int i=-N; i<=N; ++i){
		for(int j=-N; j<=N; ++j){
			for(int k=-N; k<=N; ++k){
				double norm=(i*sim.cell(t).R().col(0)+j*sim.cell(t).R().col(1)+k*sim.cell(t).R().col(2)).norm();
				energySelf+=(norm>0)?1.0/norm:0;
			}
		}
	}
	energySelf*=0.5*chargeSum;
	
	if(DEBUG_EWALD_3D>0){
		log<<"Self-Energy = "<<energySelf<<"\n";
		log<<"Interaction Energy = "<<interEnergy<<"\n";
		log<<"Total Energy = "<<interEnergy+energySelf<<"\n";
	}
	
	return interEnergy+energySelf;
}

//calculation - potential

template <class AtomT>
double Coulomb::potential(const SimAtomic<AtomT>& sim, int t, const AtomT& atom)const{
	if(DEBUG_EWALD_3D>0) log<<"potential(const SimAtomic<AtomT>&,int,const Eigen::Vector3d&)const:\n";
	double vR=0,vK=0,vP,dist;
	Eigen::Vector3d dr,r2=Eigen::Vector3d::Zero();
	
	for(unsigned int i=0; i<sim.nAtoms(); ++i){
		if(atom.specie()==sim.atom(t,i).specie() && atom.index()==sim.atom(t,i).index()) continue;
		Cell::diff(sim.atom(t,i).posn(),sim.atom(t,atom.specie(),atom.index()).posn(),dr,sim.cell(t).R(),sim.cell(t).RInv());
		for(unsigned int n=0; n<R.size(); ++n){
			dist=(dr+R[n]).norm();
			//vR+=boost::math::erfc(alpha_*dist)/dist*sim.atom(t,i).charge();
			vR+=std::erfc(alpha_*dist)/dist*sim.atom(t,i).charge();
		}
		for(unsigned int n=0; n<K.size(); ++n){
			vK+=kAmp[n]*std::cos(K[n].dot(dr))*sim.atom(t,i).charge();
		}
		r2[0]+=dr[0]*dr[0]*sim.atom(t,i).charge();
		r2[1]+=dr[1]*dr[1]*sim.atom(t,i).charge();
		r2[2]+=dr[2]*dr[2]*sim.atom(t,i).charge();
	}
	
	vK*=4*num_const::PI/sim.cell(t).vol();
	vP=-0.5/sim.cell(t).vol()*(r2[0]*Iyz+r2[1]*Izx+r2[2]*Ixy);
	
	if(DEBUG_EWALD_3D>0){
		log<<"R-space potential = "<<vR<<"\n";
		log<<"K-space potential = "<<vK<<"\n";
		log<<"Polarization potential = "<<vP<<"\n";
		log<<"Total potential = "<<vR+vK+vP<<"\n";
	}
	
	return vR+vK+vP;
}

template <class AtomT>
double Coulomb::potentialBrute(const SimAtomic<AtomT>& sim, int t, const AtomT& atom, int N)const{
	if(DEBUG_EWALD_3D>0) log<<"potentialBrute(const SimAtomic<AtomT>&,int,const Eigen::Vector3d&)const:\n";
	double v=0;
	Eigen::Vector3d dr;
	
	for(unsigned int n=0; n<sim.nAtoms(); ++n){
		if(sim.atom(t,n).specie()==atom.specie() && sim.atom(t,n).index()==atom.index()) continue;
		Cell::diff(sim.atom(t,n).posn(),sim.atom(t,atom.specie(),atom.index()).posn(),dr,sim.cell(t).R(),sim.cell(t).RInv());
		for(int i=-N; i<=N; ++i){
			for(int j=-N; j<=N; ++j){
				for(int k=-N; k<=N; ++k){
					v+=1.0/(dr+i*sim.cell(t).R().col(0)+j*sim.cell(t).R().col(1)+k*sim.cell(t).R().col(2)).norm()*sim.atom(t,n).charge();
				}
			}
		}
	}
	
	return v;
}

//**********************************************************************************************************
//Dipole Class
//**********************************************************************************************************

class Dipole: public Utility{
private:
	//summation variables
	std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > kMats;
	Eigen::Matrix3d matS;//surface term
	
	//log
	logging::DebugLogger log;
public:
	//constructors/destructors
	Dipole():log("Ewald3D::Dipole"){defaults();};
	Dipole(const Cell& cell, double prec){init(cell,prec);};
	Dipole(const Dipole& d){init(d.cell(),d.prec());};
	~Dipole(){};
	
	//operators
	Dipole& operator=(const Dipole& d);
	friend std::ostream& operator<<(std::ostream& out, const Dipole& d);
	
	//member functions
	void defaults(){};
	void clear(){defaults();};
	void init(const Cell& cell, double prec);
	void init(const Cell& cell, int nAtoms, double prec);
	Eigen::Matrix3d& interMat(const Eigen::Vector3d& r, Eigen::Matrix3d& mat)const;
	Eigen::Matrix3d& interMatBrute(const Eigen::Vector3d& r, Eigen::Matrix3d& mat, int N)const;
	Eigen::Matrix3d& interMatSelf(Eigen::Matrix3d& mat)const;
	Eigen::Matrix3d& interMatSelfR(Eigen::Matrix3d& mat, int N)const;
};

}

#endif