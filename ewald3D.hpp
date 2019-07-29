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
#include "cell.hpp"
#include "structure.hpp"

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
	
	//unit cell
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > R;//the lattice vectors to sum over
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > K;//the reciprocal lattice vectors to sum over
	
	//electrostatics
	double eps_;//the dielectric constant of the outer material
public:
	Utility(){defaults();};
	Utility(double prec){init(prec);};
	Utility(const Utility& u);
	~Utility(){};
	
	//operators
	friend std::ostream& operator<<(std::ostream& out, const Utility& u);
	
	//access
	double prec()const{return prec_;};
	double rMax()const{return rMax_;};
	double kMax()const{return kMax_;};
	double alpha()const{return alpha_;};
	const double& eps()const{return eps_;};
	double& eps(){return eps_;};
	const double& weight()const{return weight_;};
	double& weight(){return weight_;};
	
	//member functions
	void defaults();
	void clear(){defaults();};
	void init(double prec);
};

//**********************************************************************************************************
//Coulomb Class
//**********************************************************************************************************

class Coulomb: public Utility{
private:
	//calculation parameters
	double vSelfR_,vSelfK_,vSelfC_;//the self-interaction strength for an ion in a periodic lattice
	
	//unit cell
	std::vector<double> kAmp;//the reciprocal space sum amplitudes
	double Ixy,Iyz,Izx;
public:
	Coulomb(){defaults();};
	Coulomb(Structure& struc, double prec){init(struc,prec);};
	~Coulomb(){};
	
	//operators
	friend std::ostream& operator<<(std::ostream& out, const Coulomb& c);
	
	//member functions
	void defaults();
	void clear(){defaults();};
	void init(const Structure& struc, double prec);
	void init_alpha(const Structure& struc, double prec=0);
	
	//calculation - energy
	double energy(const Structure& struc)const;
	double energy_single(const Structure& struc);
	double energyBrute(const Structure& struc, int N)const;
	
	//calculation - potential
	double phi(const Structure& struc, const Eigen::Vector3d& dr)const;
	double phiSelf()const;
	double potential(const Structure& struc, unsigned int n)const;
	double potentialBrute(const Structure& struc, unsigned int n, int N)const;
	
	//calculation - electric field
	Eigen::Vector3d& efield(const Structure& struc, unsigned int n, Eigen::Vector3d&)const;
	Eigen::Vector3d& efieldBrute(const Structure& struc, unsigned int n, Eigen::Vector3d&, int N)const;
};

//**********************************************************************************************************
//Dipole Class
//**********************************************************************************************************

class Dipole: public Utility{
private:
	//summation variables
	std::vector<Eigen::Matrix3d,Eigen::aligned_allocator<Eigen::Matrix3d> > kMats;
	Eigen::Matrix3d matS;//surface term
public:
	//constructors/destructors
	Dipole(){defaults();};
	Dipole(const Cell& cell, double prec){init(cell,prec);};
	~Dipole(){};
	
	//operators
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