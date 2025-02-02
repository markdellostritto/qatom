#ifndef THOLE_HPP
#define THOLE_HPP

//c libraries
#include <cstdlib>
//c++ libraries
#include <iostream>
//electrostatics
#include "idd.hpp"
#include "ewald3D.hpp"
//math
#include "math_const.hpp"
//chem
#include "ptable.hpp"
//simulation
#include "cell.hpp"
//eigen
#include "eigen.hpp"
//units
#include "units.hpp"

#ifndef DEBUG_THOLE
#define DEBUG_THOLE 0
#endif

//****************************************************
//THOLE
//****************************************************

class Thole{
private:
	//flags
		bool inter_;//whether or not to include interaction matrices
		bool alphar_;//whether or not to use volume-dependent polarizabilities
	//parameters
		double a_;//scaling constant
	//dipole interactions
		IDD::Form idd_;
	//linear solver
		eigen::LIN_SOLVER::type linSolver_;
	//matrix utilities
		Eigen::Matrix3d interMat_,selfInterMat_;
		Eigen::MatrixXd alphaC_,identityC_;
		Eigen::MatrixXd A_;
	//data
		std::vector<double> r_,r0_;
		std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > drAtom_;
public:
	//constructors/destructors
	Thole(){defaults();};
	~Thole(){};
	
	//access
	//dipole interactions
		IDD::Form& idd(){return idd_;};
		const IDD::Form& idd()const{return idd_;};
	//linear solver
		eigen::LIN_SOLVER::type& linSolver(){return linSolver_;};
		const eigen::LIN_SOLVER::type& linSolver()const{return linSolver_;};
	//matrix utilities
		Eigen::Matrix3d& interMat(){return interMat_;};
		const Eigen::Matrix3d& interMat()const{return interMat_;};
		Eigen::Matrix3d& selfInterMat(){return selfInterMat_;};
		const Eigen::Matrix3d& selfInterMat()const{return selfInterMat_;};
		Eigen::MatrixXd& alphaC(){return alphaC_;};
		const Eigen::MatrixXd& alphaC()const{return alphaC_;};
		Eigen::MatrixXd& identityC(){return identityC_;};
		const Eigen::MatrixXd& identityC()const{return identityC_;};
		Eigen::MatrixXd& A(){return A_;};
		const Eigen::MatrixXd& A()const{return A_;};
	//polarizability parameters
		double& a(){return a_;};
		const double& a()const{return a_;};
		bool& inter(){return inter_;};
		const bool& inter()const{return inter_;};
		bool& alphar(){return alphar_;};
		const bool& alphar()const{return alphar_;};
	//data
		std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& drAtom(){return drAtom_;};
		const std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >& drAtom()const{return drAtom_;};
		std::vector<double>& r(){return r_;};
		const std::vector<double>& r()const{return r_;};
		std::vector<double>& r0(){return r0_;};
		const std::vector<double>& r0()const{return r0_;};
	
	//operators
		friend std::ostream& operator<<(std::ostream& out, const Thole& thole);
	
	//polarizability
		template <class AtomT, class MoleculeT> void alpha(MoleculeT& mol);
		void init(Structure& struc, const Ewald3D::Dipole& dipole);
		void alpha(Structure& struc, const Ewald3D::Dipole& dipole);
		void alpha_cont(Structure& struc, const Ewald3D::Dipole& dipole);
	//member functions
		void defaults();
		void read(const char* str);
};

//static functions

template <class AtomT, class MoleculeT>
void Thole::alpha(MoleculeT& mol){
	if(DEBUG_THOLE>0) std::cout<<"Thole::alpha(MoleculeT&):\n";
	//local variables
	const double ke=units::consts::ke();
	//units
	double rscale=0.0;
	if(units::consts::system()==units::System::AU) rscale=units::BOHRpANG;
	else if(units::consts::system()==units::System::METAL) rscale=1.0;
	else throw std::runtime_error("Invalid units.");
	//utility vectors
	Eigen::Vector3d rVec,r00,rNew;
	
	//set the utility matrices
	if(DEBUG_THOLE>1) std::cout<<"Setting utility matrices ("<<mol.nAtoms()<<")...\n";
	A_=Eigen::MatrixXd::Zero(3*mol.nAtoms(),3*mol.nAtoms());
	alphaC_.resize(3*mol.nAtoms(),3);
	identityC_.resize(3*mol.nAtoms(),3);
	r0_.resize(mol.nAtoms());
	r_.resize(mol.nAtoms());
	drAtom_.resize(mol.nAtoms());
	
	//set the gas-phase atomic radius
	if(DEBUG_THOLE>1) std::cout<<"Setting gas-phase atomic radii.\n";
	for(unsigned int i=0; i<mol.nAtoms(); ++i) r0_[i]=PTable::covalentRadius(mol.atom(i).an())*rscale;
	for(unsigned int i=0; i<mol.nAtoms(); ++i) r_[i]=PTable::covalentRadius(mol.atom(i).an())*rscale;
	//calculate the atom-in-molecule radius
	if(DEBUG_THOLE>1) std::cout<<"Calculating the atom-in-molecule radius...\n";
	for(unsigned int i=0; i<mol.nAtoms(); ++i){
		double rMin=100;
		for(unsigned int j=0; j<mol.nAtoms(); ++j){
			if(i==j) continue;
			//calculate the radius
			rVec.noalias()=(mol.atom(j).posn()-mol.atom(i).posn());
			double dr=rVec.norm();
			double rTemp=r0_[i]-0.5*(r0_[i]+r0_[j]-dr);
			if(rTemp<rMin){
				r_[i]=rTemp;
				rMin=rTemp;
				r00.noalias()=r0_[i]*rVec/rVec.norm();
				rNew.noalias()=r_[i]*rVec/rVec.norm();
				drAtom_[i][0]=std::fabs(rNew[0]/r00[0]);
				drAtom_[i][1]=std::fabs(rNew[1]/r00[1]);
				drAtom_[i][2]=std::fabs(rNew[2]/r00[2]);
				if(rVec[0]==0 && r00[0]==0) drAtom_[i][0]=1;
				if(rVec[1]==0 && r00[1]==0) drAtom_[i][1]=1;
				if(rVec[2]==0 && r00[2]==0) drAtom_[i][2]=1;
			}
		}
	}
	//print the radius
	if(DEBUG_THOLE>2){
		std::cout<<"RADIUS:\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			std::cout<<mol.atom(i).name()<<" "<<r0_[i]<<" "<<r_[i]<<" "<<drAtom_[i].transpose()<<"\n";
		}
	}
	
	//populate the solution vector
	if(DEBUG_THOLE>1) std::cout<<"Populating the solution vector...\n";
	for(unsigned int i=0; i<mol.nAtoms(); ++i) identityC_.block<3,3>(i*3,0)=Eigen::Matrix3d::Identity();
	//populate the diagonal of the A-matrix
	if(alphar_){
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			mol.atom(i).alpha()(0,0)*=drAtom_[i][0];
			mol.atom(i).alpha()(1,1)*=drAtom_[i][1];
			mol.atom(i).alpha()(2,2)*=drAtom_[i][2];
		}
	}
	for(unsigned int i=0; i<mol.nAtoms(); ++i){
		A_.block<3,3>(i*3,i*3).noalias()=mol.atom(i).alpha().inverse()*ke;
		/*Eigen::MatrixXd alpha=mol.atom(i).alpha();
		if(alphar_){
			alpha(0,0)*=drAtom_[i][0];
			alpha(1,1)*=drAtom_[i][1];
			alpha(2,2)*=drAtom_[i][2];
		}
		A_.block<3,3>(i*3,i*3).noalias()=alpha.inverse()*ke;*/
	}
	//print the polarizability
	if(DEBUG_THOLE>2){
		std::cout<<"ALPHA:\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			std::cout<<mol.atom(i).name()<<" "<<mol.atom(i).alpha()<<"\n";
		}
	}
	
	//calculate the off-diagonal elements (interaction matrices)
	if(DEBUG_THOLE>1) std::cout<<"Calculating off-diagonal matrices...\n";
	if(inter_){
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			for(unsigned int j=i+1; j<mol.nAtoms(); ++j){
				A_.block<3,3>(j*3,i*3).noalias()=-1*ke*IDD::iTensor[idd_](
					mol.atom(j).posn()-mol.atom(i).posn(),interMat_,
					a_*IDD::scale[idd_](mol.atom(j).alpha().trace()/3.0,mol.atom(i).alpha().trace()/3.0)
				);
			}
		}
	}
	A_=A_.selfadjointView<Eigen::Lower>();
	
	//calculate the effective polarizabilities
	if(DEBUG_THOLE>1) std::cout<<"Calculating effective polarizabilities...\n";
	switch(linSolver_){
		case eigen::LIN_SOLVER::LLT: alphaC_.noalias()=A_.llt().solve(identityC_); break;
		case eigen::LIN_SOLVER::LDLT: alphaC_.noalias()=A_.ldlt().solve(identityC_); break;
		case eigen::LIN_SOLVER::HQR: alphaC_.noalias()=A_.householderQr().solve(identityC_); break;
		case eigen::LIN_SOLVER::CPHQR: alphaC_.noalias()=A_.colPivHouseholderQr().solve(identityC_); break;
		case eigen::LIN_SOLVER::PPLU: alphaC_.noalias()=A_.partialPivLu().solve(identityC_); break;
		case eigen::LIN_SOLVER::FPLU: alphaC_.noalias()=A_.fullPivLu().solve(identityC_); break;
		default: A_.fullPivLu().solve(identityC_); break;
	};
	
	//record the polarizability
	if(DEBUG_THOLE>1) std::cout<<"Recording polarizabilities...\n";
	for(unsigned int i=0; i<mol.nAtoms(); ++i){
		mol.atom(i).alpha().noalias()=alphaC_.block<3,3>(i*3,0);
	}
}

#endif