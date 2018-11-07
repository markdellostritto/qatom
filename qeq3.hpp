#ifndef QEQ3_HPP
#define QEQ3_HPP

//c libraries
#include <cmath>
//c++ libraries
#include <vector>
#include <iostream>
//chem info
#include "ptable.hpp"
//atoms/molecules
#include "atom.hpp"
#include "molecule.hpp"
//simulation
#include "sim.hpp"
//electrostatics
#include "icc.hpp"
#include "ewald3D.hpp"
//eigen
#include "eigen.hpp"
//units
#include "units.hpp"

#ifndef DEBUG_QEQ3
#define DEBUG_QEQ3 0
#endif

//************************************************************
//QEQ3
//************************************************************

class QEQ3{
private:
	//parameters
	double k_;//scaling constant
	Eigen::VectorXd chi_;//electronegativities
	std::vector<double> r_,r0_,drAtom_;
	bool scaler_;
	
	//coulomb interactions
	ICC::Form icc_;
	Eigen::MatrixXd J_;//coulomb integral matrix
	
	//matrix utilities
	Eigen::MatrixXd A_;//solving matrix
	Eigen::VectorXd b_;//constant vector
	Eigen::VectorXd x_;//solution vector
	Eigen::VectorXd I_;//idempotential
	
	//element properties
	PTable::ElectronegativityType::type chiType_;
	
	//charge
	double qTot_;
public:
	//constructors/destructors
	QEQ3(){defaults();};
	~QEQ3(){};
	
	//operators
	friend std::ostream& operator<<(std::ostream& out, const QEQ3& qeq);
	
	//access
	//atom info - vectors
		const Eigen::VectorXd& chi()const{return chi_;};
	//atom info - atomic
		double chi(unsigned int i)const{return chi_[i];};
	//coulomb interactions
		ICC::Form& icc(){return icc_;};
		const ICC::Form& icc()const{return icc_;};
	//scaling
		double& k(){return k_;};
		const double& k()const{return k_;};
	//electronegativity
		PTable::ElectronegativityType::type& chiType(){return chiType_;};
		const PTable::ElectronegativityType::type& chiType()const{return chiType_;};
		bool& scaler(){return scaler_;};
		const bool& scaler()const{return scaler_;};
	//charge
		double& qTot(){return qTot_;};
		const double& qTot()const{return qTot_;};
	
	//member functions
	void defaults();
	void clear(){defaults();};
	void load(const char* paramFile);
	template <class MoleculeT> void qt(MoleculeT& mol);
	template <class MoleculeT> void qt_new(MoleculeT& mol);
	template <class MoleculeT> void qt_jZero(MoleculeT& mol);
	template <class AtomT> void qt(unsigned int t, SimAtomic<AtomT>& sim, const Ewald3D::Coulomb& coul);
	template <class AtomT> void qt_jZero(unsigned int t, SimAtomic<AtomT>& sim, const Ewald3D::Coulomb& coul);
	template <class AtomT> void init(SimAtomic<AtomT>& sim);
	template <class AtomT> void qt_cont(unsigned int t, SimAtomic<AtomT>& sim, const Ewald3D::Coulomb& coul);
	template <class AtomT> void qt_jZero_cont(unsigned int t, SimAtomic<AtomT>& sim, const Ewald3D::Coulomb& coul);
};

//member functions

template <class MoleculeT>
void QEQ3::qt(MoleculeT& mol){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::qt(Molecule<AtomT>&):\n";
	//utility vectors
	Eigen::Vector3d rVec,r00,rNew;
	
	if(mol.nAtoms()==0) throw std::invalid_argument("No atoms given for charge transfer.");
	else if(mol.nAtoms()==1) mol.atom(0).charge()=0;
	else {
		//resize matrices
		if(DEBUG_QEQ3>0) std::cout<<"Resizing matrices...\n";
		J_.resize(mol.nAtoms(),mol.nAtoms());
		A_.resize(mol.nAtoms(),mol.nAtoms());
		chi_.resize(mol.nAtoms());
		r_.resize(mol.nAtoms());
		r0_.resize(mol.nAtoms());
		drAtom_.resize(mol.nAtoms());
		
		//local variables
		double escale=0.0,rscale=0.0;
		if(units::consts::system()==units::System::AU){
			escale=units::HARTREEpEV;
			rscale=units::BOHRpANG;
		}
		else if(units::consts::system()==units::System::METAL){
			escale=1.0;
			rscale=1.0;
		}
		else throw std::runtime_error("Invalid units.");
		const double ke=units::consts::ke();
		
		//calculate atom info
		if(DEBUG_QEQ3>0) std::cout<<"Calculating atom info...\n";
		if(DEBUG_QEQ3>1) std::cout<<"Element AN CHI J00 C\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			//electronegativity
			if(chiType_==PTable::ElectronegativityType::PAULING) chi_[i]=PTable::electronegativityPauling(mol.atom(i).an());
			else if(chiType_==PTable::ElectronegativityType::ALLEN) chi_[i]=PTable::electronegativityAllen(mol.atom(i).an());
			else if(chiType_==PTable::ElectronegativityType::MULLIKEN) chi_[i]=PTable::electronegativityMulliken(mol.atom(i).an());
			else if(chiType_==PTable::ElectronegativityType::HINZE) chi_[i]=PTable::electronegativityHinze(mol.atom(i).an(),0);
			if(chiType_==PTable::ElectronegativityType::MULLIKEN && mol.atom(i).an()==1) chi_[i]=4.528;
			chi_[i]*=escale;
			//idempotential
			if(chiType_==PTable::ElectronegativityType::HINZE) J_(i,i)=PTable::idempotentialHinze(mol.atom(i).an(),0);
			else J_(i,i)=PTable::idempotential(mol.atom(i).an());
			J_(i,i)*=escale;
			if(DEBUG_QEQ3>1){
				std::cout<<mol.atom(i).name()<<" ";
				std::cout<<mol.atom(i).an()<<" ";
				std::cout<<chi_[i]<<" ";
				std::cout<<J_(i,i)<<"\n";
			}
		}
		//print the electronegativity
		if(DEBUG_QEQ3>0){
			std::cout<<"CHI:\n";
			for(unsigned int i=0; i<mol.nAtoms(); ++i){
				std::cout<<mol.atom(i).name()<<" "<<1.0/drAtom_[i]*chi_[i]<<"\n";
			}
		}
		
		//set the gas-phase atomic radius
		if(DEBUG_QEQ3>1) std::cout<<"Setting gas-phase atomic radii.\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i) r0_[i]=PTable::covalentRadius(mol.atom(i).an())*rscale;
		for(unsigned int i=0; i<mol.nAtoms(); ++i) r_[i]=PTable::covalentRadius(mol.atom(i).an())*rscale;
		//calculate the atom-in-molecule radius
		if(DEBUG_QEQ3>1) std::cout<<"Calculating the atom-in-molecule radius...\n";
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
					drAtom_[i]=rNew.norm()/r00.norm();
				}
			}
		}
		//print the radius
		if(DEBUG_QEQ3>0){
			std::cout<<"RADIUS:\n";
			for(unsigned int i=0; i<mol.nAtoms(); ++i){
				std::cout<<mol.atom(i).name()<<" "<<r0_[i]<<" "<<r_[i]<<" "<<drAtom_[i]<<"\n";
			}
		}
		
		//calculate coulomb integrals
		if(DEBUG_QEQ3>0) std::cout<<"Calculating Coulomb Integrals...\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			for(unsigned int j=i+1; j<mol.nAtoms(); ++j){
				J_(j,i)=ke*ICC::iTensor[icc_](
					mol.atom(j).posn()-mol.atom(i).posn(),
					k_*ICC::scale[icc_](J_(i,i),J_(j,j))
				);
			}
		}
		//calculate the operator matrix
		if(DEBUG_QEQ3>0) std::cout<<"Calculating operator matrix...\n";
		A_=J_.selfadjointView<Eigen::Lower>();
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			A_(0,i)=1;
			for(unsigned int j=1; j<mol.nAtoms(); ++j){
				A_(j,i)-=J_(i,0);
			}
		}
		//calculate the solution vector
		if(DEBUG_QEQ3>0) std::cout<<"Calculating solution vector...\n";
		if(scaler_) for(unsigned int i=0; i<mol.nAtoms(); ++i) chi_[i]*=1.0/drAtom_[i];
		if(DEBUG_QEQ3>0){
			for(unsigned int i=0; i<mol.nAtoms(); ++i){
				std::cout<<"chi["<<i<<"] = "<<chi_[i]<<"\n";
			}
		}
		b_=Eigen::VectorXd::Constant(mol.nAtoms(),chi_[0]);
		b_.noalias()-=chi_;
		b_[0]=qTot_;
		//solve the linear equations
		if(DEBUG_QEQ3>0) std::cout<<"Solving the linear equations...\n";
		x_.noalias()=A_.partialPivLu().solve(b_);
		//set the atomic charges
		if(DEBUG_QEQ3>0) std::cout<<"Setting the charge...\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i) mol.atom(i).charge()=x_[i];
	}
}

template <class MoleculeT>
void QEQ3::qt_jZero(MoleculeT& mol){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::qt_jZero(Molecule<AtomT>&):\n";
	if(mol.nAtoms()==0) throw std::invalid_argument("No atoms given for charge transfer.");
	else if(mol.nAtoms()==1) mol.atom(0).charge()=0;
	else {
		Eigen::Vector3d rVec,rNew,r00;
		
		//resize matrices
		if(DEBUG_QEQ3>0) std::cout<<"Resizing matrices...\n";
		J_.resize(mol.nAtoms(),mol.nAtoms());
		A_.resize(mol.nAtoms(),mol.nAtoms());
		chi_.resize(mol.nAtoms());
		r_.resize(mol.nAtoms());
		r0_.resize(mol.nAtoms());
		drAtom_.resize(mol.nAtoms());
		
		//local variables
		double escale=0.0,rscale=0.0;
		if(units::consts::system()==units::System::AU){
			escale=units::HARTREEpEV;
			rscale=units::ANGpBOHR;
		}
		else if(units::consts::system()==units::System::METAL){
			escale=1.0;
			rscale=1.0;
		}
		else throw std::runtime_error("Invalid units.");
		const double ke=units::consts::ke();
		
		//calculate atom info
		if(DEBUG_QEQ3>0) std::cout<<"Calculating atom info...\n";
		if(DEBUG_QEQ3>1) std::cout<<"Element AN CHI J00 C\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			//electronegativity
			if(chiType_==PTable::ElectronegativityType::PAULING) chi_[i]=PTable::electronegativityPauling(mol.atom(i).an());
			else if(chiType_==PTable::ElectronegativityType::ALLEN) chi_[i]=PTable::electronegativityAllen(mol.atom(i).an());
			else if(chiType_==PTable::ElectronegativityType::MULLIKEN) chi_[i]=PTable::electronegativityMulliken(mol.atom(i).an());
			else if(chiType_==PTable::ElectronegativityType::HINZE) chi_[i]=PTable::electronegativityHinze(mol.atom(i).an(),0);
			if(chiType_==PTable::ElectronegativityType::MULLIKEN && mol.atom(i).an()==1) chi_[i]=4.528;
			chi_[i]*=escale;
			//idempotential
			J_(i,i)=mol.atom(i).jzero();
			if(DEBUG_QEQ3>1){
				std::cout<<mol.atom(i).name()<<" ";
				std::cout<<mol.atom(i).an()<<" ";
				std::cout<<chi_[i]<<" ";
				std::cout<<J_(i,i)<<"\n";
			}
		}
		
		//set the gas-phase atomic radius
		if(DEBUG_QEQ3>1) std::cout<<"Setting gas-phase atomic radii.\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i) r0_[i]=PTable::covalentRadius(mol.atom(i).an())*rscale;
		for(unsigned int i=0; i<mol.nAtoms(); ++i) r_[i]=PTable::covalentRadius(mol.atom(i).an())*rscale;
		//calculate the atom-in-molecule radius
		if(DEBUG_QEQ3>1) std::cout<<"Calculating the atom-in-molecule radius...\n";
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
					drAtom_[i]=rNew.norm()/r00.norm();
				}
			}
		}
		//print the radius
		if(DEBUG_QEQ3>0){
			std::cout<<"RADIUS:\n";
			for(unsigned int i=0; i<mol.nAtoms(); ++i){
				std::cout<<mol.atom(i).name()<<" "<<r0_[i]<<" "<<r_[i]<<" "<<drAtom_[i]<<"\n";
			}
		}
		
		//calculate coulomb integrals
		if(DEBUG_QEQ3>0) std::cout<<"Calculating Coulomb Integrals...\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			for(unsigned int j=i+1; j<mol.nAtoms(); ++j){
				J_(j,i)=ke*ICC::iTensor[icc_](
					mol.atom(j).posn()-mol.atom(i).posn(),
					k_*ICC::scale[icc_](J_(i,i),J_(j,j))
				);
			}
		}
		//calculate the operator matrix
		if(DEBUG_QEQ3>0) std::cout<<"Calculating operator matrix...\n";
		A_=J_.selfadjointView<Eigen::Lower>();
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			A_(0,i)=1;
			for(unsigned int j=1; j<mol.nAtoms(); ++j){
				A_(j,i)-=J_(i,0);
			}
		}
		//calculate the solution vector
		if(DEBUG_QEQ3>0) std::cout<<"Calculating solution vector...\n";
		if(scaler_) for(unsigned int i=0; i<mol.nAtoms(); ++i) chi_[i]*=1.0/drAtom_[i];
		if(DEBUG_QEQ3>0){
			for(unsigned int i=0; i<mol.nAtoms(); ++i){
				std::cout<<"chi["<<i<<"] = "<<chi_[i]<<"\n";
			}
		}
		b_=Eigen::VectorXd::Constant(mol.nAtoms(),chi_[0]);
		b_.noalias()-=chi_;
		b_[0]=qTot_;
		//solve the linear equations
		if(DEBUG_QEQ3>0) std::cout<<"Solving the linear equations...\n";
		x_.noalias()=A_.partialPivLu().solve(b_);
		//set the atomic charges
		if(DEBUG_QEQ3>0) std::cout<<"Setting the charge...\n";
		for(unsigned int i=0; i<mol.nAtoms(); ++i) mol.atom(i).charge()=x_[i];
	}
}

template <class AtomT>
void QEQ3::qt(unsigned int t, SimAtomic<AtomT>& sim, const Ewald3D::Coulomb& ewald){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::qt(SimAtomic<AtomT>&,const Ewald3D::Coulomb&):\n";
	if(sim.nAtoms()==0) throw std::invalid_argument("No atoms given for charge transfer.");
	else if(sim.nAtoms()==1) sim.atom(t,0).charge()=0;
	else {
		//local variables
		Eigen::Vector3d dr;
		
		//resize matrices
		if(DEBUG_QEQ3>0) std::cout<<"Resizing matrices...\n";
		J_=Eigen::MatrixXd::Constant(sim.nAtoms(),sim.nAtoms(),0);
		A_.resize(sim.nAtoms(),sim.nAtoms());
		chi_.resize(sim.nAtoms());
		
		//local variables
		double cscale=0.0,escale=0.0;
		if(units::consts::system()==units::System::AU){
			cscale=1.0;
			escale=units::HARTREEpEV;
		}
		else if(units::consts::system()==units::System::METAL){
			cscale=units::ANGpBOHR;
			escale=1.0;
		}
		else throw std::runtime_error("Invalid units.");
		
		//calculate atom info
		if(DEBUG_QEQ3>0) std::cout<<"Calculating atom info...\n";
		if(DEBUG_QEQ3>1) std::cout<<"Element AN CHI J00 C\n";
		for(unsigned int i=0; i<sim.nAtoms(); ++i){
			//electronegativity
			if(chiType_==PTable::ElectronegativityType::PAULING) chi_[i]=PTable::electronegativityPauling(sim.atom(t,i).an());
			else if(chiType_==PTable::ElectronegativityType::ALLEN) chi_[i]=PTable::electronegativityAllen(sim.atom(t,i).an());
			else if(chiType_==PTable::ElectronegativityType::MULLIKEN) chi_[i]=PTable::electronegativityMulliken(sim.atom(t,i).an());
			else if(chiType_==PTable::ElectronegativityType::HINZE) chi_[i]=PTable::electronegativityHinze(sim.atom(t,i).an(),0);
			if(chiType_==PTable::ElectronegativityType::MULLIKEN && sim.atom(t,i).an()==1) chi_[i]=4.528;
			chi_[i]*=escale;
			//idempotential
			if(chiType_==PTable::ElectronegativityType::HINZE) J_(i,i)=PTable::idempotentialHinze(sim.atom(t,i).an(),0);
			else J_(i,i)=PTable::idempotential(sim.atom(t,i).an());
			J_(i,i)*=escale;
			if(DEBUG_QEQ3>1){
				std::cout<<sim.atom(t,i).name()<<" ";
				std::cout<<sim.atom(t,i).an()<<" ";
				std::cout<<chi_[i]<<" ";
				std::cout<<J_(i,i)<<"\n";
			}
		}
		
		//set the idempotential
		if(DEBUG_QEQ3>0) std::cout<<"Setting the idempotential...\n";
		for(unsigned int i=0; i<sim.nAtoms(); ++i) J_(i,i)+=ewald.phiSelf()*ke;
		//calculate coulomb integrals
		if(DEBUG_QEQ3>0) std::cout<<"Calculating Coulomb Integrals...\n";
		for(unsigned int i=0; i<sim.nAtoms(); ++i){
			for(unsigned int j=i+1; j<sim.nAtoms(); ++j){
				Cell::diff(sim.atom(t,i).posn(),sim.atom(t,j).posn(),dr,sim.cell(t).R(),sim.cell(t).RInv());
				double a=k_*ICC::scale[icc_](J_(i,i),J_(j,j));
				J_(j,i)=ke*(
					ICC::iTensor[icc_](dr,a)//short range
					+ewald.phi(dr)-1.0/dr.norm()//long range, short range removed
				);
			}
		}
		//calculate the operator matrix
		A_=J_.selfadjointView<Eigen::Lower>();
		for(unsigned int i=0; i<sim.nAtoms(); ++i){
			A_(0,i)=1;
			for(unsigned int j=1; j<sim.nAtoms(); ++j){
				A_(j,i)-=J_(i,0);
			}
		}
		//calculate the solution vector
		b_=Eigen::VectorXd::Constant(sim.nAtoms(),chi_[0]);
		b_.noalias()-=chi_;
		b_[0]=qTot_;
		//solve the linear equations
		if(DEBUG_QEQ3>0) std::cout<<"Solving the linear equations...\n";
		x_.noalias()=A_.partialPivLu().solve(b_);
		//set the atomic charges
		if(DEBUG_QEQ3>0) std::cout<<"Setting the charge...\n";
		for(unsigned int i=0; i<sim.nAtoms(); ++i) sim.atom(t,i).charge()=x_[i];
	}
}

template <class AtomT>
void QEQ3::qt_jZero(unsigned int t, SimAtomic<AtomT>& sim, const Ewald3D::Coulomb& ewald){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::qt_jZero(SimAtomic<AtomT>&,const Ewald3D::Coulomb&):\n";
	if(sim.nAtoms()==0) throw std::invalid_argument("No atoms given for charge transfer.");
	else if(sim.nAtoms()==1) sim.atom(t,0).charge()=0;
	else {
		//local variables
		Eigen::Vector3d dr;
		
		//resize matrices
		if(DEBUG_QEQ3>0) std::cout<<"Resizing matrices...\n";
		J_=Eigen::MatrixXd::Constant(sim.nAtoms(),sim.nAtoms(),0);
		A_.resize(sim.nAtoms(),sim.nAtoms());
		chi_.resize(sim.nAtoms());
		
		//local variables
		double cscale=0.0,escale=0.0;
		if(units::consts::system()==units::System::AU){
			cscale=1.0;
			escale=units::HARTREEpEV;
		}
		else if(units::consts::system()==units::System::METAL){
			cscale=units::ANGpBOHR;
			escale=1.0;
		}
		else throw std::runtime_error("Invalid units.");
		
		//calculate atom info
		if(DEBUG_QEQ3>0) std::cout<<"Calculating atom info...\n";
		if(DEBUG_QEQ3>1) std::cout<<"Element AN CHI J00 C\n";
		for(unsigned int i=0; i<sim.nAtoms(); ++i){
			//electronegativity
			if(chiType_==PTable::ElectronegativityType::PAULING) chi_[i]=PTable::electronegativityPauling(sim.atom(t,i).an());
			else if(chiType_==PTable::ElectronegativityType::ALLEN) chi_[i]=PTable::electronegativityAllen(sim.atom(t,i).an());
			else if(chiType_==PTable::ElectronegativityType::MULLIKEN) chi_[i]=PTable::electronegativityMulliken(sim.atom(t,i).an());
			else if(chiType_==PTable::ElectronegativityType::HINZE) chi_[i]=PTable::electronegativityHinze(sim.atom(t,i).an(),0);
			if(chiType_==PTable::ElectronegativityType::MULLIKEN && sim.atom(t,i).an()==1) chi_[i]=4.528;
			chi_[i]*=escale;
			//idempotential
			J_(i,i)=sim.atom(t,i).jzero();
			if(DEBUG_QEQ3>1){
				std::cout<<sim.atom(t,i).name()<<" ";
				std::cout<<sim.atom(t,i).an()<<" ";
				std::cout<<chi_[i]<<" ";
				std::cout<<J_(i,i)<<"\n";
			}
		}
		
		//set the idempotential
		if(DEBUG_QEQ3>0) std::cout<<"Setting the idempotential...\n";
		for(unsigned int i=0; i<sim.nAtoms(); ++i) J_(i,i)+=ewald.phiSelf()*ke;
		//calculate coulomb integrals
		if(DEBUG_QEQ3>0) std::cout<<"Calculating Coulomb Integrals...\n";
		for(unsigned int i=0; i<sim.nAtoms(); ++i){
			for(unsigned int j=i+1; j<sim.nAtoms(); ++j){
				Cell::diff(sim.atom(t,i).posn(),sim.atom(t,j).posn(),dr,sim.cell(t).R(),sim.cell(t).RInv());
				double a=k_*ICC::scale[icc_](J_(i,i),J_(j,j));
				J_(j,i)=ke*(
					ICC::iTensor[icc_](dr,a)//short range
					+ewald.phi(dr)-1.0/dr.norm()//long range, short range removed
				);
			}
		}
		//calculate the operator matrix
		A_=J_.selfadjointView<Eigen::Lower>();
		for(unsigned int i=0; i<sim.nAtoms(); ++i){
			A_(0,i)=1;
			for(unsigned int j=1; j<sim.nAtoms(); ++j){
				A_(j,i)-=J_(i,0);
			}
		}
		//calculate the solution vector
		b_=Eigen::VectorXd::Constant(sim.nAtoms(),chi_[0]);
		b_.noalias()-=chi_;
		b_[0]=qTot_;
		//solve the linear equations
		if(DEBUG_QEQ3>0) std::cout<<"Solving the linear equations...\n";
		x_.noalias()=A_.partialPivLu().solve(b_);
		//set the atomic charges
		if(DEBUG_QEQ3>0) std::cout<<"Setting the charge...\n";
		for(unsigned int i=0; i<sim.nAtoms(); ++i) sim.atom(t,i).charge()=x_[i];
	}
}

template <class AtomT>
void QEQ3::init(SimAtomic<AtomT>& sim){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::init(SimAtomic<AtomT>&):\n";
	//local function variables
	Eigen::Vector3d R;
	
	//resize matrices
	if(DEBUG_QEQ3>0) std::cout<<"Resizing matrices...\n";
	J_=Eigen::MatrixXd::Constant(sim.nAtoms(),sim.nAtoms(),0);
	A_.resize(sim.nAtoms(),sim.nAtoms());
	chi_.resize(sim.nAtoms());
	I_.resize(sim.nAtoms());
	
	//local variables
	double cscale=0.0,escale=0.0;
	if(units::consts::system()==units::System::AU){
		cscale=1.0;
		escale=units::HARTREEpEV;
	}
	else if(units::consts::system()==units::System::METAL){
		cscale=units::ANGpBOHR;
		escale=1.0;
	}
	else throw std::runtime_error("Invalid units.");
	
	//calculate atom info
	if(DEBUG_QEQ3>0) std::cout<<"Calculating atom info...\n";
	if(DEBUG_QEQ3>1) std::cout<<"Element AN CHI J00 C\n";
	unsigned int cc=0;
	for(unsigned int n=0; n<sim.nSpecies(); ++n){
		unsigned int an=PTable::an(sim.atomNames(n).c_str());
		for(unsigned int m=0; m<sim.nAtoms(n); ++m){
			//electronegativity
			if(chiType_==PTable::ElectronegativityType::PAULING) chi_[cc]=PTable::electronegativityPauling(an);
			else if(chiType_==PTable::ElectronegativityType::ALLEN) chi_[cc]=PTable::electronegativityAllen(an);
			else if(chiType_==PTable::ElectronegativityType::MULLIKEN) chi_[cc]=PTable::electronegativityMulliken(an);
			else if(chiType_==PTable::ElectronegativityType::HINZE) chi_[cc]=PTable::electronegativityHinze(an,0);
			if(chiType_==PTable::ElectronegativityType::MULLIKEN && an==1) chi_[cc]=4.528;
			chi_[cc]*=escale;
			//idempotential
			if(chiType_==PTable::ElectronegativityType::HINZE) I_(cc)=PTable::idempotentialHinze(an,0);
			else I_(cc)=PTable::idempotential(an);
			I_(cc)*=escale;
			if(DEBUG_QEQ3>1){
				std::cout<<sim.atomNames(n)<<" ";
				std::cout<<an<<" ";
				std::cout<<chi_[cc]<<" ";
				std::cout<<I_(cc)<<"\n";
			}
			++cc;
		}
	}
	
	//calculate the solution vector
	if(DEBUG_QEQ3>0) std::cout<<"Calculating the solution vector...\n";
	b_=Eigen::VectorXd::Constant(sim.nAtoms(),chi_[0]);
	b_.noalias()-=chi_;
	b_[0]=qTot_;
}

template <class AtomT> 
void QEQ3::qt_cont(unsigned int t, SimAtomic<AtomT>& sim, const Ewald3D::Coulomb& ewald){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::qt_cont(SimAtomic<AtomT>&,const Ewald3D::Coulomb&):\n";
	//local variables
	Eigen::Vector3d dr;
	const double ke=units::consts::ke();
	
	//set the idempotential
	if(DEBUG_QEQ3>0) std::cout<<"Setting the idempotential...\n";
	for(unsigned int i=0; i<sim.nAtoms(); ++i) J_(i,i)=I_[i]+ewald.phiSelf()*ke;
	//calculate coulomb integrals
	if(DEBUG_QEQ3>0) std::cout<<"Calculating Coulomb Integrals...\n";
	for(unsigned int i=0; i<sim.nAtoms(); ++i){
		for(unsigned int j=i+1; j<sim.nAtoms(); ++j){
			Cell::diff(sim.atom(t,i).posn(),sim.atom(t,j).posn(),dr,sim.cell(t).R(),sim.cell(t).RInv());
			double a=k_*ICC::scale[icc_](J_(i,i),J_(j,j));
			J_(j,i)=ke*(
				ICC::iTensor[icc_](dr,a)//short range
				+ewald.phi(dr)-1.0/dr.norm()//long range, short range removed
			);
		}
	}
	//calculate the operator matrix
	if(DEBUG_QEQ3>0) std::cout<<"Calculating the operator matrix...\n";
	A_=J_.selfadjointView<Eigen::Lower>();
	for(unsigned int i=0; i<sim.nAtoms(); ++i){
		A_(0,i)=1;
		for(unsigned int j=1; j<sim.nAtoms(); ++j){
			A_(j,i)-=J_(i,0);
		}
	}
	//solve the linear equations
	if(DEBUG_QEQ3>0) std::cout<<"Solving the linear equations...\n";
	x_.noalias()=A_.partialPivLu().solve(b_);
	//set the atomic charges
	if(DEBUG_QEQ3>0) std::cout<<"Setting the charge...\n";
	for(unsigned int i=0; i<sim.nAtoms(); ++i) sim.atom(t,i).charge()=x_[i];
}

template <class AtomT> 
void QEQ3::qt_jZero_cont(unsigned int t, SimAtomic<AtomT>& sim, const Ewald3D::Coulomb& ewald){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::qt_jZero_cont(SimAtomic<AtomT>&,const Ewald3D::Coulomb&):\n";
	//local variables
	Eigen::Vector3d dr;
	const double ke=units::consts::ke();
	
	//set the idempotential
	if(DEBUG_QEQ3>0) std::cout<<"Setting the idempotential...\n";
	for(unsigned int i=0; i<sim.nAtoms(); ++i) J_(i,i)=sim.atom(t,i).jzero()+ewald.phiSelf()*ke;
	//calculate coulomb integrals
	if(DEBUG_QEQ3>0) std::cout<<"Calculating Coulomb Integrals...\n";
	for(unsigned int i=0; i<sim.nAtoms(); ++i){
		const double Jii=J_(i,i);
		for(unsigned int j=i+1; j<sim.nAtoms(); ++j){
			const double Jjj=J_(j,j);
			Cell::diff(sim.atom(t,i).posn(),sim.atom(t,j).posn(),dr,sim.cell(t).R(),sim.cell(t).RInv());
			const double a=k_*ICC::scale[icc_](Jii,Jjj);
			J_(j,i)=ke*(
				ICC::iTensor[icc_](dr,a)//short range
				+ewald.phi(dr)-1.0/dr.norm()//long range, short range removed
			);
		}
	}
	//calculate the operator matrix
	if(DEBUG_QEQ3>0) std::cout<<"Calculating the operator matrix...\n";
	A_=J_.selfadjointView<Eigen::Lower>();
	for(unsigned int i=0; i<sim.nAtoms(); ++i){
		A_(0,i)=1;
		for(unsigned int j=1; j<sim.nAtoms(); ++j){
			A_(j,i)-=J_(i,0);
		}
	}
	//solve the linear equations
	if(DEBUG_QEQ3>0) std::cout<<"Solving the linear equations...\n";
	x_.noalias()=A_.partialPivLu().solve(b_);
	//set the atomic charges
	if(DEBUG_QEQ3>0) std::cout<<"Setting the charge...\n";
	for(unsigned int i=0; i<sim.nAtoms(); ++i) sim.atom(t,i).charge()=x_[i];
}

#endif