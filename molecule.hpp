#ifndef MOLECULE_HPP
#define MOLECULE_HPP

//c++ libraries
#include <ostream>
#include <string>
#include <vector>
//Eigen
#include <Eigen/Dense>
#include <Eigen/StdVector>
//local
#include "cell.hpp"

static const double ke=1;

template <class AtomT, class... Args>
class Molecule: public Args...{
private:
	std::vector<AtomT> atoms_;
public:
	/* constructors/destructors */
	Molecule():Args()...{};
	Molecule(const Molecule<AtomT,Args...>& m);
	~Molecule(){};
	
	/* operators */
	Molecule<AtomT,Args...>& operator=(const Molecule<AtomT,Args...>& m);
	template <class AtomT_,class... Args_> friend std::ostream& operator<<(std::ostream& out, const Molecule<AtomT_,Args_...>& m);
	
	/* access */
	//atoms
		unsigned int nAtoms()const{return atoms_.size();};
		std::vector<AtomT>& atoms(){return atoms_;};
		const std::vector<AtomT>& atoms()const{return atoms_;};
		AtomT& atom(unsigned int i){return atoms_[i];};
		const AtomT& atom(unsigned int i)const{return atoms_[i];};
	
	/* member functions */
	void clear();
	void init();
	
	/* adding/removing atoms */
	//adding atoms
		void add();
		void add(AtomT& atom);
		void add(const AtomT& atom);
	//removal
		bool remove(const AtomT& atom);
		
	/* static functions */
	//template <class AtomT_,class... Args_> Eigen::Vector3d& dipole(const Molecule<AtomT_,Args_...>& mol, Eigen::Vector3d& dipole);
};

/* constructors/destructors */

template <class AtomT, class... Args>
Molecule<AtomT,Args...>::Molecule(const Molecule<AtomT,Args...>& m):Args(m)...{
	atoms_=m.atoms();
}

/* operators */

//assignment

template <class AtomT, class... Args>
Molecule<AtomT,Args...>& Molecule<AtomT,Args...>::operator=(const Molecule<AtomT,Args...>& m){
	atoms_=m.atoms();
	int arr[sizeof...(Args)]={(Args::clear(),0)...};
}

/* member functions */

template <class AtomT, class... Args>
void Molecule<AtomT,Args...>::clear(){
	atoms_.clear();
	int arr[sizeof...(Args)]={(Args::clear(),0)...};
}

template <class AtomT, class... Args>
void Molecule<AtomT,Args...>::init(){
	int arr[sizeof...(Args)]={(Args::init(),0)...};
}

/* adding/removing atoms */

template <class AtomT, class... Args>
void Molecule<AtomT,Args...>::add(){
	atoms_.push_back(AtomT());
	atoms_.back().init();
};

template <class AtomT, class... Args>
void Molecule<AtomT,Args...>::add(AtomT& atom){
	atoms_.push_back(atom);
};

template <class AtomT, class... Args>
void Molecule<AtomT,Args...>::add(const AtomT& atom){
	atoms_.push_back(atom);
};

template <class AtomT, class... Args>
bool Molecule<AtomT,Args...>::remove(const AtomT& atom){
	for(unsigned int i=0; i<atoms_.size(); ++i){
		if(atom==atoms_[i]){
			atoms_[i]=atoms_.back();
			atoms_.pop_back();
			return true;
		}
	}
	return false;
}

/* static functions */

template <class MoleculeT>
Eigen::Vector3d& dipole(const MoleculeT& mol, Eigen::Vector3d& dipole){
	dipole.setZero();
	for(unsigned int n=0; n<mol.nAtoms(); ++n){
		dipole.noalias()+=mol.atom(n).posn()*mol.atom(n).charge();
	}
	return dipole;
}

template <class MoleculeT>
Eigen::Vector3d& dipole2(const MoleculeT& mol, Eigen::Vector3d& dipole){
	dipole.setZero();
	for(unsigned int n=0; n<mol.nAtoms(); ++n){
		dipole.noalias()+=mol.atom(n).posn()*mol.atom(n).charge();
	}
	return dipole;
}

template <class MoleculeT>
Eigen::Vector3d& dipole(const MoleculeT& mol, const Cell& cell, Eigen::Vector3d& dipole){
	Eigen::Vector3d temp;
	dipole.setZero();
	//add the difference between the ith atom and the molecule posn multiplied by charge w/0.r.t the cell
	for(unsigned int i=0; i<mol.nAtoms(); ++i){
		dipole.noalias()+=Cell::diff(mol.atom(i).posn(),mol.posn(),temp,cell.R(),cell.RInv())*mol.atom(i).charge();
	}
	return dipole;
}

template <class MoleculeT>
Eigen::Vector3d& coc(const MoleculeT& mol, Eigen::Vector3d& coc){
	coc.setZero();
	for(unsigned int n=0; n<mol.nAtoms(); ++n){
		coc.noalias()+=mol.atom(n).posn()*std::fabs(mol.atom(n).charge());
	}
	return coc;
}

template <class MoleculeT>
Eigen::Vector3d& coc(const MoleculeT& mol, const Cell& cell, Eigen::Vector3d& coc){
	Eigen::Vector3d temp;
	coc.setZero();
	double charge=std::fabs(mol.atom(0).charge());
	//find difference between atom 0 and atom i w.r.t. the cell, multiply by mass w/0.r.t the cell, and add to the center of mass
	for(unsigned int i=1; i<mol.nAtoms(); ++i){
		coc.noalias()+=Cell::diff(mol.atom(i).posn(),mol.atom(0).posn(),temp,cell.R(),cell.RInv())*std::fabs(mol.atom(i).charge());
		charge+=std::fabs(mol.atom(i).charge());
	}
	coc/=charge;
	//finally, add the coc w.r.t. atom 0 to the posn of atom 0
	coc.noalias()+=mol.atom(0).posn();
	//return the final com to within the cell
	Cell::returnToCell(coc,coc,cell.R(),cell.RInv());
	return coc;
}

template <class AtomT, class... Args>
Eigen::Vector3d& com(const Molecule<AtomT,Args...>& mol, Eigen::Vector3d& com){
	com.setZero();
	double mass=0;
	for(unsigned int n=0; n<mol.nAtoms(); ++n){
		com.noalias()+=mol.atom(n).posn()*mol.atom(n).mass();
		mass+=mol.atom(n).mass();
	}
	com.noalias()/=mass;
	return com;
}

template <class AtomT,class... Args> 
Eigen::Vector3d& com(const Molecule<AtomT,Args...>& mol, const Cell& cell, Eigen::Vector3d& com){
	Eigen::Vector3d temp;
	com.setZero();
	double mass=mol.atom(0).mass();
	//find difference between atom 0 and atom i w.r.t. the cell, multiply by mass w/0.r.t the cell, and add to the center of mass
	for(unsigned int i=1; i<mol.nAtoms(); ++i){
		com.noalias()+=Cell::diff(mol.atom(i).posn(),mol.atom(0).posn(),temp,cell.R(),cell.RInv())*mol.atom(i).mass();
		mass+=mol.atom(i).mass();
	}
	com/=mass;
	//finally, add the com w.r.t. atom 0 to the posn of atom 0
	com.noalias()+=mol.atom(0).posn();
	//return the final com to within the cell
	Cell::returnToCell(com,com,cell.R(),cell.RInv());
	return com;
}

#endif
