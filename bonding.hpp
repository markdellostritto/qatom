#ifndef BONDING_HPP
#define BONDING_HPP

//c libraries
#include <cstdlib>
#include <cmath>
//c++ libraries
#include <vector>
#include <iostream>
#include <stdexcept>
//local 
#include "log.hpp"
#include "file.hpp"
#include "cell.hpp"
#include "label.hpp"
#include "sim.hpp"
//ame
#include "math_const.hpp"

#ifndef DEBUG_BOND
#define DEBUG_BOND 0
#endif

//*************************************************************************************************
//Bond State struct
//*************************************************************************************************

template <class AtomT>
class BondState{
private:
	unsigned int nBonds_;
	std::vector<AtomT> atoms_;
public:
	//constructors/destructors
	BondState():nBonds_(0){};
	BondState(unsigned int bondOrder);
	BondState(const BondState& b);
	~BondState(){};
	
	//access
	unsigned int& nBonds(){return nBonds_;};
	const unsigned int& nBonds()const{return nBonds_;};
	std::vector<AtomT>& atoms(){return atoms_;};
	const std::vector<AtomT>& atoms()const{return atoms_;};
	AtomT& atom(unsigned int i){return atoms_[i];};
	const AtomT& atom(unsigned int i)const{return atoms_[i];};
	
	//operators
	BondState<AtomT>& operator=(const BondState<AtomT>& b);
	template <class AtomT_> friend std::ostream& operator<<(std::ostream& out, const BondState<AtomT_>& b);
	
	//members
	void clear();
	void reset();
};

//constructors/destructors

template <class AtomT>
BondState<AtomT>::BondState(unsigned int bondOrder){
	nBonds_=0;
	atoms_.resize(bondOrder);
	for(unsigned int i=0; i<atoms_.size(); ++i) atoms_[i].init();
}

template <class AtomT>
BondState<AtomT>::BondState(const BondState<AtomT>& b){
	nBonds_=b.nBonds();
	atoms_=b.atoms();
}

//operators

template <class AtomT>
BondState<AtomT>& BondState<AtomT>::operator=(const BondState<AtomT>& b){
	nBonds_=b.nBonds();
	atoms_=b.atoms();
	return *this;
}

template <class AtomT>
std::ostream& operator<<(std::ostream& out, const BondState<AtomT>& b){
	out<<"C-Bonds = (";
	for(unsigned int i=0; i<b.nBonds_-1; ++i) out<<b.atoms_[i].name()<<b.atoms_[i].index()+1<<",";
	out<<b.atoms_.back().name()<<b.atoms_.back().index()+1<<"\n";
	return out;
}

//member functions

template <class AtomT>
void BondState<AtomT>::clear(){
	nBonds_=0;
	atoms_.clear();
}

template <class AtomT>
void BondState<AtomT>::reset(){
	nBonds_=0;
	for(unsigned int i=0; i<atoms_.size(); ++i){
		atoms_[i].clear();
		atoms_[i].init();
	}
}

/*class Bonding
	This class describes bonding between atoms of a simulation.  The point of this class is to store information
relevant to both covalent and hydrogen bonds, such as bond lengths and angles used to determine whether a bond exists
between two atoms.  
	Obviously, to determine whether two atoms are bonded or not takes a great deal of information, 
including not only bond lengths and angles but the positions of all atoms for a given timestep as well as the shape
and size of the simulation cell.  Therefore, this class is most useful when inherited by a Simulation class (see below),
however, it can be used as a standalone class and can be quite useful as one (especially when parallel processing requires
separate sets of temporary objects for determining bonding states of atoms).
	The atoms bound to a given atom are determined when that atom is passed to one of the functions
of this class.  The bound atoms are then stored within the class and can be accessed directly.
	The bond parameters can be set manually, although this would be cumbersome and inflexible.  Instead, the bonding
parameters can be read from file.  The format of the parameters are as follow below.  The tags specifying paremters
are case insensitive.  The presence of extra spaces does not affect the format.  Parameters in quotation marks must be
replaced by strings, removing the quotation marks, and parameters in parentheses must be replaced by numbers, removing
the parentheses.
	Covalent Bond Lengths:
		BondLengths = {
			"specie1" - "specie2" : (bondlength1)
			"specie3" - "specie4" : (bondlength2)
			...
		}
	H-Bond Lengths:
		HBondLengths = {
			"specie1" : (hBondLength1)
			"specie2" : (hBondLength2)
		}
	H-Bond Angle:
		HBondAngle = (angle)
	Bond Orders:
		MaxBondOrders = {
			"specie1" : (bondOrder1)
			"specie2" : (bondOrder2)
			...
		}
*/

//*************************************************************************************************
//Bonding class
//*************************************************************************************************

class Bonding{
private:
	unsigned int nSpecies_;//the number of atomic species
	std::vector<std::string> atoms_;//atom names
	
	//covalent bonding parameters
	std::vector<std::vector<double> > bondLengths_;//the bond lengths for each specie
	std::vector<unsigned int> bondOrders_;//the maximum bond orders for each specie
	
	//h-bonding parameters
	double hBondAngle_;//the maximum donor-H-acceptor angle of the H-bond
	std::vector<double> hBondLengths_;//the H-bond lengths for each of the species
public:
	//constructors/destructors
	Bonding():nSpecies_(0),hBondAngle_(0){};
	Bonding(const Bonding& b);
	~Bonding(){};
	
	//operators
	Bonding& operator=(const Bonding& b);
	friend std::ostream& operator<<(std::ostream& out, const Bonding& b);
	
	//access
	unsigned int nSpecies()const{return nSpecies_;};
	const std::vector<std::vector<double> >& bondLengths()const{return bondLengths_;};
	double& bondLength(const std::string& name1, const std::string& name2);
	const double& bondLength(const std::string& name1, const std::string& name2)const;
	double& bondLength(int i, int j){return bondLengths_[i][j];};
	const double& bondLength(int i, int j)const{return bondLengths_[i][j];};
	const std::vector<unsigned int>& bondOrders()const{return bondOrders_;};
	unsigned int& bondOrder(int i){return bondOrders_[i];};
	const unsigned int& bondOrder(int i)const{return bondOrders_[i];};
	double& hBondAngle(){return hBondAngle_;};
	const double& hBondAngle()const{return hBondAngle_;};
	const std::vector<double>& hBondLengths()const{return hBondLengths_;};
	double& hBondLength(int i){return hBondLengths_[i];};
	const double& hBondLength(int i)const{return hBondLengths_[i];};
	
	//member functions
	void resize(unsigned int nSpecies);
	void loadBondParams(const char* fileName, const std::vector<std::string>& atomNames);
	void clear();
	
	//static functions
	template <class AtomT> static unsigned int bondOrder(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b);
	template <class AtomT> static unsigned int bondOrder(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b, const std::vector<AtomT>& subset);
	template <class AtomT> static std::vector<int>& bondOrder(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b, std::vector<int>& bOrder);
	template <class AtomT> static BondState<AtomT>& cBonds(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b, BondState<AtomT>& bondState);
	template <class AtomT> static BondState<AtomT>& cBonds(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b, BondState<AtomT>& bondState, const std::vector<AtomT>& subset);
	template <class AtomT> static AtomT& hBond(unsigned int ts, const AtomT& hAtom, const AtomT& donator, const SimAtomic<AtomT>& sim, const Bonding& b, AtomT& hBoundAtom);
	template <class AtomT> static AtomT& hBond(unsigned int ts, const AtomT& hAtom, const AtomT& donator, const SimAtomic<AtomT>& sim, const Bonding& b, AtomT& hBoundAtom, std::vector<AtomT>& subset);
};

//*************************************************************************************************
//Bonding Functions
//*************************************************************************************************

template <class AtomT>
unsigned int Bonding::bondOrder(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b){
	if(DEBUG_BOND>0) std::cout<<"Finding bond order of: "<<atom.name()<<atom.index()+1<<"\n";
	//local function variables
	unsigned int nBonds=0;
	double dist;
	Eigen::Vector3d temp;
	std::vector<double> minDist(b.bondOrder(atom.specie()),1000);
	
	//loop over all species
	for(unsigned int j=0; j<sim.nSpecies(); ++j){
		//skip species if instructed
		if(b.bondLength(atom.specie(),j)<=num_const::ZERO) continue;
		//loop over all atoms of the current specie
		for(unsigned int k=0; k<sim.nAtoms(j); ++k){
			//skip self-interactions
			if(j==atom.specie() && k==atom.index()) continue;
			if(DEBUG_BOND>1) std::cout<<"\tAtom: "<<sim.atom(ts,j,k).name()<<sim.atom(ts,j,k).index()+1<<"\n";
			//find the distance b/w the current H atom and the current atom of the current specie
			dist=Cell::dist(
				sim.atom(ts,atom.specie(),atom.index()).posn(),sim.atom(ts,j,k).posn(),
				temp,sim.cell(ts).R(),sim.cell(ts).RInv()
			);
			if(DEBUG_BOND>1) std::cout<<"\t\tdist="<<dist<<", bondLength="<<b.bondLength(atom.specie(),j)<<"\n";
			//check if it's within a bond length
			if(dist<=b.bondLength(atom.specie(),j)){
				if(DEBUG_BOND>1) std::cout<<"\t\tFound covalent bond.\n";
				//if this is a min dist, record the covalently bonded atom
				for(unsigned int l=0; l<b.bondOrder(atom.specie()); ++l){
					if(dist<minDist[l]){
						int m=b.bondOrder(atom.specie())-1;
						while(m>l && m>0){minDist[m]=minDist[m-1]; --m;}
						minDist[m]=dist;
						if(nBonds<b.bondOrder(atom.specie())) ++nBonds;
						break;
					}
				}
			}
		}
	}
	
	return nBonds;
}

template <class AtomT>
unsigned int Bonding::bondOrder(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b, const std::vector<AtomT>& subset){
	if(DEBUG_BOND>0) std::cout<<"Finding covalent bonds of: "<<atom<<"\n";
	//local function variables
	unsigned int nBonds=0;
	double dist;
	Eigen::Vector3d temp;
	std::vector<double> minDist(b.bondOrder(atom.specie()),1000);
	
	//loop over the subset
	for(unsigned int i=0; i<subset.size(); ++i){
		//skip species if instructed
		if(b.bondLength(subset[i].specie(),atom.specie())<=num_const::ZERO) continue;
		//skip self-interactions
		if(subset[i].specie()==atom.specie() && subset[i].index()==atom.index()) continue;
		if(DEBUG_BOND>1) std::cout<<"Atom: "<<sim.atom(ts,subset[i].specie(),subset[i].index())<<"\n";
		//find the distance b/w the current H atom and the current atom of the current specie()
		dist=Cell::dist(
			sim.atom(ts,atom.specie(),atom.index()).posn(),
			sim.atom(ts,subset[i].specie(),subset[i].index()).posn(),
			temp,sim.cell(ts).R(),sim.cell(ts).RInv()
		);
		if(DEBUG_BOND>1) std::cout<<"dist="<<dist<<"\n";
		//check if it's within a bond length
		if(dist<=b.bondLength(atom.specie(),subset[i].specie())){
			//if this is a min dist, record the covalently bonded atom
			if(DEBUG_BOND>1) std::cout<<"Found covalent bond.\n";
			for(unsigned int l=0; l<b.bondOrder(atom.specie()); ++l){
				if(dist<minDist[l]){
					int m=b.bondOrder(atom.specie())-1;
					while(m>l && m>0){minDist[m]=minDist[m-1]; --m;}
					minDist[m]=dist;
					if(nBonds<b.bondOrder(atom.specie())) ++nBonds;
					break;
				}
			}
		}
	}
	
	return nBonds;
}

template <class AtomT>
std::vector<int>& Bonding::bondOrder(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b, std::vector<int>& bOrder){
	if(DEBUG_BOND>0) std::cout<<"bondOrder(std::vector<int>&,const AtomT&,const Cell&,const StructAtomic<AtomT>&):\n";
	//local function variables
	double dist;
	Eigen::Vector3d temp;
	std::vector<double> minDist(b.bondOrder(atom.specie()),1000);
	std::vector<int> boundSpecies(b.bondOrder(atom.specie()),-1);
	
	//reset the bond order
	for(unsigned int i=0; i<bOrder.size(); ++i) bOrder[i]=0;
	
	//loop over all species
	for(unsigned int j=0; j<sim.nSpecies(); ++j){
		//skip species if instructed
		if(b.bondLength(atom.specie(),j)<=num_const::ZERO) continue;
		//loop over all atoms of the current specie
		for(unsigned int k=0; k<sim.nAtoms(j); ++k){
			//skip self-interactions
			if(j==atom.specie() && k==atom.index()) continue;
			if(DEBUG_BOND>0) std::cout<<"\tAtom: "<<sim.atom(ts,j,k)<<"\n";
			//find the distance b/w the current H atom and the current atom of the current specie
			dist=Cell::dist(
				sim.atom(ts,atom.specie(),atom.index()).posn(),sim.atom(ts,j,k).posn(),
				temp,sim.cell(ts).R(),sim.cell(ts).RInv()
			);
			if(DEBUG_BOND>1) std::cout<<"\t\tdist="<<dist<<", bondLength="<<b.bondLength(atom.specie(),j)<<"\n";
			//check if it's within a bond length
			if(dist<=b.bondLength(atom.specie(),j)){
				if(DEBUG_BOND>1) std::cout<<"\t\tFound covalent bond.\n";
				//if this is a min dist, record the covalently bonded atom
				for(unsigned int l=0; l<b.bondOrder(atom.specie()); ++l){
					if(dist<minDist[l]){
						int m=b.bondOrder(atom.specie())-1;
						while(m>l && m>0){
							boundSpecies[m]=boundSpecies[m-1];
							minDist[m]=minDist[m-1]; --m;
						}
						boundSpecies[m]=sim.atom(ts,j,k).specie();
						minDist[m]=dist;
						break;
					}
				}
			}
		}
	}
	
	for(unsigned int i=0; i<boundSpecies.size(); ++i){
		if(boundSpecies[i]>=0) ++bOrder[boundSpecies[i]];
		else break;
	}
	
	return bOrder;
}

template <class AtomT>
BondState<AtomT>& Bonding::cBonds(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b, BondState<AtomT>& bondState){
	if(DEBUG_BOND>0) std::cout<<"Finding covalent bonds of: "<<atom.name()<<atom.index()+1<<"\n";
	//local function variables
	double dist;
	Eigen::Vector3d temp;
	std::vector<double> minDist(b.bondOrder(atom.specie()),1000);
	std::vector<AtomT> boundAtoms(b.bondOrder(atom.specie()));
	for(unsigned int i=0; i<boundAtoms.size(); ++i) boundAtoms[i].init();
	for(unsigned int i=0; i<boundAtoms.size(); ++i) boundAtoms[i].name()=std::string("NULL");
	for(unsigned int i=0; i<boundAtoms.size(); ++i) boundAtoms[i].an()=0;
	
	//loop over all species
	if(DEBUG_BOND>0) std::cout<<"Finding the bound atoms...\n";
	for(unsigned int j=0; j<sim.nSpecies(); ++j){
		//skip species if instructed
		if(b.bondLength(atom.specie(),j)<=num_const::ZERO) continue;
		//loop over all atoms of the current specie
		for(unsigned int k=0; k<sim.nAtoms(j); ++k){
			//skip self-interactions
			if(j==atom.specie() && k==atom.index()) continue;
			if(DEBUG_BOND>1) std::cout<<"\tAtom: "<<sim.atom(ts,j,k).name()<<sim.atom(ts,j,k).index()+1<<"\n";
			//find the distance b/w the current H atom and the current atom of the current specie
			dist=Cell::dist(
				sim.atom(ts,atom.specie(),atom.index()).posn(),sim.atom(ts,j,k).posn(),
				temp,sim.cell(ts).R(),sim.cell(ts).RInv()
			);
			if(DEBUG_BOND>1) std::cout<<"\t\tdist="<<dist<<", bondLength="<<b.bondLength(atom.specie(),j)<<"\n";
			//check if it's within a bond length
			if(dist<=b.bondLength(atom.specie(),j)){
				if(DEBUG_BOND>1) std::cout<<"\t\tFound covalent bond.\n";
				//if this is a min dist, record the covalently bonded atom
				for(unsigned int l=0; l<b.bondOrder(atom.specie()); ++l){
					if(dist<minDist[l]){
						int m=b.bondOrder(atom.specie())-1;
						while(m>l && m>0){
							boundAtoms[m].name()=boundAtoms[m-1].name();
							boundAtoms[m].an()=boundAtoms[m-1].an();
							boundAtoms[m].specie()=boundAtoms[m-1].specie();
							boundAtoms[m].index()=boundAtoms[m-1].index();
							minDist[m]=minDist[m-1]; --m;
						}
						boundAtoms[m].name()=sim.atom(ts,j,k).name();
						boundAtoms[m].an()=sim.atom(ts,j,k).an();
						boundAtoms[m].specie()=sim.atom(ts,j,k).specie();
						boundAtoms[m].index()=sim.atom(ts,j,k).index();
						minDist[m]=dist;
						break;
					}
				}
			}
		}
	}
	
	if(DEBUG_BOND>0) std::cout<<"Setting the number of bonds in the bondState...\n";
	bool breakEarly=false;
	for(unsigned int i=0; i<boundAtoms.size(); ++i){
		if(boundAtoms[i].an()==0){bondState.nBonds()=i+1; breakEarly=true; break;}
	}
	if(!breakEarly) bondState.nBonds()=boundAtoms.size();
	if(DEBUG_BOND>0) std::cout<<"N-Bonds = "<<bondState.nBonds()<<"\n";
	if(bondState.nBonds()>bondState.atoms().size()) bondState.atoms().resize(bondState.nBonds());
	if(DEBUG_BOND>0) std::cout<<"Setting the atoms in the bondState...\n";
	for(unsigned int i=0; i<bondState.nBonds(); ++i){
		bondState.atom(i).init();
		bondState.atom(i).name()=boundAtoms[i].name();
		bondState.atom(i).specie()=boundAtoms[i].specie();
		bondState.atom(i).index()=boundAtoms[i].index();
	}
}

template <class AtomT>
BondState<AtomT>& Bonding::cBonds(unsigned int ts, const AtomT& atom, const SimAtomic<AtomT>& sim, const Bonding& b, BondState<AtomT>& bondState, const std::vector<AtomT>& subset){
	if(DEBUG_BOND>0) std::cout<<"Finding covalent bonds of: "<<atom.name()<<atom.index()+1<<"\n";
	//local function variables
	double dist;
	Eigen::Vector3d temp;
	std::vector<double> minDist(b.bondOrder(atom.specie()),1000);
	std::vector<AtomT> boundAtoms(b.bondOrder(atom.specie()));
	for(unsigned int i=0; i<boundAtoms.size(); ++i) boundAtoms[i].init();
	for(unsigned int i=0; i<boundAtoms.size(); ++i) boundAtoms[i].name()=std::string("NULL");
	
	//loop over the subset
	for(unsigned int i=0; i<subset.size(); ++i){
		//skip species if instructed
		if(b.bondLength(atom.specie(),subset[i].specie())==0) continue;
		//skip self-interactions
		if(subset[i].specie()==atom.specie() && subset[i].index()==atom.index()) continue;
		if(DEBUG_BOND>1) std::cout<<"Atom: "<<sim.atom(ts,subset[i].specie(),subset[i].index()).name()<<sim.atom(ts,subset[i].specie(),subset[i].index()).index()+1<<"\n";
		//find the distance b/w the current H atom and the current atom of the current specie()
		dist=Cell::dist(
			sim.atom(ts,atom.specie(),atom.index()).posn(),
			sim.atom(ts,subset[i].specie(),subset[i].index()).posn(),
			temp,sim.cell(ts).R(),sim.cell(ts).RInv()
		);
		if(DEBUG_BOND>1) std::cout<<"dist="<<dist<<"\n";
		//check if it's within a bond length
		if(dist<=b.bondLength(atom.specie(),subset[i].specie())){
			//if this is a min dist, record the covalently bonded atom
			if(DEBUG_BOND>1) std::cout<<"Found covalent bond.\n";
			for(unsigned int l=0; l<b.bondOrder(atom.specie()); ++l){
				if(dist<minDist[l]){
					int m=b.bondOrder(atom.specie()); --m;
					while(m>l && m>0){
						boundAtoms[m].name()=boundAtoms[m-1].name();
						boundAtoms[m].specie()=boundAtoms[m-1].specie();
						boundAtoms[m].index()=boundAtoms[m-1].index();
						minDist[m]=minDist[m-1];
						--m;
					}
					boundAtoms[m].name()=sim.atom(ts,subset[i].specie(),subset[i].index()).name();
					boundAtoms[m].specie()=sim.atom(ts,subset[i].specie(),subset[i].index()).specie();
					boundAtoms[m].index()=sim.atom(ts,subset[i].specie(),subset[i].index()).index();
					minDist[m]=dist;
					break;
				}
			}
		}
	}
	
	if(DEBUG_BOND>0) std::cout<<"Setting the number of bonds in the bondState...\n";
	for(unsigned int i=0; i<boundAtoms.size(); ++i){
		if(boundAtoms[i].name()==std::string("NULL")){bondState.nBonds()=i; break;}
	}
	if(bondState.nBonds()>bondState.atoms().size()) bondState.atoms().resize(bondState.nBonds());
	if(DEBUG_BOND>0) std::cout<<"nBonds = "<<bondState.nBonds()<<"\n";
	if(DEBUG_BOND>0) std::cout<<"Setting the atoms in the bondState...\n";
	for(unsigned int i=0; i<bondState.nBonds(); ++i){
		bondState.atom(i).init();
		bondState.atom(i).name()=boundAtoms[i].name();
		bondState.atom(i).specie()=boundAtoms[i].specie();
		bondState.atom(i).index()=boundAtoms[i].index();
	}
}

template <class AtomT>
AtomT& Bonding::hBond(unsigned int ts, const AtomT& hAtom, const AtomT& donator, const SimAtomic<AtomT>& sim, const Bonding& b, AtomT& hBoundAtom){
	//local function variables
	double dist;
	double minDist=1000;
	Eigen::Vector3d temp;
	AtomT hBoundAtom_;
	hBoundAtom_.init();
	bool hBond=false;
	
	if(DEBUG_BOND>0) std::cout<<"H Atom = "<<hAtom<<"\n";
	//find the nearest neighbor that could possibly be a H-bonded atom
	//loop over all species
	for(unsigned int j=0; j<sim.nSpecies(); ++j){
		//skip species if instructed
		if(b.hBondLength(j)<=num_const::ZERO) continue;
		//loop over all atoms of the current specie()
		for(unsigned int k=0; k<sim.nAtoms(j); ++k){
			//skip the donating atom, to which it is assumed the H-atom is already covalently bonded
			if(j==donator.specie() && k==donator.index()) continue;
			//skip self-interactions
			if(j==hAtom.specie() && k==hAtom.index()) continue;
			if(DEBUG_BOND>1) std::cout<<"Atom: "<<sim.atom(ts,j,k).name()<<sim.atom(ts,j,k).index()+1<<"\n";
			//find the distance between the H atom and the current atom
			dist=Cell::dist(
				sim.atom(ts,hAtom.specie(),hAtom.index()).posn(),sim.atom(ts,j,k).posn(),
				temp,sim.cell(ts).R(),sim.cell(ts).RInv()
			);
			if(DEBUG_BOND>1) std::cout<<"dist = "<<dist<<"\n";
			//check if it's within an H-bond length
			if(dist<=b.hBondLength(j)){
				if(DEBUG_BOND>1) std::cout<<"Found possible H-bond.\n";
				//check if it's a min dist
				if(dist<minDist){
					//hBoundAtom_=sim.atom(ts,j,k);
					hBoundAtom_.name()=sim.atom(ts,j,k).name();
					hBoundAtom_.specie()=sim.atom(ts,j,k).specie();
					hBoundAtom_.index()=sim.atom(ts,j,k).index();
					minDist=dist;
					hBond=true;
				}
			}
		}
	}
	
	//if we found a possible candidate, check the H-bond angle
	if(hBond){
		if(DEBUG_BOND>0) std::cout<<"Possibly H-bound atom: "<<sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index())<<"\n";
		Eigen::Vector3d vec1,vec2;
		double angle;
		Cell::diff(//H-atom - donator atom
			sim.atom(ts,hAtom.specie(),hAtom.index()).posn(),//H-atom
			sim.atom(ts,donator.specie(),donator.index()).posn(),//donator atom
			vec1,sim.cell(ts).R(),sim.cell(ts).RInv()
		);
		Cell::diff(//possible H-bound atom - donator atom
			sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index()).posn(),//possible H-bound atom
			sim.atom(ts,donator.specie(),donator.index()).posn(),//donator atom
			vec2,sim.cell(ts).R(),sim.cell(ts).RInv()
		);
		angle=std::acos(vec1.dot(vec2)/(vec1.norm()*vec2.norm()));
		if(DEBUG_BOND>1) std::cout<<"H-bond angle: "<<angle<<"\n";
		//check if the angle is less than the max H-bond angle
		if(angle>b.hBondAngle()){
			hBoundAtom_.clear();//if it's too large, erase the H-bonded atom
			hBond=false;
		}
	}
	
	if(hBond){
		//hBoundAtom=sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index());
		hBoundAtom.name()=sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index()).name();
		hBoundAtom.specie()=sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index()).specie();
		hBoundAtom.index()=sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index()).index();
	}
	return hBoundAtom;
}

template <class AtomT>
AtomT& Bonding::hBond(unsigned int ts, const AtomT& hAtom, const AtomT& donator, const SimAtomic<AtomT>& sim, const Bonding& b, AtomT& hBoundAtom, std::vector<AtomT>& subset){
	//local function variables
	double dist;
	double minDist=1000;
	Eigen::Vector3d temp;
	AtomT hBoundAtom_;
	
	//find the nearest neighbor that could possibly be a H-bonded atom
	//loop over the subset
	for(unsigned int i=0; i<subset.size(); ++i){
		//skip species if instructed
		if(b.hBondLength(subset[i].specie())<=num_const::ZERO) continue;
		//skip the donating atom, to which it is assumed the H-atom is already covalently bonded
		if(subset[i].specie()==donator.specie() && subset[i].index()==donator.index()) continue;
		//skip self-interactions
		if(subset[i].specie()==hAtom.specie() && subset[i].index()==hAtom.index()) continue;
		//find the distance between the H atom and the current atom
		dist=Cell::dist(
			sim.atom(ts,hAtom.specie(),hAtom.index()).posn(),
			sim.atom(ts,subset[i].specie(),subset[i].index()).posn(),
			temp,sim.cell(ts).R(),sim.cell(ts).R
		);
		if(DEBUG_BOND>1) std::cout<<"Atom: "<<sim.atom(ts,subset[i].specie(),subset[i].index())<<"\n";
		//check if it's within an H-bond length
		if(dist<=b.hBondLength(subset[i].specie())){
			if(DEBUG_BOND>1) std::cout<<"Found possible H-bond.\n";
			//check if it's a min dist
			if(dist<minDist){
				hBoundAtom_=sim.atom(ts,subset[i].specie(),subset[i].index());
				minDist=dist;
			}
		}
	}
	
	//if we found a possible candidate, check the H-bond angle
	if(hBoundAtom_.specie()>=0){
		if(DEBUG_BOND>0) std::cout<<"Possibly H-bound atom: "<<sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index())<<"\n";
		Eigen::Vector3d vec1,vec2;
		double angle;
		Cell::diff(//H-atom - donator atom
			sim.atom(ts,hAtom.specie(),hAtom.index()).posn(),//H-atom
			sim.atom(ts,donator.specie(),donator.index()).posn(),//donator atom
			vec1,sim.cell(ts).R(),sim.cell(ts).RInv()
		);
		Cell::diff(//possible H-bound atom - donator atom
			sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index()).posn(),//possible H-bound atom
			sim.atom(ts,donator.specie(),donator.index()).posn(),//donator atom
			vec2,sim.cell(ts).R(),sim.cell(ts).RInv()
		);
		angle=std::acos(vec1.dot(vec2)/(vec1.norm()*vec2.norm()));
		//check if the angle is less than the max H-bond angle
		if(DEBUG_BOND>1) std::cout<<"H-bond angle: "<<angle<<"\n";
		if(angle>b.hBondAngle()){
			//if it's too large, erase the H-bonded atom
			hBoundAtom_.clear();
		}
	}
	
	hBoundAtom=sim.atom(ts,hBoundAtom_.specie(),hBoundAtom_.index());
	return hBoundAtom;
}

#endif
