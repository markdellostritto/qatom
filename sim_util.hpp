#ifndef SIM_UTIL_HPP
#define SIM_UTIL_HPP

#include <iostream>
#include <cstdlib>
#include "sim.hpp"
#include "cell.hpp"
#include "atom.hpp"
#include "bonding.hpp"

#ifndef DEBUG_SIM_UTIL
#define DEBUG_SIM_UTIL 0
#endif

namespace sim_util{

//********************************************************************************************
//Assigning bonds - no monitoring
//********************************************************************************************

/*
Storage requires the same number of molecules in every timestep.  Thus, we do not monitor bonding
but assume that the number of bonds and the bonds between atoms do not change over the course of the
simulation.  The bonds are found only in the first timestep and then propagated for all subsequent timesteps.
*/
template <class AtomT, class MoleculeT>
void assignBonds(SimMol<AtomT,MoleculeT>& simMol, SimAtomic<AtomT>& simAtom, const Bonding& bonding){
	if(DEBUG_SIM_UTIL>0) std::cout<<"assignBonds(SimMol<AtomT,MoleculeT>&,SimAtomic<AtomT>&,const Bonding& bonding):\n";
	//local function variables
	Eigen::Vector3d vecTemp;
	MoleculeT molecule;
	std::vector<MoleculeT> moleculeList;
	std::vector<std::vector<BondState<AtomT> > > bondState(simAtom.nSpecies());
	for(unsigned int n=0; n<simAtom.nSpecies(); ++n){
		bondState[n].resize(simAtom.nAtoms(n));
	}
	
	//copy the atomic information from one simulation to the other
	if(DEBUG_SIM_UTIL>0) std::cout<<"Copying atoms from atomic simulation to molecular simulation...\n";
	if(&simMol!=&simAtom) static_cast<SimAtomicS<AtomT>&>(simMol).operator=(static_cast<SimAtomicS<AtomT>&>(simAtom));
	
	/* Calculate the bonds of the system for the first timestep */
	//calculate the bonding state of each atom
	if(DEBUG_SIM_UTIL>0) std::cout<<"Calculating the bond state of each atom...\n";
	for(unsigned int n=0; n<simMol.nSpecies(); ++n){
		for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
			//find the bond state of each atom
			if(DEBUG_SIM_UTIL>1) std::cout<<"Atom: "<<simMol.atom(0,n,m).name()<<simMol.atom(0,n,m).index()+1<<"\n";
			Bonding::cBonds(0,simMol.atom(0,n,m),simMol,bonding,bondState[n][m]);
		}
	}
	//find all bonds in the system
	if(DEBUG_SIM_UTIL>0) std::cout<<"Finding all bonds in the system...\n";
	for(unsigned int n=0; n<simMol.nSpecies(); ++n){
		for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
			if(DEBUG_SIM_UTIL>1) std::cout<<"Atom: "<<simMol.atom(0,n,m).name()<<simMol.atom(0,n,m).index()+1<<"\n";
			//loop over all bonds of the current atom
			for(unsigned int i=0; i<bondState[n][m].nBonds(); ++i){
				int bs=bondState[n][m].atom(i).specie();
				int bi=bondState[n][m].atom(i).index();
				if(DEBUG_SIM_UTIL>1) std::cout<<"Bound atom: "<<simMol.atom(0,bs,bi).name()<<simMol.atom(0,bs,bi).index()+1<<"\n";
				//check atoms with indices higher than the current index (prevents double counting of bonds)
				if(n<bs || (n==bs && m<bi)){
					if(DEBUG_SIM_UTIL>1) std::cout<<"Found bond...\n";
					for(int j=0; j<bondState[bs][bi].nBonds(); ++j){
						//if (n,m)->(bs,bi) and (n,m)<-(bs,bi), then we have a bond
						if(bondState[bs][bi].atom(j).specie()==n && bondState[bs][bi].atom(j).index()==m){
							molecule.clear();
							//make sure non-const, no copying, just transfer pointer
							molecule.add(simMol.atom(0,n,m));
							molecule.add(simMol.atom(0,bs,bi));
							moleculeList.push_back(molecule);
							if(DEBUG_SIM_UTIL>3){
								std::cout<<"\tBond=\n"<<moleculeList.back()<<"\n";
								std::cout<<"\tCOM="<<com(moleculeList.back(),simMol.cell(0),vecTemp).transpose()<<"\n";
							}
							break;
						}
					}
				}
			}
		}
	}
	
	//resize the simulation
	if(DEBUG_SIM_UTIL>0) std::cout<<"Resizing the molecular simulation...\n";
	if(DEBUG_SIM_UTIL>0) std::cout<<"N-Molecules = "<<moleculeList.size()<<"\n";
	simMol.resize(simAtom.timesteps(),moleculeList.size());
	
	//now that we have the molecules for the first timestep, we set the positions for the molecules for all subsequent timesteps
	if(DEBUG_SIM_UTIL>0) std::cout<<"Setting molecules in the simulation...\n";
	for(unsigned int t=0; t<simMol.timesteps(); ++t){
		for(unsigned int n=0; n<simMol.nMol(); ++n){
			simMol.molecule(t,n).posn().setZero();
			for(unsigned int m=0; m<moleculeList[n].nAtoms(); ++m){
				simMol.molecule(t,n).add(simMol.atom(t,moleculeList[n].atom(m).specie(),moleculeList[n].atom(m).index()));
			}
			com(simMol.molecule(t,n),simMol.cell(t),simMol.molecule(t,n).posn());
		}
	}
}

/*
Storage requires the same number of molecules in every timestep.  Thus, we do not monitor bonding
but assume that the number of bonds and the bonds between atoms do not change over the course of the
simulation.  The bonds are found only in the first timestep and then propagated for all subsequent timesteps.
*/
template <class AtomT, class MoleculeT, class AtomS>
void assignBonds(SimMol<AtomT,MoleculeT>& simMol, SimAtomic<AtomT>& simAtom, const Bonding& bonding, const std::vector<AtomS>& subset){
	if(DEBUG_SIM_UTIL>0)std::cout<<"assignBonds(SimMol<AtomT,MoleculeT>&,SimAtomic<AtomT>&,const Bonding&,const std::vector<AtomS>&):\n";
	//local function variables
	MoleculeT molecule;
	std::vector<MoleculeT> moleculeList;
	std::vector<std::vector<BondState<AtomT> > > bondState(simAtom.nSpecies());
	for(unsigned int n=0; n<simAtom.nSpecies(); ++n){
		bondState[n].resize(simAtom.nAtoms(n));
	}
	
	//copy the atomic information from one simulation to the other
	if(DEBUG_SIM_UTIL>0) std::cout<<"Copying atoms from atomic simulation to molecular simulation...\n";
	if(&simMol!=&simAtom) static_cast<SimAtomic<AtomT>&>(simMol).operator=(simAtom);
	
	/* Calculate the bonds of the system for the first timestep */
	if(DEBUG_SIM_UTIL>0) std::cout<<"Calculating the bonds for the first timestep...\n";
	//calculate the bonding state of each atom
	for(unsigned int n=0; n<simMol.nSpecies(); ++n){
		for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
			//find the bond state of each atom
			std::cout<<simMol.atom(0,n,m).name()<<simMol.atom(0,n,m).index()+1<<"\n";
			Bonding::cBonds(0,simMol.atom(0,n,m),simMol,bonding,bondState[n][m]);
		}
	}
	//find all bonds in the system for those atoms in the subset
	for(unsigned int n=0; n<subset.size(); ++n){
		//loop over all bonds of the current atom
		for(int i=0; i<bondState[subset[n].specie()][subset[n].index()].nBonds(); ++i){
			int bs=bondState[subset[n].specie()][subset[n].index()].atom(i).specie();
			int bi=bondState[subset[n].specie()][subset[n].index()].atom(i).index();
			std::cout<<"bs = "<<bs<<", bi = "<<bi<<"\n";
			//check atoms with indices higher than the current index
			//(prevents double counting of bonds)
			if(subset[n].specie()<bs || (subset[n].specie()==bs && subset[n].index()<bi)){
				for(unsigned int j=0; j<bondState[bs][bi].nBonds(); ++j){
					//if (n,m)->(bs,bi) and (n,m)<-(bs,bi), then we have a bond
					if(bondState[bs][bi].atom(j).specie()==subset[n].specie() && bondState[bs][bi].atom(j).index()==subset[n].index()){
						molecule.clear();
						//make sure non-const, no copying, just transfer pointer
						molecule.add(simMol.atom(0,subset[n].specie(),subset[n].index()));
						molecule.add(simMol.atom(0,bs,bi));
						moleculeList.push_back(molecule);
						if(DEBUG_SIM_UTIL>3){
							std::cout<<"\tBond=\n"<<moleculeList.back()<<"\n";
							std::cout<<"\tCOM="<<com(moleculeList.back(),static_cast<Cell&>(simMol)).transpose()<<"\n";
						}
						break;
					}
				}
			}
		}
	}
	
	//resize the simulation
	if(DEBUG_SIM_UTIL>0) std::cout<<"Resizing the molecular simulation...\n";
	simMol.resize(simAtom.timesteps(),moleculeList.size());
	
	//now that we have the molecules for the first timestep, we set the positions for the molecules for all subsequent timesteps
	if(DEBUG_SIM_UTIL>0) std::cout<<"Setting molecules in the simulation...\n";
	for(unsigned int t=0; t<simMol.timesteps(); ++t){
		for(unsigned int n=0; n<simMol.nMol(); ++n){
			for(unsigned int m=0; m<moleculeList[n].nAtoms(); ++m){
				//make sure non-const, no copying, just transfer pointer
				simMol.molecule(t,n).add(simMol.atom(t,moleculeList[n].atom(m).specie(),moleculeList[n].atom(m).index()));
			}
			simMol.molecule(t,n).posn()=com(simMol.molecule(t,n),static_cast<Cell&>(simMol));
		}
	}
}

//********************************************************************************************
//Assigning bonds - monitoring
//********************************************************************************************

/*
Storage requires the same number of molecules in every timestep.  We do allow bonds to break however.
The number of bonds is found in the first timestep.  For each subsequent timestep, the bonds between atoms
are found, but the total number is assumed to remain constant.  The order of the bonds in the list is
not necessarily maintained, such that the same bond may appear in different places in the list of all bonds
throughout the simulation.
*/
template <class AtomT, class MoleculeT>
void assignBondsMonitor(SimMol<AtomT,MoleculeT>& simMol, SimAtomic<AtomT>& simAtom, const Bonding& bonding){
	if(DEBUG_SIM_UTIL>0) std::cout<<"assignBondsMonitorNew(SimMol<AtomT,MoleculeT>&,SimAtomic<AtomT>&,const Bonding& bonding):\n";
	//local function variables
	MoleculeT molecule;
	std::vector<MoleculeT> moleculeList;
	std::vector<std::vector<BondState<AtomT> > > bondState(simAtom.nSpecies());
	for(unsigned int n=0; n<simAtom.nSpecies(); ++n){
		bondState[n].resize(simAtom.nAtoms(n));
	}
	Eigen::Vector3d tempVec;
	
	/* Copy atomic information from one simulation to the other */
	if(DEBUG_SIM_UTIL>0) std::cout<<"Copying atoms from atomic simulation to molecular simulation...\n";
	if(&simMol!=&simAtom) static_cast<SimAtomic<AtomT>&>(simMol).operator=(simAtom);
	
	/* Calculate the bonds of the system for the first timestep */
	//calculate the bonding state of each atom
	for(int n=0; n<simMol.nSpecies(); n++){
		for(int m=0; m<simMol.nAtoms(n); m++){
			//find the bond state of each atom
			Bonding::cBonds(0,simMol.atom(0,n,m),simMol,bonding,bondState[n][m]);
		}
	}
	//find all bonds in the system
	for(unsigned int n=0; n<simMol.nSpecies(); ++n){
		for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
			//loop over all bonds of the current atom
			for(unsigned int i=0; i<bondState[n][m].nBonds(); ++i){
				int bs=bondState[n][m].atom(i).specie();
				int bi=bondState[n][m].atom(i).index();
				//check atoms with indices higher than the current index
				//(prevents double counting of bonds)
				if(n<bs || (n==bs && m<bi)){
					for(int j=0; j<bondState[bs][bi].nBonds(); ++j){
						//if (n,m)->(bs,bi) and (n,m)<-(bs,bi), then we have a bond
						if(bondState[bs][bi].atom(j).specie()==n && bondState[bs][bi].atom(j).index()==m){
							molecule.clear();
							//make sure non-const, no copying, just transfer pointer
							molecule.add(simMol.atom(0,n,m));
							molecule.add(simMol.atom(0,bs,bi));
							moleculeList.push_back(molecule);
							if(DEBUG_SIM_UTIL>3){
								std::cout<<"\tBond=\n"<<moleculeList.back()<<"\n";
								std::cout<<"\tCOM="<<com(moleculeList.back(),simMol.cell(0),tempVec).transpose()<<"\n";
							}
							break;
						}
					}
				}
			}
		}
	}
	
	//resize the simulation
	if(DEBUG_SIM_UTIL>0) std::cout<<"Resizing the molecular simulation...\n";
	simMol.resize(simAtom.timesteps(),moleculeList.size());
	
	/*
		Calculate the bonds of the system for all timesteps
	*/
	if(DEBUG_SIM_UTIL>0) std::cout<<"Finding the bonds in the simulation...\n";
	for(unsigned int t=0; t<simMol.timesteps(); ++t){
		if(DEBUG_SIM_UTIL>1) std::cout<<"Timestep: "<<t<<"\n";
		else if(t%1000==0) std::cout<<"Timestep: "<<t<<"\n";
		int molIndex=0;
		//calculate the bonding state of each atom
		for(int n=0; n<simMol.nSpecies(); n++){
			for(int m=0; m<simMol.nAtoms(n); m++){
				//find the bond state of each atom
				Bonding::cBonds(t,simMol.atom(t,n,m),simMol,bonding,bondState[n][m]);
			}
		}
		//find all bonds in the system
		for(unsigned int n=0; n<simMol.nSpecies(); ++n){
			for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
				//loop over all bonds of the current atom
				for(unsigned int i=0; i<bondState[n][m].nBonds(); ++i){
					int bs=bondState[n][m].atom(i).specie();
					int bi=bondState[n][m].atom(i).index();
					//check atoms with indices higher than the current index
					//(prevents double counting of bonds)
					if(n<bs || (n==bs && m<bi)){
						for(int j=0; j<bondState[bs][bi].nBonds(); ++j){
							//if (n,m)->(bs,bi) and (n,m)<-(bs,bi), then we have a bond
							if(bondState[bs][bi].atom(j).specie()==n && bondState[bs][bi].atom(j).index()==m){
								//make sure non-const, no copying, just transfer pointer
								simMol.molecule(t,molIndex).add(simMol.atom(t,n,m));
								simMol.molecule(t,molIndex).add(simMol.atom(t,bs,bi));
								//set the position to the center of mass
								com(simMol.molecule(t,molIndex),simMol.cell(t),simMol.molecule(t,molIndex).posn());
								//increment the index
								++molIndex;
								if(DEBUG_SIM_UTIL>3){
									std::cout<<"\tBond=\n"<<moleculeList.back()<<"\n";
									std::cout<<"\tCOM="<<com(moleculeList.back(),simMol.cell(t),tempVec).transpose()<<"\n";
								}
								break;
							}
						}
					}
				}
			}
		}
		if(molIndex!=moleculeList.size()){
			std::cout<<"ERROR in assignBondsMonitorNew(SimMol<AtomT,MoleculeT>&,SimAtomic<AtomT>&,const Bonding& bonding):\n";
			std::cout<<"Found "<<molIndex<<" bonds in timestep "<<t<<" when number should be "<<moleculeList.size()<<"\n";
			throw std::out_of_range("Index out of range.");
		}
	}
}

/*
Storage requires the same number of molecules in every timestep.  We do allow bonds to break however.
The number of bonds is found in the first timestep.  For each subsequent timestep, the bonds between atoms
are found, but the total number is assumed to remain constant.  The order of the bonds in the list is
not necessarily maintained, such that the same bond may appear in different places in the list of all bonds
throughout the simulation.
*/
template <class AtomT, class MoleculeT, class AtomS>
void assignBondsMonitor(SimMol<AtomT,MoleculeT>& simMol, SimAtomic<AtomT>& simAtom, const Bonding& bonding, const std::vector<AtomS>& subset){
	if(DEBUG_SIM_UTIL>0)std::cout<<"assignBonds(SimMol<AtomT,MoleculeT>&,SimAtomic<AtomT>&,const Bonding&,const std::vector<AtomS>&):\n";
	//local function variables
	MoleculeT molecule;
	std::vector<MoleculeT> moleculeList;
	std::vector<std::vector<BondState<AtomT> > > bondState(simAtom.nSpecies());
	for(int n=0; n<simAtom.nSpecies(); ++n){
		bondState[n].resize(simAtom.nAtoms(n));
	}
	
	/*
		Copy atomic information from one simulation to the other
	*/
	if(DEBUG_SIM_UTIL>0) std::cout<<"Copying atoms from atomic simulation to molecular simulation...\n";
	if(&simMol!=&simAtom) static_cast<SimAtomic<AtomT>&>(simMol).operator=(simAtom);
	
	/*
		Calculate the bonds of the system for the first timestep
	*/
	//calculate the bonding state of each atom
	for(unsigned int n=0; n<simMol.nSpecies(); ++n){
		for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
			//find the bond state of each atom
			Bonding::cBonds(0,simMol.atom(0,n,m),simMol,bonding,bondState[n][m]);
		}
	}
	//find all bonds in the system for those atoms in the subset
	for(unsigned int n=0; n<subset.size(); ++n){
		//loop over all bonds of the current atom
		for(unsigned int i=0; i<bondState[subset[n].specie()][subset[n].index()].nBonds(); ++i){
			int bs=bondState[subset[n].specie()][subset[n].index()].atom(i).specie();
			int bi=bondState[subset[n].specie()][subset[n].index()].atom(i).index();
			//check atoms with indices higher than the current index
			//(prevents double counting of bonds)
			if(subset[n].specie()<bs || (subset[n].specie()==bs && subset[n].index()<bi)){
				//check if the bonded atom has a bond to the current atom
				for(int j=0; j<bondState[bs][bi].nBonds(); ++j){
					//if (n,m)->(bs,bi) and (n,m)<-(bs,bi), then we have a bond
					if(bondState[bs][bi].atom(j).specie()==subset[n].specie() && bondState[bs][bi].atom(j).index()==subset[n].index()){
						molecule.clear();
						//make sure non-const, no copying, just transfer pointer
						molecule.add(simMol.atom(0,subset[n].specie(),subset[n].index()));
						molecule.add(simMol.atom(0,bs,bi));
						moleculeList.push_back(molecule);
						if(DEBUG_SIM_UTIL>2){
							std::cout<<"\tBond="<<moleculeList.back();
							std::cout<<",("<<moleculeList.back().atom(0).name()<<moleculeList.back().atom(0).index();
							std::cout<<","<<moleculeList.back().atom(1).name()<<moleculeList.back().atom(1).index()<<")\n";
							//std::cout<<",COM="<<com(moleculeList.back(),static_cast<Cell&>(simMol)).transpose()<<"\n";
						}
						break;
					}
				}
			}
		}
	}
	
	//resize the simulation
	if(DEBUG_SIM_UTIL>0) std::cout<<"Resizing the molecular simulation...\n";
	simMol.resize(simAtom.timesteps(),moleculeList.size());
	
	/*
		Calculate the bonds of the system for all timesteps
	*/
	if(DEBUG_SIM_UTIL>0) std::cout<<"Finding the bonds in the simulation...\n";
	for(unsigned int t=0; t<simMol.timesteps(); ++t){
		if(DEBUG_SIM_UTIL>1) std::cout<<"Timestep: "<<t<<"\n";
		else if(t%1000==0) std::cout<<"Timestep: "<<t<<"\n";
		int molIndex=0;
		//calculate the bonding state of each atom
		for(unsigned int n=0; n<simMol.nSpecies(); ++n){
			for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
				//find the bond state of each atom
				Bonding::cBonds(t,simMol.atom(t,n,m),simMol,bonding,bondState[n][m]);
			}
		}
		//find all bonds in the system for those atoms in the subset
		for(unsigned int n=0; n<subset.size(); ++n){
			//loop over all bonds of the current atom
			for(unsigned int i=0; i<bondState[subset[n].specie()][subset[n].index()].nBonds(); ++i){
				int bs=bondState[subset[n].specie()][subset[n].index()].atom(i).specie();
				int bi=bondState[subset[n].specie()][subset[n].index()].atom(i).index();
				//check atoms with indices higher than the current index
				//(prevents double counting of bonds)
				if(subset[n].specie()<bs || (subset[n].specie()==bs && subset[n].index()<bi)){
					//check if the bonded atom has a bond to the current atom
					for(unsigned int j=0; j<bondState[bs][bi].nBonds(); ++j){
						//if (n,m)->(bs,bi) and (n,m)<-(bs,bi), then we have a bond
						if(bondState[bs][bi].atom(j).specie()==subset[n].specie() && bondState[bs][bi].atom(j).index()==subset[n].index()){
							//make sure non-const, no copying, just transfer pointer
							simMol.molecule(t,molIndex).add(simMol.atom(t,subset[n].specie(),subset[n].index()));
							simMol.molecule(t,molIndex).add(simMol.atom(t,bs,bi));
							//set the position to the center of mass
							com(simMol.molecule(t,molIndex),simMol.cell(t),simMol.molecule(t,molIndex).posn());
							//increment the index
							++molIndex;
							break;
						}
					}
				}
			}
		}
		if(molIndex!=moleculeList.size()){
			std::cout<<"ERROR in assignBonds(SimMol<AtomT,MoleculeT>&,SimAtomic<AtomT>&,const Bonding&,const std::vector<AtomS>&):\n";
			std::cout<<"Found "<<molIndex<<" bonds in timestep "<<t<<" when number should be "<<moleculeList.size()<<"\n";
			throw std::out_of_range("Index out of range.");
		}
	}
}

}

#endif