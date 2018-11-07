#ifndef ELECTROSTATICS_UTILITY_HPP
#define ELECTROSTATICS_UTILITY_HPP

//c libraries
#include <cmath>
//Eigen
#include <Eigen/Dense>
//string
#include "string.hpp"
#include "file.hpp"
//simulation
#include "sim.hpp"
#include "atom.hpp"
#include "molecule.hpp"
#include "bonding.hpp"
//chem
#include "ptable.hpp"

#ifndef DEBUG_ELEC_UTIL
#define DEBUG_ELEC_UTIL 0
#endif

//*****************************************************************************************
//Utility_Electrostatics
//*****************************************************************************************

struct Electrostatics{
	struct Load{
		template <class MoleculeT, class AtomT> static void alpha(const char* file, MoleculeT& mol);
		template <class AtomT> static void alpha(const char* file, SimAtomic<AtomT>& sim);
		template <class AtomT> static void charge(const char* file, SimAtomic<AtomT>& sim);
	};
	struct Print{
		template <class AtomT> static void chg_tot(const char* file, const SimAtomic<AtomT>& sim);
		template <class AtomT> static void chg(const char* file, const SimAtomic<AtomT>& sim);
		template <class AtomT> static void alpha_tot(const char* file, const SimAtomic<AtomT>& sim);
		template <class AtomT> static void dipole_tot(const char* file, const SimAtomic<AtomT>& sim, const Bonding& bonding);
		template <class MoleculeT, class AtomT> static void dipole_tot(const char* file, const SimMol<AtomT,MoleculeT>& sim);
	};
};

//Load

template <class MoleculeT, class AtomT>
void Electrostatics::Load::alpha(const char* file, MoleculeT& mol){
	if(DEBUG_ELEC_UTIL>0) std::cout<<"Electrostatics::Load::alpha<MoleculeT,AtomT>(const char*,MoleculeT&):\n";
	//local variables
	char* input=(char*)malloc(sizeof(char)*string::M);
	char* temp=(char*)malloc(sizeof(char)*string::M);
	FILE* reader=NULL;
	bool error=false;
	std::vector<std::pair<unsigned int,double> > alphas;
	
	try{
		//open the file
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O ERROR: Could not open file.");
		
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_all(string::trim_right(input,string::COMMENT));
			string::to_upper(string::copy_left(temp,input,"="));
			if(std::strcmp(temp,"ALPHA_ATOM")==0){
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strpbrk(input,"}")!=NULL) break;//break if find "}"
					std::string name=std::string(std::strtok(input,string::WS));
					unsigned int an=PTable::an(name.c_str());
					double alpha=std::atof(std::strtok(NULL,"="));
					alphas.push_back(std::pair<unsigned int,double>(an,alpha));
				}
			}
		}
		fclose(reader);
		reader==NULL;
		
		for(unsigned int i=0; i<mol.nAtoms(); ++i){
			for(unsigned int j=0; j<alphas.size(); ++j){
				if(mol.atom(i).an()==alphas[j].first){
					mol.atom(i).alpha()=Eigen::Matrix3d::Identity()*alphas[j].second;
					break;
				}
			}
		}
	}catch(std::exception& e){
		std::cout<<"ERROR in Electrostatics::Load::alpha(const char*,MoleculeT&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("Failed to load.");
}

template <class AtomT>
void Electrostatics::Load::alpha(const char* file, SimAtomic<AtomT>& sim){
	if(DEBUG_ELEC_UTIL>0) std::cout<<"Electrostatics::Load::alpha(const char*,SimAtomic<AtomT>&):\n";
	//local variables
	//typedefs
		typedef Atom<Name,Species,Index> Particle;
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
		char* atomString=(char*)malloc(sizeof(char)*string::M);
		FILE* reader=NULL;
	//atoms
		std::vector<int> atomIndices;
		std::vector<std::string> atomNames;
		std::vector<Particle> atoms;
	//misc
		bool error=false;
	
	try{
		//open the file
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O ERROR: Could not open file.");
		
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_all(string::trim_right(input,string::COMMENT));
			string::to_upper(string::copy_left(temp,input,"="));
			if(std::strcmp(temp,"ALPHA_ATOM")==0){
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strpbrk(input,"}")!=NULL) break;//break if find "}"
					//load the atom string
					string::copy_left(atomString,string::trim_all(input),"=");
					//load the atom list
					int nAtoms=atomList::loadNumAtoms(atomString);
					atomList::loadAtomIndices(atomString,atomIndices);
					atomList::loadAtomNames(atomString,atomNames);
					atoms.resize(nAtoms);
					for(unsigned int i=0; i<nAtoms; ++i){
						atoms[i].init();
						atoms[i].name()=atomNames[i];
						atoms[i].specie()=sim.speciesIndex(atomNames[i]);
						atoms[i].index()=atomIndices[i];
					}
					//load the polarizability
					double alpha=std::atof(std::strpbrk(input,"=")+1);
					//set the polarizabilities
					for(unsigned int t=0; t<sim.timesteps(); ++t){
						for(unsigned int i=0; i<nAtoms; ++i){
							sim.atom(t,atoms[i].specie(),atoms[i].index()).alpha()=Eigen::Matrix3d::Identity()*alpha;
						}
					}
				}
			}
		}
		fclose(reader);
		reader==NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Electrostatics::Load::alpha(const char*,MoleculeT&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("Failed to load.");
}

template <class AtomT>
void Electrostatics::Load::charge(const char* file, SimAtomic<AtomT>& sim){
	if(DEBUG_ELEC_UTIL>0) std::cout<<"loadChargeAtom(SimAtomic<AtomT>&,const char*):\n";
	//typedefs
		typedef Atom<Name,Species,Index> Particle;
	//local function variables
	FILE* reader=NULL;
	char* input=(char*)malloc(sizeof(char)*string::M);
	char* temp=(char*)malloc(sizeof(char)*string::M);
	char* atomString=(char*)malloc(sizeof(char)*string::M);
	int nAtoms;
	std::vector<int> atomIndices;
	std::vector<std::string> atomNames;
	std::vector<Particle> atoms;
	double charge;
	bool error=false;
	
	try{
		//open the parameter file
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open input file.");
		
		//load the charges from the parameter file
		if(DEBUG_ELEC_UTIL>0) std::cout<<"Loading charges from parameter file...\n";
		while(fgets(input, string::M, reader)!=NULL){
			string::trim_right(input,string::COMMENT);
			string::trim_all(input);
			string::copy_left(temp,input,"=");
			string::to_upper(temp);
			if(std::strcmp(temp,"CHARGES")==0){
				while(fgets(input, string::M, reader)!=NULL){
					string::trim_right(input,string::COMMENT);
					//break on "}"
					if(std::strpbrk(input,"}")!=NULL) break;
					string::trim_all(input);
					string::copy_left(atomString,input,"=");
					//load the atom list
					nAtoms=atomList::loadNumAtoms(atomString);
					atomList::loadAtomIndices(atomString,atomIndices);
					atomList::loadAtomNames(atomString,atomNames);
					atoms.resize(nAtoms);
					for(unsigned int i=0; i<nAtoms; ++i){
						atoms[i].init();
						atoms[i].name()=atomNames[i];
						atoms[i].specie()=sim.speciesIndex(atomNames[i]);
						atoms[i].index()=atomIndices[i];
						if(atoms[i].specie()<0 || atoms[i].specie()>sim.nSpecies()-1) throw std::runtime_error("Invalid atom.");
						else if(atoms[i].index()<0 || atoms[i].index()>sim.nAtoms(atoms[i].specie())-1) throw std::runtime_error("Invalid atom.");
					}
					//load the charge
					charge=std::atof(std::strpbrk(input,"=")+1);
					//print to screen
					if(DEBUG_ELEC_UTIL>0){
						int atomsPerLine=10;
						std::cout<<"Atoms:\n";
						for(int i=0; i<atoms.size()-1; i++){
							std::cout<<atoms[i].name()<<atoms[i].index()+1<<",";
							if((i+1)%atomsPerLine==0) std::cout<<"\n";
						}
						std::cout<<atoms.back().name()<<atoms.back().index()+1<<"\n";
						std::cout<<"Charge = "<<charge<<"\n";
					}
					//set the charges
					for(unsigned int t=0; t<sim.timesteps(); ++t){
						for(unsigned int i=0; i<nAtoms; ++i){
							sim.atom(t,atoms[i].specie(),atoms[i].index()).charge()=charge;
						}
					}
				}
			}
		}
	}catch(std::exception& e){
		std::cout<<"ERROR in Electrostatics::Load::Charge(const char* file, SimAtomic<AtomT>&):\n";
		std::cout<<e.what()<<"\n";
	}
	
	//free all local variables
	free(input);
	free(temp);
	free(atomString);
	if(reader!=NULL) fclose(reader);
	
	if(error) throw std::runtime_error("Data was not assigned to variables.\n");
}

//Print

template <class AtomT>
void Electrostatics::Print::chg_tot(const char* file, const SimAtomic<AtomT>& sim){
	if(DEBUG_ELEC_UTIL>0) std::cout<<"Electrostatics::Print::chg_tot(const char*,SimAtomic<AtomT>&):\n";
	//local variables
	//file i/o
		FILE* writer=NULL;
	//misc
		bool error=false;
	
	try{
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("I/O ERROR: Could not open file.");
		
		fprintf(writer, "T CHG_TOT\n");
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			double chgTot=0;
			for(unsigned int n=0; n<sim.nAtoms(); ++n) chgTot+=sim.atom(t,n).charge();
			fprintf(writer,"%i %f\n", t, chgTot);
		}
		
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Electrostatics::Print::chgTot(const char*,SimAtomic<AtomT>&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("Failed to print.");
}

template <class AtomT>
void Electrostatics::Print::chg(const char* file, const SimAtomic<AtomT>& sim){
	if(DEBUG_ELEC_UTIL>0) std::cout<<"Electrostatics::Print::chg(const char*,SimAtomic<AtomT>&):\n";
	//local variables
	//file i/o
		FILE* writer=NULL;
	//misc
		bool error=false;
		
	try{
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("I/O ERROR: Could not open file.");
		
		fprintf(writer, "T ");
		for(unsigned int n=0; n<sim.nAtoms(); ++n) fprintf(writer, "%s%i ", sim.atom(0,n).name().c_str(), sim.atom(0,n).index()+1);
		fprintf(writer, "\n");
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			fprintf(writer, "%i ", t);
			for(unsigned int n=0; n<sim.nAtoms(); ++n) fprintf(writer, "%f ", sim.atom(t,n).charge());
			fprintf(writer, "\n");
		}
		
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Electrostatics::Print::chgTot(const char*,SimAtomic<AtomT>&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("Failed to print.");
}

template <class AtomT>
void Electrostatics::Print::alpha_tot(const char* file, const SimAtomic<AtomT>& sim){
	if(DEBUG_ELEC_UTIL>0) std::cout<<"Electrostatics::Print::alpha_tot(const char*,SimAtomic<AtomT>&):\n";
	//local variables
	//file i/o
		FILE* writer=NULL;
	//misc
		bool error=false;
	
	try{
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("I/O ERROR: Could not open file.");
		
		fprintf(writer, "T XX XY XZ YX YY YZ ZX ZY ZZ\n");
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			Eigen::Matrix3d alphaTot=Eigen::Matrix3d::Zero();
			for(unsigned int n=0; n<sim.nAtoms(); ++n) alphaTot.noalias()+=sim.atom(t,n).alpha();
			fprintf(writer, "%i ", t);
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					fprintf(writer, "%f ", alphaTot(i,j));
				}
			}
			fprintf(writer,"\n");
		}
		
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Electrostatics::Print::chgTot(const char*,SimAtomic<AtomT>&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("Failed to print.");
}

template <class AtomT>
void Electrostatics::Print::dipole_tot(const char* file, const SimAtomic<AtomT>& sim, const Bonding& bonding){
	if(DEBUG_ELEC_UTIL>0) std::cout<<"Electrostatics::Print::dipole_tot(const char*,const SimAtomic<AtomT>&,const Bonding&):\n";
	//local variables
	//typedefs
		typedef Molecule<AtomT,Position,Charge,Dipole> MoleculeT;
	//file i/o
		FILE* writer=NULL;
	//bonds
		SimMol<AtomT,MoleculeT> simMol;
	//misc
		bool error=false;
	
	try{
		//open the output file
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("I/O ERROR: Could not open file.");
		
		//print the header
		fprintf(writer, "T X Y Z\n");
		
		//initialize the molecular simulation
		static_cast<SimAtomic<AtomT>&>(simMol).resize(1,sim.nAtomsVec(),sim.atomNames());
		//bond order
		std::vector<std::vector<unsigned int> > bOrder(sim.nSpecies());
		for(unsigned int n=0; n<sim.nSpecies(); ++n) bOrder[n].resize(sim.nAtoms(n));
		
		//calculate the total dipole moment
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			//initialize the cell
			simMol.cell(0).init(sim.cell(t).R());
			//local variables
			Eigen::Vector3d muTot=Eigen::Vector3d::Zero();
			//copy over the atoms
			if(DEBUG_ELEC_UTIL>1) std::cout<<"Copying over atoms...\n";
			for(unsigned int n=0; n<sim.nSpecies(); ++n){
				for(unsigned int m=0; m<sim.nAtoms(n); ++m){
					simMol.atom(0,n,m).posn().noalias()=sim.atom(t,n,m).posn();
					simMol.atom(0,n,m).charge()=sim.atom(t,n,m).charge();
				}
			}
			//calculate bonds
			if(DEBUG_ELEC_UTIL>1) std::cout<<"Calculating bonds...\n";
			sim_util::assignBonds(simMol,static_cast<SimAtomic<AtomT>&>(simMol),bonding);
			if(DEBUG_ELEC_UTIL>1) std::cout<<"NMol = "<<simMol.nMol()<<"\n";
			//find the bond order
			if(DEBUG_ELEC_UTIL>1) std::cout<<"Finding the bond order...\n";
			for(unsigned int n=0; n<sim.nSpecies(); ++n){
				for(unsigned int m=0; m<sim.nAtoms(n); ++m){
					bOrder[n][m]=Bonding::bondOrder(t,sim.atom(t,n,m),sim,bonding);
				}
			}
			if(DEBUG_ELEC_UTIL>1) std::cout<<"Divide charge by bond order...\n";
			for(unsigned int n=0; n<simMol.nSpecies(); ++n){
				for(unsigned int m=0; m<simMol.nAtoms(n); ++m){
					simMol.atom(0,n,m).charge()/=((double)bOrder[n][m]);
				}
			}
			if(DEBUG_ELEC_UTIL>1) std::cout<<"Find total molecular charge...\n";
			for(unsigned int n=0; n<simMol.nMol(); ++n){
				simMol.molecule(0,n).charge()=0;
				for(unsigned int m=0; m<simMol.molecule(0,n).nAtoms(); ++m){
					simMol.molecule(0,n).charge()+=simMol.molecule(0,n).atom(m).charge();
				}
			}
			if(DEBUG_ELEC_UTIL>1) std::cout<<"Find center of charge...\n";
			for(unsigned int n=0; n<simMol.nMol(); ++n){
				coc(simMol.molecule(0,n),simMol.cell(t),simMol.molecule(0,n).posn());
			}
			if(DEBUG_ELEC_UTIL>1) std::cout<<"Find molecular dipole moment...\n";
			for(unsigned int n=0; n<simMol.nMol(); ++n){
				dipole(simMol.molecule(0,n),simMol.cell(t),simMol.molecule(0,n).dipole());
			}
			for(unsigned int n=0; n<simMol.nMol(); ++n){
				muTot.noalias()+=simMol.molecule(0,n).dipole();
			}
			//print the total dipole moment to file
			fprintf(writer, "%i %f %f %f\n", t, muTot[0], muTot[1], muTot[2]);
		}
		
		//close the writer
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Electrostatics::Print::dipoleTot(const char*,const SimAtomic<Atomt>&,const Bonding&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
}

template <class MoleculeT, class AtomT>
void Electrostatics::Print::dipole_tot(const char* file, const SimMol<AtomT,MoleculeT>& sim){
	if(DEBUG_ELEC_UTIL>0) std::cout<<"Electrostatics::Print::dipole_tot(const char*,const SimMol<AtomT,MoleculeT>&):\n";
	//file i/o
		FILE* writer=NULL;
	//misc
		bool error=false;
	
	try{
		//open the output file
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("I/O ERROR: Could not open file.");
		
		//print the header
		fprintf(writer, "T X Y Z\n");
		
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			Eigen::Vector3d dipoleT=Eigen::Vector3d::Zero();
			for(unsigned int n=0; n<sim.nMol(); ++n){
				dipoleT.noalias()+=sim.molecule(t,n).dipole();
			}
			fprintf(writer,"%i %f %f %f\n",t,dipoleT[0],dipoleT[1],dipoleT[2]);
		}
		
		//close the file
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Electrostatics::Print::dipole_tot(const char*,const SimMol<AtomT,MoleculeT>&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
}

#endif