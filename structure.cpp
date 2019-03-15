#include "structure.hpp"

//**********************************************************************************************
//AtomType
//**********************************************************************************************

std::ostream& operator<<(std::ostream& out, const AtomType& atomT){
	if(atomT.name)		out<<"name ";
	if(atomT.an)		out<<"an ";
	if(atomT.specie)	out<<"specie ";
	if(atomT.index)		out<<"index ";
	if(atomT.mass)		out<<"mass ";
	if(atomT.charge)	out<<"charge ";
	if(atomT.posn)		out<<"posn ";
	if(atomT.velocity)	out<<"velocity ";
	if(atomT.force)		out<<"force ";
	if(atomT.dipole)	out<<"dipole ";
	if(atomT.alpha)		out<<"alpha ";
	if(atomT.jzero)		out<<"jzero ";
	if(atomT.symm)		out<<"symm ";
	if(atomT.neighlist)	out<<"neighlist ";
	if(atomT.molecule)	out<<"molecule ";
	if(atomT.frac)		out<<"frac ";
	return out;
}

void AtomType::defaults(){
	name	=false;
	an		=false;
	specie	=false;
	index	=false;
	mass	=false;
	charge	=false;
	posn	=false;
	velocity=false;
	force	=false;
	dipole	=false;
	alpha	=false;
	jzero	=false;
	symm	=false;
	neighlist	=false;
	molecule=false;
	frac	=false;
}

unsigned int AtomType::nbytes()const{
	unsigned int nbytes=0;
	nbytes+=sizeof(char)*2;
	nbytes+=sizeof(unsigned int)*an;
	nbytes+=sizeof(unsigned int)*specie;
	nbytes+=sizeof(unsigned int)*index;
	nbytes+=sizeof(double)*mass;
	nbytes+=sizeof(double)*charge;
	nbytes+=sizeof(double)*3*posn;
	nbytes+=sizeof(double)*3*velocity;
	nbytes+=sizeof(double)*3*force;
	nbytes+=sizeof(double)*3*dipole;
	nbytes+=sizeof(double)*9*alpha;
	nbytes+=sizeof(double)*jzero;
	nbytes+=sizeof(double)*16*symm;//estimate
	nbytes+=sizeof(unsigned int)*8*neighlist;//estimate
	nbytes+=sizeof(unsigned int)*molecule;
}

//**********************************************************************************************
//Interval
//**********************************************************************************************

std::ostream& operator<<(std::ostream& out, const Interval& i){
	return out<<i.beg<<":"<<i.end<<":"<<i.stride;
}

Interval Interval::read(const char* str){
	std::vector<std::string> strlist;
	string::split(str,":",strlist);
	if(strlist.size()!=2 && strlist.size()!=3) throw std::invalid_argument("Invalid interval format.");
	Interval interval;
	if(string::to_upper(strlist.at(0))=="BEG") interval.beg=1;
	else interval.beg=std::atoi(strlist.at(0).c_str());
	if(string::to_upper(strlist.at(1))=="END") interval.end=-1;
	else interval.end=std::atoi(strlist.at(1).c_str());
	if(strlist.size()==3){
		interval.stride=std::atoi(strlist.at(2).c_str());
	}
	return interval;
}

//**********************************************************************************************
//FILE_FORMAT struct
//**********************************************************************************************

FILE_FORMAT::type FILE_FORMAT::read(const std::string& str){
	if(str=="XDATCAR") return FILE_FORMAT::XDATCAR;
	else if(str=="POSCAR") return FILE_FORMAT::POSCAR;
	else if(str=="OUTCAR") return FILE_FORMAT::OUTCAR;
	else if(str=="VASP_XML") return FILE_FORMAT::VASP_XML;
	else if(str=="GAUSSIAN") return FILE_FORMAT::GAUSSIAN;
	else if(str=="DFTB") return FILE_FORMAT::DFTB;
	else if(str=="XYZ") return FILE_FORMAT::XYZ;
	else if(str=="CAR") return FILE_FORMAT::CAR;
	else if(str=="LAMMPS") return FILE_FORMAT::LAMMPS;
	else if(str=="GROMACS") return FILE_FORMAT::GROMACS;
	else if(str=="QE") return FILE_FORMAT::QE;
	else if(str=="PROPHET") return FILE_FORMAT::PROPHET;
	else if(str=="XSF") return FILE_FORMAT::XSF;
	else if(str=="AME") return FILE_FORMAT::AME;
	else return FILE_FORMAT::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, FILE_FORMAT::type& format){
	if(format==FILE_FORMAT::UNKNOWN) out<<"UNKNOWN";
	else if(format==FILE_FORMAT::XDATCAR) out<<"XDATCAR";
	else if(format==FILE_FORMAT::POSCAR) out<<"POSCAR";
	else if(format==FILE_FORMAT::OUTCAR) out<<"OUTCAR";
	else if(format==FILE_FORMAT::VASP_XML) out<<"VASP_XML";
	else if(format==FILE_FORMAT::GAUSSIAN) out<<"GAUSSIAN";
	else if(format==FILE_FORMAT::DFTB) out<<"DFTB";
	else if(format==FILE_FORMAT::XYZ) out<<"XYZ";
	else if(format==FILE_FORMAT::CAR) out<<"CAR";
	else if(format==FILE_FORMAT::LAMMPS) out<<"LAMMPS";
	else if(format==FILE_FORMAT::GROMACS) out<<"GROMACS";
	else if(format==FILE_FORMAT::QE) out<<"QE";
	else if(format==FILE_FORMAT::PROPHET) out<<"PROPHET";
	else if(format==FILE_FORMAT::XSF) out<<"XSF";
	else if(format==FILE_FORMAT::AME) out<<"AME";
	return out;
}

//**********************************************************************************************
//Structure Interface
//**********************************************************************************************

//constructors/destructors

StructureI::StructureI(const StructureI& sim){
	//cells
	cell_=sim.cell();
	//simulation info
	energy_=sim.energy();
	//atoms
	nSpecies_=sim.nSpecies();
	nAtomsT_=sim.nAtoms();
	nAtoms_.resize(nSpecies_);
	offsets_.resize(nSpecies_);
	for(unsigned int i=0; i<nSpecies_; ++i){
		nAtoms_[i]=sim.nAtoms(i);
		offsets_[i]=sim.offset(i);
	}
	atomNames_=sim.atomNames();
}

//operators

StructureI& StructureI::operator=(const StructureI& sim){
	//cells
	cell_=sim.cell();
	//simulation info
	energy_=sim.energy();
	//atoms
	nSpecies_=sim.nSpecies();
	nAtomsT_=sim.nAtoms();
	nAtoms_.resize(nSpecies_);
	offsets_.resize(nSpecies_);
	for(unsigned int i=0; i<nSpecies_; ++i){
		nAtoms_[i]=sim.nAtoms(i);
		offsets_[i]=sim.offset(i);
	}
	atomNames_=sim.atomNames();
	return *this;
}

std::ostream& operator<<(std::ostream& out, const StructureI& sim){
	out<<sim.cell()<<"\n";
	for(unsigned int i=0; i<sim.atomNames_.size(); ++i)
		out<<sim.atomNames_[i]<<" ";
	out<<"\n";
	for(unsigned int i=0; i<sim.nAtoms_.size(); ++i)
		out<<sim.nAtoms_[i]<<" ";
	return out;
}

//member functions

void StructureI::defaults(){
	if(DEBUG_STRUCTURE>0) std::cout<<"StructureI::defaults():\n";
	//cell
	cell_.clear();
	//simulation data
	energy_=0;
	//atoms
	nSpecies_=0;
	nAtomsT_=0;
	nAtoms_.clear();
	offsets_.clear();
	atomNames_.clear();
}

void StructureI::resize(const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& atomNames){
	if(DEBUG_STRUCTURE>0) std::cout<<"StructureI::resize(const std::vector<unsigned int>&,const std::vector<std::string>&):\n";
	if(nAtoms.size()!=atomNames.size()) throw std::invalid_argument("Array size mismatch.");
	nSpecies_=nAtoms.size();
	nAtoms_=nAtoms;
	atomNames_=atomNames;
	nAtomsT_=0;
	offsets_.resize(nSpecies_,0);
	for(unsigned int i=0; i<nSpecies_; ++i) nAtomsT_+=nAtoms_[i];
	for(unsigned int i=1; i<nSpecies_; ++i) offsets_[i]=offsets_[i-1]+nAtoms_[i-1];
}

//static functions

int StructureI::speciesIndex(const std::string& str, const std::vector<std::string>& atomNames){
	if(DEBUG_STRUCTURE>0) std::cout<<"StructureI::speciesIndex(const std::string&, const std::vector<std::string>&):\n";
	for(unsigned int i=0; i<atomNames.size(); ++i) if(str==atomNames[i]) return i;
	return -1;
}

int StructureI::speciesIndex(const char* str, const std::vector<std::string>& atomNames){
	if(DEBUG_STRUCTURE>0) std::cout<<"StructureI::speciesIndex(const char*,const std::vector<std::string>&):\n";
	for(unsigned int i=0; i<atomNames.size(); ++i) if(std::strcmp(str,atomNames[i].c_str())==0) return i;
	return -1;
}

std::vector<unsigned int>& StructureI::read_atoms(const char* str, std::vector<unsigned int>& ids, const StructureI& strucI){
	//local function variables
	unsigned int nAtoms=0;
	std::vector<int> atomIndices;
	std::vector<std::string> atomNames;
	std::vector<int> nameIndices;
	
	//load the number of atoms, atom indices, and atom names
	nAtoms=StructureI::read_natoms(str);
	StructureI::read_indices(str,atomIndices);
	StructureI::read_names(str,atomNames);
	nameIndices.resize(atomNames.size());
	for(unsigned int n=0; n<nameIndices.size(); ++n){
		nameIndices[n]=strucI.speciesIndex(atomNames[n]);
		if(nameIndices[n]<0) throw std::invalid_argument("ERROR: Invalid atom name in atom string.");
	}
	
	//set the atom ids
	ids.resize(nAtoms,0);
	for(unsigned int n=0; n<nAtoms; ++n) ids[n]=strucI.id(nameIndices[n],atomIndices[n]);
	
	return ids;
}

unsigned int StructureI::read_natoms(const char* str){
	const char* func_name="StructureI::read_natoms(const char*)";
	if(DEBUG_STRUCTURE>0) std::cout<<func_name<<":\n";
	//local function variables
	char* strtemp=new char[string::M];
	char* substr=new char[string::M];
	char* temp=new char[string::M];
	std::vector<std::string> substrs;
	unsigned int nStrs=0,nAtoms=0;
	bool error=false;
	
	try{
		//copy the string
		std::strcpy(strtemp,str);
		//find the number of substrings
		nStrs=string::substrN(strtemp,",");
		substrs.resize(nStrs);
		substrs[0]=std::string(std::strcpy(temp,std::strtok(strtemp,",")));
		for(unsigned int i=1; i<nStrs; ++i){
			substrs[i]=std::string(std::strcpy(temp,std::strtok(NULL,",")));
		}
		
		//parse the line by commas
		for(unsigned int i=0; i<nStrs; ++i){
			std::strcpy(substr,substrs[i].c_str());
			//find out if this substring has only one atom, or a set of atoms
			if(std::strpbrk(substr,":")==NULL) ++nAtoms; //single atom
			else {
				//set of atoms
				int beg, end;
				//find the beginning index
				std::strcpy(temp,std::strtok(substr,":"));
				if(std::strpbrk(temp,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
				else beg=std::atoi(std::strpbrk(temp,string::DIGITS))-1;
				//find the ending index
				std::strcpy(temp,std::strtok(NULL,":"));
				if(std::strpbrk(temp,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
				else end=std::atoi(std::strpbrk(temp,string::DIGITS))-1;
				//check the indices
				if(beg<0 || end<0 || end<beg) throw std::invalid_argument("Invalid atomic indices.");
				else nAtoms+=end-beg+1;
			}
		}
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	delete[] strtemp;
	delete[] substr;
	delete[] temp;
	
	if(error) throw std::invalid_argument("Invalid atom string.");
	else return nAtoms;
}
	
std::vector<int>& StructureI::read_indices(const char* str, std::vector<int>& indices){
	const char* func_name="StructureI::read_indices(const char*,std::vector<int>&)";
	if(DEBUG_STRUCTURE>0) std::cout<<func_name<<":\n";
	//local function variables
	char* strtemp=new char[string::M];
	char* substr=new char[string::M];
	char* temp=new char[string::M];
	std::vector<std::string> substrs;
	unsigned int nStrs=0;
	bool error=false;
	
	try{
		//clear the vector
		indices.clear();
		//copy the string
		std::strcpy(strtemp,str);
		//find the number of substrings
		nStrs=string::substrN(strtemp,",");
		substrs.resize(nStrs);
		substrs[0]=std::string(std::strcpy(temp,std::strtok(strtemp,",")));
		for(unsigned int i=1; i<nStrs; ++i){
			substrs[i]=std::string(std::strcpy(temp,std::strtok(NULL,",")));
		}
		
		//parse the line by commas
		for(unsigned int i=0; i<nStrs; ++i){
			std::strcpy(substr,substrs[i].c_str());
			//find out if this substring has only one atom, or a set of atoms
			if(std::strpbrk(substr,":")==NULL){
				//single atom
				if(DEBUG_STRUCTURE>1) std::cout<<"Single Atom\n";
				if(std::strpbrk(substr,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
				else indices.push_back(std::atoi(std::strpbrk(substr,string::DIGITS))-1);
			} else {
				//set of atoms
				if(DEBUG_STRUCTURE>1) std::cout<<"Set of Atoms\n";
				int beg, end;
				//find the beginning index
				std::strcpy(temp,std::strtok(substr,":"));
				if(std::strpbrk(temp,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
				else beg=std::atoi(std::strpbrk(temp,string::DIGITS))-1;
				//find the ending index
				std::strcpy(temp,std::strtok(NULL,":"));
				if(std::strpbrk(temp,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
				else end=std::atoi(std::strpbrk(temp,string::DIGITS))-1;
				//check the indices
				if(beg<0 || end<0 || end<beg) throw std::invalid_argument("Invalid atomic indices.");
				for(int j=0; j<end-beg+1; j++) indices.push_back(beg+j);
			}
		}
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	delete[] strtemp;
	delete[] substr;
	delete[] temp;
	
	if(error){
		indices.clear();
		throw std::invalid_argument("Invalid atom string.");
	} else return indices;
}

std::vector<std::string>& StructureI::read_names(const char* str, std::vector<std::string>& names){
	const char* func_name="StructureI::read_names(const char*,std::vector<std::string>&)";
	if(DEBUG_STRUCTURE>0) std::cout<<func_name<<":\n";
	//local function variables
	char* strtemp=new char[string::M];
	char* substr=new char[string::M];
	char* temp=new char[string::M];
	char* atomName=new char[string::M];
	std::vector<std::string> substrs;
	unsigned int nStrs=0;
	bool error=false;
	
	try{
		//clear the vector
		names.clear();
		//copy the string
		std::strcpy(strtemp,str);
		//find the number of substrings
		nStrs=string::substrN(strtemp,",");
		substrs.resize(nStrs);
		substrs[0]=std::string(std::strcpy(temp,std::strtok(strtemp,",")));
		for(unsigned int i=1; i<nStrs; ++i){
			substrs[i]=std::string(std::strcpy(temp,std::strtok(NULL,",")));
		}
		
		//parse the line by commas
		for(unsigned int i=0; i<nStrs; ++i){
			std::strcpy(substr,substrs[i].c_str());
			//find out if this substring has only one atom, or a set of atoms
			if(std::strpbrk(substr,":")==NULL){
				//single atom
				if(std::strpbrk(substr,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
				else names.push_back(std::string(string::trim_right(substr,string::DIGITS)));
			} else {
				//set of atoms
				int beg, end;
				//find the beginning index
				std::strcpy(temp,std::strtok(substr,":"));
				if(std::strpbrk(temp,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
				else beg=std::atoi(std::strpbrk(temp,string::DIGITS))-1;
				//find the ending index
				std::strcpy(temp,std::strtok(NULL,":"));
				if(std::strpbrk(temp,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
				else end=std::atoi(std::strpbrk(temp,string::DIGITS))-1;
				//check the indices
				if(beg<0 || end<0 || end<beg) throw std::invalid_argument("Invalid atomic indices");
				else {
					std::strcpy(atomName,string::trim_right(temp,string::DIGITS));
					for(int j=0; j<end-beg+1; ++j) names.push_back(std::string(atomName));
				}
			}
		}
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	delete[] strtemp;
	delete[] substr;
	delete[] temp;
	delete[] atomName;
	
	if(error){
		names.clear();
		throw std::invalid_argument("Invalid atom string.");
	} else return names;
}

//**********************************************************************************************
//Structure Storage
//**********************************************************************************************

//constructors/destructors

StructureS::StructureS(const StructureS& struc){
	if(DEBUG_STRUCTURE>0) std::cout<<"StructureS::StructureS(const StructureS&):\n";
	name_		=struc.name();
	an_			=struc.an();
	specie_		=struc.specie();
	index_		=struc.index();
	mass_		=struc.mass();
	charge_		=struc.charge();
	posn_		=struc.posn();
	velocity_	=struc.velocity();
	force_		=struc.force();
	dipole_		=struc.dipole();
	alpha_		=struc.alpha();
	jzero_		=struc.jzero();
	symm_		=struc.symm();
	neighlist_	=struc.neighlist();
	molecule_	=struc.molecule();
};

//operators

StructureS& StructureS::operator=(const StructureS& struc){
	if(DEBUG_STRUCTURE>0) std::cout<<"StructureS::operator=(const StructureS&):\n";
	name_		=struc.name();
	an_			=struc.an();
	specie_		=struc.specie();
	index_		=struc.index();
	mass_		=struc.mass();
	charge_		=struc.charge();
	posn_		=struc.posn();
	velocity_	=struc.velocity();
	force_		=struc.force();
	dipole_		=struc.dipole();
	alpha_		=struc.alpha();
	jzero_		=struc.jzero();
	symm_		=struc.symm();
	neighlist_	=struc.neighlist();
	molecule_	=struc.molecule();
	return *this;
}

//member functions

void StructureS::clear(){
	if(DEBUG_STRUCTURE>0) std::cout<<"StructureS::clear():\n";
	name_.clear();
	an_.clear();
	specie_.clear();
	index_.clear();
	mass_.clear();
	charge_.clear();
	posn_.clear();
	velocity_.clear();
	force_.clear();
	dipole_.clear();
	alpha_.clear();
	jzero_.clear();
	symm_.clear();
	neighlist_.clear();
	molecule_.clear();
}

AtomType StructureS::atomType(){
	AtomType atomT;
	if(name_.size()>0)		atomT.name=true;
	if(an_.size()>0)		atomT.an=true;
	if(specie_.size()>0)	atomT.specie=true;
	if(index_.size()>0)		atomT.index=true;
	if(mass_.size()>0)		atomT.mass=true;
	if(charge_.size()>0)	atomT.charge=true;
	if(posn_.size()>0)		atomT.posn=true;
	if(velocity_.size()>0)	atomT.velocity=true;
	if(force_.size()>0)		atomT.force=true;
	if(dipole_.size()>0)	atomT.dipole=true;
	if(alpha_.size()>0)		atomT.alpha=true;
	if(jzero_.size()>0)		atomT.jzero=true;
	if(symm_.size()>0)		atomT.symm=true;
	if(neighlist_.size()>0) atomT.neighlist=true;
	if(molecule_.size()>0)	atomT.molecule=true;
	return atomT;
}

//resizing

void StructureS::resize(const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& speciesNames, const AtomType& atomT){
	if(DEBUG_STRUCTURE>0) std::cout<<"StructureS::resize(const std::vector<unsigned int>&,const std::vector<std::string>&,const AtomType&):\n";
	//find the total number of atoms
	unsigned int nAtomsT=0;
	for(unsigned int i=0; i<nAtoms.size(); ++i) nAtomsT+=nAtoms[i];
	//resize the property arrays
	if(atomT.name)		name_.resize(nAtomsT);
	if(atomT.an)		an_.resize(nAtomsT);
	if(atomT.specie)	specie_.resize(nAtomsT);
	if(atomT.index)		index_.resize(nAtomsT);
	if(atomT.mass)		mass_.resize(nAtomsT);
	if(atomT.charge)	charge_.resize(nAtomsT);
	if(atomT.posn)		posn_.resize(nAtomsT);
	if(atomT.velocity)	velocity_.resize(nAtomsT);
	if(atomT.force)		force_.resize(nAtomsT);
	if(atomT.dipole)	dipole_.resize(nAtomsT);
	if(atomT.alpha)		alpha_.resize(nAtomsT);
	if(atomT.jzero)		jzero_.resize(nAtomsT);
	if(atomT.symm)		symm_.resize(nAtomsT);
	if(atomT.neighlist)	neighlist_.resize(nAtomsT);
	if(atomT.molecule)	molecule_.resize(nAtomsT);
}

//**********************************************************************************************
//Structure Atomic
//**********************************************************************************************

//operators

Structure& Structure::operator=(const Structure& sim){
	if(DEBUG_STRUCTURE>0) std::cout<<"Structure::operator=(const Structure&):\n";
	StructureI::operator=(sim);
	StructureS::operator=(sim);
	return *this;
}

std::ostream& operator<<(std::ostream& out, const Structure& sim){
	return out<<static_cast<const StructureI&>(sim);
}

//member functions

void Structure::clear(){
	if(DEBUG_STRUCTURE>0) std::cout<<"Structure::clear():\n";
	StructureI::clear();
	StructureS::clear();
}

void Structure::resize(const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& atomNames, const AtomType& atomT){
	if(DEBUG_STRUCTURE>0) std::cout<<"Structure::resize(const std::vector<unsigned int>&,const std::vector<std::string>&):\n";
	StructureI::resize(nAtoms,atomNames);
	StructureS::resize(nAtoms,atomNames,atomT);
	//set the names, species, and indices
	unsigned int count=0;
	for(unsigned int n=0; n<this->nSpecies_; ++n){
		for(unsigned int m=0; m<this->nAtoms_[n]; ++m){
			if(this->name_.size()>0) this->name_[count]=this->atomNames_[n];
			if(this->an_.size()>0) this->an_[count]=PTable::an(this->atomNames_[n].c_str());
			if(this->specie_.size()>0) this->specie_[count]=n;
			if(this->index_.size()>0) this->index_[count]=m;
			++count;
		}
	}
}

//**********************************************
// Simulation
//**********************************************

//operators

std::ostream& operator<<(std::ostream& out, const Simulation& sim){
	out<<"****************************************\n";
	out<<"************** SIMULATION **************\n";
	out<<"SIM   = "<<sim.name_<<"\n";
	out<<"TS    = "<<sim.timestep_<<"\n";
	out<<"T     = "<<sim.timesteps_<<"\n";
	out<<"CF    = "<<sim.cell_fixed_<<"\n";
	out<<"AtomT = "<<sim.atomT_<<"\n";
	out<<"************** SIMULATION **************\n";
	out<<"****************************************";
	return out;
}

//member functions

void Simulation::defaults(){
	if(DEBUG_STRUCTURE>0) std::cout<<"Simulation::defaults():\n";
	name_=std::string("SYSTEM");
	timestep_=0.5;
	timesteps_=0;
}

void Simulation::clear(){
	if(DEBUG_STRUCTURE>0) std::cout<<"Simulation::clear():\n";
	frames_.clear();
	defaults();
}

void Simulation::resize(unsigned int ts, const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& names, const AtomType& atomT){
	if(DEBUG_STRUCTURE>0) std::cout<<"Simulation::resize(unsigned int,const std::vector<unsigned int>&,const std::vector<std::string>&,const AtomType&):\n";
	timesteps_=ts;
	atomT_=atomT;
	frames_.resize(timesteps_,Structure(nAtoms,names,atomT));
}

void Simulation::resize(unsigned int ts){
	if(DEBUG_STRUCTURE>0) std::cout<<"Simulation::resize(unsigned int):\n";
	timesteps_=ts;
	frames_.resize(timesteps_);
}

//**********************************************
// Utility
//**********************************************

//read

void Utility::Read::chg_atom(const char* file, Simulation& sim){
	if(DEBUG_STRUCTURE>0) std::cout<<"Utility::Read::chg_atom(const char*,Simulation&):\n";
	char* input=new char[string::M];
	FILE* reader=fopen(file,"r");
	if(reader!=NULL){
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			fgets(input,string::M,reader);
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				sim.frame(t).charge(n)=std::atof(std::strtok(NULL,string::WS));
			}
		}
		fclose(reader);
		reader=NULL;
	} else std::cout<<"WARNING: Could not open file: \""<<file<<"\"\n";
	delete[] input;
}

void Utility::Read::alpha_atom(const char* file, Simulation& sim){
	if(DEBUG_STRUCTURE>0) std::cout<<"Utility::Read::alpha_atom(const char*,Simulation&):\n";
	char* input=new char[string::M];
	FILE* reader=fopen(file,"r");
	if(reader!=NULL){
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			fgets(input,string::M,reader);
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				for(unsigned int i=0; i<3; ++i){
					for(unsigned int j=0; j<3; ++j){
						sim.frame(t).alpha(n)(i,j)=std::atof(std::strtok(NULL,string::WS));
					}
				}
			}
		}
		fclose(reader);
		reader=NULL;
	} else std::cout<<"WARNING: Could not open file: \""<<file<<"\"\n";
	delete[] input;
}

//write

void Utility::Write::chg_tot(const char* file, const Simulation& sim){
	if(DEBUG_STRUCTURE>0) std::cout<<"Utility::Write::chg_tot(const char*,const Simulation&):\n";
	FILE* writer=fopen(file,"w");
	if(writer!=NULL){
		fprintf(writer,"#T CHG_TOT\n");
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			double chg=0;
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n) chg+=sim.frame(t).charge(n);
			fprintf(writer,"%i %f\n",t,chg);
		}
		fclose(writer);
		writer=NULL;
	} else std::cout<<"WARNING: Could not open file: \""<<file<<"\"\n";
}

void Utility::Write::chg_atom(const char* file, const Simulation& sim){
	if(DEBUG_STRUCTURE>0) std::cout<<"Utility::Write::chg_tot(const char*,const Simulation&):\n";
	FILE* writer=fopen(file,"w");
	if(writer!=NULL){
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			fprintf(writer,"%i\n",t);
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				fprintf(writer,"%s%i %f\n",sim.frame(t).name(n).c_str(),sim.frame(t).index(n)+1,sim.frame(t).charge(n));
			}
		}
		fclose(writer);
		writer=NULL;
	} else std::cout<<"WARNING: Could not open file: \""<<file<<"\"\n";
}

void Utility::Write::dipole_tot(const char* file, const Simulation& sim){
	if(DEBUG_STRUCTURE>0) std::cout<<"Utility::Write::dipole_tot(const char*,const Simulation&):\n";
	FILE* writer=fopen(file,"w");
	Eigen::Vector3d origin=Eigen::Vector3d::Zero();
	if(writer!=NULL){
		fprintf(writer,"#T DX DY DZ\n");
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			Eigen::Vector3d dipole=Eigen::Vector3d::Zero();
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				Eigen::Vector3d dr;
				Cell::diff(sim.frame(t).posn(n),origin,dr,sim.frame(t).cell().R(),sim.frame(t).cell().RInv());
				dipole.noalias()+=dr*sim.frame(t).charge(n);
			}
			dipole/=sim.frame(t).cell().vol();
			fprintf(writer,"%i ",t);
			for(unsigned int i=0; i<3; ++i){
				fprintf(writer,"%f ",dipole[i]);
			}
			fprintf(writer,"\n");
		}
		fclose(writer);
		writer=NULL;
	} else std::cout<<"WARNING: Could not open file: \""<<file<<"\"\n";
}

void Utility::Write::alpha_tot(const char* file, const Simulation& sim){
	if(DEBUG_STRUCTURE>0) std::cout<<"Utility::Write::alpha_tot(const char*,const Simulation&):\n";
	FILE* writer=fopen(file,"w");
	if(writer!=NULL){
		fprintf(writer,"#T xx xy xz yx yy yz zx zy zz\n");
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			Eigen::Matrix3d alpha_tot=Eigen::Matrix3d::Zero();
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				alpha_tot.noalias()+=sim.frame(t).alpha(n);
			}
			fprintf(writer,"%i ",t);
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					fprintf(writer,"%f ",alpha_tot(i,j));
				}
			}
			fprintf(writer,"\n");
		}
		fclose(writer);
		writer=NULL;
	} else std::cout<<"WARNING: Could not open file: \""<<file<<"\"\n";
}

void Utility::Write::alpha_atom(const char* file, const Simulation& sim){
	if(DEBUG_STRUCTURE>0) std::cout<<"Utility::Write::alpha_tot(const char*,const Simulation&):\n";
	FILE* writer=fopen(file,"w");
	if(writer!=NULL){
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			fprintf(writer,"%i\n",t);
			for(unsigned int n=0; n<sim.frame(t).nAtoms(); ++n){
				fprintf(writer,"%s%i ",sim.frame(t).name(n).c_str(),sim.frame(t).index(n)+1);
				for(unsigned int i=0; i<3; ++i){
					for(unsigned int j=0; j<3; ++j){
						fprintf(writer,"%f ",sim.frame(t).alpha(n)(i,j));
					}
				}
				fprintf(writer,"\n");
			}
		}
		fclose(writer);
		writer=NULL;
	} else std::cout<<"WARNING: Could not open file: \""<<file<<"\"\n";
}

//**********************************************
// Serialization
//**********************************************

namespace serialize{
	
	//**********************************************
	// byte measures
	//**********************************************
	
	template <> unsigned int nbytes(const ListData<unsigned int>& obj){return obj.N()?0:obj.N()*sizeof(unsigned int);};
	template <> unsigned int nbytes(const ListData<double>& obj){return obj.N()?0:obj.N()*sizeof(double);};
	template <> unsigned int nbytes(const ListData<std::string>& obj){unsigned int size=0;for(int i=0; i<obj.N(); ++i)size+=nbytes(obj[i]);return size;};
	template <> unsigned int nbytes(const ListData<Eigen::Vector3d>& obj){unsigned int size=0;for(int i=0; i<obj.N(); ++i)size+=nbytes(obj[i]);return size;};
	template <> unsigned int nbytes(const ListData<Eigen::VectorXd>& obj){unsigned int size=0;for(int i=0; i<obj.N(); ++i)size+=nbytes(obj[i]);return size;};
	
	//**********************************************
	// packing
	//**********************************************
	
	template <> void pack(const ListData<unsigned int>& obj, char* arr){unsigned int n=obj.N()?0:nbytes(obj[0]);for(int i=0; i<obj.N(); ++i)pack(obj[i],arr+i*n);};
	template <> void pack(const ListData<double>& obj, char* arr){unsigned int n=obj.N()?0:nbytes(obj[0]);for(int i=0; i<obj.N(); ++i)pack(obj[i],arr+i*n);};
	template <> void pack(const ListData<std::string>& obj, char* arr){int pos=0;if(obj.N()>0){for(int i=0;i<obj.N();++i){pack(obj[i],arr+pos);pos+=nbytes(obj[i]);}}}
	template <> void pack(const ListData<Eigen::Vector3d>& obj, char* arr){int pos=0;if(obj.N()>0){for(int i=0;i<obj.N();++i){pack(obj[i],arr+pos);pos+=nbytes(obj[i]);}}}
	template <> void pack(const ListData<Eigen::VectorXd>& obj, char* arr){int pos=0;if(obj.N()>0){for(int i=0;i<obj.N();++i){pack(obj[i],arr+pos);pos+=nbytes(obj[i]);}}}
	
	//**********************************************
	// unpacking
	//**********************************************
	
	template <> void unpack(ListData<unsigned int>& obj, const char* arr){unsigned int n=obj.N()?0:nbytes(obj[0]);for(unsigned int i=0; i<obj.N(); ++i)unpack(obj[i],arr+i*n);};
	template <> void unpack(ListData<double>& obj, const char* arr){unsigned int n=obj.N()?0:nbytes(obj[0]);for(unsigned int i=0; i<obj.N(); ++i)unpack(obj[i],arr+i*n);};
	template <> void unpack(ListData<std::string>& obj, const char* arr){int pos=0;if(obj.N()>0){for(int i=0;i<obj.N();++i){unpack(obj[i],arr+pos);pos+=nbytes(obj[i]);}}}
	template <> void unpack(ListData<Eigen::Vector3d>& obj, const char* arr){int pos=0;if(obj.N()>0){for(int i=0;i<obj.N();++i){unpack(obj[i],arr+pos);pos+=nbytes(obj[i]);}}}
	template <> void unpack(ListData<Eigen::VectorXd>& obj, const char* arr){int pos=0;if(obj.N()>0){for(int i=0;i<obj.N();++i){unpack(obj[i],arr+pos);pos+=nbytes(obj[i]);}}}
	
}