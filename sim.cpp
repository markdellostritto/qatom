#include "sim.hpp"

//**********************************************************************************************
//time_interval struct
//**********************************************************************************************

time_interval time_interval::load(const char* str){
	time_interval t;
	if(std::strcmp(str,":")==NULL) throw std::invalid_argument("Invalid time interval format: no \":\".");
	if(string::substrN(str,":")==2){
		char* temp=(char*)malloc(sizeof(char)*250);
		std::strcpy(temp,str);
		t.beg=std::atoi(std::strtok(temp,":"));
		t.end=std::atoi(std::strtok(NULL,":"));
		t.stride=1;
		free(temp); temp=NULL;
	} else if(string::substrN(str,":")==3){
		char* temp=(char*)malloc(sizeof(char)*250);
		std::strcpy(temp,str);
		t.beg=std::atoi(std::strtok(temp,":"));
		t.end=std::atoi(std::strtok(NULL,":"));
		t.stride=std::atoi(std::strtok(NULL,":"));
		free(temp); temp=NULL;
	} else throw std::invalid_argument("Invalid time interval format: too many divisions.");
	return t;
}

std::ostream& operator<<(std::ostream& out, const time_interval& t){
	return out<<"("<<t.beg<<":"<<t.end<<":"<<t.stride<<")";
}

//**********************************************************************************************
//FILE_FORMAT struct
//**********************************************************************************************

FILE_FORMAT::type FILE_FORMAT::load(const std::string& str){
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
	return out;
}

//**********************************************************************************************
//Sim Interface
//**********************************************************************************************

//constructors/destructors

SimI::SimI(const SimI& sim){
	system_=sim.system();
	timesteps_=sim.timesteps();
	beg_=sim.beg();
	end_=sim.end();
	stride_=sim.stride();
	timestep_=sim.timestep();
	periodic_=sim.periodic();
	cells_=sim.cell();
	cellFixed_=sim.cellFixed();
	energy_=sim.energy();
}

//operators

SimI& SimI::operator=(const SimI& sim){
	timesteps_=sim.timesteps();
	beg_=sim.beg();
	end_=sim.end();
	stride_=sim.stride();
	timestep_=sim.timestep();
	periodic_=sim.periodic();
	cells_=sim.cell();
	cellFixed_=sim.cellFixed();
	energy_=sim.energy();
	return *this;
}

std::ostream& operator<<(std::ostream& out, const SimI& sim){
	out<<"SYSTEM = "<<sim.system_<<"\n";
	out<<"TIMESTEP = "<<sim.timestep_<<"\n";
	out<<"TIMESTEPS= "<<sim.timesteps_<<"\n";
	out<<"INTERVAL = ("<<sim.beg_+1<<","<<sim.end_+1<<","<<sim.stride_<<")\n";
	if(sim.cellFixed_) out<<"CELL_FIXED = TRUE";
	else out<<"CELL_FIXED = FALSE";
	return out;
}

//member functions

void SimI::defaults(){
	if(DEBUG_SIM>0) std::cout<<"SimI::defaults():\n";
	system_=std::string("Simulation");
	timesteps_=0;
	beg_=0;
	end_=0;
	stride_=1;
	timestep_=0.5;
	periodic_=true;
	cellFixed_=true;
	cells_.clear();
	energy_.clear();
}

void SimI::resize(unsigned int ts){
	if(DEBUG_SIM>0) std::cout<<"SimI::resize(unsigned int):\n";
	timesteps_=ts; 
	if(beg_==0 && end_==0){
		beg_=0;
		end_=timesteps_-1;
	}
	cells_.resize(timesteps_);
	energy_.resize(timesteps_);
}

//static functions

FILE_FORMAT::type SimI::fileFormat(const char* fileName){
	//first check if there is a file extension
	if(std::strpbrk(fileName, ".")!=NULL){
		FILE_FORMAT::type format;
		char* extension=(char*)malloc(sizeof(char)*string::M);
		std::strcpy(extension,strpbrk(fileName, ".")+1);
		string::to_upper(string::trim_all(extension));
		if(std::strcmp(extension, "XDATCAR")==0) format=FILE_FORMAT::XDATCAR;
		else if(std::strcmp(extension,"OUTCAR")==0) format=FILE_FORMAT::OUTCAR;
		else if(std::strcmp(extension,"POSCAR")==0)	format=FILE_FORMAT::POSCAR;
		else if(std::strcmp(extension,"XYZ")==0) format=FILE_FORMAT::XYZ;
		else if(std::strcmp(extension,"CAR")==0) format=FILE_FORMAT::CAR;
		else if(std::strcmp(extension,"GRO")==0) format=FILE_FORMAT::GROMACS;
		else format=FILE_FORMAT::UNKNOWN;
		free(extension);
		return format;
	} else return FILE_FORMAT::UNKNOWN;
}

//**********************************************************************************************
//Sim Atomic Interface
//**********************************************************************************************

//constructors/destructors

SimAtomicI::SimAtomicI(const SimAtomicI& sim){
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

std::ostream& operator<<(std::ostream& out, const SimAtomicI& sim){
	for(unsigned int i=0; i<sim.atomNames_.size(); ++i)
		out<<sim.atomNames_[i]<<" ";
	out<<"\n";
	for(unsigned int i=0; i<sim.nAtoms_.size(); ++i)
		out<<sim.nAtoms_[i]<<" ";
	return out;
}

SimAtomicI& SimAtomicI::operator=(const SimAtomicI& sim){
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

//member functions

void SimAtomicI::defaults(){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicI::defaults():\n";
	nSpecies_=0;
	nAtomsT_=0;
	nAtoms_.clear();
	offsets_.clear();
	atomNames_.clear();
}

void SimAtomicI::resize(const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& atomNames){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicI::resize(const std::vector<unsigned int>&,const std::vector<std::string>&):\n";
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

int SimAtomicI::speciesIndex(const std::string& str, const std::vector<std::string>& atomNames){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicI::speciesIndex(const std::string&, const std::vector<std::string>&):\n";
	for(unsigned int i=0; i<atomNames.size(); ++i) if(str==atomNames[i]) return i;
	return -1;
}

int SimAtomicI::speciesIndex(const char* str, const std::vector<std::string>& atomNames){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicI::speciesIndex(const char*,const std::vector<std::string>&):\n";
	for(unsigned int i=0; i<atomNames.size(); ++i) if(std::strcmp(str,atomNames[i].c_str())==0) return i;
	return -1;
}

//**********************************************************************************************
//Sim Molecular Interface
//**********************************************************************************************

//constructors/destructors

SimMolI::SimMolI(const SimMolI& sim){
	nMol_=sim.nMol();
}

//operators

std::ostream& operator<<(std::ostream& out, const SimMolI& sim){
	return out<<"N_MOL = "<<sim.nMol_;
}

SimMolI& SimMolI::operator=(const SimMolI& sim){
	nMol_=sim.nMol();
	return *this;
}

//member functions

void SimMolI::defaults(){
	nMol_=0;
}

void SimMolI::resize(unsigned int nMol){
	nMol_=nMol;
}
