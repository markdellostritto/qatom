#ifndef MDSIMULATION_UTILITY_LAMMPS
#define MDSIMULATION_UTILITY_LAMMPS

//c libraries
#include <ctime>
//boost
#include <type_traits>
//simulation
#include "sim.hpp"
//string
#include "string.hpp"
//chem info
#include "ptable.hpp"
//list
#include "list.hpp"
//units
#include "units.hpp"

#ifndef __cplusplus
	#error A C++ compiler is required.
#endif

#ifndef DEBUG_LAMMPS
#define DEBUG_LAMMPS 2
#endif

namespace LAMMPS{

//static variables
static const char* NAMESPACE_GLOBAL="LAMMPS";

//*****************************************************
//FORMAT struct
//*****************************************************

struct Format{
	std::string in;//input file
	std::string data;//data file
	std::string dump;//dump file
	int beg,end,stride;//interval
	Format():beg(1),end(-1),stride(1){};
	static Format& load(const std::vector<std::string>& strlist, Format& format);
};

//formats
struct STYLE_ATOM{
	enum type{
		FULL,
		BOND,
		ATOMIC
	};
};
struct FORMAT_ATOM{
	int atom,mol,specie;
	int x,y,z,q;
	FORMAT_ATOM():atom(-1),mol(-1),specie(-1),x(-1),y(-1),z(-1),q(-1){};
};
struct DATA_ATOM{
	int specie,index;
	double q;
	Eigen::Vector3d posn;
};

//template functions
template <class AtomT> void read_posn(AtomT& atom, Eigen::Vector3d posn, std::false_type){};
template <class AtomT> void read_posn(AtomT& atom, Eigen::Vector3d posn, std::true_type){atom.posn()=posn;};
template <class AtomT> void read_charge(AtomT& atom, double charge, std::false_type){};
template <class AtomT> void read_charge(AtomT& atom, double charge, std::true_type){atom.charge()=charge;};

namespace DUMP{

//static variables
static const char* NAMESPACE_LOCAL="DUMP";

template <class AtomT>
void load(const char* file, SimAtomic<AtomT>& sim, int beg, int end){
	static const char* funcName="load<AtomT>(const char*,SimAtomic<AtomT>&)";
	if(DEBUG_LAMMPS>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local function variables */
	//file i/o
		FILE* reader=NULL;
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//time info
		unsigned int ts=0;//number of timesteps
		unsigned int interval=0;//requested interval
	//cell info
		Eigen::Matrix3d lv=Eigen::Matrix3d::Zero();
	//atom info
		unsigned int nSpecies=0;//the number of atomic species
		std::vector<unsigned int> type;
		std::vector<unsigned int> nAtoms;//the number of atoms in each species
		std::vector<std::string> atomNames;//the names of each species
		unsigned int N=0;
		FORMAT_ATOM formatAtom;
	//simulation info
		bool direct;//whether the coordinates are in direct or Cartesian coordinates
	//timing
		clock_t start,stop;
		double time;
	//misc
		bool error=false;
	
	try{
		//start the timer
		start=std::clock();
		
		/* open the file */
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open file.");
		
		/* clear the simulation */
		sim.clear();
		
		/* find the number of timesteps */
		if(DEBUG_LAMMPS>0) std::cout<<"Loading number of timesteps...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"TIMESTEP")!=NULL) ++ts;
		}
		
		/* read in the atom format */
		if(DEBUG_LAMMPS>0) std::cout<<"Loading atom format...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				unsigned int nTokens=string::substrN(input,string::WS)-2;
				std::strtok(input,string::WS);
				std::strtok(NULL,string::WS);
				for(unsigned int i=0; i<nTokens; ++i){
					std::strcpy(temp,std::strtok(NULL,string::WS));
					if(std::strcmp(temp,"id")==0) formatAtom.atom=i;
					else if(std::strcmp(temp,"mol")==0) formatAtom.mol=i;
					else if(std::strcmp(temp,"type")==0) formatAtom.specie=i;
					else if(std::strcmp(temp,"q")==0) formatAtom.q=i;
					else if(std::strcmp(temp,"x")==0) formatAtom.x=i;
					else if(std::strcmp(temp,"y")==0) formatAtom.y=i;
					else if(std::strcmp(temp,"z")==0) formatAtom.z=i;
				}
				break;
			}
		}
		
		/* read in the atom info */
		if(DEBUG_LAMMPS>0) std::cout<<"Loading atom info...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				unsigned int nTokens=string::substrN(input,string::WS)-2;
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strstr(input,"ITEM")!=NULL) break;
					++N;//increment the total atom number
					unsigned int specie;
					for(unsigned int i=0; i<nTokens; ++i){
						if(i==0) std::strcpy(temp,std::strtok(input,string::WS));
						else std::strcpy(temp,std::strtok(NULL,string::WS));
						if(i==formatAtom.specie){specie=std::atoi(temp); break;}
					}
					int index=-1;
					for(int i=0; i<nSpecies; ++i) if(specie==type[i]){index=i; break;}
					if(index<0){
						++nSpecies;
						type.push_back(specie);
						nAtoms.push_back(1);
					} else ++nAtoms[index];
				}
				break;
			}
		}
		atomNames.resize(nSpecies);
		for(unsigned int i=0; i<nSpecies; ++i) atomNames[i]=std::to_string(type[i]);
		
		/* print data to screen */
		if(DEBUG_LAMMPS>1){
			std::cout<<"NAME = "<<sim.system()<<"\n";
			std::cout<<"SPECIES = ";
			for(unsigned int i=0; i<atomNames.size(); ++i){
				std::cout<<atomNames[i]<<" ";
			}
			std::cout<<"\n";
			std::cout<<"NUMBERS = ";
			for(unsigned int i=0; i<nAtoms.size(); ++i){
				std::cout<<nAtoms[i]<<" ";
			}
			std::cout<<"\n";
			std::cout<<"TIMESTEPS = "<<ts<<"\n";
		}
		
		/* set the timesteps */
		sim.beg()=beg-1;
		if(end<0){
			sim.end()=ts+end;
			end=sim.end();
		} else sim.end()=end-1;
		int interval=sim.end()-sim.beg()+1;
		if(DEBUG_LAMMPS>0) std::cout<<"INTERVAL = ("<<sim.beg()<<","<<sim.end()<<")\n";
		
		/* resize the simulation */
		if(DEBUG_LAMMPS>0) std::cout<<"Allocating memory...\n";
		sim.resize(interval,nAtoms,atomNames);
		
		/* read in the atom data */
		if(DEBUG_LAMMPS>0) std::cout<<"Loading atom data...\n";
		unsigned int timestep=0;
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"BOX")!=NULL){
				if(DEBUG_LAMMPS>1) std::cout<<"Cell: "<<timestep<<"\n";
				lv.setZero();
				fgets(input,string::M,reader);
				lv(0,0)-=std::atof(std::strtok(input,string::WS));
				lv(0,0)+=std::atof(std::strtok(NULL,string::WS));
				fgets(input,string::M,reader);
				lv(1,1)-=std::atof(std::strtok(input,string::WS));
				lv(1,1)+=std::atof(std::strtok(NULL,string::WS));
				fgets(input,string::M,reader);
				lv(2,2)-=std::atof(std::strtok(input,string::WS));
				lv(2,2)+=std::atof(std::strtok(NULL,string::WS));
				sim.cell(timestep).init(lv);
			}
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				if(DEBUG_LAMMPS>1) std::cout<<"Atoms: "<<timestep<<"\n";
				if(timestep<sim.beg()){++timestep; continue;}
				else if(timestep>sim.end()) break;
				unsigned int nTokens=string::substrN(input,string::WS)-2;
				std::vector<std::string> tokens(nTokens);
				DATA_ATOM dataAtom;
				for(unsigned int i=0; i<N; ++i){
					//read in the next line
					fgets(input,string::M,reader);
					//split the line into tokens
					tokens[0]=std::string(std::strtok(input,string::WS));
					for(int j=1; j<nTokens; ++j) tokens[j]=std::string(std::strtok(NULL,string::WS));
					//read in the data
					if(formatAtom.q>=0) dataAtom.q=std::atof(tokens[formatAtom.q].c_str());
					if(formatAtom.x>=0) dataAtom.posn[0]=std::atof(tokens[formatAtom.x].c_str());
					if(formatAtom.y>=0) dataAtom.posn[1]=std::atof(tokens[formatAtom.y].c_str());
					if(formatAtom.z>=0) dataAtom.posn[2]=std::atof(tokens[formatAtom.z].c_str());
					if(formatAtom.atom>=0) dataAtom.index=std::atof(tokens[formatAtom.atom].c_str());
					if(formatAtom.specie>=0) dataAtom.specie=std::atof(tokens[formatAtom.specie].c_str());
					for(int j=dataAtom.specie-1; j>=0; --j) dataAtom.index-=nAtoms[j];
					//set the simulation data
					read_posn(sim.atom(timestep-sim.beg(),dataAtom.specie,dataAtom.index),dataAtom.posn,std::is_base_of<Position,AtomT>());
					read_charge(sim.atom(timestep-sim.beg(),dataAtom.specie,dataAtom.index),dataAtom.q,std::is_base_of<Charge,AtomT>());
				}
				++timestep;
			}
		}
		
		//stop the timer
		stop=std::clock();
		
		//print the time
		time=((double)(stop-start))/CLOCKS_PER_SEC;
		std::cout<<"Positions loaded in "<<time<<" seconds.\n";
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	if(reader!=NULL) fclose(reader);
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

template <class AtomT>
void load_posn(const char* file, SimAtomic<AtomT>& sim, int beg, int end){
	static const char* funcName="load_posn<AtomT>(const char*,SimAtomic<AtomT>&)";
	if(DEBUG_LAMMPS>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local function variables */
	//file i/o
		FILE* reader=NULL;
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//time info
		unsigned int ts=0;//number of timesteps
		unsigned int interval=0;//requested interval
	//cell
		Eigen::Matrix3d lv=Eigen::Matrix3d::Zero();
	//atom info
		unsigned int nSpecies=0;//the number of atomic species
		std::vector<unsigned int> types;
		std::vector<unsigned int> nAtoms;//the number of atoms in each species
		std::vector<std::string> atomNames;//the names of each species
		unsigned int N=0;
		FORMAT_ATOM formatAtom;
		std::vector<std::vector<unsigned int> > indices;
		std::vector<unsigned int> IDs;
	//simulation info
		bool direct;//whether the coordinates are in direct or Cartesian coordinates
	//timing
		clock_t start,stop;
		double time;
	//misc
		bool error=false;
	
	try{
		//start the timer
		start=std::clock();
		
		/* open the file */
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open file.");
		
		/* find the number of timesteps */
		if(DEBUG_LAMMPS>1) std::cout<<"Loading number of timesteps...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"TIMESTEP")!=NULL) ++ts;
		}
		
		/* read in the simulation cell */
		if(DEBUG_LAMMPS>1) std::cout<<"Loading simulation cell...\n";
		sim.periodic()=true;
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"BOX")!=NULL){
				fgets(input,string::M,reader);
				lv(0,0)-=std::atof(std::strtok(input,string::WS));
				lv(0,0)+=std::atof(std::strtok(NULL,string::WS));
				fgets(input,string::M,reader);
				lv(1,1)-=std::atof(std::strtok(input,string::WS));
				lv(1,1)+=std::atof(std::strtok(NULL,string::WS));
				fgets(input,string::M,reader);
				lv(2,2)-=std::atof(std::strtok(input,string::WS));
				lv(2,2)+=std::atof(std::strtok(NULL,string::WS));
				break;
			}
		}
		
		/* read in the atom format */
		if(DEBUG_LAMMPS>1) std::cout<<"Loading atom format...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				unsigned int nTokens=string::substrN(input,string::WS)-2;
				std::strtok(input,string::WS);
				std::strtok(NULL,string::WS);
				for(unsigned int i=0; i<nTokens; i++){
					std::strcpy(temp,std::strtok(NULL,string::WS));
					if(std::strcmp(temp,"id")==0) formatAtom.atom=i;
					else if(std::strcmp(temp,"mol")==0) formatAtom.mol=i;
					else if(std::strcmp(temp,"type")==0) formatAtom.specie=i;
					else if(std::strcmp(temp,"q")==0) formatAtom.q=i;
					else if(std::strcmp(temp,"x")==0) formatAtom.x=i;
					else if(std::strcmp(temp,"y")==0) formatAtom.y=i;
					else if(std::strcmp(temp,"z")==0) formatAtom.z=i;
				}
				break;
			}
		}
		if(DEBUG_LAMMPS>1){
			std::cout<<"FORMAT_ATOM = \n";
			std::cout<<"\tid = "<<formatAtom.atom<<"\n";
			std::cout<<"\tmol = "<<formatAtom.mol<<"\n";
			std::cout<<"\ttype = "<<formatAtom.specie<<"\n";
			std::cout<<"\tq = "<<formatAtom.q<<"\n";
			std::cout<<"\tx = "<<formatAtom.x<<"\n";
			std::cout<<"\ty = "<<formatAtom.y<<"\n";
			std::cout<<"\tz = "<<formatAtom.z<<"\n";
		}
		
		/* read in the atom info */
		if(DEBUG_LAMMPS>1) std::cout<<"Loading atom info...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				unsigned int nTokens=string::substrN(input,string::WS)-2;
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strstr(input,"ITEM")!=NULL) break;
					++N;//increment the total atom number
					//split the string into tokens
					std::vector<std::string> tokens(nTokens);
					tokens[0]=std::string(std::strtok(input,string::WS));
					for(unsigned int j=1; j<nTokens; ++j) tokens[j]=std::string(std::strtok(NULL,string::WS));
					//find the specie (may neither be zero-valued nor continuous)
					unsigned int specie=std::atoi(tokens[formatAtom.specie].c_str());
					unsigned int index=std::atoi(tokens[formatAtom.atom].c_str());
					int sIndex=-1;
					for(unsigned int i=0; i<nSpecies; ++i) if(specie==types[i]){sIndex=i; break;}
					if(sIndex<0){
						++nSpecies;
						types.push_back(specie);
						nAtoms.push_back(1);
					} else ++nAtoms[sIndex];
				}
				break;
			}
		}
		//sort the types
		list::insertionSort(types,nAtoms);
		//set the atom names
		atomNames.resize(nSpecies);
		for(unsigned int i=0; i<nSpecies; ++i) atomNames[i]=std::to_string(types[i]);
		
		/* read in the indices */
		indices.resize(types.size());
		if(DEBUG_LAMMPS>1) std::cout<<"Loading atom info...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				unsigned int nTokens=string::substrN(input,string::WS)-2;
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strstr(input,"ITEM")!=NULL) break;
					//split the string into tokens
					std::vector<std::string> tokens(nTokens);
					tokens[0]=std::string(std::strtok(input,string::WS));
					for(unsigned int j=1; j<nTokens; ++j) tokens[j]=std::string(std::strtok(NULL,string::WS));
					//find the specie (may neither be zero-valued nor continuous)
					unsigned int specie=std::atoi(tokens[formatAtom.specie].c_str());
					unsigned int index=std::atoi(tokens[formatAtom.atom].c_str());
					int sIndex=-1;
					for(unsigned int i=0; i<nSpecies; ++i) if(specie==types[i]){sIndex=i; break;}
					if(sIndex>=0) indices[sIndex].push_back(index);
				}
				break;
			}
		}
		//sort the indices
		for(unsigned int i=0; i<indices.size(); ++i) list::insertionSort(indices[i]);
		//set the IDs
		IDs.resize(N);
		for(int i=0; i<indices.size(); i++){
			for(int j=0; j<indices[i].size(); j++){
				IDs[indices[i][j]-1]=j;
			}
		}
		
		/* print data to screen */
		if(DEBUG_LAMMPS>1){
			std::cout<<"NAME = "<<sim.system()<<"\n";
			std::cout<<"R = \n"<<lv<<"\n";
			std::cout<<"SPECIES = ";
			for(unsigned int i=0; i<atomNames.size(); ++i) std::cout<<atomNames[i]<<" ";
			std::cout<<"\n";
			std::cout<<"NUMBERS = ";
			for(unsigned int i=0; i<nAtoms.size(); ++i) std::cout<<nAtoms[i]<<" ";
			std::cout<<"\n";
			std::cout<<"TYPES = ";
			for(unsigned int i=0; i<types.size(); ++i) std::cout<<types[i]<<" ";
			std::cout<<"\n";
			std::cout<<"TIMESTEPS = "<<ts<<"\n";
		}
		
		/* check the simulation data */
		if(N!=sim.nAtoms()) throw std::invalid_argument("Mismatch in number of atoms.");
		if(nSpecies!=sim.nSpecies()) throw std::invalid_argument("Mismatch in number of species.");
		
		/* set the timesteps */
		sim.beg()=beg-1;
		if(end<0){
			sim.end()=ts+end;
			end=sim.end();
		} else sim.end()=end-1;
		int interval=sim.end()-sim.beg()+1;
		if(DEBUG_LAMMPS>0) std::cout<<"Interval = ("<<sim.beg()<<","<<sim.end()<<")\n";
		if(DEBUG_LAMMPS>0) std::cout<<"Length = "<<interval<<"\n";
		
		/* resize the simulation */
		if(DEBUG_LAMMPS>0) std::cout<<"Allocating memory...\n";
		sim.resize(interval,sim.nAtomsVec(),sim.atomNames());
		
		/* read in the atom data */
		if(DEBUG_LAMMPS>0) std::cout<<"Loading atom data...\n";
		unsigned int timestep=0;
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL && timestep<sim.beg()){
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				if(DEBUG_LAMMPS>2) std::cout<<"Atoms: "<<timestep<<"\n";
				else if(timestep%1000==0) std::cout<<"Atoms: = "<<timestep<<"\n";
				for(unsigned int i=0; i<N; ++i) fgets(input,string::M,reader);
				++timestep;
			}
		}
		timestep=0;
		while(fgets(input,string::M,reader)!=NULL && timestep<interval){
			if(std::strstr(input,"BOX")!=NULL){
				if(DEBUG_LAMMPS>2) std::cout<<"Cell: "<<timestep<<"\n";
				lv.setZero();
				fgets(input,string::M,reader);
				lv(0,0)-=std::atof(std::strtok(input,string::WS));
				lv(0,0)+=std::atof(std::strtok(NULL,string::WS));
				fgets(input,string::M,reader);
				lv(1,1)-=std::atof(std::strtok(input,string::WS));
				lv(1,1)+=std::atof(std::strtok(NULL,string::WS));
				fgets(input,string::M,reader);
				lv(2,2)-=std::atof(std::strtok(input,string::WS));
				lv(2,2)+=std::atof(std::strtok(NULL,string::WS));
				sim.cell(timestep).init(lv);
			}
			if(std::strstr(input,"ITEM: ATOMS")!=NULL){
				if(DEBUG_LAMMPS>2) std::cout<<"Atoms: "<<timestep<<"\n";
				else if(timestep%1000==0) std::cout<<"Atoms: = "<<timestep<<"\n";
				unsigned int nTokens=string::substrN(input,string::WS)-2;
				std::vector<std::string> tokens(nTokens);
				DATA_ATOM dataAtom;
				for(unsigned int i=0; i<N; ++i){
					//read in the next line
					fgets(input,string::M,reader);
					//split the line into tokens
					tokens[0]=std::string(std::strtok(input,string::WS));
					for(unsigned int j=1; j<nTokens; ++j) tokens[j]=std::string(std::strtok(NULL,string::WS));
					//read in the data
					if(formatAtom.atom>=0) dataAtom.index=std::atoi(tokens[formatAtom.atom].c_str());
					if(formatAtom.specie>=0) dataAtom.specie=std::atoi(tokens[formatAtom.specie].c_str());
					if(formatAtom.q>=0) dataAtom.q=std::atof(tokens[formatAtom.q].c_str());
					if(formatAtom.x>=0) dataAtom.posn[0]=std::atof(tokens[formatAtom.x].c_str());
					if(formatAtom.y>=0) dataAtom.posn[1]=std::atof(tokens[formatAtom.y].c_str());
					if(formatAtom.z>=0) dataAtom.posn[2]=std::atof(tokens[formatAtom.z].c_str());
					for(unsigned int j=0; j<types.size(); ++j) if(dataAtom.specie==types[j]){dataAtom.specie=j; break;}
					dataAtom.index=IDs[dataAtom.index-1];
					//set the simulation data
					read_posn(sim.atom(timestep,dataAtom.specie,dataAtom.index),dataAtom.posn,std::is_base_of<Position,AtomT>());
					read_charge(sim.atom(timestep,dataAtom.specie,dataAtom.index),dataAtom.q,std::is_base_of<Charge,AtomT>());
				}
				++timestep;
			}
		}
		
		//return to cell
		for(int t=0; t<sim.timesteps(); ++t){
			for(int n=0; n<sim.nSpecies(); ++n){
				for(int m=0; m<sim.nAtoms(n); ++m){
					Cell::returnToCell(sim.atom(t,n,m).posn(),sim.atom(t,n,m).posn(),sim.cell(t).R(),sim.cell(t).RInv());
				}
			}
		}
		
		//stop the timer
		stop=std::clock();
		
		//print the time
		time=((double)(stop-start))/CLOCKS_PER_SEC;
		std::cout<<"Positions loaded in "<<time<<" seconds.\n";
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	if(reader!=NULL) fclose(reader);
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

template <class AtomT>
void write(const char* file, SimAtomic<AtomT>& sim, int beg=1, int end=-1){
	static const char* funcName="write<AtomT>(const char*,SimAtomic<AtomT>&,int,int)";
	if(DEBUG_LAMMPS>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	//local variables
	FILE* writer=NULL;
	bool error=false;
	
	try{
		//open the file
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("Unable to open file.");
		
		//set the beginning and ending timesteps
		int lbeg=beg-1;
		int lend=end;
		if(lend<0) lend=sim.timesteps()+end;
		if(lbeg<0 || lend>sim.timesteps() || lbeg>lend) throw std::invalid_argument("Invalid beginning and ending timesteps.");
		if(DEBUG_LAMMPS>0) std::cout<<"Interval = ("<<lbeg<<","<<lend<<")\n";
		
		for(unsigned int t=lbeg; t<=lend; ++t){
			fprintf(writer,"ITEM: TIMESTEP\n");
			fprintf(writer,"%i\n",t);
			fprintf(writer,"ITEM: NUMBER OF ATOMS\n");
			fprintf(writer,"%i\n",sim.nAtoms());
			fprintf(writer,"ITEM: BOX BOUNDS pp pp pp\n");
			fprintf(writer,"%f %f\n",0.0,sim.cell(t).R()(0,0));
			fprintf(writer,"%f %f\n",0.0,sim.cell(t).R()(1,1));
			fprintf(writer,"%f %f\n",0.0,sim.cell(t).R()(2,2));
			fprintf(writer,"ITEM: ATOMS id type x y z\n");
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				fprintf(writer,"%i %i %f %f %f\n",n,sim.atom(t,n).specie(),
					sim.atom(t,n).posn()[0],sim.atom(t,n).posn()[1],sim.atom(t,n).posn()[2]
				);
			}
		}
		
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

}

namespace DATA{

//static variables
static const char* NAMESPACE_LOCAL="DATA";

template <class AtomT>
void load(const char* file, SimAtomic<AtomT>& sim, STYLE_ATOM::type& styleAtom){
	static const char* funcName="load<AtomT>(const char*,SimAtomic<AtomT>&,STYLE_ATOM::type&)";
	if(DEBUG_LAMMPS>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/*
		local function variables
	*/
	//file i/o
		FILE* reader;
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//atom info
		unsigned int nSpecies=0;//the number of atomic species
		std::vector<unsigned int> nAtoms;//the number of atoms in each species
		std::vector<std::string> atomNames;//the names of each species
		unsigned int N=0;//total number of atoms
		std::vector<std::pair<int,double> > masses;
		std::vector<std::pair<int,int> > types;
	//cell info
		Eigen::Matrix3d lv=Eigen::Matrix3d::Zero();
	//miscellaneous
		bool error=false;
	
	try{
		/* open the file */
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open file.");
		
		/* clear the simulation */
		sim.clear();
		
		/* read in the number of atoms, number of species, and cell */
		if(DEBUG_LAMMPS>1) std::cout<<"Reading in number of atoms, number of species, and cell...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,"#");
			if(string::empty(input)) continue;
			//split the string into tokens
			unsigned int nTokens=string::substrN(input,string::WS);
			if(nTokens>0){
				std::vector<std::string> tokens(nTokens);
				tokens[0]=std::string(std::strtok(input,string::WS));
				for(unsigned int i=1; i<nTokens; ++i) tokens[i]=std::string(std::strtok(NULL,string::WS));
				if(nTokens==1){
					if(std::strcmp(tokens[0].c_str(),"Atoms")==0) break;
				} else if(nTokens==2){
					if(tokens[1]=="atoms") N=std::atoi(tokens[0].c_str());
				} else if(nTokens==3){
					if(tokens[1]=="atom" && tokens[2]=="types") nSpecies=std::atoi(tokens[0].c_str());
				} else if(nTokens==4){
					if(tokens[2]=="xlo" && tokens[3]=="xhi") lv(0,0)=std::atof(tokens[1].c_str())-std::atof(tokens[0].c_str());
					else if(tokens[2]=="ylo" && tokens[3]=="yhi") lv(1,1)=std::atof(tokens[1].c_str())-std::atof(tokens[0].c_str());
					else if(tokens[2]=="zlo" && tokens[3]=="zhi") lv(2,2)=std::atof(tokens[1].c_str())-std::atof(tokens[0].c_str());
				}
			}
		}
		
		/* read in the masses */
		if(DEBUG_LAMMPS>1) std::cout<<"Reading in masses...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,"#");
			if(string::empty(input)) continue;
			//split the string into tokens
			unsigned int nTokens=string::substrN(input,string::WS);
			if(nTokens>0){
				std::vector<std::string> tokens(nTokens);
				tokens[0]=std::string(std::strtok(input,string::WS));
				for(unsigned int i=1; i<nTokens; ++i) tokens[i]=std::string(std::strtok(NULL,string::WS));
				if(nTokens==1){
					if(std::strcmp(tokens[0].c_str(),"Masses")==0){
						//skip a line
						fgets(input,string::M,reader);
						//read in the atom names
						masses.resize(nSpecies);
						for(unsigned int i=0; i<nSpecies; ++i){
							fgets(input,string::M,reader);
							int index=std::atoi(std::strtok(input,string::WS));
							double mass=std::atof(std::strtok(NULL,string::WS));
							masses[i].first=index; masses[i].second=mass;
						}
						break;
					}
				}
			}
		}
		
		/* read in the atom numbers */
		if(DEBUG_LAMMPS>1) std::cout<<"Reading in atom numbers...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,"#");
			string::trim(input);
			if(std::strcmp(input,"Atoms")==0){
				types.resize(nSpecies);
				for(unsigned int i=0; i<nSpecies; ++i){
					types[i].first=masses[i].first;
					types[i].second=0;
				}
				//skip a line
				fgets(input,string::M,reader);
				//read through all the atoms
				if(styleAtom==STYLE_ATOM::FULL){
					for(int i=0; i<N; ++i){
						//read in the specie
						fgets(input,string::M,reader);
						std::strtok(input,string::WS);
						std::strtok(NULL,string::WS);
						unsigned int specie=std::atoi(std::strtok(NULL,string::WS));
						//sort the specie
						int index=-1;
						for(unsigned int j=0; j<types.size(); ++j) if(specie==types[j].first){index=j; break;}
						types[index].second++;
					}
				} else if(styleAtom==STYLE_ATOM::ATOMIC){
					for(int i=0; i<N; ++i){
						//read in the specie
						fgets(input,string::M,reader);
						std::strtok(input,string::WS);
						unsigned int specie=std::atoi(std::strtok(NULL,string::WS));
						//sort the specie
						int index=-1;
						for(unsigned int j=0; j<types.size(); ++j) if(specie==types[j].first){index=j; break;}
						types[index].second++;
					}
				}
				break;
			}
		}
		
		/* set the names and numbers */
		if(DEBUG_LAMMPS>1) std::cout<<"Setting atom names and numbers...\n";
		atomNames.resize(nSpecies);
		nAtoms.resize(nSpecies);
		for(int i=0; i<nSpecies; i++){
			atomNames[i]=std::string(PTable::elementNameMass(masses[i].second));
			nAtoms[i]=types[i].second;
		}
		
		/* resize the simulation */
		if(DEBUG_LAMMPS>1) std::cout<<"Resizing the simulation...\n";
		sim.resize(1,nAtoms,atomNames);
		sim.cell(0).init(lv);
		
		/* print data to screen */
		if(DEBUG_LAMMPS>1){
			std::cout<<"N = "<<N<<"\n";
			std::cout<<"N_SPECIES = "<<nSpecies<<"\n";
			std::cout<<"NAME = "<<sim.system()<<"\n";
			std::cout<<"CELL = \n"<<sim.cell(0)<<"\n";
			std::cout<<"SPECIES = ";
			for(unsigned int i=0; i<atomNames.size(); ++i) std::cout<<atomNames[i]<<" ";
			std::cout<<"\n";
			std::cout<<"NUMBERS = ";
			for(unsigned int i=0; i<nAtoms.size(); ++i) std::cout<<nAtoms[i]<<" ";
			std::cout<<"\n";
		}
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
	}
	
	//free local variables
	if(reader!=NULL) fclose(reader);
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

}

namespace IN{

//static variables
static const char* NAMESPACE_LOCAL="IN";

void load_style(const char* file, STYLE_ATOM::type& styleAtom);

}

template <class AtomT>
void load(const char* inFile, const char* dataFile, const char* dumpFile, SimAtomic<AtomT>& sim, int beg, int end){
	static const char* funcName="load<AtomT>(const char*,const char*,const char*,SimAtomic<AtomT>&)";
	if(DEBUG_LAMMPS>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<funcName<<":\n";
	/* local function variables */
	//atom style
		STYLE_ATOM::type styleAtom;
	//units
		double s=0.0;
		if(units::consts::system()==units::System::AU) s=units::BOHRpANG;
		else if(units::consts::system()==units::System::METAL) s=1.0;
		else throw std::runtime_error("Invalid units.");
	//miscellaneous
		bool error=false;
	//timing
		clock_t start,stop;
		double time;
		
	try{
		//start the timer
		start=std::clock();
		
		//load the atom style
		IN::load_style(inFile,styleAtom);
		
		//load the generic simulation data
		DATA::load(dataFile,sim,styleAtom);
		
		//load the positions
		DUMP::load_posn(dumpFile,sim,beg,end);
		
		//set the units
		if(s!=1.0){
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				sim.cell(t).init(sim.cell(t).R()*s);
			}
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					sim.atom(t,n).posn()*=s;
				}
			}
		}
		std::cout<<"************ CELL = \n"<<sim.cell(0)<<"\n";
		
		//stop the timer
		stop=std::clock();
		
		//print the time
		time=((double)(stop-start))/CLOCKS_PER_SEC;
		std::cout<<"Simulation loaded in "<<time<<" seconds.\n";
	}catch(std::exception& e){
		std::cout<<NAMESPACE_GLOBAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
	}
}

}

#endif
