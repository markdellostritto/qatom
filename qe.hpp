#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <chrono>
#include "sim.hpp"
#include "string.hpp"
#include "atom.hpp"
#include "cell.hpp"

#ifndef DEBUG_QE
#define DEBUG_QE 1
#endif

namespace QE{

//static variables
static const char* NAMESPACE_GLOBAL="QE";

//*****************************************************
//FORMAT struct
//*****************************************************

struct Format{
	std::string fileIn;//input
	std::string filePos;//position
	std::string fileCel;//cell
	std::string fileEvp;//energy/volume/pressure
	int beg,end,stride;
	Format():beg(1),end(-1),stride(1){};
	static Format& load(const std::vector<std::string>& strlist, Format& format);
};

//*****************************************************
//IN format
//*****************************************************

namespace IN{

//static variables
static const char* NAMESPACE_LOCAL="IN";
	
Cell& load_cell(FILE* reader, Cell& cell);
void load_atoms(FILE* reader, std::vector<std::string>& atomNames, std::vector<unsigned int>& atomNumbers);
double load_timestep(FILE* reader);

}

//*****************************************************
//CEL format
//*****************************************************

namespace CEL{

//static variables
static const char* NAMESPACE_LOCAL="CEL";

void load_cell(FILE* reader, SimI& sim);
	
}

//*****************************************************
//POS format
//*****************************************************

namespace POS{

//static variables
static const char* NAMESPACE_LOCAL="POS";

unsigned int load_timesteps(FILE* reader);

template <class AtomT>
void load_posns(FILE* reader, SimAtomic<AtomT>& sim){
	const char* func_name="load_posns(const char*,SimAtomic<AtomT>&)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//scaling
		double s=0.52917721067;//1 bohr (angstroms)
	//misc
		bool error=false;
		
	try{
		//rewind to beginning of file
		std::rewind(reader);
		
		//skip to the beginning
		for(unsigned int t=0; t<sim.beg(); ++t){
			if(DEBUG_QE>1) std::cout<<"T = "<<t<<"\n";
			else if(t%1000==0) std::cout<<"T = "<<t<<"\n";
			fgets(input,string::M,reader);//header
			for(unsigned int n=0; n<sim.nAtoms(); ++n) fgets(input,string::M,reader);
		}
		//read in the data
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			if(DEBUG_QE>1) std::cout<<"T = "<<t<<"\n";
			else if(t%1000==0) std::cout<<"T = "<<t<<"\n";
			fgets(input,string::M,reader);//header
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				fgets(input,string::M,reader);
				sim.atom_posn()(t,n)[0]=s*std::atof(std::strtok(input,string::WS));
				sim.atom_posn()(t,n)[1]=s*std::atof(std::strtok(NULL,string::WS));
				sim.atom_posn()(t,n)[2]=s*std::atof(std::strtok(NULL,string::WS));
			}
			for(unsigned int tt=1; tt<sim.stride(); ++tt){
				fgets(input,string::M,reader);//header
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					fgets(input,string::M,reader);
				}
			}
		}
		
		//return to cell
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				Cell::returnToCell(sim.atom_posn()(t,n),sim.atom_posn()(t,n),sim.cell(t).R(),sim.cell(t).RInv());
			}
		}
	}catch(std::exception& e){
		std::cout<<"Error in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("I/O Error: Failed to load.");
}

}

//*****************************************************
//EVP format
//*****************************************************

namespace EVP{
	
	//static variables
	static const char* NAMESPACE_LOCAL="EVP";
	
	void load_energy(FILE* reader, SimI& sim);
	
}

template <class AtomT>
SimAtomic<AtomT>& load(const char* infile, const char* posfile, const char* celfile, SimAtomic<AtomT>& sim, int beg, int end){
	const char* func_name="load(const char*,const char*,SimAtomic<AtomT>&,int,int)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		FILE* reader=NULL;
	//simulation parameters
		Cell cell;
		unsigned int ts=0;
		int interval=0;
		std::vector<std::string> atomNames;
		std::vector<unsigned int> atomNumbers;
	//timing
		std::chrono::high_resolution_clock::time_point start;
		std::chrono::high_resolution_clock::time_point stop;
		std::chrono::duration<double> time;
	//misc
		bool error=false;
	
	try{		
		//start the timer
		start=std::chrono::high_resolution_clock::now();
		
		//open the "in" file
		reader=fopen(infile,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"in\" file.");
		//load the cell
		IN::load_cell(reader,cell);
		//load the atom info
		IN::load_atoms(reader,atomNames,atomNumbers);
		//load the timestep
		sim.timestep()=IN::load_timestep(reader);
		//set the periodicity
		sim.periodic()=true;
		//close the "in" file
		fclose(reader);
		reader=NULL;
		
		//open the "pos" file
		reader=fopen(posfile,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"pos\" file.");
		//load the number of timesteps
		ts=POS::load_timesteps(reader);
		//set the interval
		if(beg<1) throw std::invalid_argument("Invalid beginning timestep.");
		sim.beg()=beg-1;
		if(end<0){sim.end()=ts+end; end=sim.end();} else sim.end()=end-1;
		interval=sim.end()-sim.beg()+1;
		//resize the simulation
		sim.resize(interval,atomNumbers,atomNames);
		if(DEBUG_QE>0) std::cout<<"SIM = \n"<<sim<<"\n";
		//close the "pos" file
		fclose(reader);
		reader=NULL;
		
		//open the "cel" file
		reader=fopen(celfile,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"cel\" file.");
		//load the cell
		CEL::load_cell(reader,sim);
		//close the "cel" file
		fclose(reader);
		reader=NULL;
		
		//open the "pos" file
		reader=fopen(posfile,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"pos\" file.");
		//load the positions
		POS::load_posns(reader,sim);
		//close the "pos" file
		fclose(reader);
		reader=NULL;
		
		//stop the timer
		stop=std::chrono::high_resolution_clock::now();
		
		//print the time
		time=std::chrono::duration_cast<std::chrono::duration<double> >(stop-start);
		std::cout<<"Simulation loaded in "<<time.count()<<" seconds.\n";
	}catch(std::exception& e){
		std::cout<<"Error in "<<NAMESPACE_GLOBAL<<"::"<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("I/O Error: Failed to load.");
	else return sim;
}

template <class AtomT>
SimAtomic<AtomT>& load(Format& format, SimAtomic<AtomT>& sim){
	const char* func_name="load(Format&,SimAtomic<AtomT>&)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		FILE* reader=NULL;
	//simulation parameters
		Cell cell;
		unsigned int ts=0;
		int interval=0;
		std::vector<std::string> atomNames;
		std::vector<unsigned int> atomNumbers;
	//timing
		std::chrono::high_resolution_clock::time_point start;
		std::chrono::high_resolution_clock::time_point stop;
		std::chrono::duration<double> time;
	//misc
		bool error=false;
	
	try{		
		//start the timer
		start=std::chrono::high_resolution_clock::now();
		
		if(format.fileIn.empty() && format.fileCel.empty()) throw std::runtime_error("No cell information included.");
		if(format.filePos.empty()) throw std::runtime_error("No position information included.");
		
		if(!format.fileIn.empty()){
			//open the "in" file
			reader=fopen(format.fileIn.c_str(),"r");
			if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"in\" file.");
			//load the cell
			IN::load_cell(reader,cell);
			//load the atom info
			IN::load_atoms(reader,atomNames,atomNumbers);
			//load the timestep
			sim.timestep()=IN::load_timestep(reader);
			//set the periodicity
			sim.periodic()=true;
			//close the "in" file
			fclose(reader); reader=NULL;
		}
		
		if(!format.filePos.empty()){
			//open the "pos" file
			reader=fopen(format.filePos.c_str(),"r");
			if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"pos\" file.");
			//load the number of timesteps
			ts=POS::load_timesteps(reader);
			//close the "pos" file
			fclose(reader); reader=NULL;
		} else throw std::runtime_error("Found no POS file.");
		
		//set the interval
		if(format.beg<1) throw std::invalid_argument("Invalid beginning timestep.");
		sim.beg()=format.beg-1;
		if(format.end<0){
			sim.end()=ts+format.end;
			format.end=sim.end();
		} else {
			if(format.end>ts) throw std::invalid_argument("Invalid ending timestep.");
			sim.end()=format.end-1;
		}
		interval=sim.end()-sim.beg()+1;
		
		//resize the simulation
		sim.resize(interval/format.stride,atomNumbers,atomNames);
		sim.stride()=format.stride;
		if(DEBUG_QE>0) std::cout<<"ts = "<<ts<<"\n";
		if(DEBUG_QE>0) std::cout<<"interval = "<<format.beg<<":"<<format.end<<":"<<format.stride<<"\n";
		if(DEBUG_QE>0) std::cout<<"SIM = \n"<<sim<<"\n";
		
		if(!format.fileCel.empty()){
			//open the "cel" file
			reader=fopen(format.fileCel.c_str(),"r");
			if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"cel\" file.");
			//load the cell
			CEL::load_cell(reader,sim);
			//close the "cel" file
			fclose(reader); reader=NULL;
		}
		
		if(!format.fileEvp.empty()){
			//open the "evp" file
			reader=fopen(format.fileEvp.c_str(),"r");
			if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"evp\" file.");
			//load the energy
			EVP::load_energy(reader,sim);
			//close the "evp" file
			fclose(reader); reader=NULL;
		}
		
		if(!format.filePos.empty()){
			//open the "pos" file
			reader=fopen(format.filePos.c_str(),"r");
			if(reader==NULL) throw std::runtime_error("I/O Error: Could not open \"pos\" file.");
			//load the positions
			POS::load_posns(reader,sim);
			//close the "pos" file
			fclose(reader); reader=NULL;
		}
		
		//stop the timer
		stop=std::chrono::high_resolution_clock::now();
		
		//print the time
		time=std::chrono::duration_cast<std::chrono::duration<double> >(stop-start);
		std::cout<<"Simulation loaded in "<<time.count()<<" seconds.\n";
	}catch(std::exception& e){
		std::cout<<"Error in "<<NAMESPACE_GLOBAL<<"::"<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("I/O Error: Failed to load.");
	else return sim;
}

}