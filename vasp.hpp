#ifndef VASP_HPP
#define VASP_HPP

#include <cstdio>
#include <chrono>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <Eigen/Dense>
#include "sim.hpp"
#include "cell.hpp"
#include "string.hpp"
#include "units.hpp"

#ifndef DEBUG_VASP
#define DEBUG_VASP 0
#endif

#ifndef __cplusplus
	#error A C++ compiler is required
#endif

namespace VASP{

//static variables
static const int HEADER_SIZE=7;//number of lines in the header before the atomic positions
static const char* NAMESPACE_GLOBAL="VASP";

//utilities
template <class AtomT> void setAN(AtomT& atom, std::false_type){};
template <class AtomT> void setAN(AtomT& atom, std::true_type){
	atom.an()=PTable::an(atom.name().c_str());
}

//*****************************************************
//FORMAT struct
//*****************************************************

struct Format{
	std::string xdatcar;//xdatcar
	std::string poscar;//poscar
	std::string xml;//poscar
	std::string outcar;//outcar
	std::string procar;//procar
	std::string eigenval;//eigenval
	std::string energy;//energy
	int beg,end,stride;//interval
	Format():beg(1),end(-1),stride(1){};
	static Format& load(const std::vector<std::string>& strlist, Format& format);
};

namespace XDATCAR{

//static variables
static const char* NAMESPACE_LOCAL="XDATCAR";

template <class AtomT>
void load(const char* file, SimAtomic<AtomT>& sim, int beg=1, int end=-1, int stride=1){
	static const char* funcName="load<AtomT>(const char*,SimAtomic<AtomT>&,int,int)";
	if(DEBUG_VASP>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local function variables */
	//file i/o
		FILE* reader=NULL;
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
		std::string str;
	//simulation flags
		bool direct;//whether the coordinates are in direct or Cartesian coordinates
	//time info
		unsigned int ts=0;//number of timesteps
		unsigned int interval=0;//requested interval
	//cell info
		double scale=1;
		Eigen::Matrix3d lv;
		Cell cell;
	//atom info
		unsigned int nSpecies=0;//the number of atomic species
		std::vector<unsigned int> nAtoms;//the number of atoms in each species
		std::vector<std::string> atomNames;//the names of each species
		unsigned int N=0;
	//timing
		std::chrono::high_resolution_clock::time_point start;
		std::chrono::high_resolution_clock::time_point stop;
		std::chrono::duration<double> time;
	//units
		double s=0.0;
		if(units::consts::system()==units::System::AU) s=units::BOHRpANG;
		else if(units::consts::system()==units::System::METAL) s=1.0;
		else throw std::runtime_error("Invalid units.");
	//misc
		bool error=false;
		
	try{
		//start the timer
		start=std::chrono::high_resolution_clock::now();
		
		//open the file
		if(DEBUG_VASP>1) std::cout<<"Opening file...\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Unable to open file.");
		
		//clear the simulation
		sim.clear();
		//set the periodicity (always true for VASP)
		sim.periodic()=true;
		
		/* read in the system name */
		if(DEBUG_VASP>1) std::cout<<"Reading in simulation name...\n";
		sim.system()=std::string(string::trim(fgets(input,string::M,reader)));
		
		/* load the simulation cell */
		if(DEBUG_VASP>1) std::cout<<"Loading simulation cell...\n";
		scale=std::atof(fgets(input, string::M, reader));
		if(scale<0) scale=std::pow(-scale,1.0/3.0);//negative -> volume
		/*
			Read in lattice vectors, where we take the transpose of the lattice vector matrix.
			VASP stores the matrix such that the lattice vectors are rows.
			This code store lattice vectors as columns to make coordinate transforms simpler.
		*/
		if(DEBUG_VASP>1) std::cout<<"Reading in the lattice vectors...\n";
		for(unsigned int i=0; i<sim.D; ++i){
			fgets(input, string::M, reader);
			lv(0,i)=std::atof(std::strtok(input,string::WS));
			for(unsigned int j=1; j<sim.D; ++j){
				lv(j,i)=std::atof(std::strtok(NULL,string::WS));
			}
		}
		//initialize the cell
		cell.init(s*lv,scale);
		
		/* load the atom info */
		//read in the number of species
		if(DEBUG_VASP>1) std::cout<<"Reading in the number of species...\n";
		fgets(input, string::M, reader);
		nSpecies=string::substrN(input,string::WS);
		atomNames.resize(nSpecies);
		nAtoms.resize(nSpecies);
		//read in the species names
		if(DEBUG_VASP>1) std::cout<<"Reading in the species names...\n";
		atomNames[0]=std::string(std::strtok(input,string::WS));
		for(unsigned int i=1; i<nSpecies; ++i){
			atomNames[i]=std::string(std::strtok(NULL,string::WS));
		}
		//read in the species numbers
		if(DEBUG_VASP>1) std::cout<<"Reading in the species numbers...\n";
		fgets(input, string::M, reader);
		nAtoms[0]=std::atoi(std::strtok(input,string::WS));
		for(unsigned int i=1; i<nSpecies; ++i){
			nAtoms[i]=std::atoi(std::strtok(NULL,string::WS));
		}
		//find the total number of atoms
		for(unsigned int n=0; n<nAtoms.size(); ++n) N+=nAtoms[n];
		
		/* check whether the system is direct or Cartesian*/
		fgets(input, string::M, reader);
		if(input[0]=='D') direct=true;
		else direct=false;
		
		/* check if the cell is variable or not */
		for(unsigned int n=0; n<N; ++n) fgets(input, string::M, reader);
		str=std::string(string::trim(fgets(input,string::M,reader)));
		if(str==sim.system()) sim.cellFixed()=false;
		else sim.cellFixed()=true;
		/* load the coordinate type */
		if(!sim.cellFixed()) for(unsigned int i=0; i<HEADER_SIZE; ++i) fgets(input,string::M,reader);
		
		/* find the number of timesteps */
		std::rewind(reader);
		//find the total number of lines in the file
		unsigned int nLines=0;
		for(unsigned int i=0; i<HEADER_SIZE; ++i) fgets(input, string::M, reader);
		while(fgets(input, string::M, reader)!=NULL){++nLines;};
		if(sim.cellFixed()) ts=std::floor((1.0*nLines)/(1.0*N+1.0));
		else ts=std::floor((1.0*nLines)/(1.0*N+1.0+HEADER_SIZE));
		
		//rewind
		std::rewind(reader);
		//skip the header
		if(sim.cellFixed()) for(unsigned int i=0; i<HEADER_SIZE; ++i) fgets(input,string::M,reader);
		
		//set the interval
		if(beg<1) throw std::invalid_argument("Invalid beginning timestep.");
		sim.beg()=beg-1;
		if(end<0){
			sim.end()=ts+end;
			end=sim.end();
		} else sim.end()=end-1;
		interval=sim.end()-sim.beg()+1;
		
		/* print data to screen */
		if(DEBUG_VASP>1){
			std::cout<<"NAME = "<<sim.system()<<"\n";
			std::cout<<"CELL = \n"<<cell<<"\n";
			std::cout<<"CELL_FIXED = "<<(sim.cellFixed()?"TRUE":"FALSE")<<"\n";
			std::cout<<"DIRECT = "<<(direct?"TRUE":"FALSE")<<"\n";
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
			std::cout<<"N_ATOMS = "<<N<<"\n";
			std::cout<<"TIMESTEPS = "<<ts<<"\n";
			std::cout<<"INTERVAL = ("<<sim.beg()<<","<<sim.end()<<") - "<<interval<<"\n";
			std::cout<<"STRIDE = "<<stride<<"\n";
			std::cout<<"N_STEPS = "<<interval/stride<<"\n";
		}
		
		/* resize the simulation */
		if(DEBUG_VASP>0) std::cout<<"Allocating memory...\n";
		sim.resize(interval/stride,nAtoms,atomNames);
		
		/* set the atom info */
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			for(unsigned int n=0; n<sim.nSpecies(); ++n){
				for(unsigned int m=0; m<sim.nAtoms(n); ++m){
					sim.atom(t,n,m).name()=sim.atomNames(n);
					sim.atom(t,n,m).specie()=n;
					sim.atom(t,n,m).index()=m;
					setAN(sim.atom(t,n,m),std::is_base_of<AN,AtomT>());
				}
			}
		}
		
		/* load positions */
		if(DEBUG_VASP>0) std::cout<<"Loading positions...\n";
		//skip timesteps until beg is reached
		for(unsigned int t=0; t<sim.beg(); ++t){
			if(!sim.cellFixed()) for(unsigned int i=0; i<HEADER_SIZE; ++i) fgets(input,string::M,reader); //skip header
			fgets(input,string::M,reader);//skip single line
			for(unsigned int n=0; n<N; ++n) fgets(input,string::M,reader);
		}
		//load the positions
		if(sim.cellFixed()){
			for(unsigned int t=0; t<sim.timesteps(); ++t) sim.cell(t)=cell;
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				if(DEBUG_VASP>2) std::cout<<"T = "<<t<<"\n";
				else if(t%1000==0) std::cout<<"T = "<<t<<"\n";
				fgets(input,string::M,reader);//skip line
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					fgets(input,string::M,reader);
					sim.atom_posn()(t,n)[0]=std::atof(std::strtok(input,string::WS));
					sim.atom_posn()(t,n)[1]=std::atof(std::strtok(NULL,string::WS));
					sim.atom_posn()(t,n)[2]=std::atof(std::strtok(NULL,string::WS));
				}
				//skip "stride-1" steps
				for(unsigned int tt=0; tt<stride-1; ++tt){
					fgets(input,string::M,reader);//skip line
					for(unsigned int n=0; n<sim.nAtoms(); ++n){
						fgets(input,string::M,reader);
					}
				}
			}
		} else {
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				if(DEBUG_VASP>2) std::cout<<"T = "<<t<<"\n";
				else if(t%1000==0) std::cout<<"T = "<<t<<"\n";
				//read in lattice vectors
				fgets(input,string::M,reader);//name
				scale=std::atof(fgets(input,string::M,reader));//scale
				for(unsigned int i=0; i<sim.D; ++i){
					fgets(input, string::M, reader);
					lv(0,i)=std::atof(std::strtok(input,string::WS));
					for(unsigned int j=1; j<sim.D; ++j){
						lv(j,i)=std::atof(std::strtok(NULL,string::WS));
					}
				}
				sim.cell(t).init(s*lv,scale);
				fgets(input,string::M,reader);//skip line (atom names)
				fgets(input,string::M,reader);//skip line (atom numbers)
				fgets(input,string::M,reader);//skip line (Direct or Cart)
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					fgets(input,string::M,reader);
					sim.atom_posn()(t,n)[0]=std::atof(std::strtok(input,string::WS));
					sim.atom_posn()(t,n)[1]=std::atof(std::strtok(NULL,string::WS));
					sim.atom_posn()(t,n)[2]=std::atof(std::strtok(NULL,string::WS));
				}
				//skip "stride-1" steps
				for(unsigned int tt=0; tt<stride-1; ++tt){
					fgets(input,string::M,reader);//name
					fgets(input,string::M,reader);//scale
					for(unsigned int i=0; i<sim.D; ++i){
						fgets(input,string::M,reader);//lv
					}
					fgets(input,string::M,reader);//skip line (atom names)
					fgets(input,string::M,reader);//skip line (atom numbers)
					fgets(input,string::M,reader);//skip line (Direct or Cart)
					for(unsigned int n=0; n<sim.nAtoms(); ++n){
						fgets(input,string::M,reader);
					}
				}
			}
		}
		
		/* convert to Cartesian coordinates if necessary */
		if(direct){
			if(DEBUG_VASP>0) std::cout<<"Converting to Cartesian coordinates...\n";
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					sim.atom_posn()(t,n)=sim.cell(t).R()*sim.atom_posn()(t,n);
				}
			}
		} else if(s!=1.0){
			for(unsigned int t=0; t<sim.timesteps(); ++t){
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					sim.atom_posn()(t,n)*=s;
				}
			}
		}
		
		fclose(reader);
		reader=NULL;
		
		//stop the timer
		stop=std::chrono::high_resolution_clock::now();
		
		//print the time
		time=std::chrono::duration_cast<std::chrono::duration<double> >(stop-start);
		std::cout<<"Simulation loaded in "<<time.count()<<" seconds.\n";
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free all local variables
	if(reader!=NULL) fclose(reader);
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

template <class AtomT>
void print(const char* file, const SimAtomic<AtomT>& sim, int beg=1, int end=-1, bool direct=true){
	static const char* funcName="print<AtomT>(const char*,SimAtomic<AtomT>&,int,int,bool)";
	if(DEBUG_VASP>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local function variables */
	//units
		double s=0.0;
		if(units::consts::system()==units::System::AU) s=units::ANGpBOHR;
		else if(units::consts::system()==units::System::METAL) s=1.0;
		else throw std::runtime_error("Invalid units.");
	//file i/o
		FILE* writer=NULL;
	//misc
		bool error=false;
	
	try{
		if(sim.timesteps()==0) throw std::invalid_argument("No simulation data found.");
		
		//set the timing info
		--beg;//make zero-indexed
		if(end<0) end=sim.timesteps()+end;
		else --end;//make zero-indexed
		
		//check timing info
		if(beg<0) throw std::invalid_argument("Invalid beginning timestep.");
		if(end>=sim.timesteps()) throw std::invalid_argument("Invalid ending timestep.");
		if(end<beg) throw std::invalid_argument("Invalid timestep interval.");
		
		//open the file
		if(DEBUG_VASP>1) std::cout<<"Opening file...\n";
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("I/O Error: Unable to open file.");
		
		if(sim.cellFixed()){
			//print the system name
			if(DEBUG_VASP>1) std::cout<<"Printing name...\n";
			fprintf(writer,"%s\n",sim.system().c_str());
			
			//print the simulation cell (transpose)
			if(DEBUG_VASP>1) std::cout<<"Printing cell...\n";
			fprintf(writer," %10.4f\n",sim.cell(beg).scale());
			for(unsigned int i=0; i<3; ++i){
				for(unsigned int j=0; j<3; ++j){
					fprintf(writer," %11.6f",sim.cell(0).R()(j,i)/sim.cell(0).scale());
				}
				fprintf(writer,"\n");
			}
			
			//print the atom info
			if(DEBUG_VASP>1) std::cout<<"Printing atom info...\n";
			for(unsigned int i=0; i<sim.nSpecies(); ++i){
				fprintf(writer," %4s",sim.atomNames(i).c_str());
			}
			fprintf(writer,"\n");
			for(unsigned int i=0; i<sim.nSpecies(); ++i){
				fprintf(writer," %4i",sim.nAtoms(i));
			}
			fprintf(writer,"\n");
		}
		
		//print the positions
		if(DEBUG_VASP>1) std::cout<<"Printing posns...\n";
		if(!direct){
			for(int t=beg; t<=end; ++t){
				if(!sim.cellFixed()){
					fprintf(writer,"%s\n",sim.system().c_str());
					fprintf(writer," %10.4f\n",sim.cell(t).scale());
					for(unsigned int i=0; i<3; ++i){
						for(unsigned int j=0; j<3; ++j){
							fprintf(writer," %11.6f",s*sim.cell(t).R()(j,i)/sim.cell(t).scale());
						}
						fprintf(writer,"\n");
					}
					for(unsigned int i=0; i<sim.nSpecies(); ++i) fprintf(writer," %4s",sim.atomNames(i).c_str());
					fprintf(writer,"\n");
					for(unsigned int i=0; i<sim.nSpecies(); ++i) fprintf(writer," %4i",sim.nAtoms(i));
					fprintf(writer,"\n");
				}
				fprintf(writer,"Cartesian configuration = %i\n", t+1);
				for(unsigned int n=0; n<sim.nSpecies(); ++n){
					for(unsigned int m=0; m<sim.nAtoms(n); ++m){
						fprintf(writer, "%.8f %.8f %.8f\n",
							sim.atom(t,n,m).posn()[0]*s,
							sim.atom(t,n,m).posn()[1]*s,
							sim.atom(t,n,m).posn()[2]*s
						);
					}
				}
			}
		} else {
			for(int t=beg; t<=end; ++t){
				if(!sim.cellFixed()){
					fprintf(writer,"%s\n",sim.system().c_str());
					fprintf(writer," %10.4f\n",sim.cell(t).scale());
					for(unsigned int i=0; i<3; ++i){
						for(unsigned int j=0; j<3; ++j){
							fprintf(writer," %11.6f",s*sim.cell(t).R()(j,i)/sim.cell(t).scale());
						}
						fprintf(writer,"\n");
					}
					for(unsigned int i=0; i<sim.nSpecies(); ++i) fprintf(writer," %4s",sim.atomNames(i).c_str());
					fprintf(writer,"\n");
					for(unsigned int i=0; i<sim.nSpecies(); ++i) fprintf(writer," %4i",sim.nAtoms(i));
					fprintf(writer,"\n");
				}
				fprintf(writer,"Direct configuration = %i\n", t+1);
				for(unsigned int n=0; n<sim.nSpecies(); ++n){
					for(unsigned int m=0; m<sim.nAtoms(n); ++m){
						Eigen::Vector3d posn=sim.cell(t).RInv()*sim.atom(t,n,m).posn();
						fprintf(writer, "%.8f %.8f %.8f\n",posn[0],posn[1],posn[2]);
					}
				}
			}
		}
		
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	if(writer!=NULL) fclose(writer);
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

}

namespace POSCAR{

//static variables
static const char* NAMESPACE_LOCAL="POSCAR";

template <class AtomT>
void load(const char* file, SimAtomic<AtomT>& sim){
	const char* funcName="load<AtomT>(const char*,SimAtomic<AtomT>&)";
	if(DEBUG_VASP>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<funcName<<":\n";
	/* local function variables */
	//file i/o
		FILE* reader=NULL;
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//simulation flags
		bool direct;//whether the coordinates are in direct or Cartesian coordinates
	//cell
		Cell cell;
		Eigen::Matrix3d lv;
		double scale=1;
	//atom info
		unsigned int nSpecies=0;//the number of atomic species
		std::vector<unsigned int> nAtoms;//the number of atoms in each species
		std::vector<std::string> atomNames;//the names of each species
		unsigned int N=0;
	//units
		double s=0.0;
		if(units::consts::system()==units::System::AU) s=units::BOHRpANG;
		else if(units::consts::system()==units::System::METAL) s=1.0;
		else throw std::runtime_error("Invalid units.");
	//misc
		bool error=false;
		
	try{
		//open the file
		if(DEBUG_VASP>1) std::cout<<"Opening file...\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("Unable to open file.");
		
		//clear the simulation
		sim.clear();
		//set the periodicity (always true for VASP)
		sim.periodic()=true;
		
		/* read in the system name */
		if(DEBUG_VASP>1) std::cout<<"Reading system name...\n";
		sim.system()=std::string(string::trim(fgets(input,string::M,reader)));
		
		/* load the simulation cell */
		if(DEBUG_VASP>1) std::cout<<"Loading simulation cell...\n";
		scale=std::atof(fgets(input, string::M, reader));
		if(scale<0) scale=std::pow(-1.0*scale,1.0/3.0);
		//read in the lattice vectors of the system
		/*
			Here we take the transpose of the lattice vector matrix.
			The reason is that VASP stores the lattice vector matrix such that the lattice vectors
			are the rows of the matrix.  However, it's more appropriate for operations on coordinates
			for the lattice vectors to be columns of the matrix.
		*/
		if(DEBUG_VASP>1) std::cout<<"Reading in the lattice vectors...\n";
		for(unsigned int i=0; i<sim.D; ++i){
			fgets(input, string::M, reader);
			lv(0,i)=std::atof(std::strtok(input,string::WS));
			for(unsigned int j=1; j<sim.D; ++j){
				lv(j,i)=std::atof(std::strtok(NULL,string::WS));
			}
		}
		//initialize the cell
		cell.init(s*lv,scale);
		
		/* load the atom info */
		//read in the number of species
		if(DEBUG_VASP>1) std::cout<<"Reading in the number of species...\n";
		fgets(input, string::M, reader);
		nSpecies=string::substrN(input,string::WS);
		atomNames.resize(nSpecies);
		nAtoms.resize(nSpecies);
		//read in the species names
		if(DEBUG_VASP>1) std::cout<<"Reading in the species names...\n";
		atomNames[0]=std::string(std::strtok(input,string::WS));
		for(unsigned int i=1; i<nSpecies; ++i){
			atomNames[i]=std::string(std::strtok(NULL,string::WS));
		}
		//read in the species numbers
		if(DEBUG_VASP>1) std::cout<<"Reading in the species numbers...\n";
		fgets(input, string::M, reader);
		nAtoms[0]=std::atoi(std::strtok(input,string::WS));
		for(unsigned int i=1; i<nSpecies; ++i){
			nAtoms[i]=std::atoi(std::strtok(NULL,string::WS));
		}
		
		/* load the coordinate type */
		fgets(input, string::M, reader);
		if(input[0]=='D') direct=true;
		else direct=false;
		
		/* print data to screen */
		if(DEBUG_VASP>1){
			std::cout<<"NAME = "<<sim.system()<<"\n";
			std::cout<<"DIRECT = "<<(direct?"T":"F")<<"\n";
			std::cout<<"CELL = \n"<<cell<<"\n";
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
		}
		
		/* resize the simulation */
		if(DEBUG_VASP>0) std::cout<<"Allocating memory...\n";
		sim.resize(1,nAtoms,atomNames);
		sim.cell(0)=cell;
		
		/* load positions */
		if(DEBUG_VASP>1) std::cout<<"Loading positions...\n";
		//load the positions
		for(unsigned int n=0; n<sim.nAtoms(); ++n){
			fgets(input,string::M,reader);
			sim.atom_posn()(0,n)[0]=s*std::atof(std::strtok(input,string::WS));
			sim.atom_posn()(0,n)[1]=s*std::atof(std::strtok(NULL,string::WS));
			sim.atom_posn()(0,n)[2]=s*std::atof(std::strtok(NULL,string::WS));
		}
		
		/* convert to cartesian coordinates (if necessary) */
		if(DEBUG_VASP>1) std::cout<<"Converting to Cartesian coordinates...\n";
		if(direct){
			for(unsigned int n=0; n<sim.nSpecies(); ++n){
				sim.atom_posn()(0,n)=sim.cell(0).R()*sim.atom_posn()(0,n);
			}
		}
		
		fclose(reader);
		reader=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free all local variables
	if(reader!=NULL) fclose(reader);
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

template <class AtomT>
void print(const char* file, const SimAtomic<AtomT>& sim, unsigned int t=1, bool direct=true){
	static const char* funcName="load<AtomT>(const char*,SimAtomic<AtomT>&)";
	if(DEBUG_VASP>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local function variables */
	//file i/o
		FILE* writer=NULL;
	//misc
		bool error=false;
	
	try{
		--t;//set to be zero-indexed
		if(sim.timesteps()==0) throw std::invalid_argument("No simulation data found.");
		if(t>=sim.timesteps()) throw std::invalid_argument("Invalid timestep.");
		
		//open the file
		if(DEBUG_VASP>1) std::cout<<"Opening file...\n";
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("Unable to open file.");
		
		//print the system name
		if(DEBUG_VASP>1) std::cout<<"Printing name...\n";
		fprintf(writer,"%s\n",sim.system().c_str());
		
		//print the simulation cell (transpose)
		if(DEBUG_VASP>1) std::cout<<"Printing cell...\n";
		fprintf(writer," %10.4f\n",sim.cell(0).scale());
		for(unsigned int i=0; i<3; ++i){
			for(unsigned int j=0; j<3; ++j){
				fprintf(writer," %11.6f",sim.cell(0).R()(j,i)/sim.cell(0).scale());
			}
			fprintf(writer,"\n");
		}
		
		//print the atom info
		if(DEBUG_VASP>1) std::cout<<"Printing atom info...\n";
		for(unsigned int i=0; i<sim.nSpecies(); ++i){
			fprintf(writer," %4s",sim.atomNames(i).c_str());
		}
		fprintf(writer,"\n");
		for(unsigned int i=0; i<sim.nSpecies(); ++i){
			fprintf(writer," %4i",sim.nAtoms(i));
		}
		fprintf(writer,"\n");
		
		//print the positions
		if(DEBUG_VASP>1) std::cout<<"Printing posns...\n";
		if(!direct){
			fprintf(writer,"Cart\n");
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				fprintf(writer, "%.8f %.8f %.8f\n",
					sim.atom_posn()(t,n)[0],
					sim.atom_posn()(t,n)[1],
					sim.atom_posn()(t,n)[2]
				);
			}
		} else {
			fprintf(writer,"Direct\n");
			for(unsigned int n=0; n<sim.nAtoms(); ++n){
				Eigen::Vector3d posn=sim.cell(t).RInv()*sim.atom_posn()(t,n);
				fprintf(writer, "%.8f %.8f %.8f\n",posn[0],posn[1],posn[2]);
			}
		}
		
		fclose(writer);
		writer=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	if(writer!=NULL) fclose(writer);
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

}

namespace OUTCAR{
	
//static variables
static const char* NAMESPACE_LOCAL="OUTCAR";

struct Band{
	unsigned int nBands;
	std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > k;
	std::vector<Eigen::VectorXd,Eigen::aligned_allocator<Eigen::VectorXd> > energy;
	std::vector<Eigen::VectorXi,Eigen::aligned_allocator<Eigen::VectorXi> > occ;
};

Band& load_band(const char* file, Band& band);
void print_energy(const char* file, const Band& band);

}

namespace EIGENVAL{
	
}

namespace ENERGY{
		
	//static variables
	static const char* NAMESPACE_LOCAL="ENERGY";
	
	void load(const char* file, SimI& sim);
}

namespace XML{

//static variables
static const char* NAMESPACE_LOCAL="XML";

template <class AtomT> void load_forces(std::false_type, FILE* reader, SimAtomic<AtomT>& sim){};
template <class AtomT> void load_forces(std::true_type, FILE* reader, SimAtomic<AtomT>& sim){
	static const char* funcName="load_forces(std::true_type,FILE*,SimAtomic<AtomT>&)";
	if(DEBUG_VASP>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	const char* force_str="<varray name=\"forces\" >";
	char* input=(char*)malloc(sizeof(char)*string::M);
	char* temp=(char*)malloc(sizeof(char)*string::M);
		
	rewind(reader);
	/* read in the forces */
	if(DEBUG_VASP>0) std::cout<<"Reading in forces...\n";
	std::rewind(reader);
	unsigned tt=0;
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strstr(input,"<calculation>")!=NULL){
			while(fgets(input,string::M,reader)!=NULL){
				if(std::strstr(input,"</calculation>")!=NULL){
					throw std::invalid_argument("No forces provided in calculation.");
				} else if(std::strstr(input,"forces")!=NULL){
					for(unsigned int n=0; n<sim.nAtoms(); ++n){
						fgets(input,string::M,reader); std::strtok(input,string::WS);
						sim.atom_force()(tt,n)[0]=std::atof(std::strtok(NULL,string::WS));
						sim.atom_force()(tt,n)[1]=std::atof(std::strtok(NULL,string::WS));
						sim.atom_force()(tt,n)[2]=std::atof(std::strtok(NULL,string::WS));
					}
					++tt;
					break;
				}
			}
		}
	}
	
	free(temp);
	free(input);
}

template <class AtomT>
void load(const char* file, SimAtomic<AtomT>& sim, int beg=1, int end=-1, int stride=1){
	static const char* funcName="load<AtomT>(const char*,SimAtomic<AtomT>&,int,int)";
	if(DEBUG_VASP>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local function variables */
	//file i/o
		FILE* reader=NULL;
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
		std::string str;
	//simulation flags
		bool direct;//whether the coordinates are in direct or Cartesian coordinates
	//time info
		unsigned int ts=0;//number of timesteps
		unsigned int interval=0;//requested interval
	//cell info
		double scale=1;
		Eigen::Matrix3d lv;
		Cell cell;
	//atom info
		unsigned int nSpecies=0;//the number of atomic species
		std::vector<unsigned int> nAtoms;//the number of atoms in each species
		std::vector<std::string> atomNames;//the names of each species
		unsigned int N=0;
	//timing
		std::chrono::high_resolution_clock::time_point start;
		std::chrono::high_resolution_clock::time_point stop;
		std::chrono::duration<double> time;
	//misc
		bool error=false;
		unsigned int tt;
	
	try{
		
		/* open the file */
		if(DEBUG_VASP>0) std::cout<<"Opening the file...\n";
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open file.");
		
		/* load the number of timesteps */
		if(DEBUG_VASP>0) std::cout<<"Loading the timesteps...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"<calculation>")!=NULL) ++ts;
		}
		if(DEBUG_VASP>0) std::cout<<"ts = "<<ts<<"\n";
		
		/* read in the atom info */
		if(DEBUG_VASP>0) std::cout<<"Reading in the atom info...\n";
		std::rewind(reader);
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"atominfo")!=NULL){
				fgets(input,string::M,reader);
				std::strtok(input,">");
				N=std::atoi(std::strtok(NULL,"<"));
				for(unsigned int i=0; i<6; ++i) fgets(input,string::M,reader);
				for(unsigned int i=0; i<N; ++i){
					fgets(input,string::M,reader);
					std::strtok(input,">"); std::strtok(NULL,">");
					std::string name=std::strtok(NULL,"<");
					bool match=false;
					for(unsigned int j=0; j<atomNames.size(); ++j){
						if(name==atomNames[j]){
							++nAtoms[j];
							match=true; 
							break;
						}
					}
					if(!match){
						atomNames.push_back(name);
						nAtoms.push_back(1);
					}
				}
				break;
			}
		}
		if(DEBUG_VASP>0){
			std::cout<<"ATOM_NAMES = "; for(unsigned int i=0; i<atomNames.size(); ++i) std::cout<<atomNames[i]<<" "; std::cout<<"\n";
			std::cout<<"ATOM_NUMBERS = "; for(unsigned int i=0; i<nAtoms.size(); ++i) std::cout<<nAtoms[i]<<" "; std::cout<<"\n";
		}
		if(atomNames.size()!=nAtoms.size()) throw std::runtime_error("Mismatch in atom names/numbers.");
		nSpecies=atomNames.size();
		
		/* resize the simulation */
		if(DEBUG_VASP>0) std::cout<<"Resizing the simulation...\n";
		sim.resize(ts,nAtoms,atomNames);
		
		/* load the cells */
		if(DEBUG_VASP>0) std::cout<<"Loading the cells...\n";
		tt=0;
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"<calculation>")!=NULL){
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strstr(input,"</calculation>")!=NULL){
						throw std::invalid_argument("No cell provided in calculation.");
					} else if(std::strstr(input,"basis")!=NULL){
						fgets(input,string::M,reader);
						std::strtok(input,string::WS);
						lv(0,0)=std::atof(std::strtok(NULL,string::WS));
						lv(1,0)=std::atof(std::strtok(NULL,string::WS));
						lv(2,0)=std::atof(std::strtok(NULL,string::WS));
						fgets(input,string::M,reader);
						std::strtok(input,string::WS);
						lv(0,1)=std::atof(std::strtok(NULL,string::WS));
						lv(1,1)=std::atof(std::strtok(NULL,string::WS));
						lv(2,1)=std::atof(std::strtok(NULL,string::WS));
						fgets(input,string::M,reader);
						std::strtok(input,string::WS);
						lv(0,2)=std::atof(std::strtok(NULL,string::WS));
						lv(1,2)=std::atof(std::strtok(NULL,string::WS));
						lv(2,2)=std::atof(std::strtok(NULL,string::WS));
						sim.cell(tt++).init(lv);
						break;
					}
				}
			}
		}
		
		/* read in the positions */
		if(DEBUG_VASP>0) std::cout<<"Reading in positions...\n";
		std::rewind(reader);
		tt=0;
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"<calculation>")!=NULL){
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strstr(input,"</calculation>")!=NULL){
						throw std::invalid_argument("No positions provided in calculation.");
					} else if(std::strstr(input,"positions")!=NULL){
						for(unsigned int n=0; n<sim.nAtoms(); ++n){
							fgets(input,string::M,reader); std::strtok(input,string::WS);
							sim.atom_posn()(tt,n)[0]=std::atof(std::strtok(NULL,string::WS));
							sim.atom_posn()(tt,n)[1]=std::atof(std::strtok(NULL,string::WS));
							sim.atom_posn()(tt,n)[2]=std::atof(std::strtok(NULL,string::WS));
						}
						++tt;
						break;
					}
				}
			}
		}
		
		/* read in the energies */
		if(DEBUG_VASP>0) std::cout<<"Reading in the energies...\n";
		std::rewind(reader);
		tt=0;
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"<calculation>")!=NULL){
				double energy=0;
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strstr(input,"</calculation>")!=NULL) break;
					if(std::strstr(input,"e_fr_energy")!=NULL){
						std::strtok(input,">"); energy=std::atof(std::strtok(NULL,"<"));
					}
				}
				sim.energy(tt++)=energy;
			}
		}
		
		//read in the forces
		if(DEBUG_VASP>0) std::cout<<"Reading in the forces...\n";
		load_forces(std::is_base_of<Force,AtomT>(),reader,sim);
	}catch(std::exception& e){
		std::cout<<"ERROR in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	/* transform to Cartesian coordinates */
	if(DEBUG_VASP>0) std::cout<<"Transforming to Cartesian coordinates...\n";
	for(unsigned int t=0; t<sim.timesteps(); ++t){
		for(unsigned int n=0; n<sim.nAtoms(); ++n){
			sim.atom_posn()(t,n)=sim.cell(t).R()*sim.atom_posn()(t,n);
			Cell::returnToCell(sim.atom_posn()(t,n),sim.atom_posn()(t,n),sim.cell(t).R(),sim.cell(t).RInv());
		}
	}
	
	//free all local variables
	if(reader!=NULL) fclose(reader);
	free(input);
	free(temp);
	
	if(error) throw std::runtime_error("I/O Exception Occurred.");
}

}

}

#endif
