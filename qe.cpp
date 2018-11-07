#include "qe.hpp"

namespace QE{

//*****************************************************
//FORMAT struct
//*****************************************************

Format& Format::load(const std::vector<std::string>& strlist, Format& format){
	for(unsigned int i=0; i<strlist.size(); ++i){
		if(strlist[i]=="-pos"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-pos\" option.");
			else format.filePos=strlist[i+1];
		} else if(strlist[i]=="-in"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-in\" option.");
			else format.fileIn=strlist[i+1];
		} else if(strlist[i]=="-cel"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-cel\" option.");
			else format.fileCel=strlist[i+1];
		} else if(strlist[i]=="-evp"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-evp\" option.");
			else format.fileEvp=strlist[i+1];
		} else if(strlist[i]=="-interval"){
			if(i==strlist.size()-1) throw std::invalid_argument("No interval specified for \"-int\" option.");
			else if(string::substrN(strlist[i+1].c_str(),":")<2) throw std::invalid_argument("Interval must have \":\" separating beginning and end.");
			else if(string::substrN(strlist[i+1].c_str(),":")==2){
				char* temp=(char*)malloc(sizeof(char)*250);
				std::strcpy(temp,strlist[i+1].c_str());
				format.beg=std::atoi(std::strtok(temp,":"));
				format.end=std::atoi(std::strtok(NULL,":"));
				free(temp); temp=NULL;
			} else if(string::substrN(strlist[i+1].c_str(),":")==3){
				char* temp=(char*)malloc(sizeof(char)*250);
				std::strcpy(temp,strlist[i+1].c_str());
				format.beg=std::atoi(std::strtok(temp,":"));
				format.end=std::atoi(std::strtok(NULL,":"));
				format.stride=std::atoi(std::strtok(NULL,":"));
				free(temp); temp=NULL;
			} else throw std::invalid_argument("Interval has too many \":\"");
		}
	}
	return format;
}

//*****************************************************
//CEL format
//*****************************************************

namespace CEL{

void load_cell(FILE* reader, SimI& sim){
	const char* func_name="load_cell(const char*,SimI&)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
	//cell
		Eigen::Matrix3d lv;
		double s=0.52917721067;
	//misc
		bool error=false;
		
	try{
		//rewind the reader
		std::rewind(reader);
		
		//skip to the beginning
		for(unsigned int t=0; t<sim.beg(); ++t){
			if(DEBUG_QE>1) std::cout<<"t = "<<t<<"\n";
			fgets(input,string::M,reader);//header
			for(unsigned int n=0; n<3; ++n) fgets(input,string::M,reader);
		}
		//read in the data
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			if(DEBUG_QE>1) std::cout<<"t = "<<t<<"\n";
			fgets(input,string::M,reader);//header
			//first line
			fgets(input,string::M,reader);
			lv(0,0)=std::atof(std::strtok(input,string::WS));
			lv(0,1)=std::atof(std::strtok(NULL,string::WS));
			lv(0,2)=std::atof(std::strtok(NULL,string::WS));
			//second line
			fgets(input,string::M,reader);
			lv(1,0)=std::atof(std::strtok(input,string::WS));
			lv(1,1)=std::atof(std::strtok(NULL,string::WS));
			lv(1,2)=std::atof(std::strtok(NULL,string::WS));
			//third line
			fgets(input,string::M,reader);
			lv(2,0)=std::atof(std::strtok(input,string::WS));
			lv(2,1)=std::atof(std::strtok(NULL,string::WS));
			lv(2,2)=std::atof(std::strtok(NULL,string::WS));
			//set the cell
			lv*=s;
			sim.cell(t).init(lv);
			for(unsigned int tt=1; tt<sim.stride(); ++tt){
				fgets(input,string::M,reader);//header
				for(unsigned int n=0; n<3; ++n) fgets(input,string::M,reader);
			}
		}
		
		lv=sim.cell(0).R();
		for(unsigned int t=0; t<sim.timesteps(); ++t){
			if((lv-sim.cell(t).R()).norm()>num_const::ZERO){
				sim.cellFixed()=false; break;
			}
		}
	}catch(std::exception& e){
		std::cout<<"Error in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
}
	
}

//*****************************************************
//IN format
//*****************************************************

namespace IN{

Cell& load_cell(FILE* reader, Cell& cell){
	const char* func_name="load_cell(const char*,SimAtomic<AtomT>&)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* tag=(char*)malloc(sizeof(char)*string::M);
		char* option=(char*)malloc(sizeof(char)*string::M);
	//cell
		int ibrav=0;
		double s=1.0;
		std::vector<double> lvp(6,0);
		Eigen::Matrix3d lv;
	//misc
		bool error=false;
		
	try{
		//rewind the reader
		std::rewind(reader);
		//read in ibrav
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_all(input);//remove all spaces
			string::copy_left(tag,input,"=");//copy the tag
			if(std::strcmp(tag,"ibrav")==0){
				ibrav=std::atoi(std::strpbrk(input,"=")+1);
				break;
			}
		}
		if(DEBUG_QE>1) std::cout<<"IBRAV = "<<ibrav<<"\n";
		
		//rewind the reader
		std::rewind(reader);
		if(ibrav==0){
			//read in the cell parameters
			while(fgets(input,string::M,reader)!=NULL){
				string::trim(input);//remove beg/end spaces
				string::copy_left(tag,input,string::WS);//copy the tag
				if(std::strpbrk(input,string::WS)!=NULL){//if there's an option...
					std::strtok(input,string::WS);
					string::trim_all(string::to_upper(std::strcpy(option,std::strtok(NULL,string::WS))));//copy the option
				}
				if(std::strcmp(tag,"CELL_PARAMETERS")==0){
					//read in the unit, convert to angstrom in necessary
					if(!string::empty(option)){
						if(std::strcmp(option,"ANGSTROM")==0) s=1;
						else if(std::strcmp(option,"BOHR")==0) s=0.52917721067;
					}
					if(DEBUG_QE>1) std::cout<<"option = "<<option<<"\n";
					//first lattice vector
					fgets(input,string::M,reader);
					lv(0,0)=s*std::atof(std::strtok(input,string::WS));
					lv(1,0)=s*std::atof(std::strtok(NULL,string::WS));
					lv(2,0)=s*std::atof(std::strtok(NULL,string::WS));
					//second lattice vector
					fgets(input,string::M,reader);
					lv(0,1)=s*std::atof(std::strtok(input,string::WS));
					lv(1,1)=s*std::atof(std::strtok(NULL,string::WS));
					lv(2,1)=s*std::atof(std::strtok(NULL,string::WS));
					//third lattice vector
					fgets(input,string::M,reader);
					lv(0,2)=s*std::atof(std::strtok(input,string::WS));
					lv(1,2)=s*std::atof(std::strtok(NULL,string::WS));
					lv(2,2)=s*std::atof(std::strtok(NULL,string::WS));
					//break
					break;
				}
			}
		} else {
			//read in the cell parameters
			while(fgets(input,string::M,reader)!=NULL){
				string::trim_all(input);//remove all spaces
				string::copy_left(tag,input,"(");//copy the tag
				if(std::strcmp(tag,"celldm")==0){
					std::strtok(input," \t\n()=");
					int index=std::atoi(std::strtok(NULL," \t\n()="));
					double val=std::atof(std::strtok(NULL," \t\n()="));
					lvp[index]=val;
				}
			}
			if(ibrav==0){
				lv=Eigen::Matrix3d::Identity()*lvp[0];
			} else if(ibrav==1){
				lv(0,0)=-1; lv(1,0)=0; lv(2,0)=1;
				lv(0,1)=0; lv(1,1)=1; lv(2,1)=1;
				lv(0,2)=-1; lv(1,2)=1; lv(2,2)=0;
				lv*=lvp[0]*0.5;
			} else if(ibrav==2){
				lv(0,0)=1; lv(1,0)=1; lv(2,0)=1;
				lv(0,1)=-1; lv(1,1)=1; lv(2,1)=1;
				lv(0,2)=-1; lv(1,2)=-1; lv(2,2)=1;
				lv*=lvp[0]*0.5;
			} else if(ibrav==4){
				lv(0,0)=1; lv(1,0)=0; lv(2,0)=0;
				lv(0,1)=-0.5; lv(1,1)=0.5*std::sqrt(3.0); lv(2,1)=0;
				lv(0,2)=0; lv(1,2)=0; lv(2,2)=lvp[2];
				lv*=lvp[0];
			} else if(ibrav==6){
				lv=Eigen::Matrix3d::Identity();
				lv(2,2)=lvp[2];
				lv*=lvp[0];
			} else if(ibrav==8){
				lv=Eigen::Matrix3d::Identity();
				lv(0,0)=lvp[0];
				lv(1,1)=lvp[1];
				lv(2,2)=lvp[2];
			} else if(ibrav==12){
				lv=Eigen::Matrix3d::Identity();
				lv(0,0)=lvp[0];
				lv(0,1)=lvp[0]*lvp[1]*std::cos(lvp[3]);
				lv(1,1)=lvp[0]*lvp[1]*std::cos(lvp[3]);
				lv(2,2)=lvp[0]*lvp[2];
			} else throw std::invalid_argument("Invalid IBRAV.");
			lv*=0.529;//convert to angstroms
		}
		
		if(DEBUG_QE>1) std::cout<<"lv = \n"<<lv<<"\n";
		
		//initialize the cell
		cell.init(lv);
	}catch(std::exception& e){
		std::cout<<"Error in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	//free local variables
	free(input);
	free(tag);
	free(option);
	
	if(error) throw std::runtime_error("I/O Error: Failed to load.");
	else return cell;
}

void load_atoms(FILE* reader, std::vector<std::string>& atomNames, std::vector<unsigned int>& atomNumbers){
	const char* func_name="load_atoms(FILE*,std::vector<std::string>&,std::vector<std::string>&)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
		char* tag=(char*)malloc(sizeof(char)*string::M);
		char* option=(char*)malloc(sizeof(char)*string::M);
	//atoms
		unsigned int nAtoms=0;
		unsigned int nTypes=0;
	//misc
		bool error=false;
		
	try{
		//rewind the reader
		std::rewind(reader);
		
		//read in the atom info
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_all(string::trim_right(input,","));
			string::copy_left(temp,input,"=");
			if(std::strcmp(temp,"nat")==0){
				nAtoms=std::atoi(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"ntyp")==0){
				nTypes=std::atoi(std::strpbrk(input,"=")+1);
			}
		}
		
		//rewind the reader
		std::rewind(reader);
		
		//read in the atom types
		atomNames.clear();
		atomNumbers.clear();
		while(fgets(input,string::M,reader)!=NULL){
			string::trim(input);//remove all spaces
			string::copy_left(tag,input,string::WS);//copy the tag
			if(std::strcmp(tag,"ATOMIC_POSITIONS")==0){
				for(unsigned int n=0; n<nAtoms; ++n){
					fgets(input,string::M,reader);
					if(string::substrN(input,string::WS)!=4) throw std::runtime_error("Invalid atom format.");
					//read in the name
					std::string name=std::string(std::strtok(input,string::WS));
					bool match=false;
					for(unsigned int i=0; i<atomNames.size(); ++i){
						if(name==atomNames[i]){
							match=true;
							++atomNumbers[i];
						}
					}
					if(!match){
						atomNames.push_back(name);
						atomNumbers.push_back(1);
					}
				}
				break;
			}
		}
		
		unsigned int nAtomsList=0;
		for(unsigned int i=0; i<atomNumbers.size(); ++i) nAtomsList+=atomNumbers[i];
		
		if(DEBUG_QE>0){
			std::cout<<"nTypes = "<<nTypes<<"\n";
			std::cout<<"nAtomNames = "<<atomNames.size()<<"\n";
			std::cout<<"nAtoms = "<<nAtoms<<"\n";
			std::cout<<"nAtomsList = "<<nAtomsList<<"\n";
			std::cout<<"atomNames = \n";
			for(unsigned int i=0; i<atomNames.size(); ++i){
				std::cout<<"\t"<<atomNames[i]<<"\n";
			}
		}
		
		if(atomNames.size()!=nTypes) throw std::runtime_error("Incompatible number of types and atom specification.");
		if(nAtoms!=nAtomsList) throw std::runtime_error("Incompatible number of atoms and atom specification.");
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

double load_timestep(FILE* reader){
	const char* func_name="load_timestep(FILE*)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//timestep
		double ts;
		unsigned int print=0;
	//misc
		bool error=false;
		
	try{
		//rewind the reader
		std::rewind(reader);
		
		//read in the cell
		if(DEBUG_QE>0) std::cout<<"Reading in the timestep...\n";
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_all(string::trim_right(input,","));
			string::copy_left(temp,input,"=");
			if(std::strcmp(temp,"dt")==0){
				ts=std::atof(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"iprint")==0){
				print=std::atoi(std::strpbrk(input,"=")+1);
			}
		}
		
		if(DEBUG_QE>0){
			std::cout<<"ts = "<<ts<<"\n";
			std::cout<<"print = "<<print<<"\n";
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
	else return ts*print*0.02418884326505;//convert from a.u. to fs
}

}

//*****************************************************
//POS format
//*****************************************************

namespace POS{

unsigned int load_timesteps(FILE* reader){
	const char* func_name="load_timesteps(const char*)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
	//timesteps
		unsigned int ts=0;
	//misc
		bool error=false;
		
	try{
		//rewind the reader
		std::rewind(reader);
		
		//read in the timesteps
		while(fgets(input,string::M,reader)!=NULL){
			if(string::substrN(input,string::WS)==2) ++ts;
		}
		
		if(DEBUG_QE>0){
			std::cout<<"ts = "<<ts<<"\n";
		}
	}catch(std::exception& e){
		std::cout<<"Error in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	free(input);
	
	if(error) throw std::runtime_error("I/O Error: Failed to load.");
	return ts;
}
	
}

namespace EVP{
	
void load_energy(FILE* reader, SimI& sim){
	const char* func_name="load_energy(FILE*,SimAtomic<AtomT>&)";
	if(DEBUG_QE>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
	//timesteps
		unsigned int ts=0;
	//misc
		bool error=false;
		
	try{
		//rewind the reader
		std::rewind(reader);
		
		//read in the timesteps
		fgets(input,string::M,reader);//skip the header
		//skip to beginning
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strpbrk(input,"#")!=NULL) continue;//skip any lines with "#"
			if(ts==sim.beg()){
				std::strtok(input,string::WS);
				std::strtok(NULL,string::WS);
				std::strtok(NULL,string::WS);
				std::strtok(NULL,string::WS);
				std::strtok(NULL,string::WS);
				sim.energy(ts++)=std::atof(std::strtok(NULL,string::WS))*13.605693009*2;
				for(unsigned int tt=1; tt<sim.stride(); ++tt) fgets(input,string::M,reader);
				break;
			} else ++ts;
		}
		ts=1;
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strpbrk(input,"#")!=NULL) continue;//skip any lines with "#"
			if(string::empty(input)) continue;//skip any empty lines
			if(ts==sim.timesteps()) break;
			std::strtok(input,string::WS);
			std::strtok(NULL,string::WS);
			std::strtok(NULL,string::WS);
			std::strtok(NULL,string::WS);
			std::strtok(NULL,string::WS);
			sim.energy(ts++)=std::atof(std::strtok(NULL,string::WS))*13.605693009*2;
			for(unsigned int tt=1; tt<sim.stride(); ++tt) fgets(input,string::M,reader);
		}
	}catch(std::exception& e){
		std::cout<<"Error in "<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<func_name<<":\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	free(input);
	
	if(error) throw std::runtime_error("I/O Error: Failed to load.");
}

}

}