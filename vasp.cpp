#include "vasp.hpp"

namespace VASP{

//*****************************************************
//FORMAT struct
//*****************************************************

Format& Format::load(const std::vector<std::string>& strlist, Format& format){
	for(unsigned int i=0; i<strlist.size(); ++i){
		if(strlist[i]=="-xdatcar"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-xdatcar\" option.");
			else format.xdatcar=strlist[i+1];
		} else if(strlist[i]=="-poscar"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-poscar\" option.");
			else format.poscar=strlist[i+1];
		} else if(strlist[i]=="-xml"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-poscar\" option.");
			else format.xml=strlist[i+1];
		} else if(strlist[i]=="-energy"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-energy\" option.");
			else format.energy=strlist[i+1];
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
	
namespace OUTCAR{

Band& load_band(const char* file, Band& band){
	const char* funcName="load_band(const char*,Band&)";
	if(DEBUG_VASP>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local variables*/
	//file i/o
		FILE* reader=NULL;
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//misc
		bool error=false;
		
	try{
		//open the file
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open OUTCAR file.");
		
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strstr(input,"NBANDS")!=NULL){
				std::strtok(input,"=");
				std::strtok(NULL,"=");
				std::strtok(NULL,"=");
				band.nBands=std::atoi(std::strtok(NULL,"="));
			}
			if(std::strstr(input,"E-fermi")!=NULL){
				while(fgets(input,string::M,reader)!=NULL){
					if(std::strstr(input,"k-point")!=NULL){
						//read in the k-point
						std::strtok(input,string::WS);
						std::strtok(NULL,string::WS);
						std::strtok(NULL,string::WS);
						band.k.push_back(Eigen::Vector3d::Zero());
						band.k.back()[0]=std::atof(std::strtok(NULL,string::WS));
						band.k.back()[1]=std::atof(std::strtok(NULL,string::WS));
						band.k.back()[2]=std::atof(std::strtok(NULL,string::WS));
						//read in the band energy and occupation
						band.energy.push_back(Eigen::VectorXd::Constant(band.nBands,1,0));
						band.occ.push_back(Eigen::VectorXi::Constant(band.nBands,1,0));
						fgets(input,string::M,reader);
						for(unsigned int i=0; i<band.nBands; ++i){
							fgets(input,string::M,reader);
							std::strtok(input,string::WS);
							band.energy.back()[i]=std::atof(std::strtok(NULL,string::WS));
							band.occ.back()[i]=std::atoi(std::strtok(NULL,string::WS));
						}
					}
				}
				break;
			}
		}
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
	else return band;
}

void print_energy(const char* file, const Band& band){
	const char* funcName="print_energy(const char*,const Band&)";
	if(DEBUG_VASP>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local variables*/
	//file i/o
		FILE* writer=NULL;
	//misc
		bool error=false;
		
	try{
		//open the file
		writer=fopen(file,"w");
		if(writer==NULL) throw std::runtime_error("I/O Error: Could not open OUTCAR file.");
		
		fprintf(writer, "KX KY KZ ");
		for(unsigned int n=0; n<band.nBands; ++n) fprintf(writer, "N%i ", n);
		fprintf(writer, "\n");
		
		for(unsigned int n=0; n<band.k.size(); ++n){
			fprintf(writer, "%f %f %f ",band.k[n][0],band.k[n][1],band.k[n][2]);
			for(unsigned int i=0; i<band.nBands; ++i){
				fprintf(writer, "%f ",band.energy[n][i]);
			}
			fprintf(writer, "\n");
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

namespace EIGENVAL{
	
}

namespace ENERGY{
		
	void load(const char* file, SimI& sim){
		//local variables
		FILE* reader=NULL;
		char* input=(char*)malloc(sizeof(char)*string::M);
		bool error=false;
		
		try{
			//open the file
			reader=fopen(file,"r");
			if(reader==NULL) throw std::runtime_error("Could not open file.");
			
			//read in the energy
			unsigned int t=0;
			while(fgets(input,string::M,reader)!=NULL){
				if(string::empty(input)) continue;
				if(t==sim.timesteps()) break;
				sim.energy(t++)=std::atof(input);
				//sim.energy(t++)=std::atof(input)*23.061;
				for(unsigned int tt=1; tt<sim.stride(); ++tt) fgets(input,string::M,reader);
			}
		}catch(std::exception& e){
			std::cout<<"ERROR in "<<NAMESPACE_LOCAL<<"::"<<"energy(const char*,SimI&):\n";
			std::cout<<e.what()<<"\n";
			error=true;
		}
		
		free(input);
		
		if(error) throw std::runtime_error("I/O Error.");
	}
}

}
