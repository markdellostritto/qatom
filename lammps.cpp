#include "lammps.hpp"

namespace LAMMPS{

//*****************************************************
//FORMAT struct
//*****************************************************

Format& Format::load(const std::vector<std::string>& strlist, Format& format){
	for(unsigned int i=0; i<strlist.size(); ++i){
		if(strlist[i]=="-in"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-in\" option.");
			else format.in=strlist[i+1];
		} else if(strlist[i]=="-data"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-data\" option.");
			else format.data=strlist[i+1];
		} else if(strlist[i]=="-dump"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-dump\" option.");
			else format.dump=strlist[i+1];
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
	
namespace IN{

void load_style(const char* file, STYLE_ATOM::type& styleAtom){
	static const char* funcName="load_style(const char*,STYLE_ATOM::type&)";
	if(DEBUG_LAMMPS>0) std::cout<<NAMESPACE_GLOBAL<<"::"<<NAMESPACE_LOCAL<<"::"<<funcName<<":\n";
	/* local function variables */
	//file i/o
		FILE* reader;
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//miscellaneous
		bool error=false;
	
	try{
		/* open the file */
		reader=fopen(file,"r");
		if(reader==NULL) throw std::runtime_error("I/O Error: Could not open file.");
		
		/* read in atom style */
		if(DEBUG_LAMMPS>1) std::cout<<"Reading in atom style...\n";
		while(fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,"#");
			if(string::empty(input)) continue;
			//break string into tokens
			unsigned int nTokens=string::substrN(input,string::WS);
			if(nTokens>0){
				std::vector<std::string> tokens(nTokens);
				tokens[0]=std::string(std::strtok(input,string::WS));
				for(unsigned int i=1; i<nTokens; ++i){
					tokens[i]=std::string(std::strtok(NULL,string::WS));
				}
				if(nTokens==2){
					if(std::strcmp(tokens[0].c_str(),"atom_style")==0){
						if(std::strcmp(tokens[1].c_str(),"full")==0) styleAtom=STYLE_ATOM::FULL;
						if(std::strcmp(tokens[1].c_str(),"atomic")==0) styleAtom=STYLE_ATOM::ATOMIC;
					}
				}
			}
		}
		
		/* print data to screen */
		if(DEBUG_LAMMPS>1){
			std::cout<<"ATOM_STYLE = ";
			if(styleAtom==STYLE_ATOM::FULL) std::cout<<"FULL\n";
			if(styleAtom==STYLE_ATOM::ATOMIC) std::cout<<"ATOMIC\n";
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

}
