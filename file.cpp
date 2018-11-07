#include "file.hpp"

namespace atomList{
	
	int atomIndex(const std::string& str, const std::vector<std::string>& atomNames){
		for(unsigned int i=0; i<atomNames.size(); i++)
			if(str==atomNames[i]) return i;
		return -1;
	}

	unsigned int loadNumAtoms(const char* str){
		const char* func_name="loadNumAtoms(const char*)";
		if(DEBUG_FILE_OPS>0) std::cout<<func_name<<":\n";
		//local function variables
		char* strtemp=(char*)malloc(sizeof(char)*string::M);
		char* substr=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
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
		
		free(strtemp);
		free(substr);
		free(temp);
		
		if(error) throw std::invalid_argument("Invalid atom string.");
		else return nAtoms;
	}
	
	std::vector<int>& loadAtomIndices(const char* str, std::vector<int>& atomIndices){
		const char* func_name="loadAtomIndices(const char*,std::vector<int>&)";
		if(DEBUG_FILE_OPS>0) std::cout<<func_name<<":\n";
		//local function variables
		char* strtemp=(char*)malloc(sizeof(char)*string::M);
		char* substr=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
		std::vector<std::string> substrs;
		unsigned int nStrs=0;
		bool error=false;
		
		try{
			//clear the vector
			atomIndices.clear();
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
					if(DEBUG_FILE_OPS>1) std::cout<<"Single Atom\n";
					if(std::strpbrk(substr,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
					else atomIndices.push_back(std::atoi(std::strpbrk(substr,string::DIGITS))-1);
				} else {
					//set of atoms
					if(DEBUG_FILE_OPS>1) std::cout<<"Set of Atoms\n";
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
					for(int j=0; j<end-beg+1; j++) atomIndices.push_back(beg+j);
				}
			}
		}catch(std::exception& e){
			std::cout<<"ERROR in "<<func_name<<":\n";
			std::cout<<e.what()<<"\n";
			error=true;
		}
		
		free(strtemp);
		free(substr);
		free(temp);
		
		if(error){
			atomIndices.clear();
			throw std::invalid_argument("Invalid atom string.");
		} else return atomIndices;
	}
	
	std::vector<std::string>& loadAtomNames(const char* str, std::vector<std::string>& atomNames){
		const char* func_name="loadAtomNames(const char*,std::vector<std::string>&)";
		if(DEBUG_FILE_OPS>0) std::cout<<func_name<<":\n";
		//local function variables
		char* strtemp=(char*)malloc(sizeof(char)*string::M);
		char* substr=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
		char* atomName=(char*)malloc(sizeof(char)*string::M);
		std::vector<std::string> substrs;
		unsigned int nStrs=0;
		bool error=false;
		
		try{
			//clear the vector
			atomNames.clear();
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
					else atomNames.push_back(std::string(string::trim_right(substr,string::DIGITS)));
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
						for(int j=0; j<end-beg+1; ++j) atomNames.push_back(std::string(atomName));
					}
				}
			}
		}catch(std::exception& e){
			std::cout<<"ERROR in "<<func_name<<":\n";
			std::cout<<e.what()<<"\n";
			error=true;
		}
		
		free(temp);
		free(atomName);
		
		if(error){
			atomNames.clear();
			throw std::invalid_argument("Invalid atom string.");
		} else return atomNames;
	}
	
	char* atomName(const char* str, char* atomName){
		return string::copy_left(atomName,str,string::DIGITS);
	}
	
	int atomIndex(const char* str){
		return std::atoi(std::strpbrk(str,string::DIGITS));
	}
}
