#ifndef FILE_OPS_H
#define FILE_OPS_H

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <map>
#include "string.hpp"

#ifndef DEBUG_FILE_OPS
#define DEBUG_FILE_OPS 0
#endif

namespace atomList{
	
	int atomIndex(const std::string& str, const std::vector<std::string>& atomNames);
	
	unsigned int loadNumAtoms(const char* str);
	std::vector<int>& loadAtomIndices(const char* str, std::vector<int>& atomIndices);
	std::vector<std::string>& loadAtomNames(const char* str, std::vector<std::string>& atomNames);
	
	char* atomName(const char* str, char* atomName);
	int atomIndex(const char* str);
	
	template <class AtomT>
	std::vector<AtomT>& loadAtoms(const char* atomStr, std::vector<AtomT>& atoms){
		//local function variables
		unsigned int nAtoms=0;
		std::vector<int> atomIndices;
		std::vector<std::string> atomNames;
		
		//load the number of atoms, atom indices, and atom names
		//std::cout<<"Loading nAtoms, atom indices, and atom names...\n";
		nAtoms=atomList::loadNumAtoms(atomStr);
		atomList::loadAtomIndices(atomStr,atomIndices);
		atomList::loadAtomNames(atomStr,atomNames);
		
		//build the atom list
		atoms.resize(nAtoms);
		for(unsigned int i=0; i<nAtoms; ++i){
			atoms[i].init();
			atoms[i].name()=atomNames[i];
			atoms[i].specie()=0;
			atoms[i].index()=atomIndices[i];
		}
		
		return atoms;
	}
	
	template <class AtomT>
	std::vector<AtomT>& loadAtoms(const char* atomStr, const std::vector<std::string>& speciesNames){
		//local function variables
		unsigned int nAtoms=0;
		std::vector<int> atomIndices;
		std::vector<std::string> atomNames;
		std::vector<AtomT> atoms;
		
		//load the number of atoms, atom indices, and atom names
		nAtoms=atomList::loadNumAtoms(atomStr);
		atomList::loadAtomIndices(atomStr,atomIndices);
		atomList::loadAtomNames(atomStr,atomNames);
		
		//build the atom list
		atoms.resize(nAtoms);
		for(unsigned int i=0; i<nAtoms; ++i){
			atoms[i].init();
			atoms[i]=AtomT(atomNames[i],atomIndex(atomNames[i],speciesNames),atomIndices[i]);
		}
		
		return atoms;
	}
	
	/*template <class AtomT>
	std::vector<AtomT> loadAtoms(const char* atomStr, std::vector<AtomT>& atoms){
		const char* func_name="loadAtoms(const char*,std::vector<AtomT>&,const SimAtomic<AtomT>&)";
		if(DEBUG_FILE_OPS>0) std::cout<<func_name<<":\n";
		//local function variables
		char* atomName=(char*)malloc(sizeof(char)*Constants::M);
		char* strtemp=(char*)malloc(sizeof(char)*Constants::M);
		char* substr=(char*)malloc(sizeof(char)*Constants::M);
		char* temp=(char*)malloc(sizeof(char)*Constants::M);
		std::vector<std::string> substrs;
		unsigned int nStrs=0;
		bool error=false;
		
		try{
			//clear the atom vector
			atoms.clear();
			
			//find the substrings
			std::strcpy(strtemp,atomStr);
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
					else{
						unsigned int index=std::atoi(std::strpbrk(substr,string::DIGITS))-1;
						std::strcpy(atomName,string::trim_right(substr,string::DIGITS));
						//atoms.push_back(AtomT(atomName,-1,index));
						atoms.push_back(AtomT());
						atoms.back().name()=atomName;
						atoms.back().specie()=-1;
						atoms.back().index()=index;
					}
				} else {
					//set of atoms
					if(DEBUG_FILE_OPS>1) std::cout<<"Set of Atoms\n";
					unsigned int beg, end;
					//find the beginning index
					std::strcpy(temp,std::strtok(substr,":"));
					if(std::strpbrk(temp,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
					else beg=std::atoi(std::strpbrk(temp,string::DIGITS))-1;
					//find the ending index
					std::strcpy(temp,std::strtok(NULL,":"));
					if(std::strpbrk(temp,string::DIGITS)==NULL) throw std::invalid_argument("Invalid atom specification: no index.");
					else end=std::atoi(std::strpbrk(temp,string::DIGITS))-1;
					//check the indices
					if(end<beg) throw std::invalid_argument("Invalid atomic indices.");
					//copy the atoms
					std::strcpy(atomName,string::trim_right(temp,string::DIGITS));
					for(unsigned int j=0; j<end-beg+1; ++j){
						//atoms.push_back(AtomT(atomName,-1,beg+j));
						atoms.push_back(AtomT());
						atoms.back().name()=atomName;
						atoms.back().specie()=-1;
						atoms.back().index()=beg+j;
					}
				}
			}
		}catch(std::exception& e){
			std::cout<<"ERROR in "<<func_name<<":\n";
			std::cout<<e.what()<<"\n";
			error=true;
		}
		
		free(atomName);
		free(strtemp);
		free(substr);
		free(temp);
		
		if(error){
			atoms.clear();
			throw std::invalid_argument("Invalid atom string.");
		} else return atoms;
	}*/
}

#endif
