#ifndef GAUSSIAN_HPP
#define GAUSSIAN_HPP

//c++ libraries
#include <iostream>
//boost
#include <type_traits>
//atoms
#include "property.hpp"
#include "atom.hpp"
#include "molecule.hpp"
//string
#include "string.hpp"
//chem
#include "ptable.hpp"
//simulation
#include "sim.hpp"
#include "cell.hpp"
//units
#include "units.hpp"

#ifndef __cplusplus
	#error A C++ compiler is required.
#endif

#ifndef DEBUG_MDUG
#define DEBUG_MDUG 0
#endif

namespace GAUSSIAN{

//*************************************************
//Data Formats
//*************************************************

struct FormatPosn{
	enum type{
		INPUT,
		STANDARD,
		UNKNOWN
	};
	static FormatPosn::type load(const char* str);
};
std::ostream& operator<<(std::ostream& out, const FormatPosn::type& t);

struct FormatChg{
	enum type{
		NONE,
		MULLIKEN,
		ESP,
		HIRSHFELD,
		UNKNOWN
	};
	static FormatChg::type load(const char* str);
};
std::ostream& operator<<(std::ostream& out, const FormatChg::type& t);

struct FormatCalc{
	enum type{
		SP,
		OPT,
		BOMD,
		ADMP,
		UNKNOWN
	};
	static FormatCalc::type load(const char* str);
};
std::ostream& operator<<(std::ostream& out, const FormatCalc::type& t);

struct FormatAlpha{
	enum type{
		NONE,
		DIPOLE,
		NUMERIC,
		GAMMA,
		DCSHG,
		CUBIC,
		UNKNOWN
	};
	static FormatAlpha::type load(const char* str);
};
std::ostream& operator<<(std::ostream& out, const FormatAlpha::type& t);

struct FormatVersion{
	enum type{
		NONE,
		g09,
		g16,
		UNKNOWN
	};
	static FormatVersion::type load(const char* str);
};
std::ostream& operator<<(std::ostream& out, const FormatVersion::type& t);

//*****************************************************
//FORMAT struct
//*****************************************************

struct Format{
	//files
	std::string log;//log file
	std::string com;//com file
	
	//job options
	std::string basis;//basis set
	std::string chem;//model chemistry
	std::string job;//job
	std::string verbosity;
	std::vector<std::string> options;//job options
	std::vector<std::string> optionsInt;//integral options
	std::vector<std::string> optionsSCF;//scf options
	
	//compute options
	unsigned int nProc;//number of shared processors
	unsigned int nGB;//number of gb to use
	
	//data formats
	FormatCalc::type formatCalc;//calculation format
	FormatPosn::type formatPosn;//position format
	FormatChg::type formatChg;//charge format
	FormatAlpha::type formatAlpha;//alpha format
	FormatVersion::type formatVersion;//version
	
	//molecule data
	std::string name;//name of the molecule
	double charge;//net charge
	double spinMult;//spin multiplicity
	
	//interval
	int beg,end,stride;
	
	//constructors/destructors
	Format():beg(1),end(-1),stride(1),spinMult(1),charge(0),name("MOL"),nProc(1),nGB(1),verbosity("N"){};
	
	//loading
	static Format& load(const std::vector<std::string>& strlist, Format& format);
};

//*************************************************
//Gaussian File Class
//*************************************************

struct GaussFile{
	typedef Atom<Name,Species,Index,Position> Particle;
	
	//job options
	std::string basis;//basis set
	std::string chem;//model chemistry
	std::string job;//job
	std::vector<std::string> options;//job options
	std::vector<std::string> optionsInt;//integral options
	std::vector<std::string> optionsSCF;//scf options
	
	//compute options
	unsigned int nProc;//number of shared processors
	unsigned int nGB;//number of gb to use
	
	//molecule data
	std::string name;//name of the molecule
	double charge;//net charge
	double spinMult;//spin multiplicity
	std::vector<Particle> atoms;//the atoms of the molecule
	
	//data formats
	FormatCalc::type formatCalc;
	FormatPosn::type formatPosn;
	FormatChg::type formatChg;
	FormatAlpha::type formatAlpha;
	FormatVersion::type formatVersion;
};

void print(const char* fileName, const GaussFile& gauss);

template <class AtomT>
void print(const char* fileName, const Format& gauss, const SimAtomic<AtomT>& sim, unsigned int t){
	if(DEBUG_MDUG>0) std::cout<<"print(const char*,const Format&,const SimAtomic<AtomT>&,unsigned int):\n";
	
	FILE* writer=fopen(fileName,"w");
	if(writer==NULL) throw std::runtime_error("File could not be opened.");
	
	//print checkpoint file
	fprintf(writer, "%%Chk=%s_%s\n", sim.system().c_str(), gauss.job.c_str());
	//print the compute info
	fprintf(writer, "%%NProcShared=%i\n", gauss.nProc);
	fprintf(writer, "%%Mem=%iGB\n", gauss.nGB);
	//print the job section
	fprintf(writer, "#%s %s/%s %s", gauss.verbosity.c_str(), gauss.chem.c_str(), gauss.basis.c_str(), gauss.job.c_str());
	if(gauss.options.size()>0){
		fprintf(writer, "=(");
		for(unsigned int i=0; i<gauss.options.size()-1; ++i)
			fprintf(writer,"%s,",gauss.options[i].c_str());
		fprintf(writer,"%s) ",gauss.options.back().c_str());
	}
	fprintf(writer," ");
	if(gauss.optionsInt.size()>0){
		fprintf(writer, "Int=(");
		for(unsigned int i=0; i<gauss.optionsInt.size()-1; ++i) 
				fprintf(writer,"%s,",gauss.optionsInt[i].c_str());
		fprintf(writer,"%s)",gauss.optionsInt.back().c_str());
	}
	fprintf(writer," ");
	if(gauss.optionsSCF.size()>0){
		fprintf(writer, "SCF=(");
		for(unsigned int i=0; i<gauss.optionsSCF.size()-1; ++i)
			fprintf(writer,"%s,",gauss.optionsSCF[i].c_str());
		fprintf(writer,"%s)",gauss.optionsSCF.back().c_str());
	}
	fprintf(writer,"\n\n");
	
	//print the name
	fprintf(writer,"%s\n\n", sim.system().c_str());
	
	//print the charge and multiplicity
	fprintf(writer, "%i %i\n", (int)(gauss.charge), (int)(gauss.spinMult));
	//print the atoms
	for(unsigned int i=0; i<sim.nAtoms(); ++i){
		fprintf(writer, "%s  %12.8f  %12.8f  %12.8f\n",
			sim.atom(t,i).name().c_str(),
			sim.atom(t,i).posn()[0],
			sim.atom(t,i).posn()[1],
			sim.atom(t,i).posn()[2]
		);
	}
	fprintf(writer,"\n");
	
	fclose(writer);
}

std::string loadName(FILE* reader);

void loadVersion(FILE* reader, GaussFile& gauss);
void loadVersion(FILE* reader, Format& format);

template <class AtomT, class MoleculeT>
void loadStruct(FILE* reader, const GaussFile& gauss, MoleculeT& mol){
	if(DEBUG_MDUG>0) std::cout<<"loadStruct<AtomT,MoleculeT>(const char*,const GaussFile&,MoleculeT&):\n";
	
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//property tags
		char* tag=(char*)malloc(sizeof(char)*string::M);
	
	double s=0.0;
	if(units::consts::system()==units::System::AU) s=units::BOHRpANG;
	else if(units::consts::system()==units::System::METAL) s=1.0;
	else throw std::runtime_error("Invalid units.");
	
	//set the tags
	if(gauss.formatPosn==FormatPosn::STANDARD) std::strcpy(tag,"Standard orientation:");
	else if(gauss.formatPosn==FormatPosn::INPUT) std::strcpy(tag,"Input orientation:");
	
	rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strcmp(string::trim(input),tag)==0){
			//atom data
			unsigned int an=0;
			std::string name;
			unsigned int specie=0,index=0,specieMax=0;
			Eigen::Vector3d posn;
			//skip 4 lines
			fgets(input,string::M,reader);
			fgets(input,string::M,reader);
			fgets(input,string::M,reader);
			fgets(input,string::M,reader);
			//read in the atoms
			while(fgets(input,string::M,reader)!=NULL){
				string::trim(input);
				if(string::substrN(input,string::WS)<=1) break;//break if there are no substrings
				std::strtok(input,string::WS);//break into tokens, skip the first
				//find the name
				an=(unsigned int)std::atoi(std::strtok(NULL,string::WS));
				name=std::string(PTable::elementName(an));
				//skip the next token
				std::strtok(NULL,string::WS);
				//find the coordinates
				posn[0]=std::atof(std::strtok(NULL,string::WS));
				posn[1]=std::atof(std::strtok(NULL,string::WS));
				posn[2]=std::atof(std::strtok(NULL,string::WS));
				//add the atom to the molecule
				AtomT atom;
				atom.init();
				atom.name()=name;
				atom.an()=an;
				atom.specie()=specie;
				atom.index()=index;
				atom.posn()=s*posn;
				mol.add(atom);
				//increment the index
				++index;
			}
			break;
		}
	}
	
	//free local variables
	free(input);
	free(temp);
	free(tag);
}

template <class AtomT>
void loadStruct(FILE* reader, Format& format, SimAtomic<AtomT>& sim){
	if(DEBUG_MDUG>0) std::cout<<"loadStruct<AtomT,MoleculeT>(const char*,const GaussFile&,MoleculeT&):\n";
	
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//property tags
		char* tag=(char*)malloc(sizeof(char)*string::M);
	//atom info
		std::vector<std::string> atomNames;
		std::vector<unsigned int> nAtoms;
	//timesteps
		unsigned int ts=0;
		
	//set the tags
	std::strcpy(tag,"Symbolic Z-matrix");
	
	double s=0.0;
	if(units::consts::system()==units::System::AU) s=1.0;
	else if(units::consts::system()==units::System::METAL) s=units::ANGpBOHR;
	else throw std::runtime_error("Invalid units.");
	
	//read in the timesteps
	if(DEBUG_MDUG>0) std::cout<<"Reading in timesteps...\n";
	rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strstr(input,"MaxPoints")!=NULL){
			char* tempstr=std::strstr(input,"MaxPoints");
			std::strtok(tempstr,"=");
			ts=std::atoi(std::strtok(NULL,",)"));
			break;
		}
	}
	if(ts==0) ts=50;
	++ts;
	if(DEBUG_MDUG>0) std::cout<<"ts = "<<ts<<"\n";
	
	//set the interval
	sim.beg()=format.beg-1;
	if(format.end<0){
		sim.end()=ts+format.end;
		format.end=sim.end();
	} else sim.end()=format.end-1;
	unsigned int interval=sim.end()-sim.beg()+1;
	if(DEBUG_MDUG>0) std::cout<<"(beg,end,stride) = ("<<sim.beg()<<","<<sim.end()<<","<<format.stride<<")\n";
	if(DEBUG_MDUG>0) std::cout<<"interval = "<<interval<<"\n";
	if(sim.beg()<0) throw std::invalid_argument("Invalid beginning timestep.");
	if(sim.end()>=ts) throw std::invalid_argument("Invalid ending timestep.");
	
	//read in atom info
	rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strstr(input,tag)!=NULL){
			fgets(input,string::M,reader);//skip a line
			while(!string::empty(fgets(input,string::M,reader))){
				std::string name=std::string(std::strtok(input,string::WS));
				bool match=false;
				for(unsigned int n=0; n<atomNames.size(); ++n){
					if(name==atomNames[n]){
						++nAtoms[n];
						match=true; break;
					}
				}
				if(!match){
					atomNames.push_back(name);
					nAtoms.push_back(1);
				}
			}
		}
	}
	if(DEBUG_MDUG>0){
		std::cout<<"ATOM_NAMES = "; for(unsigned int i=0; i<atomNames.size(); ++i) std::cout<<atomNames[i]<<" "; std::cout<<"\n";
		std::cout<<"ATOM_NUMBERS = "; for(unsigned int i=0; i<nAtoms.size(); ++i) std::cout<<nAtoms[i]<<" "; std::cout<<"\n";
	}
	
	//resize the simulation
	if(DEBUG_MDUG>0) std::cout<<"Resizing the simulation...\n";
	std::cout<<"TS-TOT = "<<interval/format.stride+interval%format.stride<<"\n";
	sim.resize(interval/format.stride+interval%format.stride,nAtoms,atomNames);
	
	//read in positions
	if(DEBUG_MDUG>0) std::cout<<"Loading positions...\n";
	rewind(reader);
	unsigned int t=0,tt=0;
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strstr(input,"Cartesian coordinates")){
			if(DEBUG_MDUG>0 && t%1000) std::cout<<"t = "<<t<<"\n";
			else if(DEBUG_MDUG>1) std::cout<<"t = "<<t<<"\n";
			if(t>=sim.beg() && (t-sim.beg())%format.stride==0){
				for(unsigned int n=0; n<sim.nAtoms(); ++n){
					fgets(input,string::M,reader);
					std::strtok(input,"=");
					std::strtok(NULL,"=");
					sim.atom(tt,n).posn()[0]=s*std::atof(std::strtok(NULL,"D"));
					sim.atom(tt,n).posn()[0]*=std::pow(10,std::atof(std::strtok(NULL,"DY")));
					sim.atom(tt,n).posn()[1]=s*std::atof(std::strtok(NULL,"=D"));
					sim.atom(tt,n).posn()[1]*=std::pow(10,std::atof(std::strtok(NULL,"DZ")));
					sim.atom(tt,n).posn()[2]=s*std::atof(std::strtok(NULL,"=D"));
					sim.atom(tt,n).posn()[2]*=std::pow(10,std::atof(std::strtok(NULL,"=")));
				}
				++tt;
			}
			++t;
		}
	}
	
	//free local variables
	if(DEBUG_MDUG>0) std::cout<<"Freeing local variables...\n";
	free(input);
	free(temp);
	free(tag);
}

template <class AtomT, class MoleculeT>
void loadDipole(FILE* reader, const GaussFile& gauss, MoleculeT& mol, std::false_type){}
template <class AtomT, class MoleculeT>
void loadDipole(FILE* reader, const GaussFile& gauss, MoleculeT& mol, std::true_type){
	if(DEBUG_MDUG>0) std::cout<<"loadDipole(const char*,const GaussFile&,MoleculeT&,std::true_type):\n";
	
	/* local function variables */
	//file i/o
	char* input=(char*)malloc(sizeof(char)*string::M);
	char* temp=(char*)malloc(sizeof(char)*string::M);
	
	double s=0.0;
	if(units::consts::system()==units::System::AU) s=0.393430307;
	else if(units::consts::system()==units::System::METAL) s=0.20819434;
	else throw std::runtime_error("Invalid units.");
	
	rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strpbrk(input,string::WS)!=0 && !string::empty(input)){
			std::strcpy(temp,std::strtok(input,string::WS));
			if(std::strcmp(temp,"Dipole")==0){
				fgets(input,string::M,reader);
				std::strtok(input,"=");
				mol.dipole()[0]=std::atof(std::strtok(NULL,string::WS));
				std::strtok(NULL,"=");
				mol.dipole()[1]=std::atof(std::strtok(NULL,string::WS));
				std::strtok(NULL,"=");
				mol.dipole()[2]=std::atof(std::strtok(NULL,string::WS));
				//convert
				mol.dipole()*=s;
				break;
			}
		}
	}
	
	free(input);
	free(temp);
}

template <class AtomT, class MoleculeT>
void loadChg(FILE* reader, const GaussFile& gauss, MoleculeT& mol, std::false_type){};
template <class AtomT, class MoleculeT>
void loadChg(FILE* reader, const GaussFile& gauss, MoleculeT& mol, std::true_type){
	if(DEBUG_MDUG>0) std::cout<<"loadChg(const char*,const GaussFile&,MoleculeT&,std::true_type):\n";
	
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//property tags
		char* tag=(char*)malloc(sizeof(char)*string::M);
	
	//load the charges
	rewind(reader);
	if(gauss.formatChg==FormatChg::MULLIKEN){
		if(DEBUG_MDUG>0) std::cout<<"FORMAT: MULLIKEN\n";
		//MULLIKEN
		if(gauss.formatVersion==FormatVersion::g09) std::strcpy(tag,"Mulliken atomic charges:");
		else if(gauss.formatVersion==FormatVersion::g16) std::strcpy(tag,"Mulliken charges and spin densities:");
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strcmp(string::trim(input),tag)==0){
				fgets(input,string::M,reader);//skip a line
				for(unsigned int i=0; i<mol.nAtoms(); ++i){
					std::strtok(fgets(input,string::M,reader),string::WS);
					std::strtok(NULL,string::WS);
					mol.atom(i).charge()=std::atof(std::strtok(NULL,string::WS));
				}
				break;
			}
		}
	} else if(gauss.formatChg==FormatChg::ESP){
		if(DEBUG_MDUG>0) std::cout<<"FORMAT: ESP\n";
		//ESP
		std::strcpy(tag,"Charges from ESP fit");
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strpbrk(input,",")!=NULL){
				string::copy_left(temp,input,",");
				std::strcpy(input,temp);
			}
			if(std::strcmp(string::trim(input),tag)==0){
				fgets(input,string::M,reader);//skip a line
				for(unsigned int i=0; i<mol.nAtoms(); ++i){
					std::strtok(fgets(input,string::M,reader),string::WS);
					std::strtok(NULL,string::WS);
					mol.atom(i).charge()=std::atof(std::strtok(NULL,string::WS));
				}
				break;
			}
		}
	} else if(gauss.formatChg==FormatChg::HIRSHFELD){
		if(DEBUG_MDUG>0) std::cout<<"FORMAT: HIRSHFELD\n";
		//HIRSHFELD
		if(gauss.formatVersion==FormatVersion::g09) std::strcpy(tag,"Hirshfeld spin densities");
		else if(gauss.formatVersion==FormatVersion::g16) std::strcpy(tag,"Hirshfeld charges");
		while(fgets(input,string::M,reader)!=NULL){
			if(std::strpbrk(input,",")!=NULL){
				string::copy_left(temp,input,",");
				std::strcpy(input,temp);
			}
			if(std::strcmp(string::trim(input),tag)==0){
				fgets(input,string::M,reader);//skip a line
				for(unsigned int i=0; i<mol.nAtoms(); ++i){
					if(gauss.formatVersion==FormatVersion::g09){
						std::strtok(fgets(input,string::M,reader),string::WS);
						std::strtok(NULL,string::WS);
						std::strtok(NULL,string::WS);
						mol.atom(i).charge()=std::atof(std::strtok(NULL,string::WS));
					}
					else if(gauss.formatVersion==FormatVersion::g16){
						std::strtok(fgets(input,string::M,reader),string::WS);
						std::strtok(NULL,string::WS);
						mol.atom(i).charge()=std::atof(std::strtok(NULL,string::WS));
					}
				}
				break;
			}
		}
	}
	
	//free local variables
	free(input);
	free(temp);
	free(tag);
}

template <class AtomT, class MoleculeT>
void loadAlpha(FILE* reader, const GaussFile& gauss, MoleculeT& mol, std::false_type){};
template <class AtomT, class MoleculeT>
void loadAlpha(FILE* reader, const GaussFile& gauss, MoleculeT& mol, std::true_type){
	if(DEBUG_MDUG>0) std::cout<<"loadAlpha(const char*,const GaussFile&,MoleculeT&,std::true_type):\n";
	
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//property tags
		char* tag=(char*)malloc(sizeof(char)*string::M);
		const char* delim=":,";
	
	//set the units
	double s=0;
	if(units::consts::system()==units::System::AU) s=1.0;
	else if(units::consts::system()==units::System::METAL) s=std::pow(units::ANGpBOHR,3);
	else throw std::runtime_error("Invalid units.");
	
	if(
		gauss.formatAlpha==FormatAlpha::DIPOLE ||
		gauss.formatAlpha==FormatAlpha::NUMERIC
	) std::strcpy(tag,"Exact polarizability");
	else if(
		gauss.formatAlpha==FormatAlpha::GAMMA ||
		gauss.formatAlpha==FormatAlpha::CUBIC ||
		gauss.formatAlpha==FormatAlpha::DCSHG 
	) std::strcpy(tag,"Dipole polarizability");
	
	//load the polarizability
	rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		string::copy_left(temp,input,delim);
		if(std::strcmp(string::trim(temp),tag)==0){
			if(
				gauss.formatAlpha==FormatAlpha::DIPOLE ||
				gauss.formatAlpha==FormatAlpha::NUMERIC
			){
				if(DEBUG_MDUG>1 && gauss.formatAlpha==FormatAlpha::DIPOLE) std::cout<<"FORMAT = DIPOLE\n";
				if(DEBUG_MDUG>1 && gauss.formatAlpha==FormatAlpha::NUMERIC) std::cout<<"FORMAT = DIPOLE\n";
				std::strtok(input,delim);
				//read in the polarizability
				Eigen::Matrix3d alpha;
				alpha(0,0)=std::atof(std::strtok(NULL,string::WS));//xx
				alpha(0,1)=std::atof(std::strtok(NULL,string::WS));//xy
				alpha(1,1)=std::atof(std::strtok(NULL,string::WS));//yy
				alpha(0,2)=std::atof(std::strtok(NULL,string::WS));//xz
				alpha(1,2)=std::atof(std::strtok(NULL,string::WS));//yz
				alpha(2,2)=std::atof(std::strtok(NULL,string::WS));//zz
				//reflect
				alpha(1,0)=alpha(0,1);
				alpha(2,0)=alpha(0,2);
				alpha(2,1)=alpha(1,2);
				if(DEBUG_MDUG>1) std::cout<<"alpha = \n"<<alpha<<"\n";
				mol.alpha()=s*alpha;
			}
			else if(
				gauss.formatAlpha==FormatAlpha::GAMMA ||
				gauss.formatAlpha==FormatAlpha::CUBIC ||
				gauss.formatAlpha==FormatAlpha::DCSHG 
			){
				if(DEBUG_MDUG>1) std::cout<<"Format: GAMMA-CUBIC-DCSHG\n";
				//skip 5 lines
				fgets(input,string::M,reader);//skip a line
				fgets(input,string::M,reader);//skip a line
				fgets(input,string::M,reader);//skip a line
				fgets(input,string::M,reader);//skip a line
				fgets(input,string::M,reader);//skip a line
				//read in the polarizability
				Eigen::Matrix3d alpha;
				//xx
					fgets(input,string::M,reader);
					std::strtok(input,string::WS);
					std::strcpy(temp,std::strtok(NULL,string::WS));
					alpha(0,0)=std::atof(std::strtok(temp,"D"));
					alpha(0,0)*=std::pow(10,std::atof(std::strtok(NULL,"D")));
				//yx
					fgets(input,string::M,reader);
					std::strtok(input,string::WS);
					std::strcpy(temp,std::strtok(NULL,string::WS));
					alpha(0,1)=std::atof(std::strtok(temp,"D"));
					alpha(0,1)*=std::pow(10,std::atof(std::strtok(NULL,"D")));
				//yy
					fgets(input,string::M,reader);
					std::strtok(input,string::WS);
					std::strcpy(temp,std::strtok(NULL,string::WS));
					alpha(1,1)=std::atof(std::strtok(temp,"D"));
					alpha(1,1)*=std::pow(10,std::atof(std::strtok(NULL,"D")));
				//zx
					fgets(input,string::M,reader);
					std::strtok(input,string::WS);
					std::strcpy(temp,std::strtok(NULL,string::WS));
					alpha(0,2)=std::atof(std::strtok(temp,"D"));
					alpha(0,2)*=std::pow(10,std::atof(std::strtok(NULL,"D")));
				//zy
					fgets(input,string::M,reader);
					std::strtok(input,string::WS);
					std::strcpy(temp,std::strtok(NULL,string::WS));
					alpha(2,1)=std::atof(std::strtok(temp,"D"));
					alpha(2,1)*=std::pow(10,std::atof(std::strtok(NULL,"D")));
				//zz
					fgets(input,string::M,reader);
					std::strtok(input,string::WS);
					std::strcpy(temp,std::strtok(NULL,string::WS));
					alpha(2,2)=std::atof(std::strtok(temp,"D"));
					alpha(2,2)*=std::pow(10,std::atof(std::strtok(NULL,"D")));
				//reflect
					alpha(1,0)=alpha(0,1);
					alpha(2,0)=alpha(0,2);
					alpha(2,1)=alpha(1,2);
				//assign
					mol.alpha()=s*alpha;
			}
		}
	}
	
	//free local variables
	free(input);
	free(temp);
	free(tag);
};

template <class AtomT, class MoleculeT>
void loadGamma(FILE* reader, const GaussFile& gauss, MoleculeT& mol, std::false_type){};
template <class AtomT, class MoleculeT>
void loadGamma(FILE* reader, const GaussFile& gauss, MoleculeT& mol, std::true_type){
	if(DEBUG_MDUG>0) std::cout<<"loadAlpha(const char*,const GaussFile&,MoleculeT&,std::true_type):\n";
	
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//property tags
		char* tag=(char*)malloc(sizeof(char)*string::M);
		
	if(
		gauss.formatAlpha==FormatAlpha::GAMMA ||
		gauss.formatAlpha==FormatAlpha::CUBIC ||
		gauss.formatAlpha==FormatAlpha::DCSHG 
	) std::strcpy(tag,"Second dipole hyperpolarizability");
	
	//load the second hyperpolarizability
	rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strpbrk(input,",:")!=NULL) string::copy_left(temp,input,";,");
		if(std::strcmp(string::trim(temp),tag)==0){
			//skip 7 lines
			for(unsigned int i=0; i<7; ++i) fgets(input,string::M,reader);//skip a line
			//read in the second polarizability
			//xxxx
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[0][0][0][0]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[0][0][0][0]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//xxyx
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[0][0][1][0]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[0][0][1][0]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//xxyy
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[0][0][1][1]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[0][0][1][1]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//yxyy
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[1][0][1][1]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[1][0][1][1]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//yyyy
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[1][1][1][1]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[1][1][1][1]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//xxzx
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[0][0][2][0]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[0][0][2][0]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//xxzy
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[0][0][2][1]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[0][0][2][1]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//yxzy
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[1][0][2][1]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[1][0][2][1]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//yyzy
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[1][1][2][1]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[1][1][2][1]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//xxzz
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[0][0][2][2]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[0][0][2][2]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//yxzz
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[1][0][2][2]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[1][0][2][2]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//yyzz
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[1][1][2][2]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[1][1][2][2]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//zxzz
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[2][0][2][2]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[2][0][2][2]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//zyzz
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[2][1][2][2]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[2][1][2][2]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//zzzz
				fgets(input,string::M,reader);
				std::strtok(input,string::WS);
				std::strcpy(temp,std::strtok(NULL,string::WS));
				mol.gamma()[2][2][2][2]=std::atof(std::strtok(temp,"D"));
				mol.gamma()[2][2][2][2]*=std::pow(10,std::atof(std::strtok(NULL,"D")));
			//reflect
			//xxyx
				mol.gamma()[0][0][0][1]=mol.gamma()[0][0][1][0];
				mol.gamma()[0][1][0][0]=mol.gamma()[0][0][1][0];
				mol.gamma()[1][0][0][0]=mol.gamma()[0][0][1][0];
			//xxyy
				mol.gamma()[0][1][0][1]=mol.gamma()[0][0][1][1];
				mol.gamma()[0][1][1][0]=mol.gamma()[0][0][1][1];
				mol.gamma()[1][0][0][1]=mol.gamma()[0][0][1][1];
				mol.gamma()[1][0][1][0]=mol.gamma()[0][0][1][1];
				mol.gamma()[1][1][0][0]=mol.gamma()[0][0][1][1];
			//yxyy
				mol.gamma()[0][1][1][1]=mol.gamma()[1][0][1][1];
				mol.gamma()[1][1][0][1]=mol.gamma()[1][0][1][1];
				mol.gamma()[1][1][1][0]=mol.gamma()[1][0][1][1];
			//xxzx
				mol.gamma()[0][0][0][2]=mol.gamma()[0][0][2][0];
				mol.gamma()[0][2][0][0]=mol.gamma()[0][0][2][0];
				mol.gamma()[2][0][0][0]=mol.gamma()[0][0][2][0];
			//xxzy
				mol.gamma()[0][0][1][2]=mol.gamma()[0][0][2][1];
				mol.gamma()[0][1][0][2]=mol.gamma()[0][0][2][1];
				mol.gamma()[0][1][2][0]=mol.gamma()[0][0][2][1];
				mol.gamma()[0][2][0][1]=mol.gamma()[0][0][2][1];
				mol.gamma()[0][2][1][0]=mol.gamma()[0][0][2][1];
				mol.gamma()[1][0][0][2]=mol.gamma()[0][0][2][1];
				mol.gamma()[1][0][2][0]=mol.gamma()[0][0][2][1];
				mol.gamma()[1][2][0][0]=mol.gamma()[0][0][2][1];
				mol.gamma()[2][0][0][1]=mol.gamma()[0][0][2][1];
				mol.gamma()[2][0][1][0]=mol.gamma()[0][0][2][1];
				mol.gamma()[2][1][0][0]=mol.gamma()[0][0][2][1];
			//yxzy
				mol.gamma()[0][1][1][2]=mol.gamma()[1][0][2][1];
				mol.gamma()[0][1][2][1]=mol.gamma()[1][0][2][1];
				mol.gamma()[0][2][1][1]=mol.gamma()[1][0][2][1];
				mol.gamma()[1][0][1][2]=mol.gamma()[1][0][2][1];
				mol.gamma()[1][1][0][2]=mol.gamma()[1][0][2][1];
				mol.gamma()[1][1][2][0]=mol.gamma()[1][0][2][1];
				mol.gamma()[1][2][0][1]=mol.gamma()[1][0][2][1];
				mol.gamma()[1][2][1][0]=mol.gamma()[1][0][2][1];
				mol.gamma()[2][0][1][1]=mol.gamma()[1][0][2][1];
				mol.gamma()[2][1][0][1]=mol.gamma()[1][0][2][1];
				mol.gamma()[2][1][1][0]=mol.gamma()[1][0][2][1];
			//yyzy	
				mol.gamma()[1][1][1][2]=mol.gamma()[1][1][2][1];
				mol.gamma()[1][2][1][1]=mol.gamma()[1][1][2][1];
				mol.gamma()[2][1][1][1]=mol.gamma()[1][1][2][1];
			//xxzz
				mol.gamma()[0][2][0][2]=mol.gamma()[0][0][2][2];
				mol.gamma()[0][2][2][0]=mol.gamma()[0][0][2][2];
				mol.gamma()[2][0][0][2]=mol.gamma()[0][0][2][2];
				mol.gamma()[2][0][2][0]=mol.gamma()[0][0][2][2];
				mol.gamma()[2][2][0][0]=mol.gamma()[0][0][2][2];
			//yxzz
				mol.gamma()[0][1][2][2]=mol.gamma()[1][0][2][2];
				mol.gamma()[0][2][1][2]=mol.gamma()[1][0][2][2];
				mol.gamma()[0][2][2][1]=mol.gamma()[1][0][2][2];
				mol.gamma()[1][0][2][2]=mol.gamma()[1][0][2][2];
				mol.gamma()[1][2][0][2]=mol.gamma()[1][0][2][2];
				mol.gamma()[1][2][2][0]=mol.gamma()[1][0][2][2];
				mol.gamma()[2][0][1][2]=mol.gamma()[1][0][2][2];
				mol.gamma()[2][1][0][2]=mol.gamma()[1][0][2][2];
				mol.gamma()[2][1][2][0]=mol.gamma()[1][0][2][2];
				mol.gamma()[2][2][0][1]=mol.gamma()[1][0][2][2];
				mol.gamma()[2][2][1][0]=mol.gamma()[1][0][2][2];
			//yyzz
				mol.gamma()[1][2][1][2]=mol.gamma()[1][1][2][2];
				mol.gamma()[1][2][2][1]=mol.gamma()[1][1][2][2];
				mol.gamma()[2][1][1][2]=mol.gamma()[1][1][2][2];
				mol.gamma()[2][1][2][1]=mol.gamma()[1][1][2][2];
				mol.gamma()[2][2][1][1]=mol.gamma()[1][1][2][2];
			//zxzz
				mol.gamma()[0][2][2][2]=mol.gamma()[2][0][2][2];
				mol.gamma()[2][2][0][2]=mol.gamma()[2][0][2][2];
				mol.gamma()[2][2][2][0]=mol.gamma()[2][0][2][2];
			//zyzz
				mol.gamma()[1][2][2][2]=mol.gamma()[2][1][2][2];
				mol.gamma()[2][2][1][2]=mol.gamma()[2][1][2][2];
				mol.gamma()[2][2][2][1]=mol.gamma()[2][1][2][2];
		}
	}
	
	//free local variables
	free(input);
	free(temp);
	free(tag);
}

template <class AtomT, class MoleculeT>
void load(const char* file, GaussFile& gauss, MoleculeT& mol){
	if(DEBUG_MDUG>0) std::cout<<"load(const char*,const GaussFile&,MoleculeT&):\n";
	
	FILE* reader=fopen(file,"r");
	if(reader==NULL) throw std::runtime_error("File could not be opened.");
	
	//clear the molecule
	//mol.clear();
	mol.init();
	
	//load the version
	if(DEBUG_MDUG>0) std::cout<<"Loading version...\n";
	loadVersion(reader,gauss);
	
	//load the atoms and positions
	if(DEBUG_MDUG>0) std::cout<<"Loading positions...\n";
	loadStruct<AtomT,MoleculeT>(reader,gauss,mol);
	
	//load the charges
	if(DEBUG_MDUG>0) std::cout<<"Loading charge...\n";
	loadChg<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Charge,AtomT>());
	
	//load the dipole moment
	if(DEBUG_MDUG>0) std::cout<<"Loading dipole...\n";
	loadDipole<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Dipole,MoleculeT>());
	
	//load the polarizability
	if(DEBUG_MDUG>0) std::cout<<"Loading polarizability...\n";
	loadAlpha<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Alpha,MoleculeT>());
	
	fclose(reader);
}

template <class AtomT>
void load(const char* file, GaussFile& gauss, SimAtomic<AtomT>& sim, int beg=0, int end=-1, int stride=1){
	if(DEBUG_MDUG>0) std::cout<<"load(const char*,const GaussFile&,SimAtomic<AtomT>&):\n";
	
	FILE* reader=fopen(file,"r");
	if(reader==NULL) throw std::runtime_error("File could not be opened.");
	
	//clear the simulation
	if(DEBUG_MDUG>0) std::cout<<"Clearing simulation...\n";
	sim.clear();
	
	//load the version
	if(DEBUG_MDUG>0) std::cout<<"Loading version...\n";
	loadVersion(reader,gauss);
	
	//load the atoms and positions
	if(DEBUG_MDUG>0) std::cout<<"Loading positions...\n";
	loadStruct<AtomT>(reader,gauss,sim,beg,end,stride);
	
	//set the unit cell
	if(DEBUG_MDUG>0) std::cout<<"Setting the unit cell...\n";
	for(unsigned int t=0; t<sim.timesteps(); ++t) sim.cell(t).init(Eigen::Vector3d::Identity()*-1);
	
	//load the charges
	//if(DEBUG_MDUG>0) std::cout<<"Loading charge...\n";
	//loadChg<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Charge,AtomT>());
	
	//load the dipole moment
	//if(DEBUG_MDUG>0) std::cout<<"Loading dipole...\n";
	//loadDipole<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Dipole,MoleculeT>());
	
	//load the polarizability
	//if(DEBUG_MDUG>0) std::cout<<"Loading polarizability...\n";
	//loadAlpha<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Alpha,MoleculeT>());
	
	fclose(reader);
}

template <class AtomT>
void load(Format& format, SimAtomic<AtomT>& sim){
	if(DEBUG_MDUG>0) std::cout<<"load(GaussFormat&,SimAtomic<AtomT>&):\n";
	
	FILE* reader=fopen(format.log.c_str(),"r");
	if(reader==NULL) throw std::runtime_error("File could not be opened.");
	
	//clear the simulation
	if(DEBUG_MDUG>0) std::cout<<"Clearing simulation...\n";
	sim.clear();
	
	//load the version
	if(DEBUG_MDUG>0) std::cout<<"Loading version...\n";
	loadVersion(reader,format);
	
	//load the name
	if(DEBUG_MDUG>0) std::cout<<"Loading name...\n";
	sim.system()=loadName(reader);
	
	//load the atoms and positions
	if(DEBUG_MDUG>0) std::cout<<"Loading positions...\n";
	loadStruct<AtomT>(reader,format,sim);
	
	//set the unit cell
	if(DEBUG_MDUG>0) std::cout<<"Setting the unit cell...\n";
	Eigen::Vector3d max=Eigen::Vector3d::Zero();
	Eigen::Vector3d min=Eigen::Vector3d::Zero();
	Eigen::Vector3d avg=Eigen::Vector3d::Zero();
	for(unsigned int t=0; t<sim.timesteps(); ++t){
		Eigen::Vector3d temp=Eigen::Vector3d::Zero();
		for(unsigned int n=0; n<sim.nAtoms(); ++n){
			if(sim.atom(t,n).posn()[0]>max[0]) max[0]=sim.atom(t,n).posn()[0];
			if(sim.atom(t,n).posn()[1]>max[1]) max[1]=sim.atom(t,n).posn()[1];
			if(sim.atom(t,n).posn()[2]>max[2]) max[2]=sim.atom(t,n).posn()[2];
			if(sim.atom(t,n).posn()[0]<min[0]) min[0]=sim.atom(t,n).posn()[0];
			if(sim.atom(t,n).posn()[1]<min[1]) min[1]=sim.atom(t,n).posn()[1];
			if(sim.atom(t,n).posn()[2]<min[2]) min[2]=sim.atom(t,n).posn()[2];
			temp.noalias()+=sim.atom(t,n).posn();
		}
		avg.noalias()+=temp/sim.nAtoms();
	}
	avg/=sim.timesteps();
	Eigen::Matrix3d lv=Eigen::Matrix3d::Identity()*1.5*(max.maxCoeff()-min.minCoeff());
	std::cout<<"max = "<<max.transpose()<<"\n";
	std::cout<<"min = "<<min.transpose()<<"\n";
	std::cout<<"avg = "<<avg.transpose()<<"\n";
	std::cout<<"lv = \n"<<lv<<"\n";
	for(unsigned int t=0; t<sim.timesteps(); ++t) sim.cell(t).init(lv);
	
	//move the atoms such that the average position is at the center of the cell
	if(DEBUG_MDUG>0) std::cout<<"Moving the average to the center of the cell...\n";
	Eigen::Vector3d center;
	center[0]=lv(0,0)*0.5;
	center[1]=lv(1,1)*0.5;
	center[2]=lv(2,2)*0.5;
	for(unsigned int t=0; t<sim.timesteps(); ++t){
		for(unsigned n=0; n<sim.nAtoms(); ++n){
			sim.atom(t,n).posn()+=center-avg;
			Cell::returnToCell(sim.atom(t,n).posn(),sim.atom(t,n).posn(),sim.cell(t).R(),sim.cell(t).RInv());
		}
	}
	
	//load the charges
	//if(DEBUG_MDUG>0) std::cout<<"Loading charge...\n";
	//loadChg<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Charge,AtomT>());
	
	//load the dipole moment
	//if(DEBUG_MDUG>0) std::cout<<"Loading dipole...\n";
	//loadDipole<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Dipole,MoleculeT>());
	
	//load the polarizability
	//if(DEBUG_MDUG>0) std::cout<<"Loading polarizability...\n";
	//loadAlpha<AtomT,MoleculeT>(reader,gauss,mol,std::is_base_of<Alpha,MoleculeT>());
	
	fclose(reader);
}

}

#endif
