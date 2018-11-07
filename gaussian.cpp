#include "gaussian.hpp"

namespace GAUSSIAN{

//*************************************************
//Data Formats
//*************************************************

FormatPosn::type FormatPosn::load(const char* str){
	if(std::strcmp(str,"INPUT")==0) return FormatPosn::INPUT;
	else if(std::strcmp(str,"STANDARD")==0) return FormatPosn::STANDARD;
	else return FormatPosn::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, const FormatPosn::type& t){
	if(t==FormatPosn::INPUT) out<<"INPUT";
	else if(t==FormatPosn::STANDARD) out<<"STANDARD";
	else out<<"UNKNOWN";
	return out;
}

FormatChg::type FormatChg::load(const char* str){
	if(std::strcmp(str,"NONE")==0) return FormatChg::NONE;
	else if(std::strcmp(str,"MULLIKEN")==0) return FormatChg::MULLIKEN;
	else if(std::strcmp(str,"ESP")==0) return FormatChg::ESP;
	else if(std::strcmp(str,"HIRSHFELD")==0) return FormatChg::HIRSHFELD;
	else return FormatChg::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, const FormatChg::type& t){
	if(t==FormatChg::NONE) out<<"NONE";
	else if(t==FormatChg::MULLIKEN) out<<"MULLIKEN";
	else if(t==FormatChg::ESP) out<<"ESP";
	else if(t==FormatChg::HIRSHFELD) out<<"HIRSHFELD";
	else out<<"UNKNOWN";
	return out;
}

FormatCalc::type FormatCalc::load(const char* str){
	if(std::strcmp(str,"OPT")==0) return FormatCalc::OPT;
	else if(std::strcmp(str,"SP")==0) return FormatCalc::SP;
	else if(std::strcmp(str,"BOMD")==0) return FormatCalc::BOMD;
	else if(std::strcmp(str,"ADMP")==0) return FormatCalc::ADMP;
	else return FormatCalc::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, const FormatCalc::type& t){
	if(t==FormatCalc::OPT) out<<"OPT";
	else if(t==FormatCalc::BOMD) out<<"BOMD";
	else if(t==FormatCalc::ADMP) out<<"ADMP";
	else out<<"UNKNOWN";
	return out;
}

FormatAlpha::type FormatAlpha::load(const char* str){
	if(std::strcmp(str,"NONE")==0) return FormatAlpha::NONE;
	else if(std::strcmp(str,"DIPOLE")==0) return FormatAlpha::DIPOLE;
	else if(std::strcmp(str,"NUMERIC")==0) return FormatAlpha::NUMERIC;
	else if(std::strcmp(str,"GAMMA")==0) return FormatAlpha::GAMMA;
	else if(std::strcmp(str,"DCSHG")==0) return FormatAlpha::DCSHG;
	else if(std::strcmp(str,"CUBIC")==0) return FormatAlpha::CUBIC;
	else return FormatAlpha::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, const FormatAlpha::type& t){
	if(t==FormatAlpha::NONE) out<<"NONE";
	else if(t==FormatAlpha::DIPOLE) out<<"DIPOLE";
	else if(t==FormatAlpha::NUMERIC) out<<"NUMERIC";
	else if(t==FormatAlpha::GAMMA) out<<"GAMMA";
	else if(t==FormatAlpha::DCSHG) out<<"DCSHG";
	else if(t==FormatAlpha::CUBIC) out<<"CUBIC";
	else out<<"UNKNOWN";
	return out;
}

FormatVersion::type FormatVersion::load(const char* str){
	if(std::strcmp(str,"NONE")==0) return FormatVersion::NONE;
	else if(std::strcmp(str,"g09")==0) return FormatVersion::g09;
	else if(std::strcmp(str,"g16")==0) return FormatVersion::g16;
	else return FormatVersion::UNKNOWN;
}

std::ostream& operator<<(std::ostream& out, const FormatVersion::type& t){
	if(t==FormatVersion::NONE) out<<"NONE";
	else if(t==FormatVersion::g09) out<<"g09";
	else if(t==FormatVersion::g16) out<<"g16";
	else out<<"UNKNOWN";
	return out;
}

//*****************************************************
//FORMAT struct
//*****************************************************

Format& Format::load(const std::vector<std::string>& strlist, Format& format){
	for(unsigned int i=0; i<strlist.size(); ++i){
		if(strlist[i]=="-log"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-log\" option.");
			else format.log=strlist[i+1];
		} else if(strlist[i]=="-com"){
			if(i==strlist.size()-1) throw std::invalid_argument("No file specified for \"-com\" option.");
			else format.com=strlist[i+1];
		} else if(strlist[i]=="-calc"){
			if(i==strlist.size()-1) throw std::invalid_argument("No option specified for \"-calc\" option.");
			format.formatCalc=FormatCalc::load(strlist[i+1].c_str());
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
		} else if(strlist[i]=="-basis"){
			if(i==strlist.size()-1) throw std::invalid_argument("No interval specified for \"-basis\" option.");
			format.basis=strlist[i+1];
		} else if(strlist[i]=="-chem"){
			if(i==strlist.size()-1) throw std::invalid_argument("No interval specified for \"-chem\" option.");
			format.chem=strlist[i+1];
		} else if(strlist[i]=="-job"){
			if(i==strlist.size()-1) throw std::invalid_argument("No interval specified for \"-job\" option.");
			format.job=strlist[i+1];
		} else if(strlist[i]=="-nproc"){
			if(i==strlist.size()-1) throw std::invalid_argument("No interval specified for \"-nproc\" option.");
			format.nProc=std::atoi(strlist[i+1].c_str());
		} else if(strlist[i]=="-nGB"){
			if(i==strlist.size()-1) throw std::invalid_argument("No interval specified for \"-nGB\" option.");
			format.nGB=std::atoi(strlist[i+1].c_str());
		} else if(strlist[i]=="-verbosity"){
			if(i==strlist.size()-1) throw std::invalid_argument("No interval specified for \"-nGB\" option.");
			format.verbosity=strlist[i+1];
		} 
	}
	return format;
}

//*************************************************
//Gaussian File Class
//*************************************************

std::string loadName(FILE* reader){
	char* input=(char*)malloc(sizeof(char)*string::M);
	std::string name="MOL";
	std::rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strpbrk(input,"#")!=NULL){
			fgets(input,string::M,reader);
			fgets(input,string::M,reader);
			name=string::trim_all(fgets(input,string::M,reader));
			break;
		}
	}
	free(input);
	return name;
}

void loadVersion(FILE* reader, GaussFile& gauss){
	if(DEBUG_MDUG>0) std::cout<<"loadVersion(FILE*,const GaussFile&):\n";
	
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//tags
		const char* strG09="Gaussian 09";
		const char* strG16="Gaussian 16";
	
	//read in format
	gauss.formatVersion=FormatVersion::NONE;
	rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		if(std::strstr(input,strG09)!=NULL){gauss.formatVersion=FormatVersion::g09;break;}
		else if(std::strstr(input,strG16)!=NULL){gauss.formatVersion=FormatVersion::g16;break;}
	}
	
	//free local variables
	free(input);
	free(temp);
	
	if(gauss.formatVersion==FormatVersion::NONE) throw std::runtime_error("Could not find Gaussian Version.");
}

void loadVersion(FILE* reader, Format& format){
	if(DEBUG_MDUG>0) std::cout<<"loadVersion(FILE*,const GaussFile&):\n";
	
	/* local function variables */
	//file i/o
		char* input=(char*)malloc(sizeof(char)*string::M);
		char* temp=(char*)malloc(sizeof(char)*string::M);
	//tags
		const char* strG09="Gaussian 09";
		const char* strG16="Gaussian 16";
	
	//read in format
	format.formatVersion=FormatVersion::NONE;
	rewind(reader);
	while(fgets(input,string::M,reader)!=NULL){
		string::trim(string::trim_right(input,":"));
		if(std::strcmp(input,strG09)==0) format.formatVersion=FormatVersion::g09;
		else if(std::strcmp(input,strG16)==0) format.formatVersion=FormatVersion::g16;
	}
	
	//free local variables
	free(input);
	free(temp);
	
	if(format.formatVersion==FormatVersion::NONE) throw std::runtime_error("Could not find Gaussian Version.");
}

void print(const char* fileName, const GaussFile& gauss){
	if(DEBUG_MDUG>0) std::cout<<"print(const char*,const GaussFile):\n";
	
	FILE* writer=fopen(fileName,"w");
	if(writer==NULL) throw std::runtime_error("File could not be opened.");
	
	//print checkpoint file
	fprintf(writer, "%%Chk=%s_%s\n", gauss.name.c_str(), gauss.job.c_str());
	//print the compute info
	fprintf(writer, "%%NProcShared=%i\n", gauss.nProc);
	fprintf(writer, "%%Mem=%iGB\n", gauss.nGB);
	//print the job section
	fprintf(writer, "#P %s/%s %s", gauss.chem.c_str(), gauss.basis.c_str(), gauss.job.c_str());
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
	fprintf(writer,"%s\n\n", gauss.name.c_str());
	
	//print the charge and multiplicity
	fprintf(writer, "%i %i\n", (int)(gauss.charge), (int)(gauss.spinMult));
	//print the atoms
	for(unsigned int i=0; i<gauss.atoms.size(); ++i){
		fprintf(writer, "%s  %12.8f  %12.8f  %12.8f\n",
			gauss.atoms[i].name().c_str(),
			gauss.atoms[i].posn()[0],
			gauss.atoms[i].posn()[1],
			gauss.atoms[i].posn()[2]
		);
	}
	fprintf(writer,"\n");
	
	fclose(writer);
}

}