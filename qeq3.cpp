#include "qeq3.hpp"

//************************************************************
//QEQ3
//************************************************************

//operators

std::ostream& operator<<(std::ostream& out, const QEQ3& qeq){
	out<<"******************************************************\n";
	out<<"************************ QEQ3 ************************\n";
	out<<"K       = "<<qeq.k_<<"\n";
	out<<"Q_TOT   = "<<qeq.qTot()<<"\n";
	out<<"ICC     = "<<qeq.icc()<<"\n";
	out<<"CHI     = "<<qeq.chiType_<<"\n";
	out<<"SCALE_R = "<<qeq.scaler_<<"\n";
	out<<"************************ QEQ3 ************************\n";
	out<<"******************************************************";
}

//member functions

void QEQ3::defaults(){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::defaults():\n";
	qTot_=0;
	k_=1;
	chiType_=PTable::ElectronegativityType::MULLIKEN;
	scaler_=false;
	icc_=ICC::Form::UNKNOWN;
}

void QEQ3::load(const char* paramFile){
	if(DEBUG_QEQ3>0) std::cout<<"QEQ3::load(const char*):\n";
	//local file variables
	FILE* reader=NULL;
	char* input=(char*)malloc(sizeof(char)*string::M);
	char* temp=(char*)malloc(sizeof(char)*string::M);
	
	//open the file
	reader=fopen(paramFile,"r");
	if(reader==NULL) throw std::runtime_error("I/O Error: Could not open parameter file.");
	
	//read in the parameters
	while(fgets(input,string::M,reader)!=NULL){
		string::trim_right(input,string::COMMENT);
		string::copy_left(temp,input,"=");
		string::trim_all(string::to_upper(temp));
		if(std::strcmp(temp,"K")==0){
			k_=std::atof(std::strpbrk(input,"=")+1);
		} else if(std::strcmp(temp,"Q_TOT")==0){
			qTot_=std::atof(std::strpbrk(input,"=")+1);
		} else if(std::strcmp(temp,"ICC")==0){
			icc_=ICC::load(string::trim_all(std::strpbrk(input,"=")+1));
		} else if(std::strcmp(temp,"SCALE_R")==0){
			scaler_=std::atof(std::strpbrk(input,"=")+1);
		} else if(std::strcmp(temp,"ELECTRONEGATIVITY")==0){
			std::strcpy(temp,std::strpbrk(input,"=")+1);
			string::trim_all(string::to_upper(temp));
			if(std::strcmp(temp,"MULLIKEN")==0) chiType_=PTable::ElectronegativityType::MULLIKEN;
			else if(std::strcmp(temp,"PAULING")==0) chiType_=PTable::ElectronegativityType::PAULING;
			else if(std::strcmp(temp,"ALLEN")==0) chiType_=PTable::ElectronegativityType::ALLEN;
			else if(std::strcmp(temp,"HINZE")==0) chiType_=PTable::ElectronegativityType::HINZE;
			else throw std::invalid_argument("Unrecognized electronegativity type.");
		}
	}
	
	//close the file
	fclose(reader);
	reader=NULL;
	
	//free local variables
	free(input);
	free(temp);
	
	//check the parameters
	if(icc_==ICC::Form::UNKNOWN) throw std::runtime_error("Invalid ICC type.");
	if(k_<0) throw std::runtime_error("Invalid interaction parameter.");
	if(chiType_==PTable::ElectronegativityType::UNKNOWN) throw std::runtime_error("Invalid electronegativity type.");
}
