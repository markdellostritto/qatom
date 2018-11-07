#include "thole.hpp"

//****************************************************
//THOLE
//****************************************************

//operators

std::ostream& operator<<(std::ostream& out, const Thole& thole){
	out<<"*****************************************************\n";
	out<<"*********************** THOLE ***********************\n";
	out<<"A          = "<<thole.a_<<"\n";
	out<<"ALPHA_R    = "<<thole.alphar_<<"\n";
	out<<"INTER      = "<<thole.inter_<<"\n";
	out<<"IDD        = "<<thole.idd_<<"\n";
	out<<"LIN_SOLVER = "<<thole.linSolver_<<"\n";
	out<<"*********************** THOLE ***********************\n";
	out<<"*****************************************************";
	return out;
}

//member functions

void Thole::defaults(){
	if(DEBUG_THOLE>0) std::cout<<"Thole::defaults():\n";
	a_=1;
	idd_=IDD::Form::UNKNOWN;
	linSolver_=eigen::LIN_SOLVER::FPLU;
	inter_=true;
	alphar_=true;
}

//static functions

Thole& Thole::load(const char* file, Thole& thole){
	if(DEBUG_THOLE>0) std::cout<<"Thole::load(const char*,Thole&):\n";
	//==== local variables ====
	char* input=new char[string::M];
	char* temp=new char[string::M];
	FILE* reader=NULL;
	bool error=false;
	try{
		//==== open the parameter file ====
		reader=fopen(file,"r");
		if(reader==NULL) throw std::invalid_argument("I/O Error: Could not open parameter file.");
		//==== read in the parameters ====
		while(std::fgets(input,string::M,reader)!=NULL){
			string::to_upper(string::trim_right(input,string::COMMENT));
			string::trim_all(string::copy_left(temp,input,"="));
			if(std::strcmp(temp,"A")==0){
				thole.a()=std::atof(std::strpbrk(input,"=")+1);
			} else if(std::strcmp(temp,"INTER")==0){
				thole.inter()=string::boolean(string::trim_all(std::strpbrk(input,"=")+1));
			} else if(std::strcmp(temp,"ALPHA_R")==0){
				thole.alphar()=string::boolean(string::trim_all(std::strpbrk(input,"=")+1));
			} else if(std::strcmp(temp,"IDD")==0){
				thole.idd()=IDD::load(string::trim_all(std::strpbrk(input,"=")+1));
			} else if(std::strcmp(temp,"LIN_SOLVER")==0){
				thole.linSolver()=eigen::LIN_SOLVER::load(string::trim_all(std::strpbrk(input,"=")+1));
			}
		}
		//==== close the file ====
		fclose(reader);
		reader=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Thole::load(const char*,const Thole&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	delete[] input;
	delete[] temp;
	if(error) throw std::runtime_error("Failure to load.");
	else return thole;
}