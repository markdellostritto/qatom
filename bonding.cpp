#include "bonding.hpp"

//constructors/destructors

Bonding::Bonding(const Bonding& b){
	nSpecies_=b.nSpecies();
	bondLengths_=b.bondLengths();
	bondOrders_=b.bondOrders();
	hBondAngle_=b.hBondAngle();
	hBondLengths_=b.hBondLengths();
}

//operators

Bonding& Bonding::operator=(const Bonding& b){
	nSpecies_=b.nSpecies();
	bondLengths_=b.bondLengths();
	bondOrders_=b.bondOrders();
	hBondAngle_=b.hBondAngle();
	hBondLengths_=b.hBondLengths();
	return *this;
}

bool operator==(const Bonding& b1, const Bonding b2){
	if(b1.nSpecies()!=b2.nSpecies()) return false;
	else{
		for(unsigned int i=0; i<b1.nSpecies(); ++i){
			for(unsigned int j=0; j<b2.nSpecies(); ++j){
				if(std::fabs(b1.bondLength(i,j)-b2.bondLength(i,j))>num_const::ZERO) return false;
			}
		}
		for(unsigned int i=0; i<b1.nSpecies(); ++i){
			if(b1.bondOrder(i)!=b2.bondOrder(i)) return false;
		}
		if(std::fabs(b1.hBondAngle()-b2.hBondAngle())>num_const::ZERO) return false;
		for(unsigned int i=0; i<b1.nSpecies(); ++i){
			if(std::fabs(b1.hBondLength(i)-b2.hBondLength(i))>num_const::ZERO) return false;
		}
		return true;
	}
}

std::ostream& operator<<(std::ostream& out, const Bonding& b){
	std::cout<<"****************************************************\n";
	std::cout<<"********************* BONDING *********************\n";
	if(b.nSpecies()>0){
		out<<"Max Bond Orders = ";
		for(unsigned int i=0; i<b.nSpecies(); ++i){
			out<<b.bondOrder(i)<<" ";
		}
	}
	if(b.bondLengths().size()>0){
		out<<"\n";
		out<<"Covalent Bond Lengths = \n";
		for(unsigned int i=0; i<b.nSpecies()-1; ++i){
			for(unsigned int j=0; j<b.nSpecies(); ++j){
				out<<b.bondLength(i,j)<<" ";
			}
			out<<"\n";
		}
		for(unsigned int j=0; j<b.nSpecies(); ++j){
			out<<b.bondLength(b.nSpecies()-1,j)<<" ";
		}
	}
	if(b.hBondLengths().size()>0){
		out<<"\n";
		out<<"H-Bond Lengths = ";
		for(unsigned int i=0; i<b.nSpecies(); ++i){
			out<<b.hBondLength(i)<<" ";
		}
	}
	if(b.hBondAngle()!=0) out<<"\nH-Bond Angle = "<<b.hBondAngle()*180.0/num_const::PI<<"\n";
	std::cout<<"********************* BONDING *********************\n";
	std::cout<<"****************************************************";
	return out;
}

//member functions

double& Bonding::bondLength(const std::string& name1, const std::string& name2){
	int i=-1,j=-1;
	for(unsigned int n=0; n<atoms_.size(); ++n){
		if(name1==atoms_[n]){i=n; break;}
	}
	for(unsigned int n=0; n<atoms_.size(); ++n){
		if(name2==atoms_[n]){j=n; break;}
	}
	if(i<0) throw std::runtime_error(std::string("Bonding has no atom ")+name1);
	if(j<0) throw std::runtime_error(std::string("Bonding has no atom ")+name2);
	return bondLengths_[i][j];
}

const double& Bonding::bondLength(const std::string& name1, const std::string& name2)const{
	int i=-1,j=-1;
	for(unsigned int n=0; n<atoms_.size(); ++n){
		if(name1==atoms_[n]){i=n; break;}
	}
	for(unsigned int n=0; n<atoms_.size(); ++n){
		if(name2==atoms_[n]){j=n; break;}
	}
	if(i<0) throw std::runtime_error(std::string("Bonding has no atom ")+name1);
	if(j<0) throw std::runtime_error(std::string("Bonding has no atom ")+name2);
	return bondLengths_[i][j];
}
	
void Bonding::clear(){
	nSpecies_=0;
	bondLengths_.clear();
	bondOrders_.clear();
	hBondAngle_=0;
	hBondLengths_.clear();
}

void Bonding::resize(unsigned int nSpecies){
	nSpecies_=nSpecies;
	bondLengths_.resize(nSpecies,std::vector<double>(nSpecies,0));
	bondOrders_.resize(nSpecies,0);
	hBondLengths_.resize(nSpecies,0);
}

void Bonding::loadBondParams(const char* fileName, const std::vector<std::string>& names){
	if(DEBUG_BOND>0) std::cout<<"Bonding::loadBondParams(const char*,const std::vector<std::string>&):\n";
	//local function variables
	FILE* reader=NULL;
	char* input=(char*)malloc(sizeof(char)*string::M);
	std::vector<std::string> strlist;
	bool error=false;
	
	try{
		//open the file
		if(DEBUG_BOND>0) std::cout<<"Opening the parameter file...\n";
		reader=fopen(fileName,"r");
		if(reader==NULL) throw std::runtime_error("I/O ERROR: Could not open parameter file.\n");
		
		//resize
		if(DEBUG_BOND>0) std::cout<<"Resizing object...\n";
		resize(names.size());
		
		//read in the parameters
		if(DEBUG_BOND>0) std::cout<<"Reading in parameters...\n";
		while(fgets(input,string::M,reader)!=NULL){
			//split the string into tokens
			if(string::empty(input)) continue;
			strlist.resize(string::substrN(input,string::WS));
			strlist[0]=std::string(std::strtok(input,string::WS));
			for(unsigned int i=1; i<strlist.size(); ++i) strlist[i]=std::string(std::strtok(NULL,string::WS));
			//check for the parameter type
			if(strlist[0]=="bond"){
				if(strlist[1]=="length"){
					if(strlist.size()!=5) throw std::runtime_error("Invalid bond length format.");
					int i=-1,j=-1;
					for(unsigned int n=0; n<names.size(); ++n) if(names[n]==strlist[2]) i=n;
					for(unsigned int n=0; n<names.size(); ++n) if(names[n]==strlist[3]) j=n;
					if(i<0) throw std::runtime_error(std::string("No atom \"")+strlist[2]+std::string("\" found."));
					if(j<0) throw std::runtime_error(std::string("No atom \"")+strlist[3]+std::string("\" found."));
					bondLengths_[i][j]=std::atof(strlist[4].c_str());
					bondLengths_[j][i]=bondLengths_[i][j];
				} else if(strlist[1]=="order"){
					if(strlist.size()!=4) throw std::runtime_error("Invalid bond length format.");
					int i=-1;
					for(unsigned int n=0; n<names.size(); ++n) if(names[n]==strlist[2]) i=n;
					if(i<0) throw std::runtime_error(std::string("No atom \"")+strlist[2]+std::string("\" found."));
					bondOrders_[i]=std::atoi(strlist[3].c_str());
				}
			} else if(strlist[0]=="hbond"){
				if(strlist[1]=="length"){
					if(strlist.size()!=4) throw std::runtime_error("Invalid bond length format.");
					int i=-1;
					for(unsigned int n=0; n<names.size(); ++n) if(names[n]==strlist[2]) i=n;
					if(i<0) throw std::runtime_error(std::string("No atom \"")+strlist[2]+std::string("\" found."));
					hBondLengths_[i]=std::atoi(strlist[3].c_str());
				} else if(strlist[1]=="angle") hBondAngle_=std::atof(strlist[2].c_str());
			}
		}
		
		atoms_=names;
		
		//close the file
		if(DEBUG_BOND>0) std::cout<<"Closing parameter file...\n";
		fclose(reader);
		reader=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Bonding::loadParameters(const char*,const std::vector<std::string>&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	
	if(error) throw std::runtime_error("Could not load bond parameters.");
	else if(DEBUG_BOND>0) std::cout<<"Bonding parameters loaded.\n";
}