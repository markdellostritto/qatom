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

void Thole::read(const char* file){
	if(DEBUG_THOLE>0) std::cout<<"Thole::read(const char*):\n";
	//==== local variables ====
	char* input=new char[string::M];
	std::vector<std::string> strlist;
	FILE* reader=NULL;
	bool error=false;
	try{
		//==== open the parameter file ====
		reader=fopen(file,"r");
		if(reader==NULL) throw std::invalid_argument("I/O Error: Could not open parameter file.");
		//==== read in the parameters ====
		while(std::fgets(input,string::M,reader)!=NULL){
			string::trim_right(input,string::COMMENT);
			if(string::split(input,string::WS,strlist)==0) continue;
			string::to_upper(strlist.at(0));
			if(strlist.at(0)=="A"){
				a_=std::atof(strlist.at(1).c_str());
			} else if(strlist.at(0)=="INTER"){
				inter_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="ALPHA_R"){
				alphar_=string::boolean(strlist.at(1).c_str());
			} else if(strlist.at(0)=="IDD"){
				idd_=IDD::read(string::to_upper(strlist.at(1)).c_str());
			} else if(strlist.at(0)=="LIN_SOLVER"){
				linSolver_=eigen::LIN_SOLVER::read(string::to_upper(strlist.at(1)).c_str());
			}
		}
		//==== close the file ====
		fclose(reader);
		reader=NULL;
	}catch(std::exception& e){
		std::cout<<"ERROR in Thole::read(const char*,const Thole&):\n";
		std::cout<<e.what()<<"\n";
		error=true;
	}
	delete[] input;
	if(error) throw std::runtime_error("Failure to read.");
}

//alpha - struc

void Thole::init(Structure& struc, const Ewald3D::Dipole& dipole){
	if(DEBUG_THOLE>0) std::cout<<"Thole::init(Structure&,unsigned int,const Ewald3D::Dipole&):\n";
	
	//set the utility matrices
	if(DEBUG_THOLE>0) std::cout<<"Setting utility matrices...\n";
	A_=Eigen::MatrixXd::Zero(3*struc.nAtoms(),3*struc.nAtoms());
	alphaC_.resize(3*struc.nAtoms(),3);
	identityC_.resize(3*struc.nAtoms(),3);
	
	r0_.resize(struc.nAtoms());
	r_.resize(struc.nAtoms());
	drAtom_.resize(struc.nAtoms());
	
	dipole.interMatSelfR(selfInterMat_,100);
	
	//populate the solution vector
	if(DEBUG_THOLE>1) std::cout<<"Populating solution vector...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) identityC_.block<3,3>(i*3,0)=Eigen::Matrix3d::Identity();
}

void Thole::alpha(Structure& struc, const Ewald3D::Dipole& dipole){
	if(DEBUG_THOLE>0) std::cout<<"Thole::alpha(Structure&,unsigned int,const Ewald3D::Dipole&):\n";
	//set the utility matrices
	A_=Eigen::MatrixXd::Zero(3*struc.nAtoms(),3*struc.nAtoms());
	alphaC_.resize(3*struc.nAtoms(),3);
	identityC_.resize(3*struc.nAtoms(),3);
	const double ke=units::consts::ke();
	//units
	double rscale=0.0;
	if(units::consts::system()==units::System::AU) rscale=units::BOHRpANG;
	else if(units::consts::system()==units::System::METAL) rscale=1.0;
	else throw std::runtime_error("Invalid units.");
	
	r0_.resize(struc.nAtoms());
	r_.resize(struc.nAtoms());
	drAtom_.resize(struc.nAtoms());
	
	Eigen::Vector3d rVec,r00,rNew;
	Eigen::Matrix3d selfInterMat;
	dipole.interMatSelfR(selfInterMat,50);
	
	//set the gas-phase atomic radius
	if(DEBUG_THOLE>1) std::cout<<"Setting gas-phase atomic radii.\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) r0_[i]=PTable::covalentRadius(struc.an(i))*rscale;
	for(unsigned int i=0; i<struc.nAtoms(); ++i) r_[i]=PTable::covalentRadius(struc.an(i))*rscale;
	//calculate the atom-in-molecule radius
	if(DEBUG_THOLE>1) std::cout<<"Calculating the atom-in-molecule radius...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		double rMin=100;
		for(unsigned int j=0; j<struc.nAtoms(); ++j){
			if(i==j) continue;
			Cell::diff(struc.posn(j),struc.posn(i),rVec,struc.cell().R(),struc.cell().RInv());
			double dr=rVec.norm();
			double rTemp=r0_[i]-0.5*(r0_[i]+r0_[j]-dr);
			if(rTemp<rMin){
				r_[i]=rTemp;
				rMin=rTemp;
				rVec/=rVec.norm();
				r00.noalias()=r0_[i]*rVec;
				rNew.noalias()=r_[i]*rVec;
				drAtom_[i][0]=std::fabs(rNew[0]/r00[0]);
				drAtom_[i][1]=std::fabs(rNew[1]/r00[1]);
				drAtom_[i][2]=std::fabs(rNew[2]/r00[2]);
				if(rVec[0]==0 && r00[0]==0) drAtom_[i][0]=1;
				if(rVec[1]==0 && r00[1]==0) drAtom_[i][1]=1;
				if(rVec[2]==0 && r00[2]==0) drAtom_[i][2]=1;
			}
		}
	}
	//print the radius
	if(DEBUG_THOLE>2){
		std::cout<<"RADIUS:\n";
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			std::cout<<struc.name(i)<<" "<<r0_[i]<<" "<<r_[i]<<" "<<drAtom_[i].transpose()<<"\n";
		}
	}
	//print the radius
	if(DEBUG_THOLE>2){
		std::cout<<"ALPHA:\n";
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			Eigen::Matrix3d alpha=struc.alpha(i);
			alpha(0,0)*=drAtom_[i][0];
			alpha(1,1)*=drAtom_[i][1];
			alpha(2,2)*=drAtom_[i][2];
			std::cout<<struc.name(i)<<" "<<alpha<<"\n";
		}
	}
	
	//populate the solution vector
	if(DEBUG_THOLE>1) std::cout<<"Populating solution vector...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) identityC_.block<3,3>(i*3,0)=Eigen::Matrix3d::Identity();
	
	//populate the diagonal of the A-matrix
	if(DEBUG_THOLE>1) std::cout<<"Populating diagonal of A-matrix...\n";
	if(alphar_){
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			struc.alpha(i)(0,0)*=drAtom_[i][0];
			struc.alpha(i)(1,1)*=drAtom_[i][1];
			struc.alpha(i)(2,2)*=drAtom_[i][2];
		}
	}
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		A_.block<3,3>(i*3,i*3).noalias()=struc.alpha(i).inverse()*ke;
	}
	
	//calculate the interaction matrices
	if(DEBUG_THOLE>1) std::cout<<"Calculating the off-diagonal elements...\n";
	if(inter_){
		if(DEBUG_THOLE>1) std::cout<<"Calculating self-interactions...\n";
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			A_.block<3,3>(i*3,i*3).noalias()-=ke*selfInterMat;
		}
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			for(unsigned int j=i+1; j<struc.nAtoms(); ++j){
				//calculate distance
				Cell::diff(struc.posn(j),struc.posn(i),rVec,struc.cell().R(),struc.cell().RInv());
				//dipole interaction - periodic
				A_.block<3,3>(j*3,i*3).noalias()=-1*ke*dipole.interMat(rVec,interMat_);
				//dipole interaction - short - ideal
				A_.block<3,3>(j*3,i*3).noalias()+=ke*IDD::iTensor[IDD::Form::IDEAL](rVec,interMat_,1.0);
				//dipole interaction - short - modified
				A_.block<3,3>(j*3,i*3).noalias()=-1*ke*IDD::iTensor[idd_](rVec,interMat_,
					a_*IDD::scale[idd_](struc.alpha(j).trace()/3.0,struc.alpha(i).trace()/3.0)
				);
			}
		}
	}
	
	//set the total matrix
	A_=A_.selfadjointView<Eigen::Lower>();
	
	//calculate the effective polarizabilities
	if(DEBUG_THOLE>1) std::cout<<"Calculating effective polarizabilities...\n";
	switch(linSolver_){
		case eigen::LIN_SOLVER::LLT: alphaC_.noalias()=A_.llt().solve(identityC_); break;
		case eigen::LIN_SOLVER::LDLT: alphaC_.noalias()=A_.ldlt().solve(identityC_); break;
		case eigen::LIN_SOLVER::HQR: alphaC_.noalias()=A_.householderQr().solve(identityC_); break;
		case eigen::LIN_SOLVER::CPHQR: alphaC_.noalias()=A_.colPivHouseholderQr().solve(identityC_); break;
		case eigen::LIN_SOLVER::PPLU: alphaC_.noalias()=A_.partialPivLu().solve(identityC_); break;
		case eigen::LIN_SOLVER::FPLU: alphaC_.noalias()=A_.fullPivLu().solve(identityC_); break;
		default: A_.fullPivLu().solve(identityC_); break;
	};
	
	//record the polarizability
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		struc.alpha(i).noalias()=alphaC_.block<3,3>(i*3,0);
	}
}

void Thole::alpha_cont(Structure& struc, const Ewald3D::Dipole& dipole){
	if(DEBUG_THOLE>0) std::cout<<"Thole::alpha_cont(Structure&,unsigned int,Ewald3D::Dipole&):\n";
	//local variables
	//scaling
	const double ke=units::consts::ke();
	//units
	double rscale=0.0;
	if(units::consts::system()==units::System::AU) rscale=units::BOHRpANG;
	else if(units::consts::system()==units::System::METAL) rscale=1.0;
	else throw std::runtime_error("Invalid units.");
	//utility vectors
	Eigen::Vector3d rVec,r00,rNew;
	
	//set the gas-phase atomic radius
	if(DEBUG_THOLE>1) std::cout<<"Setting gas-phase atomic radii.\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) r0_[i]=PTable::covalentRadius(struc.an(i))*rscale;
	for(unsigned int i=0; i<struc.nAtoms(); ++i) r_[i]=PTable::covalentRadius(struc.an(i))*rscale;
	//calculate the atom-in-molecule radius
	if(DEBUG_THOLE>1) std::cout<<"Calculating the atom-in-molecule radius...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		double rMin=100;
		for(unsigned int j=0; j<struc.nAtoms(); ++j){
			if(i==j) continue;
			Cell::diff(struc.posn(j),struc.posn(i),rVec,struc.cell().R(),struc.cell().RInv());
			const double dr=rVec.norm();
			const double rTemp=r0_[i]-0.5*(r0_[i]+r0_[j]-dr);
			if(rTemp<rMin){
				r_[i]=rTemp;
				rMin=rTemp;
				rVec/=rVec.norm();
				r00.noalias()=r0_[i]*rVec;
				rNew.noalias()=r_[i]*rVec;
				drAtom_[i][0]=std::fabs(rNew[0]/r00[0]);
				drAtom_[i][1]=std::fabs(rNew[1]/r00[1]);
				drAtom_[i][2]=std::fabs(rNew[2]/r00[2]);
				if(rVec[0]==0 && r00[0]==0) drAtom_[i][0]=1;
				if(rVec[1]==0 && r00[1]==0) drAtom_[i][1]=1;
				if(rVec[2]==0 && r00[2]==0) drAtom_[i][2]=1;
			}
		}
	}
	//print the radius
	if(DEBUG_THOLE>2){
		std::cout<<"RADIUS:\n";
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			std::cout<<struc.name(i)<<" "<<r0_[i]<<" "<<r_[i]<<" "<<drAtom_[i].transpose()<<"\n";
		}
	}
	//print the radius
	if(DEBUG_THOLE>2){
		std::cout<<"ALPHA:\n";
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			Eigen::Matrix3d alpha=struc.alpha(i);
			alpha(0,0)*=drAtom_[i][0];
			alpha(1,1)*=drAtom_[i][1];
			alpha(2,2)*=drAtom_[i][2];
			std::cout<<struc.name(i)<<" "<<alpha<<"\n";
		}
	}
	
	//populate the diagonal of the A-matrix
	if(DEBUG_THOLE>1) std::cout<<"Populating diagonal of A-matrix...\n";
	if(alphar_){
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			struc.alpha(i)(0,0)*=drAtom_[i][0];
			struc.alpha(i)(1,1)*=drAtom_[i][1];
			struc.alpha(i)(2,2)*=drAtom_[i][2];
		}
	}
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		A_.block<3,3>(i*3,i*3).noalias()=struc.alpha(i).inverse()*ke;
	}
	
	//calculate the interaction matrices
	if(inter_){
		if(DEBUG_THOLE>1) std::cout<<"Calculating self-interactions...\n";
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			A_.block<3,3>(i*3,i*3).noalias()-=ke*selfInterMat_;
		}
		if(DEBUG_THOLE>1) std::cout<<"Calculating the off-diagonal elements...\n";
		for(unsigned int i=0; i<struc.nAtoms(); ++i){
			for(unsigned int j=i+1; j<struc.nAtoms(); ++j){
				//calculate distance
				Cell::diff(struc.posn(j),struc.posn(i),rVec,struc.cell().R(),struc.cell().RInv());
				//dipole interaction - periodic
				A_.block<3,3>(j*3,i*3).noalias()=-1*ke*dipole.interMat(rVec,interMat_);
				//dipole interaction - short - ideal
				A_.block<3,3>(j*3,i*3).noalias()+=ke*IDD::iTensor[IDD::Form::IDEAL](rVec,interMat_,1.0);
				//dipole interaction - short - modified
				A_.block<3,3>(j*3,i*3).noalias()=-1*ke*IDD::iTensor[idd_](rVec,interMat_,
					a_*IDD::scale[idd_](struc.alpha(j).trace()/3.0,struc.alpha(i).trace()/3.0)
				);
			}
		}
	}
	
	//set the total matrix
	A_=A_.selfadjointView<Eigen::Lower>();
	
	//calculate the effective polarizabilities
	if(DEBUG_THOLE>1) std::cout<<"Calculating effective polarizabilities...\n";
	switch(linSolver_){
		case eigen::LIN_SOLVER::LLT:   alphaC_.noalias()=A_.llt().solve(identityC_); break;
		case eigen::LIN_SOLVER::LDLT:  alphaC_.noalias()=A_.ldlt().solve(identityC_); break;
		case eigen::LIN_SOLVER::HQR:   alphaC_.noalias()=A_.householderQr().solve(identityC_); break;
		case eigen::LIN_SOLVER::CPHQR: alphaC_.noalias()=A_.colPivHouseholderQr().solve(identityC_); break;
		case eigen::LIN_SOLVER::PPLU:  alphaC_.noalias()=A_.partialPivLu().solve(identityC_); break;
		case eigen::LIN_SOLVER::FPLU:  alphaC_.noalias()=A_.fullPivLu().solve(identityC_); break;
		default: A_.fullPivLu().solve(identityC_); break;
	};
	
	//record the polarizability
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		struc.alpha(i).noalias()=alphaC_.block<3,3>(i*3,0);
	}
}
