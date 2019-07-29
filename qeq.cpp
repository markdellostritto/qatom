#include "qeq.hpp"

//************************************************************
//QEQ
//************************************************************

//operators

std::ostream& operator<<(std::ostream& out, const QEQ& qeq){
	out<<"******************************************************\n";
	out<<"************************ QEQ ************************\n";
	out<<"K       = "<<qeq.k_<<"\n";
	out<<"Q_TOT   = "<<qeq.qTot()<<"\n";
	out<<"ICC     = "<<qeq.icc()<<"\n";
	out<<"CHI     = "<<qeq.chiType_<<"\n";
	out<<"SCALE_R = "<<qeq.scaler_<<"\n";
	out<<"************************ QEQ ************************\n";
	out<<"******************************************************";
}

//member functions

void QEQ::defaults(){
	if(DEBUG_QEQ>0) std::cout<<"QEQ::defaults():\n";
	qTot_=0;
	k_=1;
	chiType_=PTable::ElectronegativityType::MULLIKEN;
	scaler_=false;
	icc_=ICC::Form::UNKNOWN;
}

void QEQ::read(const char* paramFile){
	if(DEBUG_QEQ>0) std::cout<<"QEQ::read(const char*):\n";
	//local file variables
	FILE* reader=NULL;
	std::vector<std::string> strlist;
	char* input=new char[string::M];
	
	//open the file
	reader=fopen(paramFile,"r");
	if(reader==NULL) throw std::runtime_error("I/O Error: Could not open parameter file.");
	
	//read in the parameters
	while(fgets(input,string::M,reader)!=NULL){
		string::trim_right(input,string::COMMENT);
		if(string::split(input,string::WS,strlist)==0) continue;
		string::to_upper(strlist.at(0));
		if(strlist.at(0)=="K"){
			k_=std::atof(strlist.at(1).c_str());
		} else if(strlist.at(0)=="Q_TOT"){
			qTot_=std::atof(strlist.at(1).c_str());
		} else if(strlist.at(0)=="ICC"){
			icc_=ICC::read(string::to_upper(strlist.at(1)).c_str());
		} else if(strlist.at(0)=="SCALE_R"){
			scaler_=std::atof(strlist.at(1).c_str());
		} else if(strlist.at(0)=="CHI"){
			chiType_=PTable::ElectronegativityType::read(string::to_upper(strlist.at(1)).c_str());
		}
	}
	
	//close the file
	fclose(reader);
	reader=NULL;
	
	//free local variables
	delete[] input;
	
	//check the parameters
	if(icc_==ICC::Form::UNKNOWN) throw std::runtime_error("Invalid ICC type.");
	if(k_<0) throw std::runtime_error("Invalid interaction parameter.");
	if(chiType_==PTable::ElectronegativityType::UNKNOWN) throw std::runtime_error("Invalid electronegativity type.");
}

//charge transfer

void QEQ::qt(Structure& struc, const Ewald3D::Coulomb& ewald){
	if(DEBUG_QEQ>0) std::cout<<"QEQ::qt(Structure&,const Ewald3D::Coulomb&):\n";
	//local variables
	const double ke=units::consts::ke();
	Eigen::Vector3d dr;
	
	//resize matrices
	if(DEBUG_QEQ>0) std::cout<<"Resizing matrices...\n";
	J_=Eigen::MatrixXd::Constant(struc.nAtoms(),struc.nAtoms(),0);
	A_.resize(struc.nAtoms(),struc.nAtoms());
	chi_.resize(struc.nAtoms());
	
	//local variables
	double cscale=0.0,escale=0.0;
	if(units::consts::system()==units::System::AU){
		cscale=1.0;
		escale=units::HARTREEpEV;
	}
	else if(units::consts::system()==units::System::METAL){
		cscale=units::ANGpBOHR;
		escale=1.0;
	}
	else throw std::runtime_error("Invalid units.");
	
	//calculate atom info
	if(DEBUG_QEQ>0) std::cout<<"Calculating atom info...\n";
	if(DEBUG_QEQ>1) std::cout<<"Element AN CHI J00 C\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		//electronegativity
		if(chiType_==PTable::ElectronegativityType::PAULING) chi_[i]=PTable::electronegativityPauling(struc.an(i));
		else if(chiType_==PTable::ElectronegativityType::ALLEN) chi_[i]=PTable::electronegativityAllen(struc.an(i));
		else if(chiType_==PTable::ElectronegativityType::MULLIKEN) chi_[i]=PTable::electronegativityMulliken(struc.an(i));
		else if(chiType_==PTable::ElectronegativityType::HINZE) chi_[i]=PTable::electronegativityHinze(struc.an(i),0);
		if(chiType_==PTable::ElectronegativityType::MULLIKEN && struc.an(i)==1) chi_[i]=4.528;
		chi_[i]*=escale;
		//idempotential
		if(chiType_==PTable::ElectronegativityType::HINZE) J_(i,i)=PTable::idempotentialHinze(struc.an(i),0);
		else J_(i,i)=PTable::idempotential(struc.an(i));
		J_(i,i)*=escale;
		if(DEBUG_QEQ>1){
			std::cout<<struc.name(i)<<" ";
			std::cout<<struc.an(i)<<" ";
			std::cout<<chi_[i]<<" ";
			std::cout<<J_(i,i)<<"\n";
		}
	}
	
	//set the idempotential
	if(DEBUG_QEQ>0) std::cout<<"Setting the idempotential...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) J_(i,i)+=ewald.phiSelf()*ke;
	//calculate coulomb integrals
	if(DEBUG_QEQ>0) std::cout<<"Calculating Coulomb Integrals...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		for(unsigned int j=i+1; j<struc.nAtoms(); ++j){
			Cell::diff(struc.posn(i),struc.posn(j),dr,struc.cell().R(),struc.cell().RInv());
			const double a=k_*ICC::scale[icc_](J_(i,i),J_(j,j));
			const double drn=dr.norm();
			J_(j,i)=ke*(
				ICC::iTensor[icc_](drn,a)//short range
				+ewald.phi(struc,dr)-1.0/drn//long range, short range removed
			);
		}
	}
	//calculate the operator matrix
	A_=J_.selfadjointView<Eigen::Lower>();
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		A_(0,i)=1;
		for(unsigned int j=1; j<struc.nAtoms(); ++j){
			A_(j,i)-=J_(i,0);
		}
	}
	//calculate the solution vector
	b_=Eigen::VectorXd::Constant(struc.nAtoms(),chi_[0]);
	b_.noalias()-=chi_;
	b_[0]=qTot_;
	//solve the linear equations
	if(DEBUG_QEQ>0) std::cout<<"Solving the linear equations...\n";
	x_.noalias()=A_.partialPivLu().solve(b_);
	//set the atomic charges
	if(DEBUG_QEQ>0) std::cout<<"Setting the charge...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) struc.charge(i)=x_[i];
}

void QEQ::qt_jZero(Structure& struc, const Ewald3D::Coulomb& ewald){
	if(DEBUG_QEQ>0) std::cout<<"QEQ::qt_jZero(Structure&,const Ewald3D::Coulomb&):\n";
	//local variables
	const double ke=units::consts::ke();
	Eigen::Vector3d dr;
	
	//resize matrices
	if(DEBUG_QEQ>0) std::cout<<"Resizing matrices...\n";
	J_=Eigen::MatrixXd::Constant(struc.nAtoms(),struc.nAtoms(),0);
	A_.resize(struc.nAtoms(),struc.nAtoms());
	chi_.resize(struc.nAtoms());
	
	//local variables
	double cscale=0.0,escale=0.0;
	if(units::consts::system()==units::System::AU){
		cscale=1.0;
		escale=units::HARTREEpEV;
	}
	else if(units::consts::system()==units::System::METAL){
		cscale=units::ANGpBOHR;
		escale=1.0;
	}
	else throw std::runtime_error("Invalid units.");
	
	//calculate atom info
	if(DEBUG_QEQ>0) std::cout<<"Calculating atom info...\n";
	if(DEBUG_QEQ>1) std::cout<<"Element AN CHI J00 C\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		//electronegativity
		if(chiType_==PTable::ElectronegativityType::PAULING) chi_[i]=PTable::electronegativityPauling(struc.an(i));
		else if(chiType_==PTable::ElectronegativityType::ALLEN) chi_[i]=PTable::electronegativityAllen(struc.an(i));
		else if(chiType_==PTable::ElectronegativityType::MULLIKEN) chi_[i]=PTable::electronegativityMulliken(struc.an(i));
		else if(chiType_==PTable::ElectronegativityType::HINZE) chi_[i]=PTable::electronegativityHinze(struc.an(i),0);
		if(chiType_==PTable::ElectronegativityType::MULLIKEN && struc.an(i)==1) chi_[i]=4.528;
		chi_[i]*=escale;
		//idempotential
		J_(i,i)=struc.jzero(i);
		if(DEBUG_QEQ>1){
			std::cout<<struc.name(i)<<" ";
			std::cout<<struc.an(i)<<" ";
			std::cout<<chi_[i]<<" ";
			std::cout<<J_(i,i)<<"\n";
		}
	}
	
	//set the idempotential
	if(DEBUG_QEQ>0) std::cout<<"Setting the idempotential...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) J_(i,i)+=ewald.phiSelf()*ke;
	//calculate coulomb integrals
	if(DEBUG_QEQ>0) std::cout<<"Calculating Coulomb Integrals...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		for(unsigned int j=i+1; j<struc.nAtoms(); ++j){
			Cell::diff(struc.posn(i),struc.posn(j),dr,struc.cell().R(),struc.cell().RInv());
			const double a=k_*ICC::scale[icc_](J_(i,i),J_(j,j));
			const double drn=dr.norm();
			J_(j,i)=ke*(
				ICC::iTensor[icc_](drn,a)//short range
				+ewald.phi(struc,dr)-1.0/drn//long range, short range removed
			);
		}
	}
	//calculate the operator matrix
	A_=J_.selfadjointView<Eigen::Lower>();
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		A_(0,i)=1;
		for(unsigned int j=1; j<struc.nAtoms(); ++j){
			A_(j,i)-=J_(i,0);
		}
	}
	//calculate the solution vector
	b_=Eigen::VectorXd::Constant(struc.nAtoms(),chi_[0]);
	b_.noalias()-=chi_;
	b_[0]=qTot_;
	//solve the linear equations
	if(DEBUG_QEQ>0) std::cout<<"Solving the linear equations...\n";
	x_.noalias()=A_.partialPivLu().solve(b_);
	//set the atomic charges
	if(DEBUG_QEQ>0) std::cout<<"Setting the charge...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) struc.charge(i)=x_[i];

}

void QEQ::init(Structure& struc){
	if(DEBUG_QEQ>0) std::cout<<"QEQ::init(Structure&):\n";
	//local function variables
	Eigen::Vector3d R;
	
	//resize matrices
	if(DEBUG_QEQ>0) std::cout<<"Resizing matrices...\n";
	J_=Eigen::MatrixXd::Constant(struc.nAtoms(),struc.nAtoms(),0);
	A_.resize(struc.nAtoms(),struc.nAtoms());
	chi_.resize(struc.nAtoms());
	x_.resize(struc.nAtoms());
	I_.resize(struc.nAtoms());
	
	//local variables
	double cscale=0.0,escale=0.0;
	if(units::consts::system()==units::System::AU){
		cscale=1.0;
		escale=units::HARTREEpEV;
	}
	else if(units::consts::system()==units::System::METAL){
		cscale=units::ANGpBOHR;
		escale=1.0;
	}
	else throw std::runtime_error("Invalid units.");
	
	//calculate atom info
	if(DEBUG_QEQ>0) std::cout<<"Calculating atom info...\n";
	if(DEBUG_QEQ>-1) std::cout<<"Element AN CHI J00 C\n";
	unsigned int cc=0;
	for(unsigned int n=0; n<struc.nSpecies(); ++n){
		unsigned int an=PTable::an(struc.atomNames(n).c_str());
		for(unsigned int m=0; m<struc.nAtoms(n); ++m){
			//electronegativity
			if(chiType_==PTable::ElectronegativityType::PAULING) chi_[cc]=PTable::electronegativityPauling(an);
			else if(chiType_==PTable::ElectronegativityType::ALLEN) chi_[cc]=PTable::electronegativityAllen(an);
			else if(chiType_==PTable::ElectronegativityType::MULLIKEN) chi_[cc]=PTable::electronegativityMulliken(an);
			else if(chiType_==PTable::ElectronegativityType::HINZE) chi_[cc]=PTable::electronegativityHinze(an,0);
			if(chiType_==PTable::ElectronegativityType::MULLIKEN && an==1) chi_[cc]=4.528;
			chi_[cc]*=escale;
			//idempotential
			if(chiType_==PTable::ElectronegativityType::HINZE) I_(cc)=PTable::idempotentialHinze(an,0);
			if(chiType_==PTable::ElectronegativityType::MULLIKEN && an==1) I_(cc)=13.8904;
			else I_(cc)=PTable::idempotential(an);
			I_(cc)*=escale;
			if(DEBUG_QEQ>-1){
				std::cout<<struc.atomNames(n)<<" ";
				std::cout<<an<<" ";
				std::cout<<chi_[cc]<<" ";
				std::cout<<I_(cc)<<"\n";
			}
			++cc;
		}
	}
	
	//calculate the solution vector
	if(DEBUG_QEQ>0) std::cout<<"Calculating the solution vector...\n";
	b_=Eigen::VectorXd::Constant(struc.nAtoms(),chi_[0]);
	b_.noalias()-=chi_;
	b_[0]=qTot_;
}

void QEQ::qt_cont(Structure& struc, const Ewald3D::Coulomb& ewald){
	if(DEBUG_QEQ>0) std::cout<<"QEQ::qt_cont(Structure&,const Ewald3D::Coulomb&):\n";
	//local variables
	const double ke=units::consts::ke();
	Eigen::Vector3d dr;
	
	//set the idempotential
	if(DEBUG_QEQ>0) std::cout<<"Setting the idempotential...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) J_(i,i)=I_[i]+ewald.phiSelf()*ke;
	//calculate coulomb integrals
	if(DEBUG_QEQ>0) std::cout<<"Calculating Coulomb Integrals...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		const double Jii=J_(i,i);
		for(unsigned int j=i+1; j<struc.nAtoms(); ++j){
			Cell::diff(struc.posn(i),struc.posn(j),dr,struc.cell().R(),struc.cell().RInv());
			const double drn=dr.norm();
			J_(j,i)=ke*(
				ICC::iTensor[icc_](drn,k_*ICC::scale[icc_](Jii,J_(j,j)))//short range
				+ewald.phi(struc,dr)-1.0/drn//long range, short range removed
			);
		}
	}
	//calculate the operator matrix
	if(DEBUG_QEQ>0) std::cout<<"Calculating the operator matrix...\n";
	A_=J_.selfadjointView<Eigen::Lower>();
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		A_(0,i)=1;
		for(unsigned int j=1; j<struc.nAtoms(); ++j){
			A_(j,i)-=J_(i,0);
		}
	}
	//solve the linear equations
	if(DEBUG_QEQ>0) std::cout<<"Solving the linear equations...\n";
	x_.noalias()=A_.partialPivLu().solve(b_);
	//set the atomic charges
	if(DEBUG_QEQ>0) std::cout<<"Setting the charge...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) struc.charge(i)=x_[i];
}

void QEQ::qt_jZero_cont(Structure& struc, const Ewald3D::Coulomb& ewald){
	if(DEBUG_QEQ>0) std::cout<<"QEQ::qt_jZero_cont(Structure&,const Ewald3D::Coulomb&):\n";
	//local variables
	const double ke=units::consts::ke();
	Eigen::Vector3d dr;
	
	//set the idempotential
	if(DEBUG_QEQ>0) std::cout<<"Setting the idempotential...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) J_(i,i)=struc.jzero(i)+ewald.phiSelf()*ke;
	//calculate coulomb integrals
	if(DEBUG_QEQ>0) std::cout<<"Calculating Coulomb Integrals...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		const double Jii=J_(i,i);
		for(unsigned int j=i+1; j<struc.nAtoms(); ++j){
			Cell::diff(struc.posn(i),struc.posn(j),dr,struc.cell().R(),struc.cell().RInv());
			const double drn=dr.norm();
			J_(j,i)=ke*(
				ICC::iTensor[icc_](drn,k_*ICC::scale[icc_](Jii,J_(j,j)))//short range
				+ewald.phi(struc,dr)-1.0/drn//long range, short range removed
			);
		}
	}
	//calculate the operator matrix
	if(DEBUG_QEQ>0) std::cout<<"Calculating the operator matrix...\n";
	A_=J_.selfadjointView<Eigen::Lower>();
	for(unsigned int i=0; i<struc.nAtoms(); ++i){
		A_(0,i)=1;
		for(unsigned int j=1; j<struc.nAtoms(); ++j){
			A_(j,i)-=J_(i,0);
		}
	}
	//solve the linear equations
	if(DEBUG_QEQ>0) std::cout<<"Solving the linear equations...\n";
	x_.noalias()=A_.partialPivLu().solve(b_);
	//set the atomic charges
	if(DEBUG_QEQ>0) std::cout<<"Setting the charge...\n";
	for(unsigned int i=0; i<struc.nAtoms(); ++i) struc.charge(i)=x_[i];
}
