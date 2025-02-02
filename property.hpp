#ifndef PROPERTY_HPP
#define PROPERTY_HPP

//#define INCLUDE_ORBITAL

//no bounds checking for Eigen objects
#define EIGEN_NO_DEBUG

#include <string>
#include <Eigen/Dense>
#ifdef INCLUDE_ORBITAL
	#include "slater.hpp"
#endif

#ifndef __cplusplus
	#error A C++ compiler is required.
#endif


#ifndef DEBUG_PROPERTY
#define DEBUG_PROPERTY 0
#endif

//*************************************************************************************
//=================================================================================
//------------------------------- BASE CLASSES -------------------------------
//=================================================================================
//*************************************************************************************

//*********************************************************************************
//Property Class 
//*********************************************************************************

template <class T>
class Property{
protected:
	T* ptr_;
	unsigned short* count_;
public:
	//constructors/destructors
	Property();
	Property(const Property& p);
	~Property();
	//operators
	const Property<T>& operator=(const Property<T>& p);
	//access
	T* ptr()const{return ptr_;};
	unsigned short* count()const{return count_;};
	//member functions
	void clear();
	void reset(T* ptr);
	void set(T* ptr);
	void copy(T* ptr);
	void init();
};

//constructors/destructors

template <class T>
Property<T>::Property(){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::Property():\n";
	ptr_=NULL;
	count_=NULL;
}

template <class T>
Property<T>::Property(const Property<T>& p){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::Property(const Property<T>& p):\n";
	if(p.count()!=NULL){count_=p.count(); ++(*count_);}
	else count_=NULL;
	ptr_=p.ptr();
}

template <class T>
Property<T>::~Property(){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::~Property():\n";
	clear();
}

//operators

template <class T>
const Property<T>& Property<T>::operator=(const Property<T>& p){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::operator=(const Property<T>&):\n";
	clear();
	if(p.count()!=NULL){count_=p.count(); ++(*count_);}
	else count_=NULL;
	ptr_=p.ptr();
	return *this;
}

template <class T>
bool operator==(const Property<T>& p1, const Property<T>& p2){
	return (p1.ptr()==p2.ptr());
}

template <class T>
bool operator!=(const Property<T>& p1, const Property<T>& p2){
	return (p1.ptr()!=p2.ptr());
}

//member functions

template <class T>
void Property<T>::clear(){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::clear():\n";
	if(count_!=NULL){
		if(*count_==1){
			delete ptr_; delete count_; 
			ptr_=NULL; count_=NULL;
		} else {
			--(*count_);
			ptr_=NULL; count_=NULL;
		}
	}
}

template <class T>
void Property<T>::reset(T* ptr){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::reset(T*):\n";
	clear();
	ptr_=ptr;
	count_=new unsigned short(1);
}

template <class T>
void Property<T>::set(T* ptr){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::set(T*):\n";
	clear();
	ptr_=ptr;
}

template <class T>
void Property<T>::copy(T* ptr){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::copy(T*):\n";
	clear();
	ptr_=new T(*ptr);
	count_=new unsigned short(1);
}

template <class T>
void Property<T>::init(){
	if(DEBUG_PROPERTY>0) std::cout<<"Property<T>::init():\n";
	clear();
	ptr_=new T();
	count_=new unsigned short(1);
}

//*************************************************************************************
//=================================================================================
//------------------------------- DERIVED CLASSES -------------------------------
//=================================================================================
//*************************************************************************************

//*********************************************************************************
//Name Class
//*********************************************************************************

class Name: public Property<std::string>{
public:
	//constructors/destructors
	Name(){};
	Name(const Name& arg):Property<std::string>(arg){};
	~Name(){};
	//operators
	Name& operator=(const Name& arg){Property<std::string>::operator=(arg); return *this;};
	//access
	std::string& name(){return *ptr_;};
	const std::string& name()const{return *ptr_;};
};

//*********************************************************************************
//Atomic Number Class
//*********************************************************************************

class AN: public Property<unsigned int>{
public:
	//constructors/destructors
	AN(){};
	AN(const AN& arg):Property<unsigned int>(arg){};
	~AN(){};
	//operators
	AN& operator=(const AN& arg){Property<unsigned int>::operator=(arg); return *this;};
	//access
	unsigned int& an(){return *ptr_;};
	const unsigned int& an()const{return *ptr_;};
};

//*********************************************************************************
//Species Class
//*********************************************************************************

class Species: public Property<unsigned int>{
public:
	//constructors/destructors
	Species(){};
	Species(const Species& arg):Property<unsigned int>(arg){};
	~Species(){};
	//operators
	Species& operator=(const Species& arg){Property<unsigned int>::operator=(arg); return *this;};
	//access
	unsigned int& specie(){return *ptr_;};
	const unsigned int& specie()const{return *ptr_;};
};

//*********************************************************************************
//Index Class
//*********************************************************************************

class Index: public Property<unsigned int>{
public:
	//constructors/destructors
	Index(){};
	Index(const Index& arg):Property<unsigned int>(arg){};
	~Index(){};
	//operators
	Index& operator=(const Index& arg){Property<unsigned int>::operator=(arg); return *this;};
	//access
	unsigned int& index(){return *ptr_;};
	const unsigned int& index()const{return *ptr_;};
};

//*******************************************************************************
//Mass
//*******************************************************************************

class Mass: public Property<double>{
public:
	//constructors/destructors
	Mass(){};
	Mass(const Mass& arg):Property<double>(arg){};
	~Mass(){};
	//operators
	Mass& operator=(const Mass& arg){Property<double>::operator=(arg); return *this;};
	//access
	double& mass(){return *ptr_;};
	const double& mass()const{return *ptr_;};
};

//*******************************************************************************
//Charge
//*******************************************************************************

class Charge: public Property<double>{
public:
	//constructors/destructors
	Charge(){};
	Charge(const Charge& arg):Property<double>(arg){};
	~Charge(){};
	//operators
	Charge& operator=(const Charge& arg){Property<double>::operator=(arg); return *this;};
	//access
	double& charge(){return *ptr_;};
	const double& charge()const{return *ptr_;};
};

//*********************************************************************************
//Position Class
//*********************************************************************************

class Position: public Property<Eigen::Vector3d>{
public:
	//constructors/destructors
	Position(){};
	Position(const Position& arg):Property<Eigen::Vector3d>(arg){};
	~Position(){};
	//operators
	Position& operator=(const Position& arg){Property<Eigen::Vector3d>::operator=(arg); return *this;};
	//access
	Eigen::Vector3d& posn(){return *ptr_;};
	const Eigen::Vector3d& posn()const{return *ptr_;};
};

//*********************************************************************************
//Velocity Class
//*********************************************************************************

class Velocity: public Property<Eigen::Vector3d>{
public:
	//constructors/destructors
	Velocity(){};
	Velocity(const Velocity& arg):Property<Eigen::Vector3d>(arg){};
	~Velocity(){};
	//operators
	Velocity& operator=(const Velocity& arg){Property<Eigen::Vector3d>::operator=(arg); return *this;};
	//access
	Eigen::Vector3d& velocity(){return *ptr_;};
	const Eigen::Vector3d& velocity()const{return *ptr_;};
};

//*********************************************************************************
//Force Class
//*********************************************************************************

class Force: public Property<Eigen::Vector3d>{
public:
	//constructors/destructors
	Force(){};
	Force(const Force& arg):Property<Eigen::Vector3d>(arg){};
	~Force(){};
	//operators
	Force& operator=(const Force& arg){Property<Eigen::Vector3d>::operator=(arg); return *this;};
	//access
	Eigen::Vector3d& force(){return *ptr_;};
	const Eigen::Vector3d& force()const{return *ptr_;};
};

//*******************************************************************************
//Dipole
//*******************************************************************************

class Dipole: public Property<Eigen::Vector3d>{
public:
	//constructors/destructors
	Dipole(){};
	Dipole(const Dipole& arg):Property<Eigen::Vector3d>(arg){};
	~Dipole(){};
	//operators
	Dipole& operator=(const Dipole& arg){Property<Eigen::Vector3d>::operator=(arg); return *this;};
	//access
	Eigen::Vector3d& dipole(){return *ptr_;};
	const Eigen::Vector3d& dipole()const{return *ptr_;};
};

//*******************************************************************************
//Alpha
//*******************************************************************************

class Alpha: public Property<Eigen::Matrix3d>{
public:
	//constructors/destructors
	Alpha(){};
	Alpha(const Alpha& arg):Property<Eigen::Matrix3d>(arg){};
	~Alpha(){};
	//operators
	Alpha& operator=(const Alpha& arg){Property<Eigen::Matrix3d>::operator=(arg); return *this;};
	//access
	Eigen::Matrix3d& alpha(){return *ptr_;};
	const Eigen::Matrix3d& alpha()const{return *ptr_;};
};

//*******************************************************************************
//Chi
//*******************************************************************************

class Chi: public Property<double>{
public:
	//constructors/destructors
	Chi(){};
	Chi(const Chi& arg):Property<double>(arg){};
	~Chi(){};
	//operators
	Chi& operator=(const Chi& arg){Property<double>::operator=(arg); return *this;};
	//access
	double& chi(){return *ptr_;};
	const double& chi()const{return *ptr_;};
};

//*******************************************************************************
//Radius
//*******************************************************************************

class Radius: public Property<double>{
public:
	//constructors/destructors
	Radius(){};
	Radius(const Radius& arg):Property<double>(arg){};
	~Radius(){};
	//operators
	Radius& operator=(const Radius& arg){Property<double>::operator=(arg); return *this;};
	//access
	double& radius(){return *ptr_;};
	const double& radius()const{return *ptr_;};
};

//*******************************************************************************
//ZEff
//*******************************************************************************

class ZEff: public Property<double>{
public:
	//constructors/destructors
	ZEff(){};
	ZEff(const ZEff& arg):Property<double>(arg){};
	~ZEff(){};
	//operators
	ZEff& operator=(const ZEff& arg){Property<double>::operator=(arg); return *this;};
	//access
	double& zeff(){return *ptr_;};
	const double& zeff()const{return *ptr_;};
};

//*******************************************************************************
//JZero
//*******************************************************************************

class JZero: public Property<double>{
public:
	//constructors/destructors
	JZero(){};
	JZero(const JZero& arg):Property<double>(arg){};
	~JZero(){};
	//operators
	JZero& operator=(const JZero& arg){Property<double>::operator=(arg); return *this;};
	//access
	double& jzero(){return *ptr_;};
	const double& jzero()const{return *ptr_;};
};

//*******************************************************************************
//Orbital
//*******************************************************************************

#ifdef INCLUDE_ORBITAL
class Orbital: public Property<slater::SlaterS>{
public:
	//constructors/destructors
	Orbital(){};
	Orbital(const Orbital& arg):Property<slater::SlaterS>(arg){};
	~Orbital(){};
	//operators
	Orbital& operator=(const Orbital& arg){Property<slater::SlaterS>::operator=(arg); return *this;};
	//access
	slater::SlaterS& orbital(){return *ptr_;};
	const slater::SlaterS& orbital()const{return *ptr_;};
};
#endif

/*namespace serial{
	
	//**********************************************
	// byte measures
	//**********************************************

	template <> unsigned int nbytes(const Name& obj);
	template <> unsigned int nbytes(const AN& obj);
	template <> unsigned int nbytes(const Species& obj);
	template <> unsigned int nbytes(const Index& obj);
	template <> unsigned int nbytes(const Mass& obj);
	template <> unsigned int nbytes(const Charge& obj);
	template <> unsigned int nbytes(const Position& obj);
	template <> unsigned int nbytes(const Velocity& obj);
	template <> unsigned int nbytes(const Force& obj);
	template <> unsigned int nbytes(const Dipole& obj);
	template <> unsigned int nbytes(const Alpha& obj);
	template <> unsigned int nbytes(const Chi& obj);
	template <> unsigned int nbytes(const Radius& obj);
	template <> unsigned int nbytes(const ZEff& obj);
	template <> unsigned int nbytes(const Symm& obj);
	
	//**********************************************
	// packing
	//**********************************************

	template <> void pack(const Name& obj, char* arr);
	template <> void pack(const AN& obj, char* arr);
	template <> void pack(const Species& obj, char* arr);
	template <> void pack(const Index& obj, char* arr);
	template <> void pack(const Mass& obj, char* arr);
	template <> void pack(const Charge& obj, char* arr);
	template <> void pack(const Position& obj, char* arr);
	template <> void pack(const Velocity& obj, char* arr);
	template <> void pack(const Force& obj, char* arr);
	template <> void pack(const Dipole& obj, char* arr);
	template <> void pack(const Alpha& obj, char* arr);
	template <> void pack(const Chi& obj, char* arr);
	template <> void pack(const Radius& obj, char* arr);
	template <> void pack(const ZEff& obj, char* arr);
	template <> void pack(const Symm& obj, char* arr);
	
	//**********************************************
	// unpacking
	//**********************************************

	template <> void unpack(Name& obj, const char* arr);
	template <> void unpack(AN& obj, const char* arr);
	template <> void unpack(Species& obj, const char* arr);
	template <> void unpack(Index& obj, const char* arr);
	template <> void unpack(Mass& obj, const char* arr);
	template <> void unpack(Charge& obj, const char* arr);
	template <> void unpack(Position& obj, const char* arr);
	template <> void unpack(Velocity& obj, const char* arr);
	template <> void unpack(Force& obj, const char* arr);
	template <> void unpack(Dipole& obj, const char* arr);
	template <> void unpack(Alpha& obj, const char* arr);
	template <> void unpack(Chi& obj, const char* arr);
	template <> void unpack(Radius& obj, const char* arr);
	template <> void unpack(ZEff& obj, const char* arr);
	template <> void unpack(Symm& obj, const char* arr);
	
}
*/
#endif