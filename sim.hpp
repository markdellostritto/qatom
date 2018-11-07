#ifndef SIM_HPP
#define SIM_HPP

// c libraries
#include <cstdlib>
// c++ libraries
#include <iostream>
#include <stdexcept>
#include <string>
// eigen libraries
#include <Eigen/Dense>
// local libraries - cell
#include "cell.hpp"
// local libraries - string utilities
#include "string.hpp"
// local libraries - chemical constants
#include "ptable.hpp"
// local libraries - atom properties
#include "property.hpp"

#ifndef DEBUG_SIM
#define DEBUG_SIM 0
#endif

//**********************************************************************************************
//time_interval struct
//**********************************************************************************************

struct time_interval{
	int beg,end,stride;
	static time_interval load(const char*);
};
std::ostream& operator<<(std::ostream& out, const time_interval& t);

//**********************************************************************************************
//FILE_FORMAT struct
//**********************************************************************************************

struct FILE_FORMAT{
	enum type{
		UNKNOWN,//Unknown format
		XDATCAR,//VASP xdatcar file
		POSCAR,//VASP poscar file
		OUTCAR,//VASP outcar file
		VASP_XML,//VASP XML file
		GAUSSIAN,//Gaussian output file
		DFTB,//DFTB output files
		XYZ,//XYZ file
		CAR,//CAR file
		LAMMPS,//LAMMPS input,data,dump files
		GROMACS,//GROMACS trajectory files
		QE,//quantum espresso output files
		PROPHET,//PROPhet xml file
		XSF//xcrysden format
	};
	static FILE_FORMAT::type load(const std::string& str);
};

std::ostream& operator<<(std::ostream& out, FILE_FORMAT::type& format);

//**********************************************************************************************
//List Atomic
//**********************************************************************************************

template <class T>
class ListAtomic{
private:
	unsigned int ts_;//total number of timesteps
	unsigned int nAtoms_;//total number of atoms
	T* list_;//the list of atoms
public:
	//macros
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW//just in case we are storing Eigen objects
	
	//constructors/destructors
	ListAtomic(){init();};
	ListAtomic(unsigned int ts, unsigned int nAtoms){init(); resize(ts,nAtoms);};
	ListAtomic(const ListAtomic<T>& l);//deep copy
	~ListAtomic();
	
	//operators
	ListAtomic<T>& operator=(const ListAtomic<T>& l);//deep copy
	
	//access
	//primitives
		unsigned int ts()const{return ts_;};
		unsigned int nAtoms()const{return nAtoms_;};
		unsigned int size()const{return ts_*nAtoms_;};
	//list access - operators
		T& operator()(unsigned int t, unsigned int n){return list_[t*nAtoms_+n];};
		const T& operator()(unsigned int t, unsigned int n)const{return list_[t*nAtoms_+n];};
		T& operator[](unsigned int i){return list_[i];};
		const T& operator[](unsigned int i)const{return list_[i];};
	//list access - stored array
		const T* atomList()const{return list_;};
	
	//static member functions
	static unsigned int size(const ListAtomic<T>& list){return sizeof(T)*list.size()*list.nAtoms();};
	
	//member functions
	void init();
	void clear();
	bool empty()const{return (ts_==0)?true:false;};
	void resize(unsigned int ts, unsigned int nAtoms);
};

//constructors/destructors

template <class T>
ListAtomic<T>::ListAtomic(const ListAtomic<T>& l){
	if(DEBUG_SIM>0) std::cout<<"ListAtomic<T>::ListAtomic(const ListAtomic<T>&):\n";
	init(); resize(l.ts(),l.nAtoms());
	for(unsigned int i=0; i<size(); ++i) list_[i]=l[i];
}

template <class T>
ListAtomic<T>::~ListAtomic(){
	if(DEBUG_SIM>0) std::cout<<"ListAtomic<T>::~ListAtomic()\n";
	clear();
}

//operators

template <class T>
ListAtomic<T>& ListAtomic<T>::operator=(const ListAtomic<T>& l){
	if(DEBUG_SIM>0) std::cout<<"ListAtomic<T>::operator=(const ListAtomic<T>&):\n";
	ListAtomic<T> templist(l);
	resize(templist.ts(),templist.nAtoms());
	for(unsigned int i=0; i<size(); ++i) list_[i]=templist[i];
	return *this;
}

//member functions

template <class T>
void ListAtomic<T>::init(){
	if(DEBUG_SIM>0) std::cout<<"ListAtomic<T>::init():\n";
	list_=NULL;
	ts_=0; nAtoms_=0;
}

template <class T>
void ListAtomic<T>::clear(){
	if(DEBUG_SIM>0) std::cout<<"ListAtomic<T>::clear():\n";
	if(list_!=NULL) delete[] list_;
	list_=NULL;
	ts_=0; nAtoms_=0;
}

template <class T>
void ListAtomic<T>::resize(unsigned int ts, unsigned int nAtoms){
	if(DEBUG_SIM>0) std::cout<<"ListAtomic<T>::resize(unsigned int,unsigned int):\n";
	try{
		if(list_!=NULL) delete[] list_;
		ts_=ts; nAtoms_=nAtoms;
		list_=new T[ts_*nAtoms_];
	}catch(std::exception& e){
		clear();
		throw e;
	}
}

//operators

template <class T>
bool operator==(const ListAtomic<T>& l1, const ListAtomic<T>& l2){
	return (l1.ts()==l2.ts() && l1.nAtoms()==l2.nAtoms());
}

template <class T>
bool operator!=(const ListAtomic<T>& l1, const ListAtomic<T>& l2){
	return !(l1.ts()==l2.ts() && l1.nAtoms()==l2.nAtoms());
}

//**********************************************************************************************
//Sim Interface
//**********************************************************************************************

class SimI{
public:
	static const unsigned int D=3;
protected:
	std::string system_;//the name of the simulation
	
	unsigned int timesteps_;//total number of timesteps
	unsigned int beg_;//the beginning of the time interval
	unsigned int end_;//the end of the time interval
	unsigned int stride_;
	
	double timestep_;//the length in picoseconds of each timestep
	
	bool periodic_;//whether the simulation is periodic or not
	
	std::vector<Cell> cells_;//unit cells
	std::vector<double> energy_;//energy
	std::vector<double> temp_;//temperature
	bool cellFixed_;//whether the cell is static
public:
	//constructors/destructors
	SimI(){defaults();};
	SimI(const SimI& arg);
	~SimI(){};
	
	//operators
	SimI& operator=(const SimI& arg);
	friend std::ostream& operator<<(std::ostream& out, const SimI& sim);
	
	//access
	//name
		std::string& system(){return system_;};
		const std::string& system()const{return system_;};
	//interval
		unsigned int timesteps()const{return timesteps_;};
		unsigned int& beg(){return beg_;};
		const unsigned int& beg()const{return beg_;};
		unsigned int& end(){return end_;};
		const unsigned int& end()const{return end_;};
		unsigned int& stride(){return stride_;};
		const unsigned int& stride()const{return stride_;};
	//timestep
		double& timestep(){return timestep_;};
		const double& timestep()const{return timestep_;};
	//periodicity
		bool& periodic(){return periodic_;};
		const bool& periodic()const{return periodic_;};
	//cells
		const std::vector<Cell>& cell()const{return cells_;};
		Cell& cell(unsigned int i){return cells_[i];};
		const Cell& cell(unsigned int i)const{return cells_[i];};
		bool& cellFixed(){return cellFixed_;};
		const bool& cellFixed()const{return cellFixed_;};
	//energy
		const std::vector<double>& energy()const{return energy_;};
		double& energy(unsigned int i){return energy_[i];};
		const double& energy(unsigned int i)const{return energy_[i];};
	
	//member functions
	void defaults();
	void clear(){defaults();};
	bool empty() const{return (timesteps_==0)?true:false;};
	void resize(unsigned int ts);
	
	//static functions
	static FILE_FORMAT::type fileFormat(const char* fileName);
};

//**********************************************************************************************
//Sim Atomic Interface
//**********************************************************************************************

class SimAtomicI{
protected:
	unsigned int nSpecies_;//number of atomic species
	unsigned int nAtomsT_;//total number of atoms
	std::vector<unsigned int> nAtoms_;//the number of atoms of each species
	std::vector<unsigned int> offsets_;
	std::vector<std::string> atomNames_;//the names of each species
public:
	SimAtomicI(){defaults();};
	SimAtomicI(const SimAtomicI& sim);
	SimAtomicI(const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& atomNames){resize(nAtoms,atomNames);};
	~SimAtomicI(){};
	
	//operators
	friend std::ostream& operator<<(std::ostream& out, const SimAtomicI& sim);
	SimAtomicI& operator=(const SimAtomicI& sim);
	
	//access
	unsigned int nSpecies()const{return nSpecies_;};
	unsigned int nAtoms()const{return nAtomsT_;};
	const std::vector<unsigned int>& nAtomsVec()const{return nAtoms_;};
	unsigned int nAtoms(unsigned int i)const{return nAtoms_[i];};
	unsigned int offset(unsigned int i)const{return offsets_[i];};
	const std::string& atomNames(unsigned int i)const{return atomNames_[i];};
	const std::vector<std::string>& atomNames()const{return atomNames_;};
	
	//member functions
	void defaults();
	void clear(){defaults();};
	void resize(const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& atomNames);
	int speciesIndex(const std::string& str)const{return speciesIndex(str,atomNames_);};
	int speciesIndex(const char* str)const{return speciesIndex(str,atomNames_);};
	
	//static functions
	static int speciesIndex(const std::string& str, const std::vector<std::string>& atomNames);
	static int speciesIndex(const char* str, const std::vector<std::string>& atomNames);
};

//**********************************************************************************************
//Sim Molecular Interface
//**********************************************************************************************

class SimMolI{
protected:
	unsigned int nMol_;//number of molecules
public:
	SimMolI(){defaults();};
	SimMolI(const SimMolI& sim);
	SimMolI(unsigned int nMol){resize(nMol);};
	
	//operator
	friend std::ostream& operator<<(std::ostream& out, const SimMolI& sim);
	SimMolI& operator=(const SimMolI& sim);
	
	//access
	unsigned int nMol()const{return nMol_;};
	
	//member functions
	void defaults();
	void clear(){defaults();};
	void resize(unsigned int nMol);
};

//**********************************************************************************************
//Sim Atomic Storage
//**********************************************************************************************

template <class AtomT>
class SimAtomicS: public SimAtomicI{
protected:
	//the atoms of the simulation
	ListAtomic<AtomT> atoms_;
	
	//atom properties
	ListAtomic<std::string> atom_name_;//names
	ListAtomic<unsigned int> atom_an_;//atomic numbers
	ListAtomic<double> atom_mass_;//masses
	ListAtomic<unsigned int> atom_specie_;//species 
	ListAtomic<unsigned int> atom_index_;//indices
	ListAtomic<Eigen::Vector3d> atom_posn_;//positions
	ListAtomic<Eigen::Vector3d> atom_force_;//forces
	ListAtomic<double> atom_charge_;//charges
	ListAtomic<Eigen::Vector3d> atom_dipole_;//dipoles
	ListAtomic<Eigen::Matrix3d> atom_alpha_;//polarizabilities
	ListAtomic<double> atom_zeff_;//effective potential
	ListAtomic<double> atom_chi_;//electronegativity
	ListAtomic<double> atom_jzero_;//idempotential
	
	//resizing functions
	void resizeAtomArrays(unsigned int ts, unsigned int nAtomsT);
	void resizeAtomName(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomName(std::false_type, unsigned int ts, unsigned int nAtomsT){};//names
	void resizeAtomAN(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomAN(std::false_type, unsigned int ts, unsigned int nAtomsT){};//atomic numbers
	void resizeAtomMass(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomMass(std::false_type, unsigned int ts, unsigned int nAtomsT){};//masses
	void resizeAtomSpecie(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomSpecie(std::false_type, unsigned int ts, unsigned int nAtomsT){};//species
	void resizeAtomIndex(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomIndex(std::false_type, unsigned int ts, unsigned int nAtomsT){};//indices
	void resizeAtomPosn(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomPosn(std::false_type, unsigned int ts, unsigned int nAtomsT){};//positions
	void resizeAtomForce(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomForce(std::false_type, unsigned int ts, unsigned int nAtomsT){};//forces
	void resizeAtomCharge(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomCharge(std::false_type, unsigned int ts, unsigned int nAtomsT){};//charges
	void resizeAtomDipole(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomDipole(std::false_type, unsigned int ts, unsigned int nAtomsT){};//dipoles
	void resizeAtomAlpha(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomAlpha(std::false_type, unsigned int ts, unsigned int nAtomsT){};//polarizabilities
	void resizeAtomZEff(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomZEff(std::false_type, unsigned int ts, unsigned int nAtomsT){};//polarizabilities
	void resizeAtomChi(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomChi(std::false_type, unsigned int ts, unsigned int nAtomsT){};//electronegativity
	void resizeAtomJZero(std::true_type, unsigned int ts, unsigned int nAtomsT); void resizeAtomJZero(std::false_type, unsigned int ts, unsigned int nAtomsT){};//idempotential
	
	//assigment functions
	void assignAtomArrays();
	void assignAtomName(std::true_type); void assignAtomName(std::false_type){};//names
	void assignAtomAN(std::true_type); void assignAtomAN(std::false_type){};//atomic numbers
	void assignAtomMass(std::true_type); void assignAtomMass(std::false_type){};//masses
	void assignAtomSpecie(std::true_type); void assignAtomSpecie(std::false_type){};//species
	void assignAtomPosn(std::true_type); void assignAtomPosn(std::false_type){};//indices
	void assignAtomForce(std::true_type); void assignAtomForce(std::false_type){};//indices
	void assignAtomIndex(std::true_type); void assignAtomIndex(std::false_type){};//positions
	void assignAtomCharge(std::true_type); void assignAtomCharge(std::false_type){};//charges
	void assignAtomDipole(std::true_type); void assignAtomDipole(std::false_type){};//dipoles
	void assignAtomAlpha(std::true_type); void assignAtomAlpha(std::false_type){};//polarizabilities
	void assignAtomZEff(std::true_type); void assignAtomZEff(std::false_type){};//polarizabilities
	void assignAtomChi(std::true_type); void assignAtomChi(std::false_type){};//electronegativity
	void assignAtomJZero(std::true_type); void assignAtomJZero(std::false_type){};//idempotential
public:
	//constructors/destructors
	SimAtomicS(){};
	SimAtomicS(const SimAtomicS<AtomT>& sim);
	~SimAtomicS(){};
	
	//operators
	SimAtomicS<AtomT>& operator=(const SimAtomicS<AtomT>& sim);
	template <class AtomT_> friend std::ostream& operator<<(std::ostream& out, const SimAtomicS<AtomT_>& sim);
	
	//access - properties
	//names
		ListAtomic<std::string>& atom_name(){return atom_name_;};
		const ListAtomic<std::string>& atom_name()const{return atom_name_;};
	//atomic numbers
		ListAtomic<unsigned int>& atom_an(){return atom_an_;};
		const ListAtomic<unsigned int>& atom_an()const{return atom_an_;};
	//masses
		ListAtomic<double>& atom_mass(){return atom_mass_;};
		const ListAtomic<double>& atom_mass()const{return atom_mass_;};
	//species
		ListAtomic<unsigned int>& atom_specie(){return atom_specie_;};
		const ListAtomic<unsigned int>& atom_specie()const{return atom_specie_;};
	//species
		ListAtomic<unsigned int>& atom_index(){return atom_index_;};
		const ListAtomic<unsigned int>& atom_index()const{return atom_index_;};
	//positions
		ListAtomic<Eigen::Vector3d>& atom_posn(){return atom_posn_;};
		const ListAtomic<Eigen::Vector3d>& atom_posn()const{return atom_posn_;};
	//positions
		ListAtomic<Eigen::Vector3d>& atom_force(){return atom_force_;};
		const ListAtomic<Eigen::Vector3d>& atom_force()const{return atom_force_;};
	//charges
		ListAtomic<double>& atom_charge(){return atom_charge_;};
		const ListAtomic<double>& atom_charge()const{return atom_charge_;};
	//dipoles
		ListAtomic<Eigen::Vector3d>& atom_dipole(){return atom_dipole_;};
		const ListAtomic<Eigen::Vector3d>& atom_dipole()const{return atom_dipole_;};
	//polarizabilities
		ListAtomic<Eigen::Matrix3d>& atom_alpha(){return atom_alpha_;};
		const ListAtomic<Eigen::Matrix3d>& atom_alpha()const{return atom_alpha_;};
	//zeff
		ListAtomic<double>& atom_zeff(){return atom_zeff_;};
		const ListAtomic<double>& atom_zeff()const{return atom_zeff_;};
	//chi
		ListAtomic<double>& atom_chi(){return atom_chi_;};
		const ListAtomic<double>& atom_chi()const{return atom_chi_;};
	//jzero
		ListAtomic<double>& atom_jzero(){return atom_jzero_;};
		const ListAtomic<double>& atom_jzero()const{return atom_jzero_;};
	
	//access - atoms
	const ListAtomic<AtomT>& atoms()const{return atoms_;};
	AtomT& atom(unsigned int t, unsigned int n){return atoms_(t,n);};
	const AtomT& atom(unsigned int t, unsigned int n)const{return atoms_(t,n);};
	AtomT& atom(unsigned int t, unsigned int n, unsigned int m);
	const AtomT& atom(unsigned int t, unsigned int n, unsigned int m)const;
	
	//member functions
	void clear();
	void resize(unsigned int ts, const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& speciesNames);
};

//constructors/destructors

template <class AtomT>
SimAtomicS<AtomT>::SimAtomicS(const SimAtomicS<AtomT>& sim):SimAtomicI(sim){
	//set the property arrays
	atom_name_=sim.atom_name();
	atom_an_=sim.atom_an();
	atom_mass_=sim.atom_mass();
	atom_specie_=sim.atom_specie();
	atom_index_=sim.atom_index();
	atom_posn_=sim.atom_posn();
	atom_force_=sim.atom_force();
	atom_charge_=sim.atom_charge();
	atom_dipole_=sim.atom_dipole();
	atom_alpha_=sim.atom_alpha();
	atom_zeff_=sim.atom_zeff();
	atom_chi_=sim.atom_chi();
	atom_jzero_=sim.atom_jzero();
	//resize the atom array
	atoms_.resize(sim.atoms().ts(),sim.atoms().nAtoms());
	//assign the atom arrays
	assignAtomArrays();
};

//operators

template <class AtomT>
SimAtomicS<AtomT>& SimAtomicS<AtomT>::operator=(const SimAtomicS<AtomT>& sim){
	SimAtomicI::operator=(sim);
	//set the property arrays
	atom_name_=sim.atom_name();
	atom_an_=sim.atom_an();
	atom_mass_=sim.atom_mass();
	atom_specie_=sim.atom_specie();
	atom_index_=sim.atom_index();
	atom_posn_=sim.atom_posn();
	atom_force_=sim.atom_force();
	atom_charge_=sim.atom_charge();
	atom_dipole_=sim.atom_dipole();
	atom_alpha_=sim.atom_alpha();
	atom_zeff_=sim.atom_zeff();
	atom_chi_=sim.atom_chi();
	atom_jzero_=sim.atom_jzero();
	//resize the atom array
	atoms_.resize(sim.atoms().ts(),sim.atoms().nAtoms());
	//assign the atom arrays
	assignAtomArrays();
	return *this;
}

template <class AtomT>
std::ostream& operator<<(std::ostream& out, const SimAtomicS<AtomT>& sim){
	return out<<static_cast<const SimAtomicI&>(sim);
}

//member functions

template <class AtomT>
void SimAtomicS<AtomT>::clear(){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::clear():\n";
	SimAtomicI::clear();
	//clear property arrays
	atom_name_.clear();
	atom_an_.clear();
	atom_mass_.clear();
	atom_specie_.clear();
	atom_index_.clear();
	atom_posn_.clear();
	atom_force_.clear();
	atom_charge_.clear();
	atom_dipole_.clear();
	atom_alpha_.clear();
	atom_zeff_.clear();
	atom_chi_.clear();
	atom_jzero_.clear();
	//clear the atoms
	atoms_.clear();
}

template <class AtomT>
void SimAtomicS<AtomT>::resize(unsigned int ts, const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& speciesNames){
	//find the total number of atoms
	unsigned int nAtomsT=0;
	for(unsigned int i=0; i<nAtoms.size(); ++i) nAtomsT+=nAtoms[i];
	//resize SimAtomicI
	SimAtomicI::resize(nAtoms,speciesNames);
	//resize the atom/molecule arrays
	atoms_.resize(ts,nAtomsT);
	//resize the property arrays
	resizeAtomArrays(ts,nAtomsT);
	//assign atom arrays
	assignAtomArrays();
	//set the names, species, and indices
	unsigned int count;
	for(unsigned int t=0; t<ts; ++t){
		count=0;
		for(unsigned int n=0; n<nSpecies_; ++n){
			for(unsigned int m=0; m<nAtoms_[n]; ++m){
				atom_name_(t,count)=atomNames_[n];
				atom_an_(t,count)=PTable::an(atomNames_[n].c_str());
				atom_specie_(t,count)=n;
				atom_index_(t,count)=m;
				++count;
			}
		}
	}
}

//access - atoms

template <class AtomT>
AtomT& SimAtomicS<AtomT>::atom(unsigned int t, unsigned int n, unsigned int m){
	return atoms_(t,offsets_[n]+m);
}

template <class AtomT>
const AtomT& SimAtomicS<AtomT>::atom(unsigned int t, unsigned int n, unsigned int m)const{
	return atoms_(t,offsets_[n]+m);
}

//resizing atom properties

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomArrays(unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomic<AtomT>::resizeArrays():\n";
	resizeAtomName(std::is_base_of<Name,AtomT>(),ts,nAtomsT);
	resizeAtomAN(std::is_base_of<AN,AtomT>(),ts,nAtomsT);
	resizeAtomMass(std::is_base_of<Mass,AtomT>(),ts,nAtomsT);
	resizeAtomSpecie(std::is_base_of<Species,AtomT>(),ts,nAtomsT);
	resizeAtomIndex(std::is_base_of<Index,AtomT>(),ts,nAtomsT);
	resizeAtomPosn(std::is_base_of<Position,AtomT>(),ts,nAtomsT);
	resizeAtomForce(std::is_base_of<Force,AtomT>(),ts,nAtomsT);
	resizeAtomCharge(std::is_base_of<Charge,AtomT>(),ts,nAtomsT);
	resizeAtomDipole(std::is_base_of<Dipole,AtomT>(),ts,nAtomsT);
	resizeAtomAlpha(std::is_base_of<Alpha,AtomT>(),ts,nAtomsT);
	resizeAtomZEff(std::is_base_of<ZEff,AtomT>(),ts,nAtomsT);
	resizeAtomChi(std::is_base_of<Chi,AtomT>(),ts,nAtomsT);
	resizeAtomJZero(std::is_base_of<JZero,AtomT>(),ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomName(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomName(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_name_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomAN(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomAN(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_an_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomMass(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomMass(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_mass_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomSpecie(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomSpecie(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_specie_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomIndex(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomIndex(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_index_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomPosn(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomPosn(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_posn_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomForce(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomForce(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_force_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomCharge(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomCharge(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_charge_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomDipole(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomDipole(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_dipole_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomAlpha(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomAlpha(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_alpha_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomZEff(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomZEff(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_zeff_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomChi(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomChi(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_chi_.resize(ts,nAtomsT);
}

template <class AtomT>
void SimAtomicS<AtomT>::resizeAtomJZero(std::true_type, unsigned int ts, unsigned int nAtomsT){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::resizeAtomJZero(std::true_type, unsigned int ts, unsigned int nAtomsT):\n";
	atom_jzero_.resize(ts,nAtomsT);
}

//assigning atom properties

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomArrays(){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomArrays():\n";
	assignAtomName(std::is_base_of<Name,AtomT>());
	assignAtomAN(std::is_base_of<AN,AtomT>());
	assignAtomMass(std::is_base_of<Mass,AtomT>());
	assignAtomSpecie(std::is_base_of<Species,AtomT>());
	assignAtomIndex(std::is_base_of<Index,AtomT>());
	assignAtomPosn(std::is_base_of<Position,AtomT>());
	assignAtomForce(std::is_base_of<Force,AtomT>());
	assignAtomCharge(std::is_base_of<Charge,AtomT>());
	assignAtomDipole(std::is_base_of<Dipole,AtomT>());
	assignAtomAlpha(std::is_base_of<Alpha,AtomT>());
	assignAtomZEff(std::is_base_of<ZEff,AtomT>());
	assignAtomChi(std::is_base_of<Chi,AtomT>());
	assignAtomJZero(std::is_base_of<JZero,AtomT>());
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomName(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomName(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Name&>(atoms_[i]).set(&atom_name_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomAN(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomAN(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<AN&>(atoms_[i]).set(&atom_an_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomMass(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomMass(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Mass&>(atoms_[i]).set(&atom_mass_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomSpecie(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomSpecie(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Species&>(atoms_[i]).set(&atom_specie_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomIndex(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomIndex(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Index&>(atoms_[i]).set(&atom_index_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomPosn(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomPosn(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Position&>(atoms_[i]).set(&atom_posn_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomForce(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomForce(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Force&>(atoms_[i]).set(&atom_force_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomCharge(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomCharge(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Charge&>(atoms_[i]).set(&atom_charge_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomDipole(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomDipole(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Dipole&>(atoms_[i]).set(&atom_dipole_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomAlpha(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomAlpha(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Alpha&>(atoms_[i]).set(&atom_alpha_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomZEff(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomZEff(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<ZEff&>(atoms_[i]).set(&atom_zeff_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomChi(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomChi(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<Chi&>(atoms_[i]).set(&atom_chi_[i]);
	}
}

template <class AtomT>
void SimAtomicS<AtomT>::assignAtomJZero(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimAtomicS<AtomT>::assignAtomJZero(std::true_type):\n";
	for(unsigned int i=0; i<atoms_.size(); ++i){
		static_cast<JZero&>(atoms_[i]).set(&atom_jzero_[i]);
	}
}

//**********************************************************************************************
//Sim Molecular Storage
//**********************************************************************************************

template <class MoleculeT>
class SimMolS: public SimMolI{
protected:
	//molecules
	ListAtomic<MoleculeT> mols_;//names
	
	//molecule properties
	ListAtomic<std::string> mol_name_;//names
	ListAtomic<double> mol_mass_;//masses
	ListAtomic<unsigned int> mol_index_;//indices
	ListAtomic<Eigen::Vector3d> mol_posn_;//positions
	ListAtomic<double> mol_charge_;//charges
	ListAtomic<Eigen::Vector3d> mol_dipole_;//dipoles
	ListAtomic<Eigen::Matrix3d> mol_alpha_;//polarizabilities
	
	//resizing functions
	void resizeMolArrays(unsigned int ts, unsigned int nMolsT);
	void resizeMolName(std::true_type, unsigned int ts, unsigned int nMolsT); void resizeMolName(std::false_type, unsigned int ts, unsigned int nMolsT){};//names
	void resizeMolMass(std::true_type, unsigned int ts, unsigned int nMolsT); void resizeMolMass(std::false_type, unsigned int ts, unsigned int nMolsT){};//masses
	void resizeMolIndex(std::true_type, unsigned int ts, unsigned int nMolsT); void resizeMolIndex(std::false_type, unsigned int ts, unsigned int nMolsT){};//indices
	void resizeMolPosn(std::true_type, unsigned int ts, unsigned int nMolsT); void resizeMolPosn(std::false_type, unsigned int ts, unsigned int nMolsT){};//positions
	void resizeMolCharge(std::true_type, unsigned int ts, unsigned int nMolsT); void resizeMolCharge(std::false_type, unsigned int ts, unsigned int nMolsT){};//charges
	void resizeMolDipole(std::true_type, unsigned int ts, unsigned int nMolsT); void resizeMolDipole(std::false_type, unsigned int ts, unsigned int nMolsT){};//dipoles
	void resizeMolAlpha(std::true_type, unsigned int ts, unsigned int nMolsT); void resizeMolAlpha(std::false_type, unsigned int ts, unsigned int nMolsT){};//polarizabilities
	
	//assigment functions
	void assignMolArrays();
	void assignMolName(std::true_type); void assignMolName(std::false_type){};//names
	void assignMolMass(std::true_type); void assignMolMass(std::false_type){};//masses
	void assignMolIndex(std::true_type); void assignMolIndex(std::false_type){};//positions
	void assignMolPosn(std::true_type); void assignMolPosn(std::false_type){};//indices
	void assignMolCharge(std::true_type); void assignMolCharge(std::false_type){};//charges
	void assignMolDipole(std::true_type); void assignMolDipole(std::false_type){};//dipoles
	void assignMolAlpha(std::true_type); void assignMolAlpha(std::false_type){};//polarizabilities
	
public:
	//constructors/destructors
	SimMolS(){};
	SimMolS(const SimMolS<MoleculeT>& sim);
	~SimMolS(){};
	
	//operators
	SimMolS<MoleculeT>& operator=(const SimMolS<MoleculeT>& sim);
	template <class MoleculeT_> friend std::ostream& operator<<(std::ostream& out, const SimMolS<MoleculeT_>& sim);
	
	//access - properties
	//names
		ListAtomic<std::string>& mol_name(){return mol_name_;};
		const ListAtomic<std::string>& mol_name()const{return mol_name_;};
	//masses
		ListAtomic<double>& mol_mass(){return mol_mass_;};
		const ListAtomic<double>& mol_mass()const{return mol_mass_;};
	//indices
		ListAtomic<unsigned int>& mol_index(){return mol_index_;};
		const ListAtomic<unsigned int>& mol_index()const{return mol_index_;};
	//positions
		ListAtomic<Eigen::Vector3d>& mol_posn(){return mol_posn_;};
		const ListAtomic<Eigen::Vector3d>& mol_posn()const{return mol_posn_;};
	//charges
		ListAtomic<double>& mol_charge(){return mol_charge_;};
		const ListAtomic<double>& mol_charge()const{return mol_charge_;};
	//dipoles
		ListAtomic<Eigen::Vector3d>& mol_dipole(){return mol_dipole_;};
		const ListAtomic<Eigen::Vector3d>& mol_dipole()const{return mol_dipole_;};
	//polarizabilities
		ListAtomic<Eigen::Matrix3d>& mol_alpha(){return mol_alpha_;};
		const ListAtomic<Eigen::Matrix3d>& mol_alpha()const{return mol_alpha_;};
	
	//access - molecules
	const ListAtomic<MoleculeT>& molecules()const{return mols_;};
	MoleculeT& molecule(unsigned int t, unsigned int n){return mols_(t,n);};
	const MoleculeT& molecule(unsigned int t, unsigned int n)const{return mols_(t,n);};
	
	//member functions
	void clear();
	void resize(unsigned int ts, unsigned int nMols);
};

//constructors/destructors

template <class MoleculeT>
SimMolS<MoleculeT>::SimMolS(const SimMolS<MoleculeT>& sim):SimMolI(sim){
	//set the property arrays
	mol_name_=sim.mol_name();
	mol_mass_=sim.mol_mass();
	mol_index_=sim.mol_index();
	mol_posn_=sim.mol_posn();
	mol_charge_=sim.mol_charge();
	mol_dipole_=sim.mol_dipole();
	mol_alpha_=sim.mol_alpha();
	//resize the atom/molecule arrays
	mols_.resize(sim.molecules().timesteps(),sim.molecules().nAtoms());
	//assign the atom properties to the per-atom arrays
	assignMolArrays();
}

//operators

template <class MoleculeT>
SimMolS<MoleculeT>& SimMolS<MoleculeT>::operator=(const SimMolS<MoleculeT>& sim){
	SimMolI::operator=(sim);
	//set the property arrays
	mol_name_=sim.mol_name();
	mol_mass_=sim.mol_mass();
	mol_index_=sim.mol_index();
	mol_posn_=sim.mol_posn();
	mol_charge_=sim.mol_charge();
	mol_dipole_=sim.mol_dipole();
	mol_alpha_=sim.mol_alpha();
	//resize the atom/molecule arrays
	mols_.resize(sim.mols().ts(),sim.molecules().nAtoms());//"nAtoms" gives total "atoms", same as number of molecules
	//assign the atom properties to the per-atom arrays
	assignMolArrays();
	//return object
	return *this;
}

template <class MoleculeT>
std::ostream& operator<<(std::ostream& out, const SimMolS<MoleculeT>& sim){
	return out<<static_cast<const SimMolI&>(sim);
}

//member functions

template <class MoleculeT>
void SimMolS<MoleculeT>::clear(){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::clear():\n";
	SimMolI::clear();
	//clear the property arrays
	mol_name_.clear();
	mol_mass_.clear();
	mol_index_.clear();
	mol_posn_.clear();
	mol_charge_.clear();
	mol_dipole_.clear();
	mol_alpha_.clear();
	//clear the molecule array
	mols_.clear();
}

template <class MoleculeT>
void SimMolS<MoleculeT>::resize(unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resize(unsigned int,unsigned int):\n";
	//resize the atom/molecule arrays
	nMol_=nMolsT;
	mols_.resize(ts,nMolsT);
	//resize the property arrays
	this->resizeMolArrays(ts,nMolsT);
	//assign the atom properties to the per-atom arrays
	this->assignMolArrays();
	//set the molecule indices
	if(DEBUG_SIM>0) std::cout<<"Setting molecule indices...\n";
	for(unsigned int t=0; t<ts; ++t){
		for(unsigned int n=0; n<nMolsT; ++n){
			mol_index_(t,n)=n;
		}
	}
}

//resizing atom properties

template <class MoleculeT>
void SimMolS<MoleculeT>::resizeMolArrays(unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resizeMolArrays(unsigned int,unsigned int):\n";
	resizeMolName(std::is_base_of<Name,MoleculeT>(),ts,nMolsT);
	resizeMolMass(std::is_base_of<Mass,MoleculeT>(),ts,nMolsT);
	resizeMolIndex(std::is_base_of<Index,MoleculeT>(),ts,nMolsT);
	resizeMolPosn(std::is_base_of<Position,MoleculeT>(),ts,nMolsT);
	resizeMolCharge(std::is_base_of<Charge,MoleculeT>(),ts,nMolsT);
	resizeMolDipole(std::is_base_of<Dipole,MoleculeT>(),ts,nMolsT);
	resizeMolAlpha(std::is_base_of<Alpha,MoleculeT>(),ts,nMolsT);
}

template <class MoleculeT>
void SimMolS<MoleculeT>::resizeMolName(std::true_type, unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resizeMolName(std::true_type, unsigned int ts, unsigned int nMolsT):\n";
	mol_name_.resize(ts,nMolsT);
}

template <class MoleculeT>
void SimMolS<MoleculeT>::resizeMolMass(std::true_type, unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resizeMolMass(std::true_type, unsigned int ts, unsigned int nMolsT):\n";
	mol_mass_.resize(ts,nMolsT);
}

template <class MoleculeT>
void SimMolS<MoleculeT>::resizeMolIndex(std::true_type, unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resizeMolIndex(std::true_type, unsigned int ts, unsigned int nMolsT):\n";
	mol_index_.resize(ts,nMolsT);
}

template <class MoleculeT>
void SimMolS<MoleculeT>::resizeMolPosn(std::true_type, unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resizeMolPosn(std::true_type, unsigned int ts, unsigned int nMolsT):\n";
	mol_posn_.resize(ts,nMolsT);
}

template <class MoleculeT>
void SimMolS<MoleculeT>::resizeMolCharge(std::true_type, unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resizeMolCharge(std::true_type, unsigned int ts, unsigned int nMolsT):\n";
	mol_charge_.resize(ts,nMolsT);
}

template <class MoleculeT>
void SimMolS<MoleculeT>::resizeMolDipole(std::true_type, unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resizeMolDipole(std::true_type, unsigned int ts, unsigned int nMolsT):\n";
	mol_dipole_.resize(ts,nMolsT);
}

template <class MoleculeT>
void SimMolS<MoleculeT>::resizeMolAlpha(std::true_type, unsigned int ts, unsigned int nMolsT){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::resizeMolAlpha(std::true_type, unsigned int ts, unsigned int nMolsT):\n";
	mol_alpha_.resize(ts,nMolsT);
}

//assigning molecule properties

template <class MoleculeT>
void SimMolS<MoleculeT>::assignMolArrays(){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::assignMolArrays():\n";
	assignMolName(std::is_base_of<Name,MoleculeT>());
	assignMolMass(std::is_base_of<Mass,MoleculeT>());
	assignMolIndex(std::is_base_of<Index,MoleculeT>());
	assignMolPosn(std::is_base_of<Position,MoleculeT>());
	assignMolCharge(std::is_base_of<Charge,MoleculeT>());
	assignMolDipole(std::is_base_of<Dipole,MoleculeT>());
	assignMolAlpha(std::is_base_of<Alpha,MoleculeT>());
}

template <class MoleculeT>
void SimMolS<MoleculeT>::assignMolName(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::assignName(std::true_type):\n";
	for(unsigned int i=0; i<mols_.size(); ++i){
		static_cast<Name&>(mols_[i]).set(&mol_name_[i]);
	}
}

template <class MoleculeT>
void SimMolS<MoleculeT>::assignMolMass(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::assignMass(std::true_type):\n";
	for(unsigned int i=0; i<mols_.size(); ++i){
		static_cast<Mass&>(mols_[i]).set(&mol_mass_[i]);
	}
}

template <class MoleculeT>
void SimMolS<MoleculeT>::assignMolIndex(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::assignIndex(std::true_type):\n";
	for(unsigned int i=0; i<mols_.size(); ++i){
		static_cast<Index&>(mols_[i]).set(&mol_index_[i]);
	}
}

template <class MoleculeT>
void SimMolS<MoleculeT>::assignMolPosn(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::assignPosn(std::true_type):\n";
	for(unsigned int i=0; i<mols_.size(); ++i){
		static_cast<Position&>(mols_[i]).set(&mol_posn_[i]);
	}
}

template <class MoleculeT>
void SimMolS<MoleculeT>::assignMolCharge(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::assignCharge(std::true_type):\n";
	for(unsigned int i=0; i<mols_.size(); ++i){
		static_cast<Charge&>(mols_[i]).set(&mol_charge_[i]);
	}
}

template <class MoleculeT>
void SimMolS<MoleculeT>::assignMolDipole(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::assignDipole(std::true_type):\n";
	for(unsigned int i=0; i<mols_.size(); ++i){
		static_cast<Dipole&>(mols_[i]).set(&mol_dipole_[i]);
	}
}

template <class MoleculeT>
void SimMolS<MoleculeT>::assignMolAlpha(std::true_type){
	if(DEBUG_SIM>0) std::cout<<"SimMolS<MoleculeT>::assignAlpha(std::true_type):\n";
	for(unsigned int i=0; i<mols_.size(); ++i){
		static_cast<Alpha&>(mols_[i]).set(&mol_alpha_[i]);
	}
}

//**********************************************************************************************
//Sim Atomic
//**********************************************************************************************

template <class AtomT>
class SimAtomic: public SimI, public SimAtomicS<AtomT>{
public:
	//constructors/destructors
	SimAtomic(){};
	SimAtomic(const SimAtomic<AtomT>& sim):SimI(sim),SimAtomicS<AtomT>(sim){};
	~SimAtomic(){};
	
	//operators
	SimAtomic<AtomT>& operator=(const SimAtomic<AtomT>& sim);
	template <class AtomT_> friend std::ostream& operator<<(std::ostream& out, const SimAtomic<AtomT_>& sim);
	
	//member functions
	void clear();
	void resize(unsigned int ts, const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& speciesNames);
};

//operators

template <class AtomT>
SimAtomic<AtomT>& SimAtomic<AtomT>::operator=(const SimAtomic<AtomT>& sim){
	if(DEBUG_SIM>0) std::cout<<"SimAtomic<AtomT>::operator=(const SimAtomic<AtomT>&):\n";
	SimI::operator=(sim);
	SimAtomicS<AtomT>::operator=(sim);
	return *this;
}

template <class AtomT>
std::ostream& operator<<(std::ostream& out, const SimAtomic<AtomT>& sim){
	out<<static_cast<const SimI&>(sim)<<"\n";
	out<<static_cast<const SimAtomicS<AtomT>&>(sim);
	return out;
}

//member functions

template <class AtomT>
void SimAtomic<AtomT>::clear(){
	if(DEBUG_SIM>0) std::cout<<"SimAtomic<AtomT>::clear():\n";
	SimI::clear();
	SimAtomicS<AtomT>::clear();
}

template <class AtomT>
void SimAtomic<AtomT>::resize(unsigned int ts, const std::vector<unsigned int>& nAtoms, const std::vector<std::string>& speciesNames){
	if(DEBUG_SIM>0) std::cout<<"SimAtomic<AtomT>::resize(unsigned int,const std::vector<unsigned int>&,const std::vector<std::string>&):\n";
	SimI::resize(ts);
	SimAtomicS<AtomT>::resize(ts,nAtoms,speciesNames);
}

//**********************************************************************************************
//Sim Molecular
//**********************************************************************************************

template <class AtomT, class MoleculeT>
class SimMol: public SimMolS<MoleculeT>, public SimAtomic<AtomT>{
public:
	//constructors/destructors
	SimMol(){};
	SimMol(const SimMol<AtomT,MoleculeT>& sim):SimMolS<MoleculeT>(sim),SimAtomic<AtomT>(sim){};
	~SimMol(){};
	
	//operators
	SimMol<AtomT,MoleculeT>& operator=(const SimMol<AtomT,MoleculeT>& sim);
	template <class AtomT_,class MoleculeT_> friend std::ostream& operator<<(std::ostream& out, const SimMol<AtomT_,MoleculeT_>& sim);
	
	//member functions
	void clear();
	void resize(unsigned int ts, unsigned int nMols);
};

//operators

template <class AtomT,class MoleculeT>
SimMol<AtomT,MoleculeT>& SimMol<AtomT,MoleculeT>::operator=(const SimMol<AtomT,MoleculeT>& sim){
	if(DEBUG_SIM>0) std::cout<<"SimMol<AtomT,MoleculeT>::operator=(const SimMol<AtomT,MoleculeT>&):\n";
	SimAtomic<AtomT>::operator=(sim);
	SimMolS<MoleculeT>::operator=(sim);
}

template <class AtomT,class MoleculeT>
std::ostream& operator<<(std::ostream& out, const SimMol<AtomT,MoleculeT>& sim){
	out<<static_cast<const SimAtomic<AtomT>&>(sim)<<"\n";
	out<<static_cast<const SimMolS<MoleculeT>&>(sim);
	return out;
}

//member functions

template <class AtomT,class MoleculeT>
void SimMol<AtomT,MoleculeT>::clear(){
	if(DEBUG_SIM>0) std::cout<<"SimMol<AtomT,MoleculeT>::clear():\n";
	SimAtomic<AtomT>::clear();
	SimMolS<MoleculeT>::clear();
}

template <class AtomT,class MoleculeT>
void SimMol<AtomT,MoleculeT>::resize(unsigned int ts, unsigned int nMols){
	if(DEBUG_SIM>0) std::cout<<"SimMol<AtomT,MoleculeT>::resize(unsigned int,const std::vector<unsigned int>&,const std::vector<std::string>&):\n";
	if(ts!=this->timesteps_) throw std::invalid_argument("Invalid number of timesteps.");
	SimMolS<MoleculeT>::resize(ts,nMols);
}

#endif