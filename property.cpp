#include "property.hpp"

namespace serial{
	
	//**********************************************
	// byte measures
	//**********************************************

	template <> unsigned int nbytes(const Name& obj){return obj.name().size()+1;};
	template <> unsigned int nbytes(const AN& obj){return sizeof(unsigned int);};
	template <> unsigned int nbytes(const Species& obj){return sizeof(unsigned int);};
	template <> unsigned int nbytes(const Index& obj){return sizeof(unsigned int);};
	template <> unsigned int nbytes(const Mass& obj){return sizeof(double);};
	template <> unsigned int nbytes(const Charge& obj){return sizeof(unsigned int);};
	template <> unsigned int nbytes(const Position& obj){return nbytes(obj.posn());};
	template <> unsigned int nbytes(const Velocity& obj){return nbytes(obj.velocity());};
	template <> unsigned int nbytes(const Force& obj{return nbytes(obj.force());};
	template <> unsigned int nbytes(const Dipole& obj){return nbytes(obj.dipole());};
	template <> unsigned int nbytes(const Alpha& obj){return nbytes(obj.alpha());};
	template <> unsigned int nbytes(const Chi& obj){return sizeof(double);};
	template <> unsigned int nbytes(const Radius& obj){return sizeof(double);};
	template <> unsigned int nbytes(const ZEff& obj){return sizeof(double);};
	template <> unsigned int nbytes(const Symm& obj){return nbytes(obj.symm());};
	
	//**********************************************
	// packing
	//**********************************************

	template <> void pack(const Name& obj, char* arr){arr[0]=obj.name().size();std::memcpy(arr+1,obj.name().c_str(),obj.name().size());};
	template <> void pack(const AN& obj, char* arr){pack(obj.an(),arr);};
	template <> void pack(const Species& obj, char* arr){pack(obj.specie(),arr);};
	template <> void pack(const Index& obj, char* arr){pack(obj.index(),arr);};
	template <> void pack(const Mass& obj, char* arr){pack(obj.mass(),arr);};
	template <> void pack(const Charge& obj, char* arr){pack(obj.charge(),arr);};
	template <> void pack(const Position& obj, char* arr){pack(obj.posn(),arr);};
	template <> void pack(const Velocity& obj, char* arr){pack(obj.velocity(),arr);};
	template <> void pack(const Force& obj, char* arr){pack(obj.force(),arr);};
	template <> void pack(const Dipole& obj, char* arr){pack(obj.dipole(),arr);};
	template <> void pack(const Alpha& obj, char* arr){pack(obj.alpha(),arr);};
	template <> void pack(const Chi& obj, char* arr){pack(obj.chi(),arr);};
	template <> void pack(const Radius& obj, char* arr){pack(obj.radius(),arr);};
	template <> void pack(const ZEff& obj, char* arr){pack(obj.zeff(),arr);};
	template <> void pack(const Symm& obj, char* arr){pack(obj.symm(),arr);};
	
	//**********************************************
	// unpacking
	//**********************************************

	template <> void unpack(Name& obj, const char* arr){obj.name().resize(arr[0]); std::memcopy(arr+1,obj.name().c_str(),arr[0]);};
	template <> void unpack(AN& obj, const char* arr){unpack(obj.an(),arr);};
	template <> void unpack(Species& obj, const char* arr){unpack(obj.species(),arr);};
	template <> void unpack(Index& obj, const char* arr){unpack(obj.index(),arr);};
	template <> void unpack(Mass& obj, const char* arr){unpack(obj.mass(),arr);};
	template <> void unpack(Charge& obj, const char* arr){unpack(obj.charge(),arr);};
	template <> void unpack(Position& obj, const char* arr){unpack(obj.posn(),arr);};
	template <> void unpack(Velocity& obj, const char* arr){unpack(obj.velocity(),arr);};
	template <> void unpack(Force& obj, const char* arr){unpack(obj.force(),arr);};
	template <> void unpack(Dipole& obj, const char* arr){unpack(obj.dipole(),arr);};
	template <> void unpack(Alpha& obj, const char* arr){unpack(obj.alpha(),arr);};
	template <> void unpack(Chi& obj, const char* arr){unpack(obj.chi(),arr);};
	template <> void unpack(Radius& obj, const char* arr){unpack(obj.radius(),arr);};
	template <> void unpack(ZEff& obj, const char* arr){unpack(obj.zeff(),arr);};
	template <> void unpack(Symm& obj, const char* arr){unpack(obj.symm(),arr);};
	
}