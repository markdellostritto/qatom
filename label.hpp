#ifndef LABEL_HPP
#define LABEL_HPP

#include <ostream>

class Label{
protected:
	int specie_;
	int index_;
public:
	//constructors/destructors
	Label():specie_(-1),index_(-1){};
	Label(const Label& l):specie_(l.specie()),index_(l.index()){};
	Label(int specie, int index):specie_(specie),index_(index){};
	~Label(){};
	
	//operators
	Label& operator=(const Label& l){specie_=l.specie(); index_=l.index(); return *this;};
	friend std::ostream& operator<<(std::ostream& out, const Label& l){return out<<"("<<l.specie()+1<<","<<l.index()+1<<")";};
	
	//access
	int& specie(){return specie_;};
	const int& specie()const{return specie_;};
	int& index(){return index_;};
	const int& index()const{return index_;};
	
	//member functions
	void clear(){specie_=-1; index_=-1;};
};

inline bool operator==(const Label& l1, const Label& l2){
	return (l1.specie()==l2.specie() && l1.index()==l2.index());
}

inline bool operator!=(const Label& l1, const Label& l2){
	return !(l1.specie()==l2.specie() && l1.index()==l2.index());
}

#endif
