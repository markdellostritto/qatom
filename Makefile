# include directories
DIR1   = /usr/local/include/eigen-eigen-67e894c6cd8f/
INC    = $(DIR1)
# compiler 
CXXFLAGS = -fopenmp -std=gnu++11 -w -O3 $(foreach d, $(INC), -I$d) 
CXX     = g++ 

# objects list
objects_all = sim.o cell.o bonding.o ptable.o \
		idd.o icc.o ewald3D.o thole.o qeq3.o \
		fft.o signal.o interpolation.o \
		string.o file.o vasp.o lammps.o qe.o gaussian.o \
		eigen.o serialize.o math_func.o units.o \
		raman_thole.o ir_qeq.o nlo_atom.o 

objects_nlo = sim.o cell.o bonding.o ptable.o \
		idd.o icc.o ewald3D.o thole.o qeq3.o \
		fft.o signal.o interpolation.o \
		string.o file.o vasp.o lammps.o qe.o gaussian.o \
		eigen.o serialize.o math_func.o units.o

objects_ir = sim.o cell.o bonding.o ptable.o \
		icc.o ewald3D.o qeq3.o \
		fft.o signal.o interpolation.o \
		string.o file.o vasp.o lammps.o qe.o gaussian.o \
		eigen.o serialize.o math_func.o units.o

objects_raman = sim.o cell.o bonding.o ptable.o \
		idd.o ewald3D.o thole.o \
		fft.o signal.o interpolation.o \
		string.o file.o vasp.o lammps.o qe.o gaussian.o \
		eigen.o serialize.o math_func.o units.o

raman: raman_thole.cpp $(objects_raman)
	$(CXX) $(CXXFLAGS) -o raman_thole.exe raman_thole.cpp $(objects_raman) -lfftw3
ir: ir_qeq.cpp $(objects_ir)
	$(CXX) $(CXXFLAGS) -o ir_qeq.exe ir_qeq.cpp $(objects_ir) -lfftw3
nlo: nlo_atom.cpp $(objects_nlo)
	$(CXX) $(CXXFLAGS) -o nlo_atom.exe nlo_atom.cpp $(objects_nlo) -lfftw3

string.o: string.cpp
	$(CXX) $(CXXFLAGS) -c string.cpp
file.o: file.cpp string.hpp
	$(CXX) $(CXXFLAGS) -c file.cpp
math_special.o: math_special.cpp math_const.hpp
	$(CXX) $(CXXFLAGS) -c math_special.cpp
math_func.o: math_func.cpp
	$(CXX) $(CXXFLAGS) -c math_func.cpp
ptable.o: string.cpp
	$(CXX) $(CXXFLAGS) -c ptable.cpp
serialize.o: serialize.cpp
	$(CXX) $(CXXFLAGS) -c serialize.cpp
units.o: units.cpp
	$(CXX) $(CXXFLAGS) -c units.cpp
icc.o: icc.cpp math_const.hpp
	$(CXX) $(CXXFLAGS) -c icc.cpp
idd.o: idd.cpp math_const.hpp
	$(CXX) $(CXXFLAGS) -c idd.cpp
signal.o: signal.cpp math_const.hpp
	$(CXX) $(CXXFLAGS) -c signal.cpp
fft.o: fft.cpp math_const.hpp
	$(CXX) $(CXXFLAGS) -c fft.cpp 
interpolation.o: interpolation.cpp
	$(CXX) $(CXXFLAGS) -c interpolation.cpp
eigen.o: eigen.cpp string.hpp serialize.hpp
	$(CXX) $(CXXFLAGS) -c eigen.cpp
cell.o: cell.cpp math_const.hpp math_special.hpp eigen.hpp serialize.hpp
	$(CXX) $(CXXFLAGS) -c cell.cpp
sim.o: sim.cpp string.hpp ptable.hpp cell.hpp property.hpp
	$(CXX) $(CXXFLAGS) -c sim.cpp
ewald3D.o: ewald3D.cpp sim.hpp cell.hpp log.hpp
	$(CXX) $(CXXFLAGS) -c ewald3D.cpp
thole.o: thole.cpp idd.hpp ewald3D.hpp math_const.hpp ptable.hpp cell.hpp sim.hpp eigen.hpp units.hpp
	$(CXX) $(CXXFLAGS) -c thole.cpp
qeq3.o: qeq3.cpp ptable.hpp atom.hpp molecule.hpp sim.hpp icc.hpp ewald3D.hpp eigen.hpp units.hpp
	$(CXX) $(CXXFLAGS) -c qeq3.cpp
bonding.o: bonding.cpp log.hpp file.hpp cell.hpp label.hpp sim.hpp math_const.hpp
	$(CXX) $(CXXFLAGS) -c bonding.cpp
lammps.o: lammps.cpp sim.hpp string.hpp ptable.hpp list.hpp units.hpp
	$(CXX) $(CXXFLAGS) -c lammps.cpp
vasp.o: lammps.cpp sim.hpp cell.hpp string.hpp units.hpp
	$(CXX) $(CXXFLAGS) -c vasp.cpp
qe.o: lammps.cpp sim.hpp cell.hpp string.hpp units.hpp
	$(CXX) $(CXXFLAGS) -c qe.cpp
gaussian.o: gaussian.cpp property.hpp atom.hpp molecule.hpp string.hpp ptable.hpp sim.hpp cell.hpp units.hpp
	$(CXX) $(CXXFLAGS) -c gaussian.cpp

clean: 
	rm $(objects_all)

