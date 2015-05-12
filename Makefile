XDRINC = /opt/xdrfile/include/xdrfile/
XDRLIB = /opt/xdrfile/lib/
Closest-solvent : main.o string_operate.o
	g++ -O2 -I${XDRINC} -o traj-warp main.o string_operate.o -lxdrfile  -L${XDRLIB}
main.o : main.cpp  read_ndx.h string_operate.h 
	g++ -O2 -I${XDRINC} -c main.cpp
string_operate.o : string_operate.h string_operate.cpp
	g++ -O2 -c string_operate.cpp
install:
	cp traj-warp ~/bin/
clean :
	rm *.o
