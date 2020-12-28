CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: equilibration.exe MolDyn_NVE.exe
equilibration.exe : equilibration.o
	$(CC) equilibration.o -o equilibration.exe
MolDyn_NVE.exe : MolDyn_NVE.o
	$(CC) MolDyn_NVE.o -o MolDyn_NVE.exe
equilibration.o : equilibration.cpp
	$(CC) -c equilibration.cpp -o equilibration.o
MolDyn_NVE.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
clean :
	rm *.o *.exe
