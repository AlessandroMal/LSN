CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: MolDyn_NVE.exe

MolDyn_NVE.exe : MolDyn_NVE.o stattools.o
	$(CC) stattools.o MolDyn_NVE.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
stattools.o : stattools.cpp stattools.h
	$(CC) -c stattools.cpp -o stattools.o $(CFLAGS)

clean :
	rm *.o *.exe
