CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: ex1.exe qmc1d

qmc1d: qmc1d.o random.o
	$(CC) random.o qmc1d.o -o qmc1d
ex1.exe : ex1.o random.o stattools.o
	$(CC) random.o stattools.o ex1.o -o ex1.exe
qmc1d.o : qmc1d.cpp
	$(CC) -c qmc1d.cpp -o qmc1d.o $(CFLAGS)
ex1.o : ex1.cpp
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
stattools.o : stattools.cpp stattools.h
	$(CC) -c stattools.cpp -o stattools.o $(CFLAGS)

clean :
	rm *.o *.exe qmc1d potential.dat kinetic.dat probability.dat

#ex2.exe : ex2.o random.o stattools.o
#	$(CC) random.o stattools.o ex2.o -o ex2.exe -O2 -larmadillo -llapack -lblas
#ex2.o : ex2.cpp
#	$(CC) -c ex2.cpp -o ex2.o $(CFLAGS)
