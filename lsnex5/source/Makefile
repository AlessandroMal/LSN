CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: ex1.exe

ex1.exe : ex1.o random.o stattools.o metropolis.o
	$(CC) random.o stattools.o metropolis.o ex1.o -o ex1.exe
ex1.o : ex1.cpp
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
stattools.o : stattools.cpp stattools.h
	$(CC) -c stattools.cpp -o stattools.o $(CFLAGS)
metropolis.o : metropolis.cpp metropolis.h
	$(CC) -c metropolis.cpp -o metropolis.o $(CFLAGS)

clean :
	rm *.o *.exe
