CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: equilibration.exe ex1.exe
equilibration.exe : equilibration.o random.o
	$(CC) random.o equilibration.o -o equilibration.exe
ex1.exe : ex1.o random.o
	$(CC) random.o ex1.o -o ex1.exe
equilibration.o : equilibration.cpp
	$(CC) -c equilibration.cpp -o equilibration.o
ex1.o : ex1.cpp
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe
