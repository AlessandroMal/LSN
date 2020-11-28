CC = g++
CFLAGS = -Wall -O3 --std=c++11

all: ex1.exe

ex1.exe : ex1.o random.o stattools.o
	$(CC) random.o stattools.o ex1.o -o ex1.exe
ex1.o : ex1.cpp
	$(CC) -c ex1.cpp -o ex1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
stattools.o : stattools.cpp stattools.h
	$(CC) -c stattools.cpp -o stattools.o $(CFLAGS)

clean :
	rm *.o *.exe
