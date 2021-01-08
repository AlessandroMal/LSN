CC = g++
CFLAGS = -Wall -O3 --std=c++11

qmc1d.o : qmc1d.cpp
	$(CC) -c qmc1d.cpp -o qmc1d.o $(CFLAGS)
qmc1d: qmc1d.o
	$(CC) qmc1d.o -o qmc1d
clean:
	rm *.o qmc1d potential.dat kinetic.dat probability.dat
