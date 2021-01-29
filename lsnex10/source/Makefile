CC = g++
CFLAGS = -Wall -O3 --std=c++11

ex10.2.exe : ex10.2.cpp random.o
	mpicxx ex10.2.cpp random.o -o ex10.2.exe $(CFLAGS)
ex10.1.exe : ex10.1.o random.o
	$(CC) random.o ex10.1.o -o ex10.1.exe
ex10.1.o : ex10.1.cpp
	$(CC) -c ex10.1.cpp -o ex10.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

clean :
	rm *.o *.exe
