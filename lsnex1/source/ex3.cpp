#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include "random.h"
#include "stattools.h"

//domande per galli:
//ma non potevo generare direttamente rng tra -1 e 1 per avere cos uniforme?

double rigettocostheta(){   //funzione per ottenere un coseno casuale (uniforme)
	double x;           //col metodo della generazione di punti nel quadrato
	double y;           //di lato 2 e centrato nell'origine
	Random rnd;
	while(1>0){
		x=rnd.Rannyu(-1,1);
		y=rnd.Rannyu(-1,1);
		if (x*x+y*y<1 && x*x+y*y!=0) break; //tengo il punto solo se
	}					    //all'interno del cerchio inscritto
	return x/sqrt(x*x+y*y); //restituisco direttamente il coseno invece dell'angolo
}

using namespace std;
 
int main (int argc, char *argv[]){

   const int M=100000;       //numero tot rn
   const int N=100;          //numero di blocchi
   const int L=M/N;          //numero di rn per blocco (implement. per M divisibile per N)

   vector<double> av(N);     //media che questa volta e' PI

   Random rnd=generaterng();

   double somma;	//numero di intersezioniper blocco
   const double a=1;	//lunghezza ago
   const double d=1.4;	//distanza righe

   for(int i=0; i<N; i++){ //ciclo su blocchi
	   somma=0;
	   for(int j=0; j<L; j++)
		   if(abs(rnd.Rannyu(-d/2,d/2))<=abs(a/2*rigettocostheta()))
			   somma+=1;
  	   av[i]=2*a*L/somma/d;
   }
   //ho posto la linea nell'origine e fatto cadere il centro tra -d/2 e d/2

   rnd.SaveSeed();
   stampadati(L,av,"../data/buffon.dat");
   return 0;
}
