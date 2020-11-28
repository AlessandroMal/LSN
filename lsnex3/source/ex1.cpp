#include<iostream>
#include<cmath>
#include<vector>
#include"random.h"
#include"stattools.h"
#include"optioncost.h"

int main (int argc, char *argv[]){

   const int M=4E6;	//numero tot rn
   const int N=400;	//numero di blocchi
   const int L=M/N;	//numero di rn per blocco (implement. per M divisibile per N)

   double S0=100;	//asset price at t=0
   double K=100;	//strike price
   double T=1;		//delivery time
   double r=0.1;	//tasso d'interesse 'risk free'
   double sigma=0.25;	//volatilita'
   int Nsteps=100;	//numero di intervalli in [0,T] per il cammino del GBM

   vector<double> callav(N),putav(N),callavstep(N),putavstep(N);	//risultati dei blocchi
   double sumc,sump,sumcstep,sumpstep,s,sst;				//accumulatori

   for(int i=0;i<N;i++){ //ciclo sui blocchi
	   sumc=0;
	   sump=0;
	   sumcstep=0;
	   sumpstep=0;

	   for(int j=0;j<L;j++){//ciclo su elementi del blocco
		   s=S(T,r,sigma,S0);	//salvo il calcolo dell'asset price
		   sumc+=C(r,T,K,s);	//e faccio il calcolo diretto per call
		   sump+=P(r,T,K,s);	//e put

		   sst=S0;		//inizializzo il cammino al prezzo di partenza
		   for(int k=0;k<Nsteps;k++)	//e effettuo i passi per arrivare a T
			   sst=Sstep(T/Nsteps, r, sigma, sst);

		   sumcstep+=C(r,T,K,sst);	//cosi' calcolo call e put sviluppando
		   sumpstep+=P(r,T,K,sst);	//il cammino del GBM discretizzato
	   }

	   callav[i]=sumc/L;	//calcolo le medie nei blocchi
	   putav[i]=sump/L;
	   callavstep[i]=sumcstep/L;
	   putavstep[i]=sumpstep/L;
   }
   stampadati(L, callav, "../data/directcall.dat");//stampo medie progressive per i blocchi
   stampadati(L, putav, "../data/directput.dat");
   stampadati(L, callavstep, "../data/stepcall.dat");
   stampadati(L, putavstep, "../data/stepput.dat");
   
   rnd.SaveSeed();
   return 0;
}
