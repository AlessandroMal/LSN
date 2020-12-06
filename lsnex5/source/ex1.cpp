#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <string>
#include "random.h"
#include "stattools.h"
#include "metropolis.h"

//domande: T in A e' davvero simmetrico quindi?
//devo implementare anche psi210 con accept rate 1/2 o posso lasciare mostrato cosi'?

vector<double> blockaverage(vector<double>, int, int);
void stampacoord(vector<double>, unsigned int, unsigned int, string);

int main (int argc, char *argv[]){
   const int M=1E6;	//numero tot di lanci per caloclare r medio
   const int N=100;	//numero di blocchi
   const int d=3;	//dimensionalita' dello spazio

   vector<double> start(d,0.);	//partenza del RW
   vector<double> r(M);		//M punti generati dal Metropolis
   vector<double> av_r(N);	//medie sui valori dei blocchi

//restituisco gli M valori distribuiti con p(x) tramite algoritmo di Metropolis
   r=metropolis(start, M, pGS, Tunif, Asimm);	//r medio del GS con T uniforme

//nell'ultima posizione ho salvato la frequenza di accettazione
   cout<<"acceptance rate for GS and unform T: "<<r.back()<<endl;
   r.pop_back();
//faccio le medie dei valori nello stesso blocco
   av_r=blockaverage(r,N,d);
   stampadati(M/N,av_r, "../data/av_rGSunif.dat");	//media progressiva
   stampacoord(r,d,20,"../data/av_rGSunif.xyz");	//scrivo punti accettati
   r.resize(0);						//svuoto i vector e distruggo elementi
   av_r.resize(0);

//ripeto il procedimento per le altre simulazioni: stato eccitato e T gaussiano
   r=metropolis(start, M, p210, Tunif, Asimm);
   cout<<"acceptance rate for (210) and unform T: "<<r.back()<<endl;
   r.pop_back();
   av_r=blockaverage(r,N,d);
   stampadati(M/N,av_r, "../data/av_r210unif.dat");
   stampacoord(r,d,20,"../data/av_r210unif.xyz");
   r.resize(0);	
   av_r.resize(0);

   r=metropolis(start, M, pGS, Tgauss, Asimm);
   cout<<"acceptance rate for GS and normal T: "<<r.back()<<endl;
   r.pop_back();
   av_r=blockaverage(r,N,d);
   stampadati(M/N,av_r, "../data/av_rGSgauss.dat");
   stampacoord(r,d,20,"../data/av_rGSgauss.xyz");
   r.resize(0);	
   av_r.resize(0);

   r=metropolis(start, M, p210, Tgauss, Asimm);
   cout<<"acceptance rate for (210) and normal T: "<<r.back()<<endl;
   r.pop_back();
   av_r=blockaverage(r,N,d);
   stampadati(M/N,av_r, "../data/av_r210gauss.dat");
   stampacoord(r,d,20,"../data/av_r210gauss.xyz");
   r.resize(0);	
   av_r.resize(0);

   start[0]=20;	//provo a vedere cosa succede con partenza in zona a p bassa
   r=metropolis(start, M, pGS, Tunif, Asimm);
   cout<<"acceptance rate for GS and uniform T with start far from origin: "<<r.back()<<endl;
   r.pop_back();
   av_r=blockaverage(r,N,d);
   stampadati(M/N,av_r, "../data/av_rGSfarunif.dat");
   stampacoord(r,d,20,"../data/av_rGSfarunif.xyz");
   r.resize(0);	
   av_r.resize(0);

   rnd.SaveSeed();
   return 0;
}

vector<double> blockaverage(vector<double> datas, int N, int d){
//argomenti: dati da mediare, numero di blocchi, dimensione dello spazio
//(l' ultimo mi serve se voglio generalizzare con cammini in altre dimensioni)
	int L=datas.size()/3/N;
	double sum;
	vector<double> r(d);
	vector<double> av(N);

	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
			for(int k=0; k<d; k++)
				r[k]=datas[d*i*L+d*j+k];
			//ho assegnato posizione x,y,z per calcolarne il modulo
			sum += sqrt(norm2(r));
		}
		av[i]=sum/L;
	}
	return av;
}

void stampacoord(vector<double> r, unsigned int d, unsigned int stamprate, string path){
	//input: coordinate punti accettati, dimensionalita',
	//ogni quanti punti ne stampo uno, path del file
	ofstream Output;
	Output.open(path);
	cout<<"printing "<<r.size()/d/stamprate<<" xyz coordinates in "<<path<<endl;
	for(unsigned int i=0; i<r.size()/d/stamprate; i++){
		for(unsigned int j=0; j<d; j++){
			Output<<r[i*d*stamprate+j];
			if(j!=d-1)
				Output<<" ";
		}
		Output<<endl;
	}
	Output.close();
	return;
}
