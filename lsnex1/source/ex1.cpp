#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include"random.h"
#include"stattools.h"

using namespace std;
 
int main (int argc, char *argv[]){

   const int M=100000;       //numero tot rn
   const int N=100;          //numero di blocchi. l'ho usato anche per dire su quanti
   			     //sottost di M calcolo il chi2

   const int L=M/N;          //numero di rn per blocco (implement. per M divisibile per N)
   const int subinterv=100;  //sottoinvervalli in cui divido [0,1) per il calcolo di chi2

   vector<double> av(N);     //media
   vector<double> var(N);    //media del quadrato
   vector<double> chi2(N);   //vettore di chi2
   vector<int> contatore(N); //contatore per gli intervalli di appartenenza di r

   double somma,sommavar,r;
   Random rnd=generaterng();

   for(int i=0;i<N;i++){ //ciclo sui blocchi
	   somma=0;
	   sommavar=0;
	   fill(contatore.begin(), contatore.end(), 0); //memset slightly faster ma vabe

	   for(int j=0;j<L;j++){//ciclo su elementi del blocco
		   r=rnd.Rannyu();
		   somma+=r;
		   sommavar+=pow(r-0.5,2);

		   for(int k=0;k<subinterv;k++) //ricerca seriale intervallo di appartenenza r
			   if(r<(k+1.)/subinterv){
				   contatore[k]+=1;
				   break;
			   }
	   }

	   av[i]=somma/L;
	   var[i]=sommavar/L;
	   chi2[i]=chi_quadro(contatore, (double)L/subinterv);
   }
   stampadati(L, av, "../data/meanr.dat");
   stampadati(L, var, "../data/varr.dat");
   
   ofstream Chi;
   Chi.open("../data/chi2.dat");
   for(int i=0; i<N; i++)
	   Chi<<i+1<<" "<<chi2[i]<<endl;
   
   rnd.SaveSeed();
   return 0;
}
