#include <iostream>
#include <vector>
#include "random.h"
#include "stattools.h"
#include "header.h"
using namespace std;

int main (/*int argc, char *argv[]*/){
   initialize();
   vector<double> meas;			//M punti generati dal Metropolis
   vector<double> av_meas(nblk);	//medie sui valori dei blocchi
   int acc=0;

   //start
   double alfa;
   double x=0;
   double xtrial;

   double current= hamilt_psi(x)/psiTrial(x,0);
   double trial;
//   cout<<current<<endl;

   unsigned int bin;
   for(int i=1; i<=M; i++){
	   xtrial=x + rnd.Rannyu(-delta, delta);
	   trial=hamilt_psi(xtrial)/psiTrial(xtrial,0);
	   alfa=pdf(xtrial)/pdf(x);
	   if( rnd.Rannyu()<alfa ){	//Acceptance
		   current=trial;
		   x=xtrial;
		   acc+=1;
	   }
//	   if(i%M/200==0) cout<<current<<endl;
	   meas.push_back(current);
	   //histo
	   if(printpsi==1){
		   bin=(x-x0)/bin_size;
	//	   cout<<x<<endl;
		   if(bin<psihist.size() && x>=x0) psihist[bin]+=1;
	   }
   }

//   cout<<"acceptance rate for energy simulation: "<<double(acc)/M<<endl;
//   cout<<"----------------------------------"<<endl<<endl;
//   cout<<endl;
   if(double(acc)/M >0.15 && double(acc)/M <0.85){
	   av_meas=blockaverage(meas,nblk);
	   //   for(auto el : meas) cout<<el<<endl;
	   stampadati(M/nblk,av_meas, "../data/H_av.dat");	//media progressiva
	   if(printpsi==1) printhisto();
   }else{
//	   cout<<"Warning for this run: acceptance rate out of bounds, no data printed"<<endl;
   }

   rnd.SaveSeed();
   return 0;
}
