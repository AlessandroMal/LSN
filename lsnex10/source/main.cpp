#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "random.h"
#include "header.h"
using namespace std;

int main (/*int argc, char *argv[]*/){
   rnd=generaterng();
   initialize();
   population.resize(2);

//   place_on_circ(1,Ncity);		//piazzo le città
//   place_regular_polygon(1,Ncity);
   place_in_square(1,Ncity);

   ofstream outfile("../data/SA_pathlength.dat");
   outfile<<"# T, L_min"<<endl;

   population[0]=gen_path();
   unsigned int opt_path[Ncity-1];	//permutazione ottimale (indici città)
   double minl=compute_L(population[0]);
   for(unsigned int i=0; i<Ncity-1; i++) opt_path[i]=population[0][i];

   double T=Tup;
   double deltaT=(Tup-Tlow)/Tstep;
   double r;
   while(T>Tlow){
	   for(unsigned int i=0; i<Nstep; i++){
		   population[1]=duplicate_path(population[0]);
		   r=rnd.Rannyu();
		   if(r<0.25) pair_perm_mutation(population[1]); //mutazioni trail population
		   if(r>=0.25 && r<0.5) shift_mutation(population[1]);
		   if(r>=0.5 && r<0.75) block_perm_mutation(population[1]);
		   if(r>=0.75) invert_mutation(population[1]);

		   if(rnd.Rannyu()<p_Boltzmann(population[0],population[1],T)){
			   delete population[0];
			   population[0]=population[1];
		   }else{
			   delete population[1];
		   }

		   r=compute_L(population[0]);
		   if(r<minl){
		   	for(unsigned int j=0; j<Ncity-1; j++) opt_path[j]=population[0][j];
		   	minl=r;
		   }
	   }
	   outfile<<T<<" "<<minl<<endl;
	   T-=deltaT;
   }

   cout<<"minimum population length vs T printed in ../data/SA_pathlength.dat"<<endl;
   outfile.close();
   outfile.clear();
   outfile.open("../data/SA_optpath.dat");
   for(unsigned int i=0; i<cities.size(); i++) outfile<<cities[i]<<" ";
   outfile<<endl;
   for(unsigned int i=0; i<Ncity-1; i++) outfile<<opt_path[i]<<" ";
   outfile.close();
   cout<<"optimized configuration and city coordinates printed in ../data/SA_optpath.dat"<<endl;

   rnd.SaveSeed();
   return 0;
}
