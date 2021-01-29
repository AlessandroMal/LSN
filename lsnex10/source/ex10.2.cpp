#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "random.h"
#include "mpi.h"
#include "header10.2.h"
using namespace std;

int main (int argc, char* argv[]){
   int rk;
   MPI_Init(& argc, & argv);
   MPI_Comm_rank(MPI_COMM_WORLD, & rk);
   initialize();

   if(rk==0){
   	rnd=generaterng();
   	place_in_square(1,Ncity);
   }

   if(rk!=0) cities.resize(Ncity*2);
   MPI_Bcast(&cities.front(), 2*Ncity, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD);
   if(rk!=0) rnd=generaterng(rk);

   for(unsigned int i=0; i<Nstart_pop; i++) population.push_back(gen_path());//popolo gen 0
   for(unsigned int i=0; i<population.size(); i++){
	  l.push_back(compute_L(population[i]));
	  l2.push_back(compute_L2(population[i]));
   }

   unsigned int opt_path[Ncity-1];	//permutazione ottimale (indici città)
   double minl=l[0];
   for(unsigned int i=0; i<Ncity-1; i++) opt_path[i]=population[0][i];

   for(unsigned int i=1; i<population.size(); i++)
	   if(l[i]<minl){
		   for(unsigned int j=0; j<Ncity-1; j++) opt_path[j]=population[i][j];
		   minl=l[i];
	   }
   unsigned int survivor, chosen;
   double p;
   vector<double> sortedL, sortedL2;
   ofstream outfile("../data/path_length"+to_string(rk)+".dat");
   for(unsigned int g=1; g<=Ngen; g++){	//generazioni successive

	   survivor=population.size();
	   for(unsigned int i=0; i<extinction_coeff*survivor; i++) selection_p_propto_L2();
	   selection_meanL();
	   survivor=population.size();

	   if(g%Nmerge==0 && g!=Ngen) merge_continent();

	   for(unsigned int i=survivor; i<Nstart_pop; i++){	//mutazione sopravvissuti e ripopolazione
		   chosen=rnd.Rannyu(0,survivor);
		   population.push_back(duplicate_path(population[chosen]));

		   p=rnd.Rannyu();
		   if(p<0.1) pair_perm_mutation(population[i]);
		   if(p>=0.1 && p<0.2) shift_mutation(population[i]);
		   if(p>=0.2 && p<0.3) block_perm_mutation(population[i]);
		   if(p>=0.3 && p<0.4) invert_mutation(population[i]);
		   if(rnd.Rannyu()>0.5) crossover(survivor);
		   l.push_back(compute_L(population[i]));
		   l2.push_back(compute_L2(population[i]));
		   if(l[i]<minl){
		   	for(unsigned int j=0; j<Ncity-1; j++) opt_path[j]=population[i][j];
		   	minl=l[i];
		   }
	   }

	   sortedL=l;
	   sort(sortedL.begin(), sortedL.end());
	   sortedL2=l2;
	   sort(sortedL2.begin(), sortedL2.end());

	   p=0;
	   for(unsigned int i=0; i<sortedL.size()/2; i++) p+=sortedL[i];
	   p/=sortedL.size()/2;
	   outfile<<p<<" ";
	   p=0;
	   for(unsigned int i=0; i<sortedL2.size()/2; i++) p+=sortedL2[i];
	   p/=sortedL2.size()/2;
	   outfile<<p<<" "<<minl<<endl;
   }
   outfile.close();
   outfile.clear();
   rnd.SaveSeed();

   outfile.open("../data/opt_path"+to_string(rk)+".dat");
   for(unsigned int i=0; i<cities.size(); i++) outfile<<cities[i]<<" ";
   outfile<<endl;
   for(unsigned int i=0; i<Ncity-1; i++) outfile<<opt_path[i]<<" ";
   outfile.close();

   if(rk==0){
	   cout<<"best half average l and l^2 printed in ../data/path_length.dat for each node"<<endl;
	   cout<<"optimized configuration and city coordinates printed in ../data/opt_path.dat for each node"<<endl;
   }
   MPI_Finalize();
   return 0;
}
