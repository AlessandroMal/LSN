#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include "random.h"
#include "header.h"
using namespace std;

int main (/*int argc, char *argv[]*/){
   rnd=generaterng();

//   place_on_circ(1,Ncity);		//piazzo le città
//   place_regular_polygon(1,Ncity);
   place_in_square(1,Ncity);

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
   ofstream outfile("../data/path_length.dat");
   for(unsigned int g=0; g<Ngen; g++){	//generazioni successive

	   survivor=population.size();
	   for(unsigned int i=0; i<extinction_coeff*survivor; i++) selection_p_propto_L2();
	   selection_meanL();
	   survivor=population.size();

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

	   //stampo le media della meta migliore di popolazione
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

   cout<<"best half average l and l^2 printed in ../data/path_length.dat"<<endl;

   outfile.close();
   outfile.clear();
   outfile.open("../data/opt_path.dat");
   for(unsigned int i=0; i<cities.size(); i++) outfile<<cities[i]<<" ";
   outfile<<endl;
   for(unsigned int i=0; i<Ncity-1; i++) outfile<<opt_path[i]<<" ";
   outfile.close();

   cout<<"optimized configuration and city coordinates printed in ../data/opt_path.dat"<<endl;

   rnd.SaveSeed();
   return 0;
}
