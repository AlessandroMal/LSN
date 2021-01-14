#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include "random.h"
using namespace std;

const unsigned int Ncity=30;		//numero città
const unsigned int Nstart_pop=6;	//popolazione iniziale
const unsigned int Ngen=200;		//numero generazioni
vector<double> cities;			//coordinate xy città
vector<unsigned int*> population;	//popolazione di percorsi
vector<double> l(Nstart_pop);		//misure dei percorsi
vector<double> l2(Nstart_pop);
Random rnd;

unsigned int BC_distance(unsigned int, unsigned int);
unsigned int BC(unsigned int);
	

void place_on_circ(const double r, const unsigned int n){
	double angle;
	for(unsigned int i=0; i<n; i++){
		angle=rnd.Rannyu(0,2*M_PI);
		cities.push_back(r*cos(angle));
		cities.push_back(r*sin(angle));
	}
}

void place_regular_polygon(const double r, const unsigned int n){
	double step=2*M_PI/double(n);
	//cout<<n<<endl;
	for(unsigned int i=0; i<n; i++){
		cities.push_back(r*cos(i*step));
		cities.push_back(r*sin(i*step));
	}
//	for(auto el : cities) cout<<el<<" ";
//	cout<<endl;
}

void place_in_square(const double l, const unsigned int n){
	for(unsigned int i=0; i<n; i++){
		cities.push_back(rnd.Rannyu(-l/2,l/2));
		cities.push_back(rnd.Rannyu(-l/2,l/2));
	}
}

double norm2(const unsigned int cit1, const unsigned int cit2){
	return pow(cities[cit1*2]-cities[cit2*2], 2) + pow(cities[cit1*2+1]-cities[cit2*2+1], 2);
}

unsigned int* gen_path(){
	vector<unsigned int> not_visited(Ncity-1);
	for(unsigned int i=0; i<not_visited.size(); i++) not_visited[i]=i+1;

	unsigned int* path=new unsigned int[Ncity-1];
	unsigned int visiting_ind;
	for(unsigned int i=0; i<Ncity-1; i++){
		visiting_ind=rnd.Rannyu(0,not_visited.size());
		path[i]=not_visited[visiting_ind];
		not_visited.erase(not_visited.begin()+visiting_ind);
	}
	return path;
}

unsigned int* duplicate_path(const unsigned int* path){
	unsigned int* copypath=new unsigned int[Ncity-1];
	for(unsigned int i=0; i<Ncity-1; i++) copypath[i]=path[i];
	return copypath;
}

double compute_L(const unsigned int* path){
	double l=sqrt(norm2(0,path[0]));
	for(unsigned int i=0; i<Ncity-2; i++) l+=sqrt(norm2(path[i],path[i+1]));
	l+=sqrt(norm2(path[Ncity-2],0));

	return l;
}

double compute_L2(const unsigned int* path){
	double l2=norm2(0,path[0]);
	for(unsigned int i=0; i<Ncity-2; i++) l2+=norm2(path[i],path[i+1]);
	l2+=norm2(path[Ncity-2],0);

	return l2;
}

void selection_meanL(){//dovrei sortare?? forse struct
	double av=0;
	for(auto el : l) av+=el;
	av/=l.size();	//calcolato media della generazione

	for(int i=0; i<int(l.size()); i++) //estinguo quelli sopra la media
		if(l[i]>av && l.size()!=1){
			delete population[i];
			population.erase(population.begin()+i);
			l.erase(l.begin()+i);
			l2.erase(l2.begin()+i);
			i--;
		}
}

void selection_p_propto_L2(){
	double sum=0;
	for(auto el : l2) sum+=el;
	double r=rnd.Rannyu(0,sum);

	sum=0;
	unsigned int counter=0;
	while(sum<r){
		sum+=l2[counter];
		counter++;
	}
//	cout<<"a"<<counter-1;
	delete population[counter-1];
//	population[counter-1]=NULL;
	population.erase(population.begin()+counter-1);
	l.erase(l.begin()+counter-1);
	l2.erase(l2.begin()+counter-1);
}

void pair_perm_mutation(unsigned int* path){
	unsigned int pos1=rnd.Rannyu(0,Ncity-1);
	unsigned int pos2=rnd.Rannyu(0,Ncity-1);
	//possibile anche non scambio se le pos sono uguali

	unsigned int c1=path[pos1];
	path[pos1]=path[pos2];
	path[pos2]=c1;
}


void shift_mutation(unsigned int* path){
	unsigned int begin=rnd.Rannyu(0,Ncity-1);	//posizione da cui shifto
	unsigned int m=rnd.Rannyu(1,Ncity-1);		//numero città che shifto
	unsigned int n=rnd.Rannyu(1,Ncity-1);		//#shifting position

	unsigned int shift_pos;
	unsigned int shifted[m];
	vector<unsigned int> rest;
	for(unsigned int i=begin; i<begin+m; i++) shifted[i-begin]=path[BC(i)];
	if(begin < BC(begin+m)){
		shift_pos=begin;
		for(unsigned int i=0; i<begin; i++) rest.push_back(path[i]);
		for(unsigned int i=begin+m; i<Ncity-1; i++) rest.push_back(path[i]);
	}else{
		shift_pos=0;
		for(unsigned int i=BC(begin+m); i<begin; i++) rest.push_back(path[i]);
	}

	for(unsigned int i=0; i<m; i++) path[BC(begin+n+i)]=shifted[i];	//metto il sottovettore shiftato

	unsigned int save;
	for(unsigned int i=0; i<shift_pos+n; i++){
		save=rest[0];
		rest.erase(rest.begin());
		rest.push_back(save);
	}

	for(unsigned int i=0; i<rest.size(); i++) path[BC(begin+n+m+i)]=rest[i];	//incollo tutto il resto in path
}

void block_perm_mutation(unsigned int* path){
	unsigned int pos1=rnd.Rannyu(0,Ncity-1);
	unsigned int m=rnd.Rannyu(1,int(Ncity/2));	//block length
	unsigned int pos2=rnd.Rannyu(0,Ncity-1);
	unsigned int subvec1[m];
	unsigned int save;
	if(BC_distance(pos1,pos2) < m && pos1!=pos2){
		if(pos2<pos1){
			save=pos1;
			pos1=pos2;
			pos2=save;
		}else{	//metto pos1 al vettore con indice più basso tanto è l'algoritmo non cambia
			if(BC(pos2+m) > pos1 && BC(pos2+m) < pos1+m){ //QUI! LINEA 54 IN FILE
				save=pos1;
				pos1=pos2;
				pos2=save;
			}
		}
		for(unsigned int i=0; i<m; i++) subvec1[i]=path[BC(pos1+i)]; //copio sottovettore uno
		for(unsigned int i=0; i<BC_distance(pos1,pos2); i++) path[BC(pos1+i)]=path[BC(pos1+m+i)];
	}else{
		for(unsigned int i=0; i<m; i++) subvec1[i]=path[BC(pos1+i)]; //copio sottovettore uno
		for(unsigned int i=0; i<m; i++) path[BC(pos1+i)]=path[BC(pos2+i)];
	}
	for(unsigned int i=0; i<m; i++) path[BC(pos2+i)]=subvec1[i];
}

void invert_mutation(unsigned int* path){
	unsigned int begin=rnd.Rannyu(0,Ncity-1);	//posizione da cui parte il sottvettore
	unsigned int m=rnd.Rannyu(0,Ncity-1);		//numero città che inverto (lunghezza sottovettore)

//	cout<<begin<<" "<<m<<endl;
	unsigned int inv[Ncity-1];
	for(int i=begin+m-1; i>=int(begin); i--){
	//	cout<<begin+m-1-i<<" "<<i<<endl;
	       	inv[begin+m-1 -i]=path[BC(i)];	//copio sottovettore invertito
	}
	for(unsigned int i=begin; i<begin+m; i++) path[BC(i)]=inv[i-begin];		//incollo in path
}

unsigned int BC_distance(unsigned int ind1, unsigned int ind2){
	if(ind1<Ncity/2 && ind2>=Ncity/2) ind1+=Ncity-1;
	if(ind2<Ncity/2 && ind1>=Ncity/2) ind2+=Ncity-1;
	return abs(int(ind1-ind2));
}

unsigned int BC(unsigned int ind){
	if(ind>=Ncity-1) ind=ind%(Ncity-1);
	return ind;
}

void crossover(unsigned int newfromhere){
	if(population.size()-newfromhere >=2){
		unsigned int p1_ind=rnd.Rannyu(newfromhere, population.size());
		unsigned int p2_ind=p1_ind;
		while(p1_ind==p2_ind) p2_ind=rnd.Rannyu(newfromhere, population.size());	//scelti due genitori diversi
		unsigned int* path1=population[p1_ind];
		unsigned int* path2=population[p2_ind]; //pointer ai genitori

		unsigned int start1[Ncity-1], start2[Ncity-1];			//parti terminali tagliate
		for(unsigned int i=0; i<Ncity-1; i++){
			start1[i]=path1[i];
			start2[i]=path2[i];
		}

		unsigned int cut_ind=rnd.Rannyu(1,Ncity-1);			//posizione in cui taglio i percorsi
		unsigned int fill=0;
		for(unsigned int i=0; i<Ncity-1; i++){	//riempio path1
			if(fill < Ncity-1-cut_ind){
				for(unsigned int j=cut_ind; j<Ncity-1; j++)
					if(start1[j]==start2[i]){
						path1[cut_ind+fill]=start1[j];
						fill++;
						break;
					}
			}else{
				break;
			}
		}

		fill=0;
		for(unsigned int i=0; i<Ncity-1; i++){	//riempio path2
			if(fill < Ncity-1-cut_ind){
				for(unsigned int j=cut_ind; j<Ncity-1; j++)
					if(start2[j]==start1[i]){
						path2[cut_ind+fill]=start2[j];
						fill++;
						break;
					}
			}else{
				break;
			}
		}

	}
}

int main (/*int argc, char *argv[]*/){
   rnd=generaterng();
   unsigned int opt_path[Ncity-1];	//permutazione ottimale (indici città)

//   place_on_circ(1,Ncity);		//piazzo le città
   place_regular_polygon(1,Ncity);
//   place_in_square(1,Ncity);

   for(unsigned int i=0; i<Nstart_pop; i++) population.push_back(gen_path());	//popolo gen 0
   for(unsigned int i=0; i<population.size(); i++){
	  l[i]=compute_L(population[i]);
	  l2[i]=compute_L2(population[i]);
   }

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
   for(unsigned int g=0; g<Ngen; g++){

	   survivor=population.size();
	   for(unsigned int i=0; i<0.8*survivor; i++) selection_p_propto_L2();
//	   selection_meanL();
	   survivor=population.size();

	   for(unsigned int i=survivor; i<Nstart_pop; i++){	//mutazione sopravvissuti e ripopolazione
		   chosen=rnd.Rannyu(0,survivor);
		   population.push_back(duplicate_path(population[chosen]));

	//	   p=rnd.Rannyu();
	//	   if(p<0.1) pair_perm_mutation(population[i]);
	//	   if(p>=0.1 && p<0.2) shift_mutation(population[i]);
	//	   if(p>=0.2 && p<0.3) block_perm_mutation(population[i]);
	//	   cout<<"ala"<<endl;
		   block_perm_mutation(population[i]);
	//	   if(p>=0.3 && p<0.4) invert_mutation(population[i]);
	//	   if(rnd.Rannyu()>0.5) crossover(survivor);
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

	   for(auto el : population){
	   	for(unsigned int i=0; i<Ncity-1; i++) outfile<<el[i]<<" ";
	  	 outfile<<endl;
	   }
	   p=0;
	   for(unsigned int i=0; i<sortedL.size()/2; i++) p+=sortedL[i];
	   p/=sortedL.size()/2;
	   outfile<<p<<" ";
	   p=0;
	   for(unsigned int i=0; i<sortedL2.size()/2; i++) p+=sortedL2[i];
	   p/=sortedL2.size()/2;
	   outfile<<p<<endl;
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
