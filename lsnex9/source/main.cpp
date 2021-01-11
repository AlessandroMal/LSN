#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include "random.h"
using namespace std;

const unsigned int Ncity=10;		//numero città
vector<double> cities;			//coordinate xy città
vector<unsigned int*> population;	//popolazione di percorsi
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
	unsigned int tot_to_visit=cities.size()/2 -1;
	vector<unsigned int> not_visited(tot_to_visit);
	for(unsigned int i=0; i<not_visited.size(); i++) not_visited[i]=i+1;

	unsigned int* path=new unsigned int[tot_to_visit];
	unsigned int visiting_ind;
	for(unsigned int i=0; i<tot_to_visit; i++){
		visiting_ind=rnd.Rannyu(0,not_visited.size());
		path[i]=not_visited[visiting_ind];
		not_visited.erase(not_visited.begin()+visiting_ind);
	}
	return path;
}

double compute_L(const unsigned int* path){
	double l=0;
	for(unsigned int i=0; i<Ncity-1; i++) l+=sqrt(norm2(path[i],path[i+1]));
	l+=sqrt(norm2(0,path[0]));
	l+=sqrt(norm2(path[Ncity-1],0));

	return l;
}

double compute_L2(const unsigned int* path){
	double l2=0;
	for(unsigned int i=0; i<Ncity-1; i++) l2+=norm2(path[i],path[i+1]);
	l2+=norm2(0,path[0]);
	l2+=norm2(path[Ncity-1],0);

	return l2;
}

unsigned int selection_meanL(vector<double> l){//dovrei sortare?? forse struct
	double av=0;
	unsigned int start=l.size()-population.size();
	unsigned int end=l.size();
	for(unsigned int i=start; i<end; i++) av+=l[i];
	av/=population.size();	//calcolato media della generazione

	unsigned int dead=0;
	for(unsigned int i=0; i<population.size(); i++)
		if(l[start+i]>av){
			population.erase(population.begin()+i);
			dead++;
		}
	return dead;
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
	unsigned int m=rnd.Rannyu(0,Ncity-1);		//numero città che shifto
	unsigned int n=rnd.Rannyu(0,Ncity-1);		//#shifting position

//	cout<<"robe"<<begin<<m<<n<<endl;

	unsigned int shift_pos;
	unsigned int shifted[m];
	unsigned int rest[Ncity-1-m];
	for(unsigned int i=begin; i<begin+m; i++) shifted[i-begin]=path[BC(i)];
	if(begin < BC(begin+m)){
		shift_pos=begin;
		for(unsigned int i=0; i<begin; i++) rest[i]=path[i];
		for(unsigned int i=begin+m; i<Ncity-1; i++) rest[i-m]=path[i];
	}else{
		shift_pos=0;
		for(unsigned int i=BC(begin+m); i<begin; i++) rest[i-BC(begin+m)]=path[i];
	}

	for(unsigned int i=0; i<m; i++) path[BC(begin+n)]=shifted[i];	//metto il sottovettore shiftato

	shift_pos=(shift_pos+n)%(Ncity-1-m);
	vector<unsigned int> r;
	for(unsigned int i=0; i<n; i++)


	for(unsigned int i=0; i<begin; i++) rest[i]=path[i];				//copio quello che c'è all'inizio
	for(unsigned int i=begin; i<m; i++) shifted[i]=path[begin+i];			//copio quello che c'è all'inizio

	for(unsigned int i=0; i<begin; i++) cout<<shifted[i]<<endl;

	for(unsigned int i=begin; i<begin+n; i++) shifted[BC(i)]=path[BC(m+i)];		//shifto all'indietro quello che c'è dopo, con BC

	for(unsigned int i=begin; i<begin+n; i++) cout<<shifted[BC(i)]<<endl;


	for(unsigned int i=0; i<Ncity-1; i++) path[i]=shifted[i];
}

void block_perm_mutation(unsigned int* path){
	unsigned int pos1=rnd.Rannyu(0,Ncity-1);
	unsigned int m=rnd.Rannyu(0,(Ncity-1)/2);	//block length
	unsigned int pos2=rnd.Rannyu(0,Ncity-1);

	unsigned int subvec1[m];
	for(unsigned int i=0; i<m; i++) subvec1[i]=path[BC(pos1+i)]; //copio sottovettore uno

	if(BC_distance(pos1,pos2<m)){
		if(BC(pos1+m)>pos2){
			for(unsigned int i=0; i<BC_distance(pos1,pos2); i++) path[BC(pos1+i)]=path[BC(pos1+m+i)];
			for(unsigned int i=0; i<m; i++) path[BC(pos2+i)]=subvec1[i];
		}
		if(BC(pos2+m)>pos1){
			for(unsigned int i=0; i<BC_distance(pos1,pos2); i++) path[BC(pos1+m-BC_distance(pos1,pos2)+i)]=path[BC(pos2+i)];
			for(unsigned int i=0; i<m; i++) path[BC(pos2+i)]=subvec1[i];
		}
	}else{
		for(unsigned int i=0; i<m; i++) path[BC(pos1+i)]=path[BC(pos2+i)];
		for(unsigned int i=0; i<m; i++) path[BC(pos2+i)]=subvec1[i];
	}
}

void invert_mutation(unsigned int* path){
	unsigned int begin=rnd.Rannyu(0,Ncity-1);	//posizione da cui parte il sottvettore
	unsigned int m=rnd.Rannyu(0,Ncity-1);		//numero città che inverto (lunghezza sottovettore)

	unsigned int inv[Ncity-1];
	for(unsigned int i=begin+m-1; i>=begin; i--) inv[begin+m-1 -i]=path[BC(i)];	//copio sottovettore invertito
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

void crossover(){
	unsigned int p1_ind=rnd.Rannyu(0, population.size());
	unsigned int p2_ind=p1_ind;
	while(p1_ind==p2_ind) p2_ind=rnd.Rannyu(0, population.size());	//scelti due genitori diversi
	unsigned int* path1=population[p1_ind];
	unsigned int* path2=population[p2_ind]; //pointer ai genitori

	unsigned int start1[Ncity-1], start2[Ncity-1];			//parti terminali tagliate
	for(unsigned int i=0; i<Ncity-1; i++){
		start1[i]=path1[i];
		start2[i]=path2[i];
	}

	unsigned int cut_ind=rnd.Rannyu(1,Ncity-1);			//posizione in cui taglio i percorsi
	cout<<cut_ind<<endl;
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

int main (/*int argc, char *argv[]*/){
   const unsigned int Nstart_pop=6;	//popolazione iniziale
   const unsigned int Ngen=10;		//numero generazioni
   rnd=generaterng();
   vector<double> l;			//misure dei percorsi
   vector<double> l2;
   vector<unsigned int> opt_path(Ncity-1);	//permutazione ottimale (indici città)

   place_on_circ(1,Ncity);	//piazzo le città
/*
   for(auto el : populations){
	   for(unsigned int i=0; i<Ncity-1; i++) cout<<el[i];
	   cout<<endl;
   }
*/

   unsigned int dead=Nstart_pop;
   double p;
   for(unsigned int g=0; g<Ngen; g++){
	   if(g==0) for(unsigned int i=0; i<dead; i++) population.push_back(gen_path());	//ripopolo generazione nuova
//	   cout<<"no"<<endl;
//	   for(auto el : population){ //computo lunghezze percorsi
//		   l.push_back(compute_L(el));
//		   l2.push_back(compute_L2(el));
//	   }
//	   cout<<"si"<<endl;

	   for(auto el : population){
		   for(unsigned int i=0; i<Ncity-1; i++) cout<<el[i];
		   cout<<" ";
	   }
	   cout<<endl;

//	   dead=selection_meanL(l); //selezione (restituisce numero di elementi che non l'hanno superata
	   for(unsigned int i=0; i<population.size(); i++){	//mutazioni elementi sopravvissuti
//		   p=rnd.Rannyu();
//		   if(p<0.1) pair_perm_mutation(population[i]);
//		   if(p>=0.1 && p<0.2) shift_mutation(population[i]);
//		   shift_mutation(population[i]);
//		   if(p>=0.2 && p<0.3) block_perm_mutation(population[i]);
//		   block_perm_mutation(population[i]);
//		   if(p>=0.3 && p<0.4) invert_mutation(population[i]);
	   }
//	   if(rnd.Rannyu()>0.5) crossover();

	   for(auto el : population){
		   for(unsigned int i=0; i<Ncity-1; i++) cout<<el[i];
		   cout<<" ";
	   }
	   cout<<endl<<endl;

   }

   rnd.SaveSeed();
   return 0;
}
