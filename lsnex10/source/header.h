#ifndef __Header__
#define __Header__

#include<vector>
#include <cmath>
#include"random.h"
using namespace std;

unsigned int Ncity;			//numero città
unsigned int Nstep;			//step per ogni temperatura
unsigned int Tstep;			//numero di temperature diverse nell'intervallo
double Tup;				//lim sup T
double Tlow;				//lim inf T
vector<double> cities;			//coordinate xy città
vector<unsigned int*> population;	//percorso
Random rnd;

vector<double> l;			//lunghezze path
vector<double> l2;
unsigned int Nstart_pop;		//numero di individui in popolazione

void initialize();
void place_on_circ(const double, const unsigned int);
void place_regular_polygon(const double, const unsigned int);
void place_in_square(const double, const unsigned int);
double norm2(const unsigned int, const unsigned int);
unsigned int* gen_path();
unsigned int* duplicate_path(const unsigned int*);
double compute_L(const unsigned int*);
double compute_L2(const unsigned int*);
void selection_meanL();
void selection_p_propto_L2();
void pair_perm_mutation(unsigned int*);
void shift_mutation(unsigned int*);
void block_perm_mutation(unsigned int*);
void invert_mutation(unsigned int*);
unsigned int BC_distance(unsigned int, unsigned int);
unsigned int BC(unsigned int);
void crossover(unsigned int);

double p_Boltzmann(const unsigned int*, const unsigned int*, double);

void initialize(){
	ifstream init("../data/initialize.dat");
	init >> Ncity >> Tlow >> Tup >> Tstep >> Nstep ;
	init.close();
}


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
	for(unsigned int i=0; i<n; i++){
		cities.push_back(r*cos(i*step));
		cities.push_back(r*sin(i*step));
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

void selection_meanL(){
	double av=0;
	double elite=0.05*Nstart_pop;
	for(auto el : l) av+=el;
	av/=l.size();	//calcolato media della generazione

	for(int i=0; i<int(l.size()); i++) //estinguo quelli sopra la media
		if(l[i]>av && l.size()>elite){
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
	while(BC_distance(pos1,pos2) < m){		//solo sottovettori separati
		pos1=rnd.Rannyu(0,Ncity-1);
		pos2=rnd.Rannyu(0,Ncity-1);
		m=rnd.Rannyu(1,int(Ncity/2));
	}
	unsigned int subvec1[m];
	for(unsigned int i=0; i<m; i++) subvec1[i]=path[BC(pos1+i)]; //copio sottovettore uno
	for(unsigned int i=0; i<m; i++) path[BC(pos1+i)]=path[BC(pos2+i)];
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
	unsigned int d=abs(int(ind1)-int(ind2));
	if(d > (Ncity-1)/2){
		d= Ncity-1 - d;
		if(Ncity%2 == 1) d++;
	}
	return d;
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

double p_Boltzmann(const unsigned int* current, const unsigned int* trial, double T){
	double beta=1/T;
	double b=exp( -beta*(compute_L(trial) - compute_L(current)) );
	if(b>1) b=1;
	return b;
}
#endif // __Header__
