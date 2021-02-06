#include <vector>
#include <cmath>
#include "random.h"
using namespace std;

Random rnd=generaterng();

//incertezze statistiche
double error(double sum, double sum2, int N){
	double e = sqrt( (sum2-pow(sum,2))/(N-1) );
	return e;// /(2*sqrt(sum));	
}
//modulo quadro di un vettore
double norm2(vector<double> v){
	double s=0;
	for(unsigned int i=0; i<v.size(); i++)
		s += v[i]*v[i];
	return s;
}

vector<double> rndstepD(vector<double> pos, double a){
	int r = rnd.Rannyu(0,pos.size());	//direzione di spostamento casuale
	if( rnd.Rannyu()<0.5 ) 
		pos[r] += a;
	else
		pos[r] -= a;

	return pos;
}

vector<double> rndstepC(vector<double> pos, double a){
	double * p=rnd.vec1();

	for(unsigned int i=0;i<pos.size();i++)
		pos[i]=pos[i]+a*p[i];

	return pos;
}
