#include <vector>
#include <cmath>
#include "random.h"
using namespace std;
Random rnd=generaterng();
const double a0=1;
//const double a0=5.29177772E-11;
//lascio commentato se voglio le grandezze in unita' di raggio di bohr

vector<double> metropolis(vector<double> start, int iterations, double (*p)(vector<double>), vector<double> (*T)(vector<double>, double), bool (*A)(vector<double>, vector<double>, double (*parg)(vector<double>)) ){
//la funzione ha come argomenti nell'ordine:
//start: il vettore con le coordinate di partenza
//iterations: il numero di iterazioni e quindi il numero di triplette di coordinate restituite
//p: la funzione pdf con cui voglio distribuire i punti
//T: probabilita' di transizione di prova (con argomenti punto attuale e lunghezza caratter.)
//A: acceptance rule (con argomenti punto attuale, di prova e p(x))
	vector<double> x=start;
	vector<double> xnew;
	vector<double> acc;
	int accrate=0;

	for(int i=0; i<iterations; i++){
		xnew=(*T)(x,a0);	//Transizione di prova
		if((*A)(x,xnew,p)){	//Acceptance
			x=xnew;
			accrate+=1;
		}

		for(auto el : x)
			acc.push_back(el);
	}

	acc.push_back(double(accrate)/iterations);
	return acc;
}

bool Asimm(vector<double> x, vector<double> xnew, double (*p)(vector<double>)){
//acceptance rule con T simmetrica
	if((*p)(xnew)>=(*p)(x)){
		return true;
	}else{
		if((*p)(xnew)/(*p)(x)>rnd.Rannyu()){
			return true;
		}else{
			return false;
		}
	}
}

vector<double> Tunif(vector<double> pos, double a){
	for(unsigned int i=0; i<pos.size(); i++)
		pos[i]+=rnd.Rannyu(-1.2*a,1.2*a);
	return pos;
}

vector<double> Tgauss(vector<double> pos, double a){
	for(unsigned int i=0; i<pos.size(); i++)
		pos[i]+=rnd.Gauss(0,a*3/4);
	return pos;
}

//modulo quadro di un vettore
double norm2(vector<double> v){
	double s=0;
	for(unsigned int i=0; i<v.size(); i++)
		s += v[i]*v[i];
	return s;
}

double pGS(vector<double> x){ //pdf che voglio campionare
	return 1/(M_PI*pow(a0,3))*exp(-2/a0*sqrt(norm2(x)));
}

double p210(vector<double> x){ //pdf che voglio campionare
	return pow(x[2],2)/(32*M_PI*pow(a0,5))*exp(-sqrt(norm2(x))/a0);
}
