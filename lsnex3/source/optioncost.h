#include <cmath>
#include "random.h"
Random rnd=generaterng();

//stima diretta del prezzo del bene con European option
double S(double t, double mu, double sigma, double S0){
	return S0*exp((mu-pow(sigma,2)/2)*t+sigma*rnd.Gauss(0,t));
}

//Call: calcolo del prezzo dell'opzione in prospettiva
double C(double r, double T, double K, double St){
	double max=St-K;
	if(max<0)
		max=0;
	return exp(-r*T)*max;
}

//Put
double P(double r, double T, double K, double St){
	double max=K-St;
	if(max<0)
		max=0;
	return exp(-r*T)*max;
}

//stima indiretta dell'asset price attraverso il cammino del GBM
double Sstep(double deltat, double mu, double sigma, double St){
	return St*exp((mu-pow(sigma,2)/2)*deltat+sigma*rnd.Gauss(0,1)*sqrt(deltat));
}
