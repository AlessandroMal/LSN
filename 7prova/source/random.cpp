#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;

Random :: Random(){}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("../data/seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t); //Box Muller
   return mean + x * sigma;  //trasf lineare
}

double Random :: Exppdf(double lambda){ //uso l'invertibilita' della funzione cumulativa
	double y=Rannyu();
	return -1/lambda*log(1-y);
}

double Random :: Lorentz(double mean, double gamma){
	double y=Rannyu();
	double x=gamma*tan(M_PI*(y-0.5));
	return mean+x;
}

double * Random :: sph1(){
	double * p=new double[2];
	p[0]=acos(2*Rannyu()-1);	//theta tra 0 e pi
	p[1]=Rannyu(0,2*M_PI);		//phi tra 0 e 2pi
	return p;
}

double * Random :: vec1(){
	double * versore=new double[3];
	double * punto=sph1();
	versore[0]=sin(punto[0])*cos(punto[1]);
	versore[1]=sin(punto[0])*sin(punto[1]);
	versore[2]=cos(punto[0]);

	return versore;
}


double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}

double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

Random generaterng(){

	Random rnd;

	int seed[4];
 	int p1, p2;
  	ifstream Primes("../data/Primes");
  	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
  	} else cerr << "PROBLEM: Unable to open Primes" << endl;
  	Primes.close();

  	ifstream input("../data/seed.in");
  	string property;
  	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	return rnd;
}

Random generaterng(bool restart){
//implemento la ripartenza dalla configurazione precedente
	Random rnd;

	int seed[4];
 	int p1, p2;
  	ifstream Primes("../data/Primes");
  	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
  	} else cerr << "PROBLEM: Unable to open Primes" << endl;
  	Primes.close();

  	ifstream input;

	if(restart==0)	input.open("../data/seed.in");
	else input.open("../data/seed.out");

  	if (input.is_open()){
		while ( !input.eof() ){
			input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
			rnd.SetRandom(seed,p1,p2);
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	return rnd;
}
