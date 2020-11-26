#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include"random.h"
#include"stattools.h"

using namespace std;

//funzione integranda
double integranda(double x){
	return M_PI/2*cos( M_PI*x/2);
}
//distribuzione per importance sampling: voglio gli rn distribuiti come 2*(1-r)
//che sicuramente assomiglia di piu' alla mia funzione integranda
double pdf(double y){
	return 1-sqrt(1-y);
}

int main (int argc, char *argv[]){

   const int M=1E5;	//numero tot rn
   const int N=100;	//numero di blocchi

   const int L=M/N;	//numero di rn per blocco (implement. per M divisibile per N)
   double a = 0, b = 1;	//estremi di integrazione

   double r,r_imp,sum,sum_imp;
   vector<double> I(N);			//integrale uniform sampling
   vector<double> I_imp(N);		//integrale importance sampling
   Random rnd=generaterng();
	
	for(int i=0; i<N; i++){
		sum = 0;
		sum_imp = 0;
		for(int j=0; j<L; j++){
			r = rnd.Rannyu();		//uniform sampling
			sum += integranda(r);	
			r_imp = pdf(rnd.Rannyu());	//importance sampling
			sum_imp += integranda(r_imp)/(2*(1-r_imp));
		}
		I[i] = (b-a)*sum/L;
		I_imp[i] = (b-a)*sum_imp/L;
	}

	stampadati(L,I, "../data/Iuniform.dat");
	stampadati(L,I_imp, "../data/IIS.dat");

	rnd.SaveSeed();
	return 0;
}

