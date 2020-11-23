#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include"random.h"
using namespace std;
 
int main (int argc, char *argv[]){
	
	Random rnd=generaterng();		//imposto il generatore di numeri casuali

	const int n=1E4;
	const vector<int> N = {0,1,2,10,100};	//DA METTERE IN ORDINE CRESCENTE
						//LO ZERO E' NECESSARIO

	vector<double> stand(n,0.);		//dati per S_N uniforme
	vector<double> expo(n,0.);		//dati per S_N esponenziale
	vector<double> lorentz(n,0.);		//dati per S_N lorentziani
	//settarli a zero mi e` utile per scrivere in maniera compatta il codice dopo

	ofstream st("../data/standard.dat");
	ofstream ex("../data/exponential.dat");
	ofstream lor("../data/lorentzian.dat");

	for(unsigned int k=1;k<N.size();k++){
		//con questo metodo sommo sugli rn gia` generati in precedenza
		for(int j=0;j<N[k]-N[k-1];j++)
			for(int i=0;i<n;i++){
				//ciclo per nvalues valori casuali di ciascun istogramma
				stand[i]+=rnd.Rannyu();
				expo[i]+=rnd.Exppdf(1);
				lorentz[i]+=rnd.Lorentz(0,1);
			}

		for(int j=0;j<n;j++){ //scrivo i dati di ogni istogramma sulle righe
			st<<stand[j]<<" ";
			ex<<expo[j]<<" ";
			lor<<lorentz[j]<<" ";
		}
		st<<endl;
		ex<<endl;
		lor<<endl;
	}
	st.close();
	ex.close();
	lor.close();
	rnd.SaveSeed();

   return 0;
}
