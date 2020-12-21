#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>
using namespace std;

//funzione che calcola le  medie progressive
vector<double> avprogr( vector<double> v ){

	int N = v.size();
	vector<double> meanprog(N);
	for(int i=0; i<N; i++){
		meanprog[i] = 0;
		for(int j=0; j<=i; j++)
			meanprog[i] += v[j];
		meanprog[i] /= (i+1);
	}
	return meanprog;
}
//funzione che calcola le medie progressive dei quadrati
vector<double> av2progr( vector<double> v ){

	int N = v.size();
	vector<double> mpq(N);
	for(int i=0; i<N; i++){
		mpq[i] = 0;
		for(int j=0; j<=i; j++)
			mpq[i] += pow(v[j],2);
		mpq[i] /= (i+1);
	}
	return mpq;
}

//calcolo delle incertezze statistiche sulle medie progressive
vector<double> error( vector<double> v ){

	vector<double> mp(v.size());
	mp = avprogr(v);	//medie progressive
	vector<double> mpq(v.size());
	mpq = av2progr(v);  // medie progressive dei quadrati
	vector<double> err(v.size());

	for(unsigned int i=0; i<v.size(); i++){
		if(i==0){
			err[i]=0;
		}else{
			err[i] = sqrt( mpq[i]-pow(mp[i],2) );
			err[i] /= sqrt(i);
	}
	return err;
}
//risultati data blocking su file
void stampadati(int L, vector<double> results, string path){
//passo i risultati dell'esperimento, ne faccio medie e errore e li stampo
	int nblock = results.size();
	vector<double> ave(nblock), err(nblock);
	ave = avprogr(results);
	err = error(results);
	
	ofstream Output;
	Output.open(path);
	for(int i=0; i<nblock; i++){
		Output<<(i+1)*L<<" "<<setprecision(6)<<ave[i]<<" "<<err[i]<<endl;
	}
	Output.close();

	return;
}
//calcolo del chi quadrato
double chi_quadro( vector<int> contatore, double E ){

	double chi = 0;
	for(unsigned int i=0; i<contatore.size(); i++)
		chi += pow( (contatore[i]-E),2 )/E;
	return chi;
}
