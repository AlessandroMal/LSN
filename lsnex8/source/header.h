#include <iostream>
#include <fstream>
//#include <ostream>
#include <cmath>
#include <vector>
#include "random.h"
#include "stattools.h"
using namespace std;

//parameters
double m,hbar;
double mu,sigma,delta;

//simulation
int M,nblk;
bool printpsi;
Random rnd;

//histo
double x0,x1,bin_size;
vector<int> psihist;

//functions
void initialize();
vector<double> blockaverage(vector<double>, int);
double pdf(double);
double psiTrial(double,int);
double V(double);
double hamilt_psi(double);
void printhisto();

void initialize(){
	rnd=generaterng();
	hbar=1;
	m=1;

	ifstream params;
	params.open("../data/inputVMC.dat");
	params>>mu;
	params>>sigma;
	params>>delta;
	params>>M;
	params>>nblk;
	params>>printpsi;
	params.close();
//	cout<<"Gaussians with mu="<<mu<<"    sigma="<<sigma<<endl;
//	cout<<"total steps: "<<M<<"    nblocks: "<<nblk<<endl;
//	mu=abs(mu); //operazione già fatta in annealing.py
//	sigma=abs(sigma);
//	delta=abs(delta);

	if(printpsi==1){ //preparo istogramma
		int nbins=200;
		psihist.resize(nbins);
		fill(psihist.begin(), psihist.end(), 0.);
		x0=-mu-3*sigma;
		x1= mu+3*sigma;
		bin_size=(x1-x0)/nbins;
	}
}

vector<double> blockaverage(vector<double> datas, int N){
//argomenti: dati da mediare, numero di blocchi
	int L=datas.size()/N;
	double sum;
	vector<double> av(N);
//	cout<<L<<endl;
//	cout<<N<<endl;

	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
//			if(j==L-1) cout<<datas[i*L+j]<<endl;
			sum += datas[i*L+j];
		}
		av[i]=sum/L;
	//	cout<<av[i]<<endl;//d
	}
	//cout<<endl;
	return av;
}

double pdf(double x){//if real
	double f=psiTrial(x,0);
	return f*f;
}

double psiTrial(double x, int der){
	double e_p=exp( -0.5*pow( (x+mu)/sigma, 2) );
	double e_m=exp( -0.5*pow( (x-mu)/sigma, 2) );
	double psi = e_m+e_p;
//	if(der==2) psi=pow( (x-mu)/(sigma*sigma), -2)*e_m + pow( (x+mu)/(sigma*sigma), -2)*e_p -(e_m+e_p)/(sigma*sigma);
	if(der==2){
		double psi1= pow(sigma, -2) *( (mu-x)*psi-2*mu*e_p );
		psi=pow(sigma, -2) *( (mu-x)*psi1 - psi + 2*mu*pow(sigma, -2)*(x+mu)*e_p );
	}
	return psi;
}

double V(double x){
	return x*x*(x*x - 2.5);
}

double hamilt_psi(double x){
	return -pow(hbar,2)/(2*m)*psiTrial(x,2) + V(x)*psiTrial(x,0);
}

void printhisto(){
	unsigned int tot=0;
	for(auto el : psihist) tot+=el;

	ofstream out;
	out.open("../data/psi.out");
	for(unsigned int i=0; i<psihist.size(); i++)
		out<<x0+ (i+0.5)*bin_size<<" "<<double(psihist[i])/(tot*bin_size)<<endl;
}
