#include <iostream>
#include <cmath>
#include <vector>
#include <armadillo>
#include "random.h"
#include "stattools.h"
using namespace arma;

const double mu=2;
const double sigma=1;
const double hbar=1;
const double m=1;
Random rnd=generaterng();

vector<double> blockaverage(vector<double>, int);
double p(double);
double psiTrial(double);
double V(double);
vec cdf(double, int);
bool Asimm(double, double, double (*p)(double));

int main (/*int argc, char *argv[]*/){
   const int M=1E5;	//numero tot di lanci per caloclare r medio
   const int N=100;	//numero di blocchi
   const double dx=sigma/1000;

   vector<double> meas(M);		//M punti generati dal Metropolis
   vector<double> av_meas(N);		//medie sui valori dei blocchi
   int accrate=0;

   double x=mu;
   if(abs(x)>range/2){ //controllo per x di partenza
	   std::cout<<"start resized"<<endl;
	   if(x<0) x=-range/2;
	   else x=range/2;
   }
   int xbin=(range/2 + x)/bin_size + 0.5;
//   vec psiTrial=cdf(range,nbins);
//   vec eigv=cdf(range,nbins);
//   double current=HpsiTrial(xbin)/psiTrial(-range/2 + xbin*bin_size);
//   double current=eigv(xbin);

   double current= hamil_psi(x,bin_size)/psiTrial(x);
   double trial;
   
   for(int i=0; i<M; i++){
	   x =rnd.Rannyu(x-range/20, x+range/20);
	   if(xbin>nbins-1) xbin=nbins-1;
	   if(xbin<0) xbin=0;
//	   trial=HpsiTrial(xbin)/psiTrial(-range/2 + xbin*bin_size);
	   trial=eigv(xbin);
	   if( Asimm(current,trial,p) ){	//Acceptance
		   current=trial;
		   accrate+=1;
	   }
	   meas.push_back(current);
   }


/*   for(int i=0; i<M; i++){
//	   xbin =rnd.Rannyu(xbin-double(nbins)/100, xbin+double(nbins)/100);
	   if(xbin>nbins-1) xbin=nbins-1;
	   if(xbin<0) xbin=0;
//	   trial=HpsiTrial(xbin)/psiTrial(-range/2 + xbin*bin_size);
	   trial=eigv(xbin);
	   if( Asimm(current,trial,p) ){	//Acceptance
		   current=trial;
		   accrate+=1;
	   }
	   meas.push_back(current);
   }*/

   std::cout<<"acceptance rate: "<<double(accrate)/M<<endl;
   av_meas=blockaverage(meas,N);
   stampadati(M/N,av_meas, "../data/H_av.dat");	//media progressiva

   rnd.SaveSeed();
   return 0;
}

vector<double> blockaverage(vector<double> datas, int N){
//argomenti: dati da mediare, numero di blocchi
	int L=datas.size()/N;
	double sum;
	vector<double> av(N);

	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++)
			sum += datas[i*L+j];
		av[i]=sum/L;
	}
	return av;
}

double p(double x){//if real
	return psiTrial(x)*psiTrial(x);
}

double psiTrial(double x){
	return exp( -0.5*pow( (x-mu)/sigma, 2) ) + exp( -0.5*pow( (x+mu)/sigma, 2) );
}

double V(double x){
	return x*x*(x*x - 2.5);
}

double hamilt_psi(double x, double bin_size){
	return hbar*habr/2/m * (psiTrial(x+bin_size)+psiTrial(x-bin_size)-2*psiTrial(x))/(bin_size*bin_size) + V(x)*psiTrial(x);

vec cdf(double a, int N){
	double bin_size=a/N;
//	vec k=(pow(hbar/bin_size, 2)/m)*ones(N);
//	vec v=zeros(N);

	vec k(N,fill::ones);
	k= k*(pow(hbar/bin_size, 2)/m);
	vec v(N);
	for(int i=0; i<N; i++)
		v(i)=V(-a/2 + i*bin_size);

	mat H(N,N, fill::zeros);
	H.diag()  = k+v;
	k.resize(N-1);
	H.diag(1) = -0.5*k;
	H.diag(-1)= -0.5*k;

//	vec e[2];
	vec eigvalHpsi(N);
	mat Hpsi(N,N);
	eig_sym(eigvalHpsi,Hpsi,H);
//	eig_sym(e[0],e[1],H);
//	return Hpsi;

	for(int i=0; i<N; i++)
		cout<<eigvalHpsi(i)<<endl;
	return eigvalHpsi;
}

bool Asimm(double x, double xnew, double (*p)(double)){
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
