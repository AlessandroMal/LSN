#include <iostream>
#include <fstream>
#include "RW.h"
#include "random.h"
#include "stattools.h"

//domande galli: qui non uso blocchi, dovrei riprogrammare e fare prima tutti i primi passi, quindi seguire gli L randomwalk insieme

int main (int argc, char *argv[]){
	int nstep = 100;	//numero di passi nel RW
	int M = 1E4;		//numero di RW generati
	double a=1;		//passo reticolare

	//Random walk cubic lattice
	vector<double> RW3Ddisc(3);	//posizione RW
	vector<double> avr2(nstep), mom2r2(nstep);
	//Random walk continuum
	vector<double> RW3Dcont(3);	//posizione RW
	vector<double> avr2cont(nstep), mom2r2cont(nstep);

	for(int i=0; i<M; i++){
		fill(RW3Ddisc.begin(), RW3Ddisc.end(), 0.); //memset slightly faster ma vabe
		fill(RW3Dcont.begin(), RW3Dcont.end(), 0.); //memset slightly faster ma vabe
		
		for(int j=0; j<nstep; j++){	//faccio partire i cammini
			//Random walk cubic lattice
			RW3Ddisc=rndstepD(RW3Ddisc,a);
			avr2[j] += norm2(RW3Ddisc);
			mom2r2[j] += pow(norm2(RW3Ddisc),2);

			//Random walk continuum
			RW3Dcont=rndstepC(RW3Dcont,a);
			avr2cont[j] += norm2(RW3Dcont);
			mom2r2cont[j] += pow( norm2(RW3Dcont), 2 );
		}
	}

	for(int i=0; i<nstep; i++){
		//considero tutti i cammini come esperimenti differenti
		//Random walk cubic lattice
		avr2[i] = avr2[i]/M;			//media degli r^2
		mom2r2[i] = mom2r2[i]/M;		//media dei quadrati degli r^2
		//Random walk continuum
		avr2cont[i] = avr2cont[i]/M;
		mom2r2cont[i] = mom2r2cont[i]/M;	//stessa cosa nel continuo
	}

	ofstream Cube("../data/RWlattice.dat");
	ofstream Cont("../data/RW3Dcont.dat");
	for(int i=0; i<nstep; i++){
		Cube<<sqrt(avr2[i])<<" "<<rmsclt(avr2[i], mom2r2[i], M)<<endl;
		Cont<<sqrt(avr2cont[i])<<" "<<rmsclt(avr2cont[i], mom2r2cont[i], M)<<endl;
	}

	Cube.close();
	Cont.close();
	rnd.SaveSeed();
	return 0;
}
