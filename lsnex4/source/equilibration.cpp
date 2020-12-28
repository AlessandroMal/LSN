#include <iostream>
#include "MolDyn_NVE.h"
using namespace std;

int main(){
  Input(); //Inizialization
  ofstream out;
  out.open("../data/measures/equil_temp.out");
  double t;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0){
	 t=0.;
	 for (unsigned int i=0; i<x.size(); ++i) t += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	 out<<0.5*t *2/3/x.size()<<endl;
      }
    }
  }
  out.close();
  ConfFinal(); //Write final configuration
  return 0;
}

