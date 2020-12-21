#include <iostream>
#include <fstream>
#include <string>
#include "Monte_Carlo_NVT.h"
using namespace std;

int main(int argc, char* argv[]){
  if(argc<2){
	  cout<<"need to specify the phase (solid, liquid, gas) in argv for the file"<<endl;
	  return 0;
  }

  Input(); //Inizialization
  //int nconf = 1;
  ofstream p,u;
  string s(argv[1]);
  cout<<"printing "<<nstep<<" instant values of U(t)/N, P(t) in ../data/measures/ ..."<<endl;
  u.open("../data/measures/U_inst."+s);
  p.open("../data/measures/P_inst."+s);
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep){
      Move();
      Measure();
      u<<walker[iv]/(double)npart + vtail<<endl;
      p<<rho * temp + (walker[iw] + ptail * (double)npart) / vol<<endl;
      if(istep%(nstep/10)==0) cout<<istep*100./nstep<<"%"<<endl;

//      if(istep%10 == 0){
//       ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
//        nconf += 1;
//      }
    }
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
  }
  u.close();
  p.close();
  ConfFinal(); //Write final configuration

  return 0;
}



