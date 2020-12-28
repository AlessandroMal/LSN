#include "MolDyn_NVE.h"
using namespace std;

int main(){ 
  Input(); //Inizialization
//  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
      if(istep%10 == 0){
         Measure();     	//Properties measurement
         Accumulate();		//Update block averages
//         ConfXYZ(nconf);	//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
//         nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  return 0;
}
