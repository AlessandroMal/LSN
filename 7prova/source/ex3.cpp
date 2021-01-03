#include "Monte_Carlo_NVT.h"

int main()
{ 
  Input(); //Inizialization
//  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
  {
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      Measure();
      Accumulate(); //Update block averages
//      if(istep%(nstep/10)==0) cout<<istep*100./nstep<<"%"<<endl;
//      if(istep%10 == 0){
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
//        nconf += 1;
//      }
    }
    Averages(iblk);   //Print results for current block
  }
  ConfFinal(); //Write final configuration
  return 0;
}
