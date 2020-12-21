#include <iostream>
#include "Monte_Carlo_NVT.h"
using namespace std;

int main(){ //muovo e basta per equilibrare
  Input();  //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
    for(int istep=1; istep <= nstep; ++istep){
      Move();
    }
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
  }
  ConfFinal(); //Write final configuration

  return 0;
}



