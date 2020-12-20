#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_NVT.h"

using namespace std;
//muovo e basta per equilibrare
int main(){ 
  Input(); //Inizialization
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
    for(int istep=1; istep <= nstep; ++istep){
      Move();
    }
    cout << "Acceptance rate " << accepted/attempted << endl << endl;
  }
  ConfFinal(); //Write final configuration

  return 0;
}



