#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <vector>	//have to define a vector to use functions in stattools.h
#include "MolDyn_NVE.h"
#include "stattools.h"

//domande: devo aggiungere tail correction?
//vedi il problema della pressione
//le ripartenze vanno implementate? io le ho lanciate manualmente nel numero che mi serviva

using namespace std;

int main(){ 
  Input();             //Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  datablocking(100);	//funzione che prende in enttrata il numero di blocchi
  			//e effettua datablocking sui dati dei file .dat
  ConfFinal();		//Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("../data/configurations/input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> old;
  ReadInput >> rescaling;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  if(old==1)
	  cout << "The program will return also the configuration before the final one" << endl<<endl;
  ReadInput.close();

//Prepare array for measurements
//  iv = 0; //Potential energy
//  ik = 1; //Kinetic energy
//  ie = 2; //Total energy
//  it = 3; //Temperature
//  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl;
  ReadConf.open("../data/configurations/config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Read pre-initial configuration
  ReadConf.open("../data/configurations/old.0");
  if (ReadConf.is_open()){
	  cout << "Read pre-initial configuration from file old.0 " << endl << endl;
	  for (int i=0; i<npart; ++i){
		  ReadConf >> xold[i] >> yold[i] >> zold[i];
		  xold[i] = xold[i] * box;
		  yold[i] = yold[i] * box;
		  zold[i] = zold[i] * box;
	  }
	  ReadConf.close();
  }else{
	  ReadConf.close();
//Prepare initial velocities
	  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	  double sumv[3] = {0.0, 0.0, 0.0};
	  for (int i=0; i<npart; ++i){
		  vx[i] = rand()/double(RAND_MAX) - 0.5;
		  vy[i] = rand()/double(RAND_MAX) - 0.5;
		  vz[i] = rand()/double(RAND_MAX) - 0.5;
		  sumv[0] += vx[i];
		  sumv[1] += vy[i];
		  sumv[2] += vz[i];
	  }
	  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	  double sumv2 = 0.0, fs;
	  for (int i=0; i<npart; ++i){	//somma delle v^2
		  vx[i] = vx[i] - sumv[0];
		  vy[i] = vy[i] - sumv[1];
		  vz[i] = vz[i] - sumv[2];
		  
		  sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	  }
	  sumv2 /= (double)npart;
	  
	  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	  for (int i=0; i<npart; ++i){
		  vx[i] *= fs;
		  vy[i] *= fs;
		  vz[i] *= fs;
		  
		  xold[i] = Pbc(x[i] - vx[i] * delta);
		  yold[i] = Pbc(y[i] - vy[i] * delta);
		  zold[i] = Pbc(z[i] - vz[i] * delta);
	  }
  }

  if(rescaling==1){
	  double sumv2 = 0.0, fs;
	  Move();
	  for (int i=0; i<npart; ++i){
		  vx[i] *= (x[i]-xold[i])/(delta/2);
		  vy[i] *= (y[i]-yold[i])/(delta/2);
		  vz[i] *= (z[i]-zold[i])/(delta/2);
	
		  sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	  }
	  sumv2 /= (double)npart;
	  
	  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	  for (int i=0; i<npart; ++i){
		  vx[i] *= fs;
		  vy[i] *= fs;
		  vz[i] *= fs;
		  
		  xold[i] = Pbc(x[i] - vx[i] * delta);
		  yold[i] = Pbc(y[i] - vy[i] * delta);
		  zold[i] = Pbc(z[i] - vz[i] * delta);
	  }
  }
  
  return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
 // int bin;
  double v, t, vij, w, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Pres;

  Epot.open("../data/measures/output_epot.dat",ios::app);
  Ekin.open("../data/measures/output_ekin.dat",ios::app);
  Temp.open("../data/measures/output_temp.dat",ios::app);
  Etot.open("../data/measures/output_etot.dat",ios::app);
  Pres.open("../data/measures/output_pres.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       //potential energy
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       v += vij;
       //virial for pressure
       wij = 48*(1./pow(dr,12)+0.5/pow(dr,6));
       w += wij;
     }
    }
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart;			//Potential energy per particle
    stima_kin = t/(double)npart;			//Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart;		//Temperature
    stima_etot = (t+v)/(double)npart;			//Total energy per particle
    stima_pres= rho*temp + w/3/vol;			//Pressure

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;
    Pres << stima_pres << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();
    Pres.close();

    return;
}

void datablocking(int nblocks){
	int L=nstep/nblocks; //occhio che nstep sia divisibile per nblocks
	double sum,d;
	vector<double> avinblock(nblocks);
	string obs[m_props]={"epot","ekin","etot","temp","pres"};

	for(auto el : obs){
		ifstream data;
		data.open("../data/measures/output_"+el+".dat");

		for(int i=0;i<nblocks;i++){
			sum=0;
			for(int j=0;j<L;j++){
				data>>d;
				sum+=d;
			}
			sum/=L;
			avinblock[i]=sum;
		}
		data.close();
		stampadati(L,avinblock,"../data/measures/ave_"+el+".out");
	}
}

void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final "<< endl;
  WriteConf.open("../data/configurations/config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  if(old==1){
	  cout<< "Print configuration before the final one in old.final" <<endl<<endl;
	  ofstream oldconf;
	  oldconf.open("../data/configurations/old.final");
	  for (int i=0; i<npart; ++i)
		  oldconf <<xold[i]/box << "   " <<yold[i]/box << "   " <<zold[i]/box << endl;
	  oldconf.close();
  }

  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("../data/measures/frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
