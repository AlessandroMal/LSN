#ifndef __NVT__
#define __NVT__

//Random numbers
#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include "random.h"
using namespace std;
int seed[4];
Random rnd;

//parameters, observables
int n_props, iv,iw,nbins;
double bin_size,vtail,ptail;
vector<double> walker;

// averages
vector<double> blk_av,glob_av,glob_av2;
int blk_norm,accepted,attempted;

//configuration
vector<double> x,y,z;

// thermodynamical state
double beta,temp,vol,rho,box,rcut;

// simulation
int nstep, nblk;
double delta;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, unsigned int);
double Pbc(double);
double Error(double,double,int);

void Input(void)
{
  ifstream ReadInput,ReadConf;
  int npart;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

//Read seed for random numbers
  int p1, p2;
  ifstream Primes("../data/Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("../data/seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();
  
//Read input informations
  ReadInput.open("../data/configurations/input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  box = pow(vol,1.0/3.0);
  cout << "Volume of the simulation box = " << vol << endl;
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
    
  //Tail corrections for potential energy and pressure
  vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl; 

  ReadInput >> delta;

  ReadInput >> nblk;

  ReadInput >> nstep;

  cout << "The program perform Metropolis moves with uniform translations" << endl;
  cout << "Moves parameter = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iv = 0; //Potential energy
  iw = 1; //Virial
  n_props = 2; //Number of observables

//measurement of g(r)
  nbins = 100;
  bin_size = (box/2.0)/(double)nbins;
  walker.resize(n_props+nbins);
  blk_av.resize(n_props+nbins);
  glob_av.resize(n_props+nbins);
  glob_av2.resize(n_props+nbins);

//Read initial configuration
  double xread, yread, zread;
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("../data/configurations/config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> xread >> yread >> zread;
    x.push_back(xread * box);
    y.push_back(yread * box);
    z.push_back(zread * box);
  }
  ReadConf.close();
  
//Evaluate potential energy and virial of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
  cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
  cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
}


void Move(void)
{
  unsigned int o;
  double p, energy_old, energy_new;
  double xold, yold, zold, xnew, ynew, znew;


  for(unsigned int i=0; i<x.size(); ++i)
  {
  //Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
    o = (unsigned int)(rnd.Rannyu()*x.size());

  //Old
    xold = x[o];
    yold = y[o];
    zold = z[o];

    energy_old = Boltzmann(xold,yold,zold,o);

  //New
    xnew = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
    ynew = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
    znew = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

    energy_new = Boltzmann(xnew,ynew,znew,o);

  //Metropolis test
    p = exp(beta*(energy_old-energy_new));
    if(p >= rnd.Rannyu())  
    {
    //Update
       x[o] = xnew;
       y[o] = ynew;
       z[o] = znew;
    
       accepted = accepted + 1.0;
    }
    attempted = attempted + 1.0;
  }
}

double Boltzmann(double xx, double yy, double zz, unsigned int ip)
{
  double ene=0.0;
  double dx, dy, dz, dr;

  for (unsigned int i=0; i<x.size(); ++i)
  {
    if(i != ip)
    {
// distance ip-i in pbc
      dx = Pbc(xx - x[i]);
      dy = Pbc(yy - y[i]);
      dz = Pbc(zz - z[i]);

      dr = dx*dx + dy*dy + dz*dz;
      dr = sqrt(dr);

      if(dr < rcut)
      {
        ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
      }
    }
  }

  return 4.0*ene;
}

void Measure()
{
  double v = 0.0, w = 0.0;
  double vij, wij;
  double dx, dy, dz, dr;

//reset the hystogram of g(r)
  for (int k=n_props; k<n_props+nbins; ++k) walker[k]=0.0;

//cycle over pairs of particles
  for (unsigned int i=0; i<x.size()-1; ++i)
  {
    for (unsigned int j=i+1; j<x.size(); ++j)
    {

// distance i-j in pbc
     dx = Pbc(x[i] - x[j]);
     dy = Pbc(y[i] - y[j]);
     dz = Pbc(z[i] - z[j]);

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

//update of the histogram of g(r)
     for(int k=n_props; k<n_props+nbins; ++k)
	     if(dr<(k-n_props+1)*bin_size){
		     walker[k]+=2;
		     break;
	     }

     if(dr < rcut)
     {
       vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
       wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

// contribution to energy and virial
       v += vij;
       w += wij;
     }
    }
  }
  v= 4*v;
  w= 16*w;

  walker[iv] = v;
  walker[iw] = w;
  for(int k=n_props; k<n_props+nbins; ++k)
	  walker[k]/= rho*x.size()*4./3*M_PI*(  pow( (k-n_props+1)*bin_size, 3 ) - pow( (k-n_props)*bin_size, 3 ) );

}


void Reset(int iblk){ //Reset block averages
   if(iblk == 1){
	   fill(glob_av.begin(), glob_av.end(), 0.);
	   fill(glob_av2.begin(), glob_av2.end(), 0.);
   }
   fill(blk_av.begin(), blk_av.end(), 0.);
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}

void Accumulate(void){ //Update block averages
   for(int i=0; i<n_props+nbins; ++i){
	   blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm++;
}

void Averages(int iblk) //Print results for current block
{
    ofstream out,outg;
    double stima,err;
    const int wd=12;
    cout << "Block number " << iblk << endl;
    cout << "Acceptance rate " << double(accepted)/attempted << endl << endl;
    
    out.open("../data/measures/output_epot.0",ios::app);
    stima= blk_av[iv]/blk_norm/x.size() + vtail; //Potential energy
    glob_av[iv] += stima;
    glob_av2[iv] += stima*stima;
    err= Error(glob_av[iv],glob_av2[iv],iblk);
    out << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err << endl;
    out.close();

    out.open("../data/measures/output_pres.0",ios::app);
    stima= rho * temp + (blk_av[iw]/blk_norm + ptail * x.size()) / vol; //Pressure
    glob_av[iw] += stima;
    glob_av2[iw] += stima*stima;
    err= Error(glob_av[iw],glob_av2[iw],iblk);
    out << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err << endl;
    out.close();

    out.open("../data/measures/output_gofr.0",ios::app);
    if(iblk==nblk)
	    outg.open("../data/measures/output_gave.0",ios::app);
    for(int k=n_props; k<n_props+nbins; ++k){ //g(r)
	    stima= blk_av[k]/blk_norm;
	    glob_av[k] += stima;
	    glob_av2[k] += stima*stima;
	    err= Error(glob_av[k],glob_av2[k],iblk);
	    out <<  setw(wd) << bin_size*(k-n_props+0.5) << setw(wd) << glob_av[k]/(double)iblk << setw(wd) << err << endl;
	    if(iblk==nblk)
		    outg <<  setw(wd) << bin_size*(k-n_props+0.5) << setw(wd) << glob_av[k]/(double)iblk << setw(wd) << err << endl;
    }
    out.close();
    if(iblk==nblk)
	    outg.close();
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("../data/configurations/config.final");
  for (unsigned int i=0; i<x.size(); ++i)
  {
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("../data/measures/frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << x.size() << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (unsigned int i=0; i<x.size(); ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk){
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}

#endif
