#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include <vector>
using namespace std;

//parameters, observables
int n_props, iv,ik,it,ih,iw,nbins;
double bin_size;//,vtail,ptail;
vector<double> walker;

// averages
vector<double> blk_av,glob_av,glob_av2;
int blk_norm;

//configuration
vector<double> x,y,z,xold,yold,zold,vx,vy,vz;

// thermodynamical state
double temp,vol,rho,box,rcut;

// simulation
int nstep, nblk, iprint, seed;
double delta;
bool old; 	//aggiungo il booleano per poter fare il restart con configurazione vecchia

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Force(unsigned int, int);
double Pbc(double);
double Error(double,double,int);

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  int npart;

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
/*  //Tail corrections for potential energy and pressure
  vtail = (8.0*M_PI*rho)/(9.0*pow(rcut,9)) - (8.0*M_PI*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*M_PI*rho)/(9.0*pow(rcut,9)) - (16.0*M_PI*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Tail correction for the virial           = " << ptail << endl; 
*/
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> old;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  if(old==1) cout << "The program will return also the configuration before the final one" << endl<<endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ih = 2; //Total energy
  it = 3; //Temperature
  iw = 4; //Virial
  n_props = 5; //Number of observables

//measurement of g(r)
  nbins = 100;
  bin_size = (box/2.0)/(double)nbins;
  walker.resize(n_props+nbins);
  blk_av.resize(n_props+nbins);
  glob_av.resize(n_props+nbins);
  glob_av2.resize(n_props+nbins);

//Read initial configuration
  double xread, yread, zread;
  cout << "Read initial configuration from file config.0 " << endl;
  ReadConf.open("../data/configurations/config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> xread >> yread >> zread;
    x.push_back(xread * box);
    y.push_back(yread * box);
    z.push_back(zread * box);
  }
  ReadConf.close();
  xold.resize(x.size());
  yold.resize(x.size());
  zold.resize(x.size());
  vx.resize(x.size());
  vy.resize(x.size());
  vz.resize(x.size());

  double sumv2=0. ,fs;
//Read pre-initial configuration
  if(old==1){
	  ReadConf.open("../data/configurations/old.0");
	  if (ReadConf.is_open()){
		  cout << "Read pre-initial configuration from file old.0 " << endl << endl;
		  for (int i=0; i<npart; ++i){
			  ReadConf >> xread >> yread >> zread;
			  xold[i]= xread * box;
			  yold[i]= yread * box;
			  zold[i]= zread * box;
		  }
		  Move();
		  for (int i=0; i<npart; ++i){
			  vx[i]= (x[i]-xold[i])/(delta/2);
			  vy[i]= (y[i]-yold[i])/(delta/2);
			  vz[i]= (z[i]-zold[i])/(delta/2);
			  sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		  }

	  }else{
		  cout<<"FILE old.0 NOT FOUND: PREPARING WITH RANDOM VELOCITIES."<<endl;
		  old=0;
	  }
	  ReadConf.close();
  }
//Prepare initial random velocities
  if(old==0){
	  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
	  double sumv[3] = {0.0, 0.0, 0.0};
	  for (int i=0; i<npart; ++i){
		  vx[i]= rand()/double(RAND_MAX) - 0.5 ;
		  vy[i]= rand()/double(RAND_MAX) - 0.5 ;
		  vz[i]= rand()/double(RAND_MAX) - 0.5 ;
		  sumv[0] += vx[i];
		  sumv[1] += vy[i];
		  sumv[2] += vz[i];
	  }
	  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	  for (int i=0; i<npart; ++i){	//somma delle v^2
		  vx[i] = vx[i] - sumv[0];
		  vy[i] = vy[i] - sumv[1];
		  vz[i] = vz[i] - sumv[2];	  
		  sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	  }
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

  return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[x.size()], fy[x.size()], fz[x.size()];

  for(unsigned int i=0; i<x.size(); ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(unsigned int i=0; i<x.size(); ++i){ //Verlet integration scheme

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

double Force(unsigned int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (unsigned int i=0; i<x.size(); ++i){
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

void Measure()
{
  double v = 0.0, w = 0.0, t = 0.0;
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
  for (unsigned int i=0; i<x.size(); ++i) t += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  v= 4*v;
  w= 16*w;
  t= 0.5*t;

  walker[iv] = v;
  walker[iw] = w;
  walker[ik] = t;
  walker[it] = 2./3 * t;
  walker[ih] = t+v;
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
}


void Accumulate(void){ //Update block averages
   for(int i=0; i<n_props+nbins; ++i){
	   blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm++;
}


void Averages(int iblk) //Print results for current block
{
    ofstream out;
    double stima,err;
    const int wd=12;
    cout << "Block number " << iblk << endl;
    
    out.open("../data/measures/ave_epot.out",ios::app);
    stima= blk_av[iv]/blk_norm/x.size() /* + vtail*/; //Potential energy
    glob_av[iv] += stima;
    glob_av2[iv] += stima*stima;
    err= Error(glob_av[iv],glob_av2[iv],iblk);
    out << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err << endl;
    out.close();

    out.open("../data/measures/ave_pres.out",ios::app);
    stima= rho * temp + (blk_av[iw]/blk_norm/* + ptail * (double)npart*/) / vol; //Pressure
    glob_av[iw] += stima;
    glob_av2[iw] += stima*stima;
    err= Error(glob_av[iw],glob_av2[iw],iblk);
    out << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err << endl;
    out.close();

    out.open("../data/measures/ave_ekin.out",ios::app);
    stima= blk_av[ik]/blk_norm/x.size(); //Kinetic energy
    glob_av[ik] += stima;
    glob_av2[ik] += stima*stima;
    err= Error(glob_av[ik],glob_av2[ik],iblk);
    out << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err << endl;
    out.close();

    out.open("../data/measures/ave_etot.out",ios::app);
    stima= blk_av[ih]/blk_norm/x.size(); //Kinetic energy
    glob_av[ih] += stima;
    glob_av2[ih] += stima*stima;
    err= Error(glob_av[ih],glob_av2[ih],iblk);
    out << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[ih]/(double)iblk << setw(wd) << err << endl;
    out.close();

    out.open("../data/measures/ave_temp.out",ios::app);
    stima= blk_av[it]/blk_norm/x.size(); //Kinetic energy
    glob_av[it] += stima;
    glob_av2[it] += stima*stima;
    err= Error(glob_av[it],glob_av2[it],iblk);
    out << setw(wd) << iblk <<  setw(wd) << stima << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err << endl;
    out.close();

    out.open("../data/measures/ave_gcor.out",ios::app);
    for(int k=n_props; k<n_props+nbins; ++k){ //g(r)
	    stima= blk_av[k]/blk_norm;
	    glob_av[k] += stima;
	    glob_av2[k] += stima*stima;
	    if(iblk==nblk){
		    err= Error(glob_av[k],glob_av2[k],iblk);
		    out <<  setw(wd) << bin_size*(k-n_props+0.5) << setw(wd) << glob_av[k]/(double)iblk << setw(wd) << err << endl;
	    }
    }
    out.close();
    cout << "----------------------------" << endl << endl;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final "<< endl;
  WriteConf.open("../data/configurations/config.final");

  for (unsigned int i=0; i<x.size(); ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  cout<< "Print configuration before the final one in old.final" <<endl<<endl;
  ofstream oldconf;
  oldconf.open("../data/configurations/old.final");
  for (unsigned int i=0; i<x.size(); ++i)
	  oldconf <<xold[i]/box << "   " <<yold[i]/box << "   " <<zold[i]/box << endl;
  oldconf.close();

  return;
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

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
