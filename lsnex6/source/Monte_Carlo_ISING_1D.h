#ifndef __Ising__
#define __Ising__

//Random numbers
#include "random.h"
int seed[4];
Random rnd;

//parameters, observables
const int m_props=100;
int n_props,iu,ic,im,ix,ig;
double nbins;
double walker[m_props];

// averages
double blk_av[m_props],blk_norm,accepted,attempted;
double glob_av[m_props],glob_av2[m_props];
double stima_u,stima_c,stima_m,stima_x;
double err_u,err_c,err_m,err_x,err_g;

//configuration
const int m_spin=50;
double s[m_spin];

// thermodynamical state
int nspin;
double beta,temp,J,h;

// simulation
int nstep, nblk;
bool metro, restart, instant, vsTemp;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(int);
void ConfFinal(void);
void Measure(void);
int Pbc(int);
double Gibbs(int);
double Boltzmann(int);
double Error(double,double,int);
double Min(double, double);

#endif
