//parameters, observables
const int m_props=5;
//int n_props;
//int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;
bool old; 	//aggiungo il booleano per poter fare il restart con configurazione vecchia
bool rescaling; //aggiungo booleano per operazione di rescaling

//functions
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
void datablocking(int);
double Force(int, int);
double Pbc(double);

