#ifndef __Metropolis__
#define __Metropolis__

#include <vector>
#include "random.h"
using namespace std;
extern Random rnd;
extern const double a0;

vector<double> metropolis(vector<double>, int, double (*p)(vector<double>), vector<double> (*T)(vector<double>, double), bool (*A)(vector<double>, vector<double>, double (*parg)(vector<double>)) );
bool Asimm(vector<double>, vector<double>, double (*p)(vector<double>));
vector<double> Tunif(vector<double>, double);
vector<double> Tgauss(vector<double>, double);
double pGS(vector<double>);
double p210(vector<double>);
double norm2(vector<double>);

#endif // __Metropolis__
