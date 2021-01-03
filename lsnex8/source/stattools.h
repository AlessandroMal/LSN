#ifndef __Stattools__
#define __Stattools__

#include <vector>
#include <string>
using namespace std;

vector<double> avprogr( vector<double> );
vector<double> av2progr( vector<double> );
vector<double> error( vector<double> );
void stampadati(int, vector<double>, string);
double chi_quadro( vector<int>, double );
double rmsclt(double, double, int );

#endif // __Stattools__
