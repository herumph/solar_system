#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ctime>
#include "boost/numeric/odeint.hpp"

using namespace std;
using namespace boost::numeric::odeint;

int bodies = 3;
double omega, rcom[3];

typedef std::vector<double> state_type;

//Energy calculation
double calcenergy(const state_type &x, const state_type &m)
{
	//COM for each direction
	rcom[0] = (m[0]*x[0]+m[1]*x[6]+m[2]*x[12])/(m[0]+m[1]+m[2]);
	rcom[1] = (m[0]*x[1]+m[1]*x[7]+m[2]*x[13])/(m[0]+m[1]+m[2]);
	rcom[2] = (m[0]*x[2]+m[1]*x[8]+m[2]*x[14])/(m[0]+m[1]+m[2]);

	double r1,r2,G = 6.67384*1e-8;
	double jacobi = 0;

	r1 = sqrt(pow(x[0]-x[12],2.0)+pow(x[1]-x[13],2.0));
	r2 = sqrt(pow(x[6]-x[12],2.0)+pow(x[7]-x[13],2.0));
	jacobi = 1.0/2.0*(pow(x[12],2.0)+pow(x[13],2.0)) + G*(m[0]/r1 + m[1]/r2) - (pow(x[15],2.0)+pow(x[16],2.0));

	return jacobi;	
}

int main()
{
	state_type x(6*bodies), m(bodies);
	double G = 6.67384*1e-8, Msun = 1.9891*1e33, E1;

	//Initial masses
	m[0] = 1.20*Msun, m[1] = 0.40*Msun, m[2] = 1e-10;
	//Orbit stuff
	double Porb = 0.177*86400.0;
	omega = 2.0*M_PI/Porb;
 	double semiaxis = pow(G*(m[0]+m[1])/pow(omega,2.0),1.0/3.0);

	//first body
	x[0] = 0.0, x[1] = 0.0, x[2] = 0.0, x[3] = 0.0, x[4] = 0.0, x[5] = 0.0;
	//second
	x[6] = semiaxis, x[7] = 0.0, x[8] = 0.0, x[9] = 0.0, x[10] = 0.0, x[11] = 0.0;
	//third
	x[12] = x[0] + 0.1*(x[0]-x[6]), x[13] = 0.0, x[14] = 0.0, x[15] = 0.0, x[16] = 1.0, x[17] = 0.0;

	//Initial energy
	double Ei = calcenergy(x,m);
	double E0 = Ei;

	x[0] = 0.0 - rcom[0];
	x[6] = x[0] + semiaxis;

	E1 = calcenergy(x,m);

	int nx = 100, ny = 100;
	double xmin = -2*1e11, xmax = 2*1e11, dx = (xmax-xmin)/nx;
	double ymin = -2*1e11, ymax = 2*1e11, dy = (ymax-ymin)/ny;
	double y3,x3,psi,r1,r2;

	ofstream myfile;
	myfile.open("results.txt");

	for(int j = 0; j <= ny; j++){
		y3 = ymin + (j-0.5)*dy;
		for(int i = 0; i <= nx; i++){
			x3 = xmin + (i-0.5)*dx;
			r1 = sqrt(pow(x[0]-x3,2.0) + pow(x[1]-y3,2.0));
			r2 = sqrt(pow(x[6]-x3,2.0) + pow(x[7]-y3,2.0));

			//x3,y3 are correct, psi is off, not rcom
			psi = -G*(m[0]/r1 + m[1]/r2) - 0.5*(pow(x3-rcom[0],2.0) + pow(y3-rcom[1],2.0))*pow(omega,2.0);
			myfile << setprecision(15) << x3 << " " << y3 << " " << psi << "\n";
		}
	}

	return 0;
}
