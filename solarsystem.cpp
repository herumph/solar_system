#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "boost/numeric/odeint.hpp"

using namespace boost::numeric::odeint;
using namespace std;

//Global force array. The only way to get it into the ODE function.
//F[bodies][3]
double F[11][3];
int bodies = 11;

typedef std::vector<double> state_type;

//ODEs to solve.
void ODE(const state_type &x, state_type &dxdt, const double t)
{
	double G = 6.67e-8;//4*M_PI/pow(365.0,2.0);
        //System of 3*bodies equations
        for(int i = 0; i < bodies; i++){
                //x
                dxdt[6*i] = x[6*i+3];
                dxdt[6*i+3] = G*F[i][0];
                //y
                dxdt[6*i+1] = x[6*i+4];
                dxdt[6*i+4] = G*F[i][1];
                //z
                dxdt[6*i+2] = x[6*i+5];
                dxdt[6*i+5] = G*F[i][2];
        }
}

//Energy calculation
double calcenergy(const state_type &x, const state_type &m)
{
	double energy = 0, kinetic = 0, potential = 0, mag;
	double G = 6.67e-8;//4*M_PI/pow(365.0,2.0);

	for(int i = 0; i < bodies; i++){
		kinetic = kinetic + m[i]*(pow(x[6*i+3],2.0) + pow(x[6*i+4],2.0) + pow(x[6*i+5],2.0));
		for(int j = 0; j < bodies; j++){
			if(j != i){
				mag = sqrt(pow(x[6*i]-x[6*j],2.0) + pow(x[6*i+1]-x[6*j+1],2.0) + pow(x[6*i+2]-x[6*j+2],2.0));
				potential = potential - G*m[i]*m[j]/mag;
			}
		}
	}

	energy = 0.5*(kinetic + potential);
	return energy;
}

//Angular momentum calculation
double calcmom(const state_type &x, const state_type &m)
{
	double angx = 0, angy = 0, angz = 0, ang;

	for(int i = 0; i < bodies; i++){
		angx = angx + m[i]*(x[6*i+1]*x[6*i+5]-x[6*i+2]*x[6*i+4]);
		angy = angy + m[i]*(x[6*i+2]*x[6*i+3]-x[6*i]*x[6*i+5]);
		angz = angz + m[i]*(x[6*i]*x[6*i+4]-x[6*i+1]*x[6*i+3]);
	}

	ang = sqrt(pow(angx,2.0) + pow(angy,2.0) + pow(angz,2.0));
	return ang;
}

int main()
{
	int steps = 0;
	//Errors
	double energyerr = 1e-10, interr = 1e-12, minstep = 1e-5;
	//Force softening
	double epsilon = 0.0;
	//Time step
	double h = 1e3, t = 0.0, tmax = 86400.0*365.0, tmin = 0;
	state_type x(6*bodies), xold(6*bodies), m(bodies), tempx(6*bodies);
	double energy, mag, h1=h, E1, ang;

	//Initial conditions
	m[0] = 1.99e33, m[1] = 3.303e26, m[2] = 4.870e27, m[3] = 5.976e27;
	m[4] = 6.418e26, m[5] = 1.899e30, m[6] = 5.686e29, m[7] = 8.66e28;
	m[8] = 1.030e29, m[9] = 1.0e25, m[10] = 7.348e25;
	
	//Sun initial conditions
	x[0] = -1.921e10, x[1] = -3.673e10, x[2] = -6.294e8, x[3] = 1064.112, x[4] = -228.215; x[5] = -22.618;
	//x[0] = 1.476e10, x[1] = -3.408e10, x[2] = 1.392e9, x[3] = 1049.592, x[4] = 403.244, x[5] = -24.314;
	//Mercury ICs
	x[6] = -2.542e12, x[7] = -6.518e12, x[8] = -2.987e11, x[9] = 3.561e6, x[10] = -1.527e6; x[11] = -4.514e5;
	//x[6] = 1.806e12, x[7] = -6.531e12, x[8] = -6.965e11, x[9] = 3.721e6, x[10] = 1.543e6, x[11] = -2.152e5;
	//Venus ICs
	x[12] = -6.906e12, x[13] = -8.392e12, x[14] = 2.823e11, x[15] = 2.679e6, x[16] = -2.244e6, x[17] = -1.853e5;
	//x[12] = -7.315e11, x[13] = 1.070e13, x[14] = 1.888e11, x[15] = -3.505e6, x[16] = -2.613e5, x[17] = 1.987e5;
	//Earth ICs
	x[18] = -2.712e12, x[19] = 1.442e13, x[20] = -1.025e9, x[21] = -2.975e6, x[22] = -5.558e5, x[23] = 1.166;
	//x[18] = -2.612e12, x[19] = 1.444e13, x[20] = -1.871e9, x[21] = -2.980e6, x[22] = -5.43e5, x[23] = -64.782;
	//Mars ICs
	x[24] = 1.614e13, x[25] = -1.300e13, x[26] = -6.689e11, x[27] = 1.609e6, x[28] = 2.097e6, x[29] = 4420.935;
	//x[24] = -2.261e13, x[25] = 1.039e13, x[26] = 7.724e11, x[27] = -9.220e5, x[28] = -1.993e6, x[29] = -1.914e4;
	//Jupiter ICs
	x[30] = 2.130e13, x[31] = 7.265e13, x[32] = -7.796e11, x[33] = -1.270e6, x[34] = 4.30e5, x[35] = 2.662e4;
	//x[30] = -1.989e13, x[31] = 7.504e13, x[32] = 1.323e11, x[33] = -1.279e6, x[34] = -2.727e5, x[35] = 2.976e4;
	//Saturn ICs
	x[36] = -1.209e14, x[37] = -8.254e13, x[38] = 6.247e12, x[39] = 4.923e5, x[40] = -8.001e5, x[41] = -5660.41;
	//x[36] = -1.030e14, x[37] = -1.059e14, x[38] = 5.940e12, x[39] = 6.399e5, x[40] = -6.761e5, x[41] = -1.366e4;
	//Uranus ICs
	x[42] = 2.975e14, x[43] = 3.846e13, x[44] = -3.712e12, x[45] = -9.228e4, x[46] = 6.434e5, x[47] = 3602.15;
	//x[42] = 2.939e14, x[43] = 5.865e13, x[44] = -3.590e12, x[45] = -1.382e5, x[46] = 6.361e5, x[47] = 4155.692;
	//Neptune ICs
	x[48] = 3.973e14, x[49] = -2.083e14, x[50] = -4.867e12, x[51] = 2.487e5, x[52] = 4.845e5, x[53] = -15728.10;
	//x[48] = 4.049e14, x[49] = -1.929e14, x[50] = -5.359e12, x[51] = 2.301e5, x[52] = 4.939e5, x[53] = -15443.827;
	//Pluto ICs
	x[54] = 7.641e13, x[55] = -4.772e14, x[56] = 2.896e13, x[57] = 5.453e5, x[58] = -2.387e4, x[59] = -1.547e5;
	//x[54] = 9.361e13, x[55] = -4.776e14, x[56] = 2.403e13, x[57] = 5.430e5, x[58] = -4791.767, x[59] = -1.546e5;
	//Moon ICs
	x[60] = -2.743e12, x[61] = 1.445e13, x[62] = -4.551e9, x[63] = -3.035e6, x[64] = -6.352e5, x[65] = 306.72;
	//x[60] = -2.610e12, x[61] = 1.440e13, x[62] = 8.099e8, x[63] = 2.870e6, x[64] = -5.335e5, x[65] = 4868.09;

	//Initial energy
	double Ei = calcenergy(x,m);
	double E0 = Ei;
	cout << setprecision(15) << "Initial energy: " << E0 << "\n";

	//Initial angular momentum
	double angi = calcmom(x,m);
	cout << setprecision(15) << "Initial angular momentum: " << angi << "\n";

	typedef runge_kutta_cash_karp54<state_type> error_stepper_type;

	ofstream myfile;
	myfile.open("results.txt");

	while(t <= tmax){
		if(steps % 100000 == 0){
			cout << "Time reached: " << t << "\n";
		}

		energy = calcenergy(x,m);
		ang = calcmom(x,m);
		E1 = energy/E0 - 1.0;

		if(abs(E1) > energyerr && h1 > minstep){
			//Error is too large and h1 has not hit minimum -> revert to last good step
			h1 = max(minstep, 0.5*h1);
			for(int i = 0; i < 6*bodies; i++){
				x[i] = xold[i];
			}
		}

		else{
			//Good step
			for(int i = 0; i < 6*bodies; i++){
				xold[i] = x[i];
				myfile << setprecision(15) << x[i] << " ";
			}
			myfile << t << " " << energy/Ei - 1.0 << " " << ang/angi - 1.0 << "\n";

			//Finding new forces
			for(int i = 0; i < bodies; i++){
				F[i][0] = 0;
				F[i][1] = 0;
				F[i][2] = 0;
				for(int j = 0; j < bodies; j++){
					if(j != i){
						mag = sqrt(pow(x[6*i]-x[6*j],2.0) + pow(x[6*i+1]-x[6*j+1],2.0) + pow(x[6*i+2]-x[6*j+2],2.0));
						F[i][0] = F[i][0] + m[j]*(x[6*j]-x[6*i])/pow(mag+epsilon,3.0);
						F[i][1] = F[i][1] + m[j]*(x[6*j+1]-x[6*i+1])/pow(mag+epsilon,3.0);
						F[i][2] = F[i][2] + m[j]*(x[6*j+2]-x[6*i+2])/pow(mag+epsilon,3.0);
					}
				}
			}

			t = t + h1;
			h1 = min(h,1.20*h1);
			E0 = energy;
		}

		//Integration
                integrate_adaptive(make_controlled<error_stepper_type>(interr, interr), ODE, x, t, t+h1, h1);
		steps += 1;
	}

        cout << "Final energy: " << energy << "\n";
	cout << "Final angular momentum: " << ang << "\n";
	cout << "Sun final position: " << x[0] << " " << x[1] << " " << x[2] << "\n";
	cout << "Mercury final position: " << x[6] << " " << x[7] << " " << x[8] << "\n";
	cout << "Venus final position: " << x[12] << " " << x[13] << " " << x[14] << "\n";
	cout << "Earth final position: " << x[18] << " " << x[19] << " " << x[20] << "\n";
	cout << "Moon final position: " << x[60] << " " << x[61] << " " << x[62] << "\n";
	cout << "Mars final position: " << x[24] << " " << x[25] << " " << x[26] << "\n";
	cout << "Jupiter final position: " << x[30] << " " << x[31] << " " << x[33] << "\n";
	cout << "Saturn final position: " << x[36] << " " << x[37] << " " << x[38] << "\n";
	cout << "Uranus final position: " << x[42] << " " << x[43] << " " << x[44] << "\n";
	cout << "Neptune final position: " << x[48] << " " << x[49] << " " << x[50] << "\n";
	cout << "Pluto final position: " << x[54] << " " << x[55] << " " << x[56] << "\n";
}	
