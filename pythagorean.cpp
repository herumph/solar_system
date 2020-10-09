/*************************************************************************************

Pythagorean three body problem. Can be generalized to n bodies with minor changes.

**************************************************************************************/

#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <chrono>
#include "boost/numeric/odeint.hpp"

using namespace boost::numeric::odeint;
using namespace std;

//Global force array. The only way to get it into the ODE function.
double F[3][3];
int bodies = 3;

typedef std::vector<double> state_type;

//ODEs to solve.
void ODE(const state_type &x, state_type &dxdt, const double t)
{
        //System of 3*bodies equations
        for(int i = 0; i < bodies; i++){
                //x
                dxdt[6*i] = x[6*i+3];
                dxdt[6*i+3] = F[i][0];
                //y
                dxdt[6*i+1] = x[6*i+4];
                dxdt[6*i+4] = F[i][1];
                //z
                dxdt[6*i+2] = x[6*i+5];
                dxdt[6*i+5] = F[i][2];
        }
}

//Energy calculation
double calcenergy(const state_type &x, const state_type &m)
{
	double energy = 0, kinetic = 0, potential = 0, mag;

	for(int i = 0; i < bodies; i++){
		kinetic = kinetic + m[i]*(pow(x[6*i+3],2.0) + pow(x[6*i+4],2.0) + pow(x[6*i+5],2.0));
		for(int j = 0; j < bodies; j++){
			if(j != i){
				mag = sqrt(pow(x[6*i]-x[6*j],2.0) + pow(x[6*i+1]-x[6*j+1],2.0) + pow(x[6*i+2]-x[6*j+2],2.0));
				potential = potential - m[i]*m[j]/mag;
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
	std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
	int steps = 0;
	//Errors
	double energyerr = 1e-8, interr = 1e-10, minstep = 1e-10;
	//Force softening
	double epsilon = 0.0;
	//Time step
	double h = 1e-3, t = 0.0, tmax = 10.0, tmin = 0;
	state_type x(6*bodies), xold(6*bodies), m(bodies), tempx(6*bodies);
	double energy, mag, h1=h, E1, ang;

	//Initial conditions
	m[0] = 3, m[1] = 4, m[2] = 5;
	//First body x and v.
	x[0] = 1, x[1] = 3, x[2] = 0, x[3] = 0, x[4] = 0; x[5] = 0;
	//Second body x and v.
	x[6] = -2, x[7] = -1, x[8] = 0, x[9] = 0, x[10] = 0; x[11] = 0;
	//Third body x and v.
	x[12] = 1, x[13] = -1, x[14] = 0, x[15] = 0, x[16] = 0, x[17] = 0;

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
				if(i % 6 == 0){
					myfile << setprecision(15) << t << " " << h1 << " " << energy  << " " << ang << " " << x[i] << " " << x[i+1] << " " << x[i+2] << "\n";
				}
			}

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

	std:chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
        cout << "Program ran for: " << std::chrono::duration_cast<std::chrono::minutes>(end-start).count() << " minutes, ";
        cout << std::chrono::duration_cast<std::chrono::seconds>(end-start).count()%60 << " seconds.\n";
        cout << "Max time reached: " << t << "\n";
        cout << "Final energy: " << energy << "\n";
	cout << "Final angular momentum: " << ang << "\n";
}	
