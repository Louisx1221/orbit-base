#include <iostream>
#include <math.h>

using namespace std;
#include "constants.h"
#include "math_plus.h"
#include "orbit.h"

int main() 
{
	double hp = 1000; 
	double ha = 1800;
	double jd = JulianDate(2021);

	Orbit orb;
	//Orbit orb = Orbit(oev, jd);
	orb.r_a = ha + RE;
	orb.r_p = hp + RE;
	orb.i = 90 * D2R;
	orb.Omega = 0;
	orb.omega = 60 * D2R;
	orb.f = 90 * D2R;
	orb.jd = jd;

	orb.RaRp2SmaEcc();
	orb.Coe2Oev();
	orb.Coe2Rv();
	orb.MeanMotion();
	orb.Coe2Eci();
	orb.Eci2Ecef();
	orb.Ecef2Geod();

	double oev_t[6] = {};
	double dt = 86400.0;
	J2long(orb.oev, dt, oev_t);

	system("pause");
	return 0;
}