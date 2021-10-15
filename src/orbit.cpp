//@file       : orbit.cpp
//@autor      : github.com/louisx1221
//@date       : 2021/09/25

#include <iostream>
#include "orbit.h"
#include "math_plus.h"
#include "constants.h"

//#include <cmath>
//#include <ctime>
//#include <algorithm>

//using namespace std;

// default constructor
Orbit::Orbit()
{
	for (int i = 0; i < 6; i++)
		oev[i] = {};
	for (int i = 0; i < 3; i++)
	{
		r_eci[i] = {};
		v_eci[i] = {};
		r_ecef[i] = {};
		v_ecef[i] = {};
		lla[i] = {};
	}
	a = {};
	e = {};
	i = {};
	Omega = {};
	omega = {};
	f = {};
	E = {};
	M = {};
	T = {};
	n = {};
	r = {};
	v = {};
	r_a = {};
	r_p = {};
	lon = {};
	lat = {};
	alt = {};
	jd = {};
}

// constructor
Orbit::Orbit(double oev_[], double jd_)
{
	for (int i = 0; i < 6; i++)
		oev[i] = oev_[i];
	jd = jd_;
	Oev2Coe();
}

// class destructor
Orbit::~Orbit()
{

}

void Orbit::Oev2Coe()
{
	// classical orbital elements
	a = oev[0];
	e = oev[1];
	i = oev[2];
	Omega = oev[3];
	omega = oev[4];
	f = oev[5];

	// radius and velocity
	Coe2Rv();

	// radius of apogee and perigee
	SmaEcc2RaRp();

	// mean motion
	MeanMotion();
}

void Orbit::Coe2Oev()
{
	oev[0] = a;
	oev[1] = e;
	oev[2] = i;
	oev[3] = Omega;
	oev[4] = omega;
	oev[5] = f;
}

void Orbit::Coe2Rv()
{
	r = a * (1 - pow(e, 2)) / (1 + e * cos(f));
	v = sqrt(MU * (2 / r - 1 / a));
}

void Orbit::RaRp2SmaEcc()
{
	a = (r_a + r_p) / 2;
	e = (r_a - r_p) / (r_a + r_p);
}

void Orbit::SmaEcc2RaRp()
{
	r_a = a * (1 + e);
	r_p = a * (1 - e);
}

void Orbit::MeanMotion()
{
	// eccentric anomaly
	E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(f / 2));
	E = Mod(E, TWO_PI);
	// mean anomaly
	M = E - e * sin(E);
	M = Mod(M, TWO_PI);
	T = TWO_PI * sqrt(a * a * a / MU);
	n = sqrt(MU / a / a / a);
}

void Orbit::TrueMotion()
{
	E = SolveKepler(e, M);
	E = Mod(E, TWO_PI);
	f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
	f = Mod(f, TWO_PI);
}

void Orbit::Eci2Coe()
{
	double rv[6] = { r_eci[0], r_eci[1], r_eci[2], v_eci[0], v_eci[1], v_eci[2] };
	Eci2Oev(rv, oev);
	Oev2Coe();
}

void Orbit::Coe2Eci()
{
	Coe2Oev();
	double rv[6] = {};
	Oev2Eci(oev, rv);
	for (int i = 0; i < 3; i++)
	{
		r_eci[i] = rv[i];
		v_eci[i] = rv[i + 3];
	}
}

void Orbit::Geod2Ecef()
{
	lla[0] = lon;
	lla[1] = lat;
	lla[2] = (alt < 0.0) ? 0.0 : alt;
	Lla2R(lla, r_ecef);
}

void Orbit::Ecef2Geod()
{
	R2Lla(r_ecef, lla);
	lon = lla[0];
	lat = lla[1];
	alt = lla[2];
}

void Orbit::Eci2Ecef()
{
	J2k2Ecef(jd, r_eci, v_eci, r_ecef, v_ecef);
}

void Orbit::Ecef2Eci()
{
	J2k2Ecef(jd, r_ecef, v_ecef, r_eci, v_eci);
}

void Orbit::OrbProp(double dt)
{
	double gamma_2 = -J2 / 2 * (RE / a) * (RE / a);
	double eta = sqrt(1 - e * e);
	double a_r = a / r;
	double a_mean = a + a * gamma_2 * ((3 * cos(i) * cos(i) - 1) * (a_r * a_r * a_r - eta * eta * eta) + 3 * (1 - cos(i) * cos(i)) * (a_r * a_r * a_r) * cos(2 * omega + 2 * f));
	double C_J2 = 1.5 * J2 * RE * RE * sqrt(MU) / sqrt(pow(a_mean, 7));

	double dot_Omega = -C_J2 / ((1 - e * e) * (1 - e * e)) * cos(i);
	double dot_omega = -C_J2 / ((1 - e * e) * (1 - e * e)) * (2.5 * sin(i) * sin(i) - 2);
	double dot_M = -C_J2 / sqrt((1 - e * e) * (1 - e * e) * (1 - e * e)) * (1.5 * sin(i) * sin(i) - 1);

	double n = sqrt(MU / (a_mean * a_mean * a_mean));
	omega += dot_omega * dt;
	Omega += dot_Omega * dt;
	M += (n + dot_M) * dt;
	M = Mod(M, TWO_PI);
	E = SolveKepler(e, M);
	E = Mod(E, TWO_PI);
	f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
	f = Mod(f, TWO_PI);

	Coe2Oev();
	jd += (dt / 86400.0);
}


//本函数将轨道六根数（oev）转化为地心惯性坐标系（ECI）下的位置速度矢量
//输入 轨道六根数 [半长轴，离心率，倾角，升交点赤经，近地点角距，真近点角]
//输出 ECI系下的位置、速度
void Oev2Eci(double oev[], double rv[])
{
	double sma = oev[0];
	double ecc = oev[1];
	double inc = oev[2];
	double raan = oev[3];
	double argper = oev[4];
	double tanom = oev[5];

	double slr = sma * (1 - ecc * ecc);
	double rm = slr / (1 + ecc * cos(tanom));

	double arglat = argper + tanom;
	double sarglat = sin(arglat);
	double carglat = cos(arglat);

	double c4 = sqrt(MU / slr);
	double c5 = ecc * cos(argper) + carglat;
	double c6 = ecc * sin(argper) + sarglat;

	double sinc = sin(inc);
	double cinc = cos(inc);
	double sraan = sin(raan);
	double craan = cos(raan);

	//Position vector
	rv[0] = rm * (craan * carglat - sraan * sarglat * cinc);
	rv[1] = rm * (sraan * carglat + craan * sarglat * cinc);
	rv[2] = rm * sarglat * sinc;

	// Velocity vector
	rv[3] = -c4 * (c6 * craan + c5 * sraan * cinc);
	rv[4] = -c4 * (c6 * sraan - c5 * craan * cinc);
	rv[5] = c4 * c5 * sinc;
}

//本函数将地心惯性坐标系（ECI）下的位置速度矢量转化为轨道六根数（oev）
//输入 ECI系下的位置、速度
//输出 轨道六根数 [半长轴，离心率，倾角，升交点赤经，近地点角距，真近点角]
void Eci2Oev(double rv[], double* oev)
{
	// Position and velocity magnitude
	double r[3] = { rv[0], rv[1], rv[2] };
	double rmag = sqrt(Dot(r, r));
	double v[3] = { rv[3], rv[4], rv[5] };
	double vmag = sqrt(Dot(v, v));

	// position and velocity unit vectors
	double rhat[3] = { r[0] / rmag, r[1] / rmag, r[2] / rmag, };

	// angular momentum vectors
	double hv[3] = { 0.0 };
	hv[0] = r[1] * v[2] - r[2] * v[1];
	hv[1] = r[2] * v[0] - r[0] * v[2];
	hv[2] = r[0] * v[1] - r[1] * v[0];
	double hvmag = sqrt(Dot(hv, hv));
	double hhat[3] = { hv[0] / hvmag, hv[1] / hvmag, hv[2] / hvmag };

	// eccentricity vector
	double vtmp[3] = { v[0] / MU, v[1] / MU, v[2] / MU };
	double ecc0[3] = { 0.0 };
	ecc0[0] = vtmp[1] * hv[2] - vtmp[2] * hv[1];
	ecc0[1] = vtmp[2] * hv[0] - vtmp[0] * hv[2];
	ecc0[2] = vtmp[0] * hv[1] - vtmp[1] * hv[0];
	double ecc[3] = { ecc0[0] - rhat[0], ecc0[1] - rhat[1], ecc0[2] - rhat[2] };

	// semimajor axis
	double sma = 1 / (2 / rmag - vmag * vmag / MU);
	double p = hhat[0] / (1 + hhat[2]);
	double q = -hhat[1] / (1 + hhat[2]);
	if (hhat[0] == 0)
	{
		p = 0;
	}
	if (hhat[1] == 0)
	{
		q = 0;
	}
	double const1 = 1 / (1 + p * p + q * q);
	double fhat[3] = { 0.0 };
	fhat[0] = const1 * (1 - p * p + q * q);
	fhat[1] = const1 * 2 * p * q;
	fhat[2] = -const1 * 2 * p;

	double ghat[3] = { 0.0 };
	ghat[0] = const1 * 2 * p * q;
	ghat[1] = const1 * (1 + p * p - q * q);
	ghat[2] = const1 * 2 * q;

	double h = Dot(ecc, ghat);
	double xk = Dot(ecc, fhat);
	double x1 = Dot(r, fhat);
	double y1 = Dot(r, ghat);

	// orbital eccentricity
	double eccm = sqrt(h * h + xk * xk);

	// orbital inclination
	double inc = 2 * atan(sqrt(p * p + q * q));

	// true longitude
	double xlambdat = Mod(atan2(y1, x1), 2 * PI);
	// check for equatorial orbit
	double raan = Mod(atan2(p, q), 2 * PI);
	if (inc < 0.00000001)
	{
		raan = 0;
	}

	// check for circular orbit
	double argper = Mod(Mod(atan2(h, xk), 2 * PI) - raan, 2 * PI);
	if (eccm < 0.00000001)
	{
		argper = 0;
	}

	// true anomaly
	double tanom = Mod(xlambdat - raan - argper, 2 * PI);
	if (hhat[2] == -1)
	{
		tanom = Mod(xlambdat - raan - argper, 2 * PI);
	}

	// singular value when hhat(3) == -1
	if (hhat[2] == -1)
	{
		inc = PI;
		raan = PI;
		argper = Mod(-Mod(atan2(h, xk), 2 * PI) + raan, 2 * PI);
		tanom = Mod(xlambdat - raan - argper, 2 * PI);
	}

	// orbital element vector
	oev[0] = sma;
	oev[1] = eccm;
	oev[2] = inc;
	oev[3] = raan;
	oev[4] = argper;
	oev[5] = tanom;
}

double SolveKepler(double e, double M)
{
	double e1 = -1;
	double e2 = 1e-9;
	double e3 = 1e-6;
	int k1 = 1;
	int N = 50;

	if (M > TWO_PI)
	{
		M = Mod(M, TWO_PI);
	}

	double x0 = M;
	double v1 = x0 - e * sin(x0) - M;
	double v2 = 1 - e * cos(x0);

	// iteration for x1
	double x1 = -1;
	if (abs(v1) >= e1 && abs(v2) >= e1)
	{
		x1 = x0 - v1 / v2;
		while (abs(x1 - x0) > e2 && abs(v1) > e3 && k1 < N)
		{
			k1 = k1 + 1;
			x0 = x1;
			v1 = x0 - e * sin(x0) - M;
			v2 = 1 - e * cos(x0);
			x1 = x0 - v1 / v2;
		}
		if (k1 > N)
		{
			x1 = -1; // no solution
		}
	}
	double E = x1;
	double f = Mod(2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2)), TWO_PI);
	return f;
}

void Lla2R(double lla[], double r[])
{
	// 大地坐标系经纬高转ECEF坐标系位置
	// 输入：
	// 经纬高 lla(3)
	// lla(1) = longitude(rad)
	// lla(2) = latitude(rad)
	// lla(3) = altitude(km)
	// 输出：
	// ECEF坐标系位置 r(3)   (km)

	// 地球扁率平方
	//e2 = 0.00669;
	//e = 0.08181921805
	double e2 = ECC_EARTH * ECC_EARTH;

	double v = RE / sqrt(1 - e2 * sin(lla[1]) * sin(lla[1]));
	r[0] = (v + lla[2]) * cos(lla[1]) * cos(lla[0]);
	r[1] = (v + lla[2]) * cos(lla[1]) * sin(lla[0]);
	r[2] = ((1 - e2) * v + lla[2]) * sin(lla[1]);
}

void R2Lla(double r[], double lla[])
{
	//function: From positions in ECEF Coordinate to LLA
	//input : Positions in ECEF Coordinateguan 地球固连坐标系
	//output : Geodetic latitude, Geodetic longitudeand altitude 大地纬度、大地经度和高度

	double rx = r[0];
	double ry = r[1];
	double rz = r[2];

	double lon, lat_c, lat_d, v, alt;
	if (ry >= 0.0)
		lon = acos(rx / sqrt(rx * rx + ry * ry)); // geodetic longitude
	else
		lon = -acos(rx / sqrt(rx * rx + ry * ry));
	lat_c = atan(rz / sqrt(rx * rx + ry * ry)); // geocentric latitude

	// 地球扁率平方
	double e2 = ECC_EARTH * ECC_EARTH;

	lat_d = lat_c; // geodetic latitude : initial value
	int iter_max = 5;
	for (int i = 0; i < iter_max; i++)
	{
		v = RE / sqrt(1 - e2 * sin(lat_d) * sin(lat_d));
		lat_d = atan((rz + e2 * v * sin(lat_d)) / sqrt(rx * rx + ry * ry));
	}
	alt = rx / cos(lon) / cos(lat_d) - v;

	lla[0] = lon;
	lla[1] = lat_d;
	lla[2] = alt;
}

double JulianDate(int yr, int mon, int day, int hr, int minute, double sec)
{
	/*
	功能：
	计算儒略日
	输入：
	yr, mon, day, hr, minute, sec：年月日时分秒
	输出：
	double jd：儒略日
	**/
	double jd = 367.0 * yr - floor((7.0 * (yr + floor((mon + 9.0) / 12.0))) * 0.25)
		+ floor(275.0 * mon / 9.0) + day + 1721013.5 + ((sec / 60.0 + minute) / 60.0 + hr) / 24.0;
	return jd;
}

void Ecef2J2k(double jd, double r_ECEF[], double v_ECEF[], double r_J2000[], double v_J2000[])
{
	//  ECEF坐标系位置速度转J2000坐标系位置速度
	//	输入：
	//	儒略日            jd(days)
	//	ECEF坐标系位置   r_ECEF(3)    (km)
	//	ECEF坐标系速度   v_ECEF(3)    (km / s)
	//	输出：
	//	J2000坐标系位置   r_J2000(3)    (km)
	//	J2000坐标系速度   v_J2000(3)    (km / s)

	double T = (jd - JD2000) / 36525.0;  // 儒略世纪数
	
	// 岁差矩阵
	double zeta = 0.011180860 * T + 1.464e-6 * T * T + 8.7e-8 * T * T * T;
	double z = 0.011180860 * T + 5.308e-6 * T * T + 8.9e-8 * T * T * T;
	double theta = 0.009717173 * T - 2.068e-6 * T * T - 2.02e-7 * T * T * T;

	double R1[3][3] = { {cos(-zeta), sin(-zeta), 0.0}, {sin(-zeta), cos(-zeta), 0.0}, {0.0, 0.0, 1.0} };
	double R2[3][3] = { {cos(theta), 0.0, -sin(theta)}, {0.0, 1.0, 0.0}, {sin(theta), 0.0, cos(theta)} };
	double R3[3][3] = { {cos(-z), sin(-z), 0.0}, {-sin(-z), cos(-z), 0.0}, {0.0, 0.0, 1.0} };

	double P[3][3] = { 0.0 }, P_pre[3][3] = { 0.0 };
	M33T33(R3, R2, P_pre);
	M33T33(P_pre, R1, P);

	// 地球自转旋转矩阵
	double thast = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T + 0.0093104 * T * T - 6.2e-6 * T * T * T;
	thast = Mod(thast, 86400.0);
	thast = thast / 43200 * PI;
	double R[3][3] = { {cos(thast), sin(thast), 0.0}, {-sin(thast), cos(thast), 0.0}, {0.0, 0.0, 1.0} };

	// 逆矩阵
	double W[3][3] = { 0.0 }, W_T[3][3] = { 0.0 };
	M33T33(R, P, W);
	MatTrans(W, W_T);

	// 位置
	M33T31(W_T, r_ECEF, r_J2000);

	// 速度
	double omega_earth[3] = { 0.0, 0.0, -OMEGA_EARTH }, add[3] = { 0.0 };
	Cross(omega_earth, r_J2000, add);
	M33T31(W_T, v_ECEF, v_J2000);
	for (int i = 0; i < 3; i++)
		v_J2000[i] -= add[i];
}

void J2k2Ecef(double jd, double r_j2k[], double v_j2k[], double r_ecef[], double v_ecef[])
{
	//J2000坐标系位置速度转ECEF坐标系位置速度
	//输入：
	//儒略日            jd(days)
	//J2000坐标系位置   r_J2000(3)    (km)
	//J2000坐标系速度   v_J2000(3)    (km / s)
	//输出：
	//ECEF坐标系位置   r_j2k(3)    (km)
	//ECEF坐标系速度   v_j2k(3)    (km / s)

	double T = (jd - JD2000) / 36525.0;  // 儒略世纪数

	// 岁差矩阵
	double zeta = 0.011180860 * T + 1.464e-6 * T * T + 8.7e-8 * T * T * T;
	double z = 0.011180860 * T + 5.308e-6 * T * T + 8.9e-8 * T * T * T;
	double theta = 0.009717173 * T - 2.068e-6 * T * T - 2.02e-7 * T * T * T;

	double R1[3][3] = { {cos(-zeta), sin(-zeta), 0.0}, {sin(-zeta), cos(-zeta), 0.0}, {0.0, 0.0, 1.0} };
	double R2[3][3] = { {cos(theta), 0.0, -sin(theta)}, {0.0, 1.0, 0.0}, {sin(theta), 0.0, cos(theta)} };
	double R3[3][3] = { {cos(-z), sin(-z), 0.0}, {-sin(-z), cos(-z), 0.0}, {0.0, 0.0, 1.0} };

	double P[3][3] = { 0.0 }, P_pre[3][3] = { 0.0 };
	M33T33(R3, R2, P_pre);
	M33T33(P_pre, R1, P);

	// 地球自转旋转矩阵
	double thast = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T + 0.0093104 * T * T - 6.2e-6 * T * T * T;
	thast = Mod(thast, 86400.0);
	thast = thast / 43200 * PI;
	double R[3][3] = { {cos(thast), sin(thast), 0.0}, {-sin(thast), cos(thast), 0.0}, {0.0, 0.0, 1.0} };

	// 转换矩阵
	double W[3][3] = { 0.0 }, W_T[3][3] = { 0.0 };
	M33T33(R, P, W);

	// 位置
	M33T31(W, r_j2k, r_ecef);

	// 速度
	double omega_earth[3] = { 0.0, 0.0, -OMEGA_EARTH }, add[3] = { 0.0 };
	Cross(omega_earth, r_j2k, add);
	for (int i = 0; i < 3; i++)
		add[i] += v_j2k[i];
	M33T31(W, add, v_ecef);
}

void Hohmann(double R_init, double R_fin, double* dv)
{
	//霍曼变轨
	//输入 初始轨道半径 终端轨道半径 
	//输出 两次脉冲幅值
	//Initial circular velocity
	double v_init = sqrt(MU / R_init); //km / s

	//Final circular velocity
	double v_fin = sqrt(MU / R_fin); //km / s

	//Semi - major axis of transfer orbit(km)
	double a_tx = (R_init + R_fin) / 2;

	// Required velocities at periapse and apoapse of the transfer orbit
	double V_trans_a = sqrt(2 * MU / R_init - MU / a_tx); //Initial transfer vel.needed(km / s)
	double V_trans_b = sqrt(2 * MU / R_fin - MU / a_tx);  //Final transfer vel. (km / s)

	// Change in velocities needed
	dv[0] = V_trans_a - v_init; //km / s
	dv[1] = v_fin - V_trans_b;  //km / s
}

void J2Long(double oev[], double dt, double oev_t[])
{
	double a = oev[0];
	double e = oev[1];
	double i = oev[2];
	double Omega = oev[3];
	double omega = oev[4];
	double f = oev[5];
	double E = 2 * atan(sqrt((1 - e) / (1 + e)) * tan(f / 2));
	double M = E - e * sin(E);

	double r_norm = a * (1 - e * e) / (1 + e * cos(f));
	double gamma_2 = -J2 / 2 * (RE / a) * (RE / a);
	double eta = sqrt(1 - e * e);
	double a_r = a / r_norm;
	double a_mean = a + a * gamma_2 * ((3 * cos(i) * cos(i) - 1) * (a_r * a_r * a_r - eta * eta * eta) + 3 * (1 - cos(i) * cos(i)) * (a_r * a_r * a_r) * cos(2 * omega + 2 * f));
	double C_J2 = 1.5 * J2 * RE * RE * sqrt(MU) / sqrt(pow(a_mean, 7));

	double dot_Omega = -C_J2 / ((1 - e * e) * (1 - e * e)) * cos(i);
	double dot_omega = -C_J2 / ((1 - e * e) * (1 - e * e)) * (2.5 * sin(i) * sin(i) - 2);
	double dot_M = -C_J2 / sqrt((1 - e * e) * (1 - e * e) * (1 - e * e)) * (1.5 * sin(i) * sin(i) - 1);

	double n = sqrt(MU / (a_mean * a_mean * a_mean));
	omega += dot_omega * dt;
	Omega += dot_Omega * dt;
	M += (n + dot_M) * dt;
	M = Mod(M, TWO_PI);
	E = SolveKepler(e, M);
	E = Mod(E, TWO_PI);
	f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
	f = Mod(f, TWO_PI);

	for (int i = 0; i < 3; i++)
		oev_t[i] = oev[i];
	oev_t[3] = Omega;
	oev_t[4] = omega;
	oev_t[5] = f;
}