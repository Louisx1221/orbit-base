//@file       : orbit.h
//@autor      : github.com/louisx1221
//@date       : 2021/09/25

#ifndef ORBIT_H_
#define ORBIT_H_

class Orbit
{
private:
	
public:
	Orbit();
	Orbit(double oev_[], double jd_ = 0.0);
	~Orbit();
	void Oev2Coe();
	void Coe2Oev();
	void Coe2Rv();
	void RaRp2SmaEcc();
	void SmaEcc2RaRp();
	void MeanMotion();
	void TrueMotion();
	void Coe2Eci();
	void Eci2Coe();
	void Geod2Ecef();
	void Ecef2Geod();
	void Eci2Ecef();
	void Ecef2Eci();

	double a, e, i, Omega, omega, f, E, M, T, n;
	double r, v, r_a, r_p;
	double lon, lat, alt;
	double jd;
	double oev[6], r_eci[3], v_eci[3], r_ecef[3], v_ecef[3], lla[3];
	
};

//轨道六根数转化为位置速度
void Oev2Eci(double oev[], double rv[]);

//位置速度转换为轨道六根数
void Eci2Oev(double rv[], double oev[]);

//大地坐标系经纬高转WGS84坐标系位置
void Lla2R(double lla[], double r[]);

//WGS84坐标系位置转大地坐标系经纬高
void R2Lla(double r[], double lla[]);

//儒略日
double JulianDate(int yr, int mon = 1, int day = 1, int hr = 0, int minute = 0, double sec = 0.0);

//J2000坐标系位置速度转ECEF坐标系位置速度
void J2k2Ecef(double jd, double r_j2k[], double v_j2k[], double r_ecef[], double v_ecef[]);

//ECEF坐标系位置速度转J2000坐标系位置速度
void Ecef2J2k(double jd, double r_ECEF[], double v_ECEF[], double r_J2000[], double v_J2000[]);

//求解开普勒方程
double SolveKepler(double e, double M);

//霍曼变轨
void Hohmann(double r_init, double r_fin, double * dv);

//线性J2模型
void J2long(double oev[], double dt, double oev_t[]);

//lambert轨道机动
double ** solve_lambert(double r1[], double r2[], double tof, double cw, int multi_revs);

//单脉冲星下点轨迹调整
double ground_track_adjustment(double oev_0[], double long_1, double lat_1, double gmst0, int flag_ascending_descending);

//j2模型轨道积分
void orbJ2_rkf78(double x[], double h, double ti, double tf, double *xoutadr, int neq);

#endif