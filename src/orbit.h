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

//���������ת��Ϊλ���ٶ�
void Oev2Eci(double oev[], double rv[]);

//λ���ٶ�ת��Ϊ���������
void Eci2Oev(double rv[], double oev[]);

//�������ϵ��γ��תWGS84����ϵλ��
void Lla2R(double lla[], double r[]);

//WGS84����ϵλ��ת�������ϵ��γ��
void R2Lla(double r[], double lla[]);

//������
double JulianDate(int yr, int mon = 1, int day = 1, int hr = 0, int minute = 0, double sec = 0.0);

//J2000����ϵλ���ٶ�תECEF����ϵλ���ٶ�
void J2k2Ecef(double jd, double r_j2k[], double v_j2k[], double r_ecef[], double v_ecef[]);

//ECEF����ϵλ���ٶ�תJ2000����ϵλ���ٶ�
void Ecef2J2k(double jd, double r_ECEF[], double v_ECEF[], double r_J2000[], double v_J2000[]);

//��⿪���շ���
double SolveKepler(double e, double M);

//�������
void Hohmann(double r_init, double r_fin, double * dv);

//����J2ģ��
void J2long(double oev[], double dt, double oev_t[]);

//lambert�������
double ** solve_lambert(double r1[], double r2[], double tof, double cw, int multi_revs);

//���������µ�켣����
double ground_track_adjustment(double oev_0[], double long_1, double lat_1, double gmst0, int flag_ascending_descending);

//j2ģ�͹������
void orbJ2_rkf78(double x[], double h, double ti, double tf, double *xoutadr, int neq);

#endif