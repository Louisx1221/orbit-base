//此处定义轨道力学中的一些常用数值

const double PI = 3.14159265358979324;
const double TWO_PI = 2.0 * PI;
const double HALF_PI = PI / 2.0;
const double R2D = 180.0 / PI;
const double D2R = PI / 180.0;
const double G = 6.673e-20; //km ^ 3 / km / s ^ 2
const double RE = 6378.14; //km 地球半径
const double MU = 398600.4415; //km ^ 3 / s ^ 2 引力常数
const double MASS_EARTH = 5.9742e24; //kg
const double OMEGA_EARTH = 7.2921158553e-5; //rad / s 地球自转角速度
const double ECC_EARTH = 0.08181921805; // 地球扁率
const double J2 = 1.0826269e-3; //J2摄动

const double JD2000 = 2451545.0; // 2000年1月1日UT12:00儒略日

const double A_EARTH = 149598023; //km
const double E_EARTH = 0.016708617;
const double E_LONG_PERH = 102.93734808; //deg
const double E_PER = 0.99997862; //yrs

const double R_MOON = 1737.5; //km
const double MU_MOON = 4902.801; //km ^ 3 / s ^ 2 per JPL
const double MASS_MOON = 7.3483e22; //kg

const double AU = 149597870; //km
const double MU_SUN = 1.32712428e11; //km ^ 3 / s ^ 2

const double A_MARS = 1.52367934*au; //km
const double R_MARS = 3397; //km
const double MU_MARS = 4.305e4; //km ^ 3 / s ^ 2

const double A_VENUS = 108208601; //km