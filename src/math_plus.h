//@file       : math_plus.h
//@autor      : github.com/louisx1221
//@date       : 2021/09/27

#ifndef MATH_PLUS_H_
#define MATH_PLUS_H_

double Mod(double x, double y);
double Atan3(double y, double x);
double Dot(double R1[], double R2[], int n = 3);
void Cross(double a[], double b[], double* axb);
void M33T31(double a[][3], double b[], double c[]);
void M33T33(double a[][3], double b[][3], double c[][3]);
void MatTrans(double m1[][3], double m2[][3]);

//坐标旋转
void R1(double x[], double theta, double* x1);
void R2(double x[], double theta, double* x1);
void R3(double x[], double theta, double* x1);

//盛金公式求解三次方程的解 ax ^ 3 + bx ^ 2 + cx + d = 0
void ShengJin(double a, double b, double c, double d, double * result);

#endif