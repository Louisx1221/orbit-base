//@file       : math_plus.h
//@autor      : github.com/louisx1221
//@date       : 2021/09/27

#include <iostream>
#include "constants.h"
#include "math_plus.h"

double Atan3(double y, double x)
{
	/**
	功能：
	反正切
	输入：
	double x：x坐标
	double y：y坐标
	输出：
	out：反正切
	*/
	double out = fmod(atan2(y, x), TWO_PI);
	return out;
}

double DotProduct(double R1[], double R2[], int n)
{
	/*
	功能：向量点积,out=R1.R2
	输入：
		向量：double R1[]
		向量：double R2[]
		维数：unsigned long n
	输出：返回点积
	*/
	double out = 0.0;
	for (int j = 0; j < n; j++)
		out += R1[j] * R2[j];
	return out;
}

void CrossProduct(double a[], double b[], double *c)
{
	/*
	功能：向量叉乘
	输入：
	向量：double a[]
	向量：double a[]
	输出：
	叉乘向量：double c=cross(a，b)
	*/
	c[0] = a[1] * b[2] - a[2] * b[1];
	c[1] = a[2] * b[0] - a[0] * b[2];
	c[2] = a[0] * b[1] - a[1] * b[0];
}

void M33T31(double a[][3], double b[], double c[])
{
	/* 矩阵乘以向量 */
	for (int j = 0; j < 3; j++)
		c[j] = a[0][j] * b[0] + a[1][j] * b[1] + a[2][j] * b[2];
}

void M33T33(double a[][3], double b[][3], double c[][3])
{
	/* 矩阵乘以矩阵 */
	for (int i = 0; i < 3; i++)
	{
		c[i][0] = a[i][0] * b[0][0] + a[i][1] * b[1][0] + a[i][2] * b[2][0];
		c[i][1] = a[i][0] * b[0][1] + a[i][1] * b[1][1] + a[i][2] * b[2][1];
		c[i][2] = a[i][0] * b[0][2] + a[i][1] * b[1][2] + a[i][2] * b[2][2];
	}
}

void MatTrans(double m1[][3], double m2[][3])
{
	/* 矩阵转置 */
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
			m2[i][j] = m1[j][i];
	}
}

void R1(double x[], double theta, double* x1)
{
	x1[0] = x[0];
	x1[1] = x[1] * cos(theta) + x[2] * sin(theta);
	x1[2] = -x[1] * sin(theta) + x[2] * cos(theta);
}

void R2(double x[], double theta, double* x1)
{
	x1[0] = x[0] * cos(theta) - x[2] * sin(theta);
	x1[1] = x[1];
	x1[2] = x[0] * sin(theta) + x[2] * cos(theta);
}

void R3(double x[], double theta, double* x1)
{
	x1[0] = x[0] * cos(theta) + x[1] * sin(theta);
	x1[1] = -x[0] * sin(theta) + x[1] * cos(theta);
	x1[2] = x[2];
}

void ShengJin(double a, double b, double c, double d, double* result)
{
	/************************************************************************/
	/* 盛金公式求解三次方程的解 ax^3+bx^2+cx+d=0
	德尔塔f=B^2-4AC
	这里只要了实根，虚根需要自己再整理下拿出来
	输出 [flag, x1,x2,x3] 若flag为实根个数，若flag==1，则仅x1有效，若flag==2，仅x1、x2有效，以此类推；
	*/
	/************************************************************************/
	result[1] = 0;
	result[2] = 0;
	result[3] = 0;//初始化

	double A = b * b - 3 * a * c;
	double B = b * c - 9 * a * d;
	double C = c * c - 3 * b * d;
	double f = B * B - 4 * A * C;
	double i_value;
	double Y1, Y2;
	if (fabs(A) < 1e-6 && fabs(B) < 1e-6)//公式1 A=B=0
	{
		result[0] = 1;
		result[1] = (-b / (3 * a));
	}
	else if (f > 1e-6)      //公式2 f>0
	{
		Y1 = A * b + 3 * a * (-B + sqrt(f)) / 2;
		Y2 = A * b + 3 * a * (-B - sqrt(f)) / 2;
		double Y1_value = (Y1 / fabs(Y1)) * pow((double)fabs(Y1), 1.0 / 3);
		double Y2_value = (Y2 / fabs(Y2)) * pow((double)fabs(Y2), 1.0 / 3);
		result[0] = 1;
		result[1] = (-b - Y1_value - Y2_value) / (3 * a);
	}
	else if (fabs(f) < 1e-6)   //公式3 f=0
	{
		double K = B / A;
		result[0] = 2;
		result[1] = (-b / a + K);
		result[2] = (-K / 2);
	}
	else if (f < -1e-6)   //公式4 f<0
	{
		double T = (2 * A * b - 3 * a * B) / (2 * A * sqrt(A));
		double S = acos(T);
		result[0] = 3;
		result[1] = (-b - 2 * sqrt(A) * cos(S / 3)) / (3 * a);
		result[2] = (-b + sqrt(A) * (cos(S / 3) + sqrt(3.0) * sin(S / 3))) / (3 * a);
		result[3] = (-b + sqrt(A) * (cos(S / 3) - sqrt(3.0) * sin(S / 3))) / (3 * a);
	}
}