//常用测量计算函数
#pragma once

// 求平面上两点间距离
double Dist(double X1, double Y1,double X2, double Y2);

//重载，求空间上两点间距离
double Dist(double X1, double Y1, double Z1, 
			double X2, double Y2, double Z2);

//求两点的方位角
double  Azimuth(const double &X1, const double &Y1,
				const double &X2,const double &Y2);
//符号函数
int sgn(double x);