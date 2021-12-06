#include "pch.h"//预编译头
#include "CommonSurveyFunctions.h"
#include "math.h"
#include "Angle.h"


// 求平面上两点间距离
double Dist(double X1, double Y1,double X2, double Y2)
{
	double d;
	d=sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1));
	return d;
}

//重载，求空间上两点间距离
double Dist(double X1, double Y1, double Z1, 
			double X2, double Y2, double Z2)
{
	double d;
	d=sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1));
	return d;
}

//求两点的方位角
double  Azimuth(const double &X1, const double &Y1,
				const double &X2,const double &Y2)
{
		double dx, dy;
        dx = X2 - X1;
        dy = Y2 - Y1 + EPSILON;
        return PI - sgn(dy) * PI / 2 - atan(dx / dy);
}
//符号函数
int sgn(double x)
{
    if (x>=0) 
		return(1);//非负则返回1
	else
		return(-1);
}
