#include "pch.h"
#include "Angle.h"
#include "math.h"

//重载构造函数，有缺省值
CAngle::CAngle(double value,AngleStyle style)
{
	dValue=value;
	nCurStyle=style;
}

CAngle::~CAngle(void)
{
}
//重载（）函数
double& CAngle::operator() (AngleStyle style) //指定的类型获取角度值
{
	//double dAngleValue;
    if(style==DMS)
	{
		if(nCurStyle==DEG)
		{
            dValue=Dms(dValue);
		}
		else if(nCurStyle==RAD)
		{
			dValue=Dms(dValue*180.0/PI);
		}
		nCurStyle=DMS;
        		
	}
	else if(style==DEG)
	{
	    if(nCurStyle==DMS)
		{
            dValue=Deg(dValue);
		}
		else if(nCurStyle==RAD)
		{
			dValue=dValue*180.0/PI;
		}
       nCurStyle=DEG;
	}
	else
	{
       if(nCurStyle==DMS)
		{
            dValue=Deg(dValue)*PI/180;
		}
		else if(nCurStyle==DEG)
		{
			dValue=dValue*PI/180;
		}
        nCurStyle=RAD;
	}
	return dValue;
}
//重载（）函数，该函数是常函数，只能被常CAngle对象使用
double CAngle::operator() (AngleStyle style) const //指定的类型获取角度值
{
	double dAngleValue;
    if(style==DMS)
	{
		if(nCurStyle==DEG)
		{
            dAngleValue=Dms(dValue);
		}
		else if(nCurStyle==RAD)
		{
			dAngleValue=Dms(dValue*180.0/PI);
		}
		else
		{
			dAngleValue=dValue;
		}
        		
	}
	else if(style==DEG)
	{
	    if(nCurStyle==DMS)
		{
            dAngleValue=Deg(dValue);
		}
		else if(nCurStyle==RAD)
		{
			dAngleValue=dValue*180.0/PI;
		}
       else
		{
			dAngleValue=dValue;
		}
	}
	else
	{
       if(nCurStyle==DMS)
		{
            dAngleValue=Deg(dValue)*PI/180;
		}
		else if(nCurStyle==DEG)
		{
			dAngleValue=dValue*PI/180;
		}
        else
		{
			dAngleValue=dValue;
		}
	}
	return dAngleValue;
}


//私有成员，度分秒向十进制度转换
double CAngle::Deg(double dDms) const
{
	int iDeg,iMin;
	double dSec;

	iDeg = int(dDms + EPSILON);//度//加一个很小的数，以防止取整时的出错
	iMin = int((dDms - iDeg) * 100+ EPSILON);//分
	dSec = ((dDms - iDeg) * 100 - iMin) * 100 ;//秒

	return iDeg + (double)iMin / 60 + dSec / 3600;
}

//私有成员，十进制度向度分秒转换
double CAngle::Dms(double dDeg) const
{
	int iDeg,iMin;
	double dSec;
	double dTmp;

	iDeg = int(dDeg + EPSILON);//整数部分度
	dTmp = (dDeg - iDeg) * 60;//小数部分转换成分
	iMin = int(dTmp+ EPSILON);//取分的整数部分
	dSec = (dTmp - iMin) * 60;//截取秒

	return iDeg + (double)iMin / 100 + dSec / 10000;
}

//友元重载+函数
CAngle operator + (const CAngle& m1,const CAngle& m2)
{
       CAngle addAngle(0,RAD);
	   addAngle(RAD)=m1(RAD)+m2(RAD);
	   return addAngle;
}
//友元重载-函数
CAngle operator - (const CAngle& m1,const CAngle& m2)
{
       CAngle subAngle(0,RAD);
	   subAngle(RAD)=m1(RAD)-m2(RAD);
	   return subAngle;
}