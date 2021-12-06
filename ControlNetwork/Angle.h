#pragma once

//枚举数据类型，用于代表角度形式
enum AngleStyle
{
	DEG,//十进制度
	DMS,//度分秒
	RAD//弧度制
};

/*******************************************************************
* 类名：CAngle													   *
*																   *
* 描述：角度类													   *
*																   *
* 使用方法：1、初始化：CAngle(double value=0,AngleStyle style=DMS) *
*			即默认值为0，角度为度分秒，如果value为a，AngleStyle为  *
*			Deg或RAD，则类对象的值便为a，类型为度或弧度。          *
*			2、获取指定类型的角度值的函数有两个					   *
*			第一个为double& operator() (AngleStyle style)，使用该  *
*			函数，返回引用类对象，类对象的值和类型改变（指定类型） *
*			第二个为double operator() (AngleStyle style) const，   *
*			const CAngle类型变量调用，该函数返回的类对象不能引用， *
*			且类对象的值和类型不变（原来类型）					   * 
*			3、角度加减，即CAngle operator + (const CAngle& m1,	   *
			const CAngle& m2)，CAngle operator - (const CAngle& m1,*
			const CAngle& m2)，该函数返回的类对象类型为RAD。       *
*																   *
* 历史：**日期** **理由** **签名**								   *
*       2021/7/6    无     苑艺博								   *
*																   *
* 外部类：无													   *
*******************************************************************/
class CAngle
{
public:
	CAngle(double value=0,AngleStyle style=DMS);
	~CAngle(void);
private:
	double dValue;//角度值
	AngleStyle  nCurStyle;//当前角度值类型
private:
	//设置常成员函数的作用：1.类成员不会被改变
	//2.可以被常类对象调用
	double Deg(double dDms) const;
	double Dms(double dDeg) const;
	
public:
	//获取指定的类型获取角度值，
    //由于返回的是dValue的引用，所以该值大小可以改变，即可以进行赋值
	double& operator() (AngleStyle style);   //CAanle ang1,   ang1.operator(DEG); angl1(DEG)

	//重载，获取指定的类型获取角度值，该值不可改变，const CAngle类型变量调用
	double operator() (AngleStyle style) const;
	//重载运算符+/-
    friend CAngle operator + (const CAngle& m1,const CAngle& m2);
	friend CAngle operator - (const CAngle& m1,const CAngle& m2);
};
