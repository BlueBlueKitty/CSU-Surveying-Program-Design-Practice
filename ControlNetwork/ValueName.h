#pragma once
#include "Angle.h"

/*********************************************************
* 文件名：ValueName.h                                    *
*                                                        *
* 描述：该头文件主要定义一些控制网平差类所用到的结构体   *
		变量或枚举变量。包括：控制点结构体CrtlPoint、角度*
		观测值结构体ObsAngle、距离观测结构体ObaDist和    *
		左右角枚举变量AngleStyle					     *
**
* 历史：** 日期**    理由**  签名**                      *
*		  2021/7/6    无     苑艺博                      *
**
* 外部过程：无                                           *
*********************************************************/

//控制点类型
enum PointStyle
{
	KnownPoint,//已知点
	UnknownPoint//未知点
};

/**********************************************************
* 结构体名：CtrlPoint								      *
*														  *
* 描述：控制点结构体									  *
*														  *
* 使用方法：成员有坐标dx，dy和是否求出的状态bState        *
*														  *
* 历史：**日期** **理由** **签名**						  *
*       2021/7/6    无     苑艺博						  *
*********************************************************/
struct CtrlPoint
{
	CString strID;
	double dx;//坐标x
	double dy;//坐标y；
	bool bState;//是否求出的状态，TRUE求出，FALSE未求出
	int iNum;//编号
	PointStyle ePointStyle;
	double dE, dF;//误差椭圆的长半轴和短半轴
	CAngle dAlfa;//误差椭圆长半轴的方位角
	double dMx, dMy, dMk;//点位误差
};


/********************************************************************
* 结构体名：ObsAngle												*
*																	*
* 描述：角度观测值结构体											*
*																	*
* 使用方法：成员有站点StationPoint和照准点ObjPoint（CtrlPoint指针） *
*			和角度值AngleValue（CAngle类型）						*
*																	*
* 历史：**日期** **理由** **签名**									*
*       2021/7/6    无     苑艺博									*
********************************************************************/
struct ObsAngle
{
	CtrlPoint* StationPoint;//站点
	CtrlPoint* ObjPoint;//照准点
	CAngle AngleValue;//角度值
	double AngleError;//角度改正数
};


/********************************************************************
* 结构体名：ObsDist  												*
*																	*
* 描述：距离观测值结构体											*
*																	*
* 使用方法：成员有起点StartPoint和终点EndPoint（CtrlPoint指针）     *
*			和距离值DistValue（double 类型）						*
*																	*
* 历史：**日期** **理由** **签名**									*
*       2021/7/6    无     苑艺博									*
********************************************************************/
struct ObsDist
{
	CtrlPoint* StartPoint;//起点
	CtrlPoint* EndPoint;//终点
	double DistValue;//距离值
	double DistError;//距离改正数
};
