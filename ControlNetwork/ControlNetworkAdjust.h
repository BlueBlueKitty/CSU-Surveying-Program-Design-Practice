#pragma once
#include "ValueName.h"
#include "Matrix.h"

const double EPSILON = 1.0E-12;
const double PI = 4.0 * atan(1.0);


/*******************************************************************
* 类名：ControlNetworkAdjust								       *
*																   *
* 描述：控制网平差类											   *
*																   *
* 使用方法：													   *
*																   * 
* 历史：**日期** **理由** **签名**								   *
*       2021/7/6    无     苑艺博								   *
*																   *
* 外部类：无													   *
*******************************************************************/
class ControlNetworkAdjust
{
private:
	int m_KnownPointCount;//已知点数
	CtrlPoint* m_KnownPoint;//已知点指针
	int m_UnknownPointCount;//未知点数
	CtrlPoint* m_UnknownPoint;//未知点指针
	int m_ObsAngleCount;//观测角度数
	ObsAngle* m_ObsAngle;//观测角度指针
	int m_ObsDistCount;//观测边长数
	ObsDist* m_ObsDist;//观测边长指针

	//以下是绘图用到的变量
	double dXmin;//坐标x的最小值
	double dYmin;//坐标y的最小值
	double dDetX;//坐标x的最大值和最小值的差值
	double dDetY;//坐标y的最大值和最小值的差值
	double dScale;//绘图比例
public:
	ControlNetworkAdjust();//构造函数
	ControlNetworkAdjust(ControlNetworkAdjust& A);//拷贝构造函数
	~ControlNetworkAdjust();//析构函数
private:
	friend void ReadFileContent(CStdioFile& sf, CString& strFileContent);//从文件中读取数据到CString字符串中

	CAngle Azi(const CtrlPoint& P1, const CtrlPoint& P2);//已知两个控制点，求P1->P2的方位角
	double DIST(const CtrlPoint& P1, const CtrlPoint& P2);//已知两个控制点，求P1和P2的距离
	int SplitStringArray(CString str, char split, CStringArray& aStr);// 分割字符串，得到字符串数组


	CtrlPoint* SearchPoint(CString strName);//根据点名寻找符合的已知点或未知点，返回CtrlPoint指针变量
	bool JudgeUnknownPoint();//判断未知点是否全部求出，若全部求出返回真，反之返回假   
	ObsAngle* SearchObsAngle_StationPtIsLookedPt(CString strName);//寻找测站点为所找点且照准点坐标已求出的角度观测值
	ObsAngle* SearchObsAngle(CString strStationName, CString strObjName);//寻找测站点和照准点的点号符合的角度观测值
	ObsDist* SearchObsDist_EndPtIsLookedPt(CString strName);//寻找终点为所找点且起点坐标已求出的边长观测值
	ObsDist* SearchObsDist_StartPtIsLookedPt(CString strName);//寻找起点为所找点且终点坐标已求出的边长观测值
	ObsAngle* SearchZeroDirection(int iNum);//寻找任一角度观测值中测站的零方向所对应的角度观测值


	void ObsAngleErrorMatrix(CMatrix& B, CMatrix& L, CMatrix& P);//构建观测角度误差方程和权矩阵
	void ObsDistErrorMatrix(CMatrix& B, CMatrix& L, CMatrix& P);//构建观测边长误差方程和权矩阵
	void CalBLP(CMatrix& B, CMatrix& L, CMatrix& P);//计算B、L、P
	void CalAdjustedAngleDist(CMatrix& V);//计算角度和边长改正数以及改正后边长和角度
	void PrecisionEvaluate(CMatrix& Qx,double m0);//精度评定，包括点位误差和误差椭圆参数


	void XYCal();//计算dXmin,dYmin,dDetX,dDetY
	void SetScale(CRect& rect);//设置自适应窗口的绘图比例
	POINT LP2CP(const CtrlPoint& Point, CRect& rect);//从逻辑坐标转换到客户区坐标
	void TrianglePointDrawing(CDC* pDC, const POINT& CentralPoint);//以指定点为中心画小三角形来表示已知点


	void PrintCMatrix(CMatrix& A, CString strFileName);//输出矩阵
public:
	void ReadFileData(CStdioFile& sf);//从文件中读取数据
	bool CalUnknownPoint();//计算未知点近似坐标
	void PrintRoughCalResult(CString& strOut);//向文本框中输出概算坐标
	bool AdjustCal();//平差计算
	void PrintAdjustCalResult(CString& strOut);//向文本框中输出平差坐标，并存入文件中

	bool JudgeDrawConditon();//判断是否达到绘图条件
	void ControlNetworkDrawing(CDC* pDC, CRect& rect);//绘制三角网
};
