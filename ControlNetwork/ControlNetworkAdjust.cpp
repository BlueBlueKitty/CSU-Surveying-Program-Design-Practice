#include "pch.h"
#include "ControlNetworkAdjust.h"
#include "CommonSurveyFunctions.h"
#include <locale>

//构造函数
ControlNetworkAdjust::ControlNetworkAdjust()
{
	m_KnownPointCount = 0;
	m_KnownPoint = NULL;

	m_UnknownPointCount = 0;
	m_UnknownPoint = NULL;

	m_ObsAngleCount = 0;
	m_ObsAngle = NULL;

	m_ObsDistCount = 0;
	m_ObsDist = NULL;
}

//拷贝构造函数
ControlNetworkAdjust::ControlNetworkAdjust(ControlNetworkAdjust& A)
{
	m_KnownPointCount = A.m_KnownPointCount;
	m_KnownPoint = new CtrlPoint[m_KnownPointCount];
	m_KnownPoint = A.m_KnownPoint;

	m_UnknownPointCount = A.m_UnknownPointCount;
	m_UnknownPoint = new CtrlPoint[m_UnknownPointCount];
	m_UnknownPoint = A.m_UnknownPoint;

	m_ObsAngleCount = A.m_ObsAngleCount;
	m_ObsAngle = new ObsAngle[m_ObsAngleCount];
	m_ObsAngle = A.m_ObsAngle;

	m_ObsDistCount = A.m_ObsDistCount;
	m_ObsDist = new ObsDist[m_ObsDistCount];
	m_ObsDist = A.m_ObsDist;
}

//析构函数
ControlNetworkAdjust::~ControlNetworkAdjust()
{
	if (m_KnownPoint != NULL)
	{
		delete[] m_KnownPoint;
		m_KnownPoint = NULL;
	}

	if (m_UnknownPoint != NULL)
	{
		delete[] m_UnknownPoint;
		m_UnknownPoint = NULL;
	}

	if (m_ObsAngle != NULL)
	{
		delete[] m_ObsAngle;
		m_ObsAngle = NULL;
	}

	if (m_ObsDist != NULL)
	{
		delete[] m_ObsDist;
		m_ObsDist = NULL;
	}
}

//已知两个控制点，求P1->P2的方位角
CAngle ControlNetworkAdjust::Azi(const CtrlPoint& P1, const CtrlPoint& P2)
{
	CAngle angAzi;

	angAzi(RAD) = Azimuth(P1.dx, P1.dy, P2.dx, P2.dy);

	return angAzi;
}

//已知两个控制点，求P1和P2的距离
double ControlNetworkAdjust::DIST(const CtrlPoint& P1, const CtrlPoint& P2)
{
	return Dist(P1.dx, P1.dy, P2.dx, P2.dy);
}

//分割字符串
int ControlNetworkAdjust::SplitStringArray(CString str, char split, CStringArray& aStr)
{
	int startIdx = 0;
	int idx = str.Find(split, startIdx);
	aStr.RemoveAll();//先清空
	while (-1 != idx)
	{
		CString sTmp = str.Mid(startIdx, idx - startIdx);
		aStr.Add(sTmp);
		startIdx = idx + 1;
		idx = str.Find(split, startIdx);
	}
	CString sTmp = str.Right(str.GetLength() - startIdx);
	if (!sTmp.IsEmpty())
		aStr.Add(sTmp);
	return aStr.GetSize();
}



//根据点名寻找符合的已知点或未知点，返回CtrlPoint指针变量
CtrlPoint* ControlNetworkAdjust::SearchPoint(CString strName)
{
	strName.Trim();
	CString strPointID;//修整后的点名

	//遍历已知点
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		strPointID = m_KnownPoint[i].strID.Trim();
		if (strName == strPointID)
		{
			return &m_KnownPoint[i];
		}
	}

	//遍历未知点
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		strPointID = m_UnknownPoint[i].strID.Trim();
		if (strPointID == strName)
		{
			return &m_UnknownPoint[i];
		}
	}
}


//判断未知点是否全部求出，若全部求出返回真，反之返回假
bool ControlNetworkAdjust::JudgeUnknownPoint()
{
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		if (m_UnknownPoint[i].bState == FALSE)
		{
			return FALSE;
		}
	}

	return TRUE;
}

//寻找测站点为所找点且照准点坐标已求出的角度观测值
ObsAngle* ControlNetworkAdjust::SearchObsAngle_StationPtIsLookedPt(CString strName)
{
	strName.Trim();
	CString strStationPtID;//修整后角度观测值测站点点名

	//遍历角度观测值
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		strStationPtID = m_ObsAngle[i].StationPoint->strID.Trim();
		if (strName == strStationPtID && m_ObsAngle[i].ObjPoint->bState == TRUE)
		{
			return &m_ObsAngle[i];
		}
	}

	return NULL;
}

//寻找测站点和照准点的点号符合的角度观测值
ObsAngle* ControlNetworkAdjust::SearchObsAngle(CString strStationName, CString strObjName)
{
	strStationName.Trim();
	strObjName.Trim();

	CString strStationID;//修整后的测站点的点号
	CString strObjID;//修整后的照准点的点号

	//遍历角度观测值
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		strStationID = m_ObsAngle[i].StationPoint->strID.Trim();
		strObjID = m_ObsAngle[i].ObjPoint->strID.Trim();

		if (strStationID == strStationName && strObjID == strObjName)
		{
			return &m_ObsAngle[i];
		}
	}

	return NULL;
}

//寻找终点为所找点且起点坐标已求出的边长观测值
//输入：strName为所找点点名
//输出：满足条件的边长观测值的指针变量
ObsDist* ControlNetworkAdjust::SearchObsDist_EndPtIsLookedPt(CString strName)
{
	CString strEndID;//修整后边长终点点名

	strName.Trim();

	for (int i = 0; i < m_ObsDistCount; i++)
	{
		strEndID = m_ObsDist[i].EndPoint->strID.Trim();

		if (strName == strEndID && m_ObsDist[i].StartPoint->bState == TRUE)
		{
			return &m_ObsDist[i];
		}
	}

	return NULL;
}

//寻找起点为所找点且终点坐标已求出的边长观测值
//输入：strName为所找点点名
//输出：满足条件的边长观测值的指针变量
ObsDist* ControlNetworkAdjust::SearchObsDist_StartPtIsLookedPt(CString strName)
{
	CString strStartID;//修整后边长终点点名

	strName.Trim();

	for (int i = 0; i < m_ObsDistCount; i++)
	{
		strStartID = m_ObsDist[i].StartPoint->strID.Trim();

		if (strName == strStartID && m_ObsDist[i].EndPoint->bState == TRUE)
		{
			return &m_ObsDist[i];
		}
	}

	return NULL;
}

//寻找任一角度观测值中测站的零方向所对应的角度观测值
//输入：任一角度观测值观测值在数组中的位置iNum
//输出：该角度观测值中测站的零方向所对应的角度观测值，类型为ObsAngle指针变量
ObsAngle* ControlNetworkAdjust::SearchZeroDirection(int iNum)
{
	//向后遍历角度观测值
	for (int i = iNum; i >= 0; i--)
	{
		//如果角度观测值为0，说明该角度观测值所对应的方向即为零方向
		if (m_ObsAngle[i].AngleValue(DMS) == 0)
		{
			return &m_ObsAngle[i];
		}
	}

	return NULL;
}




//构建观测角度误差方程和权矩阵
//输入：无
//输出：B、L均为误差方程V=BX+L中角度和边长合矩阵，P为角度和边长合权阵,且在B、L、和P中角度在前，边长在后
void ControlNetworkAdjust::ObsAngleErrorMatrix(CMatrix& B, CMatrix& L, CMatrix& P)
{
	double p = 180 / PI * 3600;//使弧度制转为秒需乘以的常数

	//遍历角度观测值，构建观测角度误差方程和权矩阵
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		//计算照准方向的方位角近似值
		CAngle SightAzi = Azi(*m_ObsAngle[i].StationPoint, *m_ObsAngle[i].ObjPoint);

		//计算零方向的方位角近似值
		ObsAngle* ZeroDirection = SearchZeroDirection(i);
		CAngle ZeroAzi;
		if (ZeroDirection != NULL)
		{
			ZeroAzi = Azi(*ZeroDirection->StationPoint, *ZeroDirection->ObjPoint);
		}

		//计算照准方向的边长
		double dDist = DIST(*m_ObsAngle[i].StationPoint, *m_ObsAngle[i].ObjPoint);

		//计算误差方程B的系数
		double dA = sin(SightAzi(RAD)) / dDist * p / 1000;
		double dB = cos(SightAzi(RAD)) / dDist * (-p) / 1000;

		//计算常数项
		CAngle Angle = SightAzi - ZeroAzi;//计算夹角
		//如果计算结果小于0
		if (Angle(DEG) < 0)
		{
			Angle(DEG) = Angle(DEG) + 360;
		}
		CAngle AngleError = Angle - m_ObsAngle[i].AngleValue;
		double dAngleErrorSec = AngleError(DEG) * 3600;//转为秒

		//构建误差方程和权阵
		{
			//构建B
			//如果测站点是未知点
            if (m_ObsAngle[i].StationPoint->ePointStyle == UnknownPoint)
			{
				//计算各参数在误差方程中B中的列数
				int iNumx = 2 * m_ObsAngle[i].StationPoint->iNum;//坐标x在误差方程中B中的列数
				int iNumy = iNumx + 1;//坐标y在误差方程中B中的列数

				B(i, iNumx) = dA;
				B(i, iNumy) = dB;
			}
			//如果照准点是未知点
			if (m_ObsAngle[i].ObjPoint->ePointStyle == UnknownPoint)
			{
				//计算各参数在误差方程中B中的列数
				int iNumx = 2 * m_ObsAngle[i].ObjPoint->iNum;//坐标x在误差方程中B中的列数
				int iNumy = iNumx + 1;//坐标y在误差方程中B中的列数

				B(i, iNumx) = -dA;
				B(i, iNumy) = -dB;
			}

			//构建L
			L(i, 0) = dAngleErrorSec;

			//构建P
			P(i, i) = 1;
		}
	}
}

//构建观测边长误差方程和权矩阵
//输入：无
//输出：B、L均为误差方程V=BX+L中角度和边长合矩阵，P为角度和边长合权阵,且在B、L、和P中角度在前，边长在后
void ControlNetworkAdjust::ObsDistErrorMatrix(CMatrix& B, CMatrix& L, CMatrix& P)
{
	//遍历角度观测值，构建观测角度误差方程和权矩阵
	for (int i = 0; i < m_ObsDistCount; i++)
	{
		//计算边长的方位角近似值
		CAngle DistAzi = Azi(*m_ObsDist[i].StartPoint, *m_ObsDist[i].EndPoint);

		//计算边长近似值
		double dDist = DIST(*m_ObsDist[i].StartPoint, *m_ObsDist[i].EndPoint);

		//计算常数项
		double dDistErrorMm = (dDist - m_ObsDist[i].DistValue) * 1000;//转为毫米

		//构建观测边长误差方程和权矩阵
		{
			//构建B
			//如果起点为未知点
			if (m_ObsDist[i].StartPoint->ePointStyle == UnknownPoint)
			{
				//计算各参数在误差方程中B中的列数
				int iNumx = 2 * m_ObsDist[i].StartPoint->iNum;//坐标x在误差方程中B中的列数
				int iNumy = iNumx + 1;//坐标y在误差方程中B中的列数

				B(m_ObsAngleCount + i, iNumx) = -cos(DistAzi(RAD));
				B(m_ObsAngleCount + i, iNumy) = -sin(DistAzi(RAD));
			}

			//如果终点为未知点
			if (m_ObsDist[i].EndPoint->ePointStyle == UnknownPoint)
			{
				//计算各参数在误差方程中B中的列数
				int iNumx = 2 * m_ObsDist[i].EndPoint->iNum;//坐标x在误差方程中B中的列数
				int iNumy = iNumx + 1;//坐标y在误差方程中B中的列数

				B(m_ObsAngleCount + i, iNumx) = cos(DistAzi(RAD));
				B(m_ObsAngleCount + i, iNumy) = sin(DistAzi(RAD));
			}

			//构建L
			L(m_ObsAngleCount + i, 0) = dDistErrorMm;

			//构建P
			P(m_ObsAngleCount + i, m_ObsAngleCount + i) = 1000 / m_ObsDist[i].DistValue;
		}
	}
}

//计算B、L、P
void ControlNetworkAdjust::CalBLP(CMatrix& B, CMatrix& L, CMatrix& P)
{
	//计算B中角度观测值零方向的系统误差Z的系数
	int iNumz = 2 * m_UnknownPointCount;//Z在B中的列数

	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		//如果是零方向
		if (m_ObsAngle[i].AngleValue(DEG) == 0)
		{
			B(i, iNumz) = -1;

			//把该测站上其他观测值的Z的系数赋值
			for (int j = i + 1; j < m_ObsAngleCount; j++)
			{
				if (m_ObsAngle[j].AngleValue(DEG) != 0)
				{
					B(j, iNumz) = -1;
				}
				else
				{
					break;
				}
			}

			iNumz++;
		}
	}

	//计算误差方程和权阵
	ObsAngleErrorMatrix(B, L, P);
	ObsDistErrorMatrix(B, L, P);
}

//计算角度和边长改正数以及改正后的值
void ControlNetworkAdjust::CalAdjustedAngleDist(CMatrix& V)
{
	//计算改正后的角度观测值
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		//计算角度改正数,单位秒，保留两位小数
		double dAngleError = V(i, 0);
		CString strAngleError;
		strAngleError.Format(_T("%.2f"), dAngleError);
		m_ObsAngle[i].AngleError = _tstof(strAngleError);

		//计算改正后的角度值
		//由于角度改正数一般不大于60秒，如果大于，则说明观测值存在问题
		if (m_ObsAngle[i].AngleError < 60)
		{
			m_ObsAngle[i].AngleValue(DMS) = m_ObsAngle[i].AngleValue(DMS) + m_ObsAngle[i].AngleError / 10000;//计算改正后角度值
		}
		else
		{
			MessageBox(NULL, _T("角度改正数大于60秒，角度观测值存在问题，请检查！"), _T("警告"), MB_OK);
		}
	}

	//计算改正后的边长观测值
	for (int i = 0; i < m_ObsDistCount; i++)
	{
		//计算边长改正数，单位m，保留四位小数
		double dDistError= V(m_ObsAngleCount + i, 0) / 1000;
		CString strDistError;
		strDistError.Format(_T("%.4f"), dDistError);
		m_ObsDist[i].DistError = _tstof(strDistError);

		//计算改正后的边长
		m_ObsDist[i].DistValue = m_ObsDist[i].DistValue + m_ObsDist[i].DistError;
	}
}

//精度评定，包括点位误差和误差椭圆参数
void ControlNetworkAdjust::PrecisionEvaluate(CMatrix& Qx, double m0)
{
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		//计算未知点坐标x,y在矩阵Qx中的行数
		int iNumx = 2 * m_UnknownPoint[i].iNum;
		int iNumy = iNumx + 1;

		double Qxx = Qx(iNumx, iNumx);
		double Qyy = Qx(iNumy, iNumy);
		double Qxy = Qx(iNumx, iNumy);

		//计算点位误差
		m_UnknownPoint[i].dMx = m0 * sqrt(Qxx);
		m_UnknownPoint[i].dMy = m0 * sqrt(Qyy);
		m_UnknownPoint[i].dMk = sqrt(m_UnknownPoint[i].dMx * m_UnknownPoint[i].dMx + m_UnknownPoint[i].dMy * m_UnknownPoint[i].dMy);

		//计算误差椭圆参数
		double dAlfa;//长轴方位角
		dAlfa = 0.5 * atan(2 * Qxy / (Qxx - Qyy));
		if (dAlfa < 0)
		{
			dAlfa += PI;
		}
		m_UnknownPoint[i].dAlfa(RAD) = dAlfa;
		m_UnknownPoint[i].dE = m0 * sqrt(Qxx + Qxy * tan(m_UnknownPoint[i].dAlfa(RAD)));
		m_UnknownPoint[i].dF = m0 * sqrt(Qxx + Qxy * tan(m_UnknownPoint[i].dAlfa(RAD) + PI / 2));
	}
}





//计算dXmin,dYmin,dDetX,dDetY
void ControlNetworkAdjust::XYCal()
{
	double dXmax, dYmax;

	dXmin = m_KnownPoint[0].dx;
	dXmax = m_KnownPoint[0].dx;
	dYmin = m_KnownPoint[0].dy;
	dYmax = m_KnownPoint[0].dy;

	//遍历已知点
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		dXmin = m_KnownPoint[i].dx < dXmin ? m_KnownPoint[i].dx : dXmin;
		dXmax = m_KnownPoint[i].dx > dXmax ? m_KnownPoint[i].dx : dXmax;
		dYmin = m_KnownPoint[i].dy < dYmin ? m_KnownPoint[i].dy : dYmin;
		dYmax = m_KnownPoint[i].dy > dYmax ? m_KnownPoint[i].dy : dYmax;
	}

	//遍历未知点
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		dXmin = m_UnknownPoint[i].dx < dXmin ? m_UnknownPoint[i].dx : dXmin;
		dXmax = m_UnknownPoint[i].dx > dXmax ? m_UnknownPoint[i].dx : dXmax;
		dYmin = m_UnknownPoint[i].dy < dYmin ? m_UnknownPoint[i].dy : dYmin;
		dYmax = m_UnknownPoint[i].dy > dYmax ? m_UnknownPoint[i].dy : dYmax;
	}

	dDetX = dXmax - dXmin;
	dDetY = dYmax - dYmin;
}

//设置自适应窗口的绘图比例
void ControlNetworkAdjust::SetScale(CRect& rect)
{
	double a = rect.Width();
	double b = rect.Height();
	double ry = double(rect.Width()) * 3 / 4 / dDetY;//宽度扩大到2/3，方便在右边画坐标轴和比例尺
	double rx = double(rect.Height()) * 3 / 4 / dDetX;
	dScale = rx < ry ? rx : ry;
}

//从逻辑坐标转换到客户区坐标
POINT ControlNetworkAdjust::LP2CP(const CtrlPoint& Point, CRect& rect)
{
	POINT P;
	P.y = rect.Height() - (Point.dx - dXmin) * dScale;
	P.x = (Point.dy - dYmin) * dScale;
	return P;
}

//以指定点为中心画小三角形来表示已知点
void ControlNetworkAdjust::TrianglePointDrawing(CDC* pDC, const POINT& CentralPoint)
{
	POINT Vertice[3];//三角形三个顶点

	Vertice[0].x = CentralPoint.x;
	Vertice[0].y = CentralPoint.y - 10;

	Vertice[1].x = CentralPoint.x + 8;
	Vertice[1].y = CentralPoint.y + 6;

	Vertice[2].x = CentralPoint.x - 8;
	Vertice[2].y = CentralPoint.y + 6;

	for (int i = 0; i < 2; i++)
	{
		pDC->MoveTo(Vertice[i].x, Vertice[i].y);
		pDC->LineTo(Vertice[i + 1].x, Vertice[i + 1].y);
	}

	pDC->MoveTo(Vertice[0].x, Vertice[0].y);
	pDC->LineTo(Vertice[2].x, Vertice[2].y);
}

//绘制控制网
void ControlNetworkAdjust::ControlNetworkDrawing(CDC* pDC, CRect& rect)
{
	XYCal();//计算dXmin,dYmin,dDetX,dDetY
	SetScale(rect);//设置自适应窗口的绘图比例

	//三角网的原点，可以控制三角网的整体平移
	double dOrgNetX = rect.Width() / 4;
	double dOrgNetY = -rect.Height() / 5;

	CPen pen(PS_SOLID, 1, RGB(0, 0, 0));
	CPen* pOldPen = pDC->SelectObject(&pen);

	POINT StartPoint;//起点在图上的坐标
	POINT EndPoint;//终点在图上的坐标

	//在已知点上画三角形
	StartPoint = LP2CP(m_KnownPoint[0], rect);
	EndPoint = LP2CP(m_KnownPoint[1], rect);
	StartPoint.x = StartPoint.x + dOrgNetX;
	StartPoint.y = StartPoint.y + dOrgNetY;
	EndPoint.x = EndPoint.x + dOrgNetX;
	EndPoint.y = EndPoint.y + dOrgNetY;
	TrianglePointDrawing(pDC, StartPoint);
	TrianglePointDrawing(pDC, EndPoint);

	//根据已知点之间的斜率来画另一条平行线
	int dx = StartPoint.x - EndPoint.x;
	int dy = StartPoint.y - EndPoint.y;
	double dAngle;//垂直于已知点之间直线的直线与X轴的夹角
	if (dx == 0)
	{
		dAngle = PI / 2;
		StartPoint.x = StartPoint.x + 3;
		StartPoint.y = StartPoint.y;
		EndPoint.x = EndPoint.x + 3;
		EndPoint.y = EndPoint.y;
	}
	else
	{
		dAngle = PI / 2 - atan(dy / dx);
		StartPoint.x = StartPoint.x + 3 * cos(dAngle);
		StartPoint.y = StartPoint.y - 3 * sin(dAngle);
		EndPoint.x = EndPoint.x + 3 * cos(dAngle);
		EndPoint.y = EndPoint.y - 3 * sin(dAngle);
	}

	pDC->MoveTo(StartPoint);
	pDC->LineTo(EndPoint);

	//遍历角度观测值，画出照准边
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		StartPoint = LP2CP(*m_ObsAngle[i].StationPoint, rect);
		EndPoint = LP2CP(*m_ObsAngle[i].ObjPoint, rect);
		StartPoint.x = StartPoint.x + dOrgNetX;
		StartPoint.y = StartPoint.y + dOrgNetY;
		EndPoint.x = EndPoint.x + dOrgNetX;
		EndPoint.y = EndPoint.y + dOrgNetY;

		pDC->MoveTo(StartPoint);
		pDC->LineTo(EndPoint);
	}
	

	//画误差椭圆
	
	//误差椭圆的中心
	double dOrgX;
	double dOrgY;

	//遍历未知点，画误差椭圆
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		double dStartX, dStartY, dEndX, dEndY;

		//误差椭圆的中心坐标
		POINT Org;
		Org = LP2CP(m_UnknownPoint[i], rect);//坐标转换
		dOrgX = Org.x + dOrgNetX;
		dOrgY = Org.y + dOrgNetY;


		double dF = m_UnknownPoint[i].dF * 5000;
		double dE = m_UnknownPoint[i].dE * 5000;
		double dAlfa = m_UnknownPoint[i].dAlfa(RAD);

		//绘制短半轴
		dStartX = (dF * sin(dAlfa)) * dScale + dOrgX;
		dStartY = (-dF * cos(dAlfa)) * dScale + dOrgY;
		dEndX = (-dF * sin(dAlfa)) * dScale + dOrgX;
		dEndY = (dF * cos(dAlfa)) * dScale + dOrgY;
		pDC->MoveTo(dStartX, dStartY);
		pDC->LineTo(dEndX, dEndY);

		//绘制长半轴
		dStartX = (-dE * cos(dAlfa)) * dScale + dOrgX;
		dStartY = (-dE * sin(dAlfa)) * dScale + dOrgY;
		dEndX = (dE * cos(dAlfa)) * dScale + dOrgX;
		dEndY = (dE * sin(dAlfa)) * dScale + dOrgY;
		pDC->MoveTo(dStartX, dStartY);
		pDC->LineTo(dEndX, dEndY);

		double ex, fy;
		ex = dE;
		fy = 0;
		//转换到长半轴方向上
		dStartX = (ex * cos(dAlfa)
			- fy * sin(dAlfa)
			) * dScale + dOrgX;
		dStartY = (fy * cos(dAlfa)
			+ ex * sin(dAlfa)
			) * dScale + dOrgY;
		pDC->MoveTo(dStartX, dStartY);

		for (int i = 6; i <= 360; i += 6)
		{
			//在坐标轴方向的坐标
			ex = dE * cos((i / 180.0) * PI);
			fy = dF * sin((i / 180.0) * PI);

			//转换到长半轴方向上
			dEndX = (ex * cos(dAlfa)
				- fy * sin(dAlfa)
				) * dScale + dOrgX;
			dEndY = (fy * cos(dAlfa)
				+ ex * sin(dAlfa)
				) * dScale + dOrgY;
			pDC->LineTo(dEndX, dEndY);
		}
	}
	



	//画坐标轴
	double dLenX = double(rect.Height()) / 10;//坐标轴X轴的长度
	double dLenY = double(rect.Width())  / 10;//坐标轴Y轴的长度

	//坐标轴的原点坐标
	double dOrgCoordX =  rect.left + double(rect.Width()) / 10;
	double dOrgCoordY = rect.bottom - double(rect.Height()) / 10;

	//画Y轴
	pDC->MoveTo(dOrgCoordX, dOrgCoordY);
	pDC->LineTo(dOrgCoordX + dLenY, dOrgCoordY);
	pDC->LineTo(dOrgCoordX + dLenY - 5, dOrgCoordY - 5);
	pDC->MoveTo(dOrgCoordX + dLenY, dOrgCoordY);
	pDC->LineTo(dOrgCoordX + dLenY - 5, dOrgCoordY + 5);

	//画X轴
	pDC->MoveTo(dOrgCoordX, dOrgCoordY);
	pDC->LineTo(dOrgCoordX, dOrgCoordY - dLenX);
	pDC->LineTo(dOrgCoordX - 5, dOrgCoordY - dLenX + 5);
	pDC->MoveTo(dOrgCoordX, dOrgCoordY - dLenX);
	pDC->LineTo(dOrgCoordX + 5, dOrgCoordY - dLenX + 5);



	//画比例尺
	//比例尺的起点
	double dOrgScaleX = dOrgCoordX + dLenX + double(rect.Width()) / 10;
	double dOrgScaleY = dOrgCoordY;
	pDC->MoveTo(dOrgScaleX, dOrgScaleY);
	pDC->LineTo(dOrgScaleX + 200 * dScale, dOrgScaleY);
	pDC->LineTo(dOrgScaleX + 200 * dScale, dOrgScaleY - 2);
	pDC->MoveTo(dOrgScaleX, dOrgScaleY);
	pDC->LineTo(dOrgScaleX, dOrgScaleY - 2);


	pDC->SelectObject(pOldPen);
	pen.DeleteObject();




	CFont Font;//创建字体
	VERIFY(Font.CreateFont(
		10,                        // nHeight
		0,                         // nWidth
		0,                         // nEscapement
		0,                         // nOrientation*
		FW_NORMAL,                 // nWeight
		FALSE,                     // bItalic
		FALSE,                     // bUnderline
		0,                         // cStrikeOut
		ANSI_CHARSET,              // nCharSet
		OUT_DEFAULT_PRECIS,        // nOutPrecision
		CLIP_DEFAULT_PRECIS,       // nClipPrecision
		DEFAULT_QUALITY,           // nQuality
		DEFAULT_PITCH | FF_SWISS,  // nPitchAndFamily
		_T("宋体")));                 // lpszFacename

	//pDC->SetTextColor(RGB(0, 0, 0));///改变字体颜色

	CFont* pOldFont = pDC->SelectObject(&Font);

	//写点名
	POINT Point;
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		Point = LP2CP(m_KnownPoint[i], rect);
		Point.x = Point.x + dOrgNetX;
		Point.y = Point.y + dOrgNetY;

		pDC->TextOutW(Point.x, Point.y + 10, m_KnownPoint[i].strID);
	}

	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		Point = LP2CP(m_UnknownPoint[i], rect);
		Point.x = Point.x + dOrgNetX;
		Point.y = Point.y + dOrgNetY;

		pDC->TextOutW(Point.x, Point.y + 10, m_UnknownPoint[i].strID);
	}

	//写坐标轴名
	pDC->TextOutW(dOrgCoordX - 12, dOrgCoordY - dLenX - 2, _T("X"));
	pDC->TextOutW(dOrgCoordX + dLenY - 3, dOrgCoordY + 5, _T("Y"));

	pDC->SelectObject(pOldPen);
	Font.DeleteObject();


	//写比例尺名
	pDC->TextOutW(dOrgScaleX + 70 * dScale, dOrgScaleY + 6, _T("200m"));
}




//从文件中读取数据
void ControlNetworkAdjust::ReadFileData(CStdioFile& sf)
{
	CString strLine;
	CStringArray aStrTmp;

	//读取已知点
	sf.ReadString(strLine);
	m_KnownPointCount = _ttoi(strLine);
	m_KnownPoint = new CtrlPoint[m_KnownPointCount];
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		sf.ReadString(strLine);
		SplitStringArray(strLine, ',', aStrTmp);
		m_KnownPoint[i].strID = aStrTmp[0];
		m_KnownPoint[i].dx = _tstof(aStrTmp[1]);
		m_KnownPoint[i].dy = _tstof(aStrTmp[2]);
		m_KnownPoint[i].bState = TRUE;
		m_KnownPoint[i].ePointStyle = KnownPoint;
	}

	//读取未知点
	sf.ReadString(strLine);
	m_UnknownPointCount = _ttoi(strLine);
	m_UnknownPoint = new CtrlPoint[m_UnknownPointCount];
	sf.ReadString(strLine);
	SplitStringArray(strLine, ',', aStrTmp);
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		m_UnknownPoint[i].strID = aStrTmp[i];
		m_UnknownPoint[i].dx = -1;
		m_UnknownPoint[i].dy = -1;//坐标x、y赋值为-1，方便后面判断是否已经进行概算
		m_UnknownPoint[i].bState = FALSE;
		m_UnknownPoint[i].iNum = i;
		m_UnknownPoint[i].ePointStyle = UnknownPoint;
		m_UnknownPoint[i].dE = -1;//长轴赋值为-1，方便后面判断是否已经平差计算
	}

	//读取边长观测值
	sf.ReadString(strLine);
	m_ObsDistCount = _ttoi(strLine);
	m_ObsDist = new ObsDist[m_ObsDistCount];
	for (int i = 0; i < m_ObsDistCount; i++)
	{
		sf.ReadString(strLine);
		SplitStringArray(strLine, ',', aStrTmp);
		m_ObsDist[i].StartPoint = SearchPoint(aStrTmp[0]);
		m_ObsDist[i].EndPoint = SearchPoint(aStrTmp[1]);
		m_ObsDist[i].DistValue = _tstof(aStrTmp[2]);
	}

	//读取角度观测值
	sf.ReadString(strLine);
	m_ObsAngleCount = _ttoi(strLine);
	m_ObsAngle = new ObsAngle[m_ObsAngleCount];
	for (int i = 0; i < m_ObsAngleCount; i++) 
	{
		sf.ReadString(strLine);
		SplitStringArray(strLine, ',', aStrTmp);
		m_ObsAngle[i].StationPoint = SearchPoint(aStrTmp[0]);
		m_ObsAngle[i].ObjPoint = SearchPoint(aStrTmp[1]);
		m_ObsAngle[i].AngleValue(DMS)= _tstof(aStrTmp[2]);
	}
}

//计算未知点近似坐标
bool ControlNetworkAdjust::CalUnknownPoint()
{
	if (m_UnknownPointCount == 0)
	{
		MessageBox(NULL, _T("未读取数据！"), _T("警告"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dx >= 0)
	{
		MessageBox(NULL, _T("已经完成坐标概算！"), _T("警告"), MB_OK);
		return FALSE;
	}
	else
	{
		do
		{
			//遍历未知点，计算未知点的近似坐标
			for (int i = 0; i < m_UnknownPointCount; i++)
			{
				//如果该未知点的坐标未求出
				if (m_UnknownPoint[i].bState == FALSE)
				{
					ObsAngle* LookedAngle2 = NULL;//测站点和照准点点均已求出的角度观测值
					ObsAngle* LookedAngle1 = NULL;//测站点坐标已求出且照准点为所求未知点的角度观测值

					//寻找终点为未知点且起点坐标已求出的边长观测值
					ObsDist* LookedDist = SearchObsDist_EndPtIsLookedPt(m_UnknownPoint[i].strID);
					
					//如果边长观测值存在
					if (LookedDist != NULL)
					{
						//寻找测站点为Dist起点且照准点坐标已求出的角度观测值
						LookedAngle2 = SearchObsAngle_StationPtIsLookedPt(LookedDist->StartPoint->strID);

						//寻找照准边为Dist所在边的角度观测值
						LookedAngle1 = SearchObsAngle(LookedDist->StartPoint->strID, LookedDist->EndPoint->strID);
					}
					else
					{
						//寻找起点为未知点且终点坐标已求出的边长观测值
						LookedDist = SearchObsDist_StartPtIsLookedPt(m_UnknownPoint[i].strID);

						if (LookedDist != NULL)
						{
							//寻找测站点为Dist起点且照准点坐标已求出的角度观测值
							LookedAngle2 = SearchObsAngle_StationPtIsLookedPt(LookedDist->EndPoint->strID);

							//寻找照准边为Dist所在边的角度观测值
							LookedAngle1 = SearchObsAngle(LookedDist->EndPoint->strID, LookedDist->StartPoint->strID);
						}
					}
						
					//如果两个角度观测值存在，证明该未知点坐标可求出
					if (LookedAngle1 != NULL && LookedAngle2 != NULL)
					{
						//计算LookedAngle1和LookedAngle2之间的夹角
						CAngle Angle = LookedAngle1->AngleValue - LookedAngle2->AngleValue;

						//计算LookedAngle2所在照准边的方位角
						CAngle LookedAngle2Azi = Azi(*LookedAngle2->StationPoint, *LookedAngle2->ObjPoint);

						//计算LookedAngle1所在照准边的方位角
						CAngle LookedAngle1Azi = LookedAngle2Azi + Angle;

						//如果求得LookedAngle1所在照准边的方位角大于等于360,减去360
						if (LookedAngle1Azi(DEG) > 360)
						{
							LookedAngle1Azi(DEG) = LookedAngle1Azi(DEG) - 360;
						}

						double dx, dy;//坐标增量

						//计算坐标增量
						dx = LookedDist->DistValue * cos(LookedAngle1Azi(RAD));
						dy = LookedDist->DistValue * sin(LookedAngle1Azi(RAD));

						//计算照准点的近似坐标
						LookedAngle1->ObjPoint->dx = LookedAngle1->StationPoint->dx + dx;
						LookedAngle1->ObjPoint->dy = LookedAngle1->StationPoint->dy + dy;
						LookedAngle1->ObjPoint->bState = TRUE;
					}
				}
			}
		} while (!JudgeUnknownPoint());

		return TRUE;
	}
}

//输出矩阵到文件
void ControlNetworkAdjust::PrintCMatrix(CMatrix& A, CString strFileName)
{
	CStdioFile sf;
	if (!sf.Open(strFileName, CFile::modeWrite | CFile::modeCreate)) return;

	CString strResult;
	CString strContent;

	for (int i = 0; i < A.Row(); i++)
	{
		for (int j = 0; j < A.Col(); j++)
		{
			strResult.Format(_T("%-15.4f"), A(i, j));
			strContent += strResult;
		}
		strContent += _T("\n");
	}

	sf.WriteString(strContent);

	sf.Close();
}

//向文本框中输出概算坐标
void ControlNetworkAdjust::PrintRoughCalResult(CString& strOut)
{
	strOut.Empty();

	CString strResult;

	strOut += (_T("*****************坐标概算*****************\r\n\r\n"));

	strResult.Format(_T("%-20s%-22s%-22s\r\n\r\n"), _T("点号"), _T("X坐标"), _T("Y坐标"));
	strOut += strResult;

	//输出已知点坐标
	strOut += (_T("*****************已知点*******************\r\n"));
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-20.6f%-20.6f\r\n"), m_KnownPoint[i].strID, m_KnownPoint[i].dx, m_KnownPoint[i].dy);
		strOut += strResult;
	}

	//输出未知点坐标
	strOut += (_T("\r\n*****************未知点*******************\r\n"));
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-20.6f%-20.6f\r\n"), m_UnknownPoint[i].strID, m_UnknownPoint[i].dx, m_UnknownPoint[i].dy);
		strOut += strResult;
	}

	return;
}

//向文本框中输出平差坐标，并存入文件中
void ControlNetworkAdjust::PrintAdjustCalResult(CString& strOut)
{
	CString strFileData;
	CString strResult;

	strFileData += (_T("******************************平差计算************************************\n"));

	//输出方向观测成果
	strFileData += (_T("\n******************************方向观测成果********************************\n"));
	strResult.Format(_T("%-13s%-13s%-13s%-10s%-17s%-10s\n"), _T("测站"), _T("照准"), _T("方向值(dms)"), _T("改正数(s)"),
		_T("平差后值(dms)"), _T("备注"));
	strFileData += strResult;
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		//如果是零方向
		if (m_ObsAngle[i].AngleValue(DEG) == 0)
		{
			strResult.Format(_T("%-15s%-15s%-18.6f\n"), m_ObsAngle[i].StationPoint->strID,
				m_ObsAngle[i].ObjPoint->strID, m_ObsAngle[i].AngleValue(DMS));
			strFileData += strResult;
		}
		else//如果不是零方向
		{
			strResult.Format(_T("%-15s%-15s%-18.6f%-15.2f%-18.6f\n"), m_ObsAngle[i].StationPoint->strID, m_ObsAngle[i].ObjPoint->strID,
				m_ObsAngle[i].AngleValue(DMS) - m_ObsAngle[i].AngleError / 10000, m_ObsAngle[i].AngleError, m_ObsAngle[i].AngleValue(DMS));
			strFileData += strResult;
		}
	}

	//输出距离观测成果
	strFileData += (_T("\n******************************距离观测成果********************************\n"));
	strResult.Format(_T("%-13s%-15s%-14s%-10s%-13s%-20s\n"), _T("测站"), _T("照准"), _T("距离(m)"), _T("改正数(m)"), _T("平差后值(m)"),
		_T("方位角(dms)"));
	strFileData += strResult;
	for (int i = 0; i < m_ObsDistCount; i++)
	{
		CAngle DistAzi = Azi(*m_ObsDist[i].StartPoint, *m_ObsDist[i].EndPoint);//计算方位角
		strResult.Format(_T("%-15s%-15s%-18.4f%-15.4f%-18.4f%-18.6f\n"), m_ObsDist[i].StartPoint->strID, m_ObsDist[i].EndPoint->strID,
			m_ObsDist[i].DistValue - m_ObsDist[i].DistError, m_ObsDist[i].DistError, m_ObsDist[i].DistValue, DistAzi(DMS));
		strFileData += strResult;
	}

	//输出平面点位误差
	strFileData += (_T("\n******************************平面点位误差********************************\n"));
	strResult.Format(_T("%-13s%-12s%-10s%-15s%-16s%-20s\n"), _T("点名"), _T("长轴(m)"), _T("短轴(m)"), _T("长轴方位角(dms)"),
		_T("点位中误差(m)"), _T("备注"));
	strFileData += strResult;
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-15.4f%-15.4f%-23.6f%-18.4f\n"), m_UnknownPoint[i].strID, m_UnknownPoint[i].dE,
			m_UnknownPoint[i].dF, m_UnknownPoint[i].dAlfa(DMS), m_UnknownPoint[i].dMk);
		strFileData += strResult;
	}

	//输出控制点成果
	strFileData += (_T("\n******************************控制点成果*********************************\n"));
	strResult.Format(_T("%-19s%-20s%-20s%-13s%-20s\n"), _T("点名"), _T("X(m)"), _T("Y(m)"), _T("H(m)"), _T("备注"));
	strFileData += strResult;
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-18.4f%-18.4f%-17s%-20s\n"), m_KnownPoint[i].strID, m_KnownPoint[i].dx, m_KnownPoint[i].dy, _T(" "), _T("未知点"));
		strFileData += strResult;
	}
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-18.4f%-18.4f%-17s%-20s\n"), m_UnknownPoint[i].strID, m_UnknownPoint[i].dx, m_UnknownPoint[i].dy, _T(" "), _T(" "));
		strFileData += strResult;
	}

	//输出到文件中
	CString strFileName = _T("AdjustReult.dat");
	CStdioFile sf;
	if (!sf.Open(strFileName, CFile::modeReadWrite | CFile::modeCreate)) return;

	sf.WriteString(strFileData);
	MessageBox(NULL, _T("平差结果已保存到“AdjustReult.dat”文件中！"), _T("提示"), MB_OK);

	sf.SeekToBegin();//文件指针定位到开头

	ReadFileContent(sf,strOut);//从文件读取数据

	sf.Close();
}


//判断是否达到绘图条件
bool ControlNetworkAdjust::JudgeDrawConditon()
{
	if (m_UnknownPointCount == 0)
	{
		MessageBox(NULL, _T("未读取数据！"), _T("警告"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dx < 0)
	{
		MessageBox(NULL, _T("未坐标概算！"), _T("警告"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dE < 0)
	{
		MessageBox(NULL, _T("未平差计算！"), _T("警告"), MB_OK);
		return FALSE;
	}
	else
	{
		return TRUE;
	}
}

//平差计算
bool ControlNetworkAdjust::AdjustCal()
{
	if (m_UnknownPointCount == 0)
	{
		MessageBox(NULL, _T("未读取数据！"), _T("警告"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dx < 0)
	{
		MessageBox(NULL, _T("未坐标概算！"), _T("警告"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dE >= 0)
	{
		MessageBox(NULL, _T("已经完成平差计算！"), _T("警告"), MB_OK);
		return FALSE;
	}
	else
	{
		//定义误差方程和权阵
		CMatrix X(m_UnknownPointCount * 3 + m_KnownPointCount, 1);
		CMatrix B(m_ObsAngleCount + m_ObsDistCount, m_UnknownPointCount * 3 + m_KnownPointCount);
		CMatrix L(m_ObsAngleCount + m_ObsDistCount, 1);
		CMatrix P(m_ObsAngleCount + m_ObsDistCount, m_ObsAngleCount + m_ObsDistCount);
		CMatrix V(m_ObsAngleCount + m_ObsDistCount, 1);


		CMatrix VT;
		CMatrix BT;//B的转置矩阵
		CMatrix NBB;
		CMatrix NBB1;//NBB的逆矩阵

		//平差计算直到满足精度要求
		double dXmax;//X中坐标x或y增量的最大值
		do
		{
			//计算B、L、P
			CalBLP(B, L, P);

			BT = ~B;//求B的转置矩阵
			NBB = BT * P * B;//求NBB
			NBB1 = NBB.Inv();//求NBB的逆矩阵
			X = -1 * NBB1 * BT * P * L;//求X
			V = B * X + L;//计算改正数V

			//求未知点坐标的近似平差值
			for (int i = 0; i < m_UnknownPointCount; i++)
			{
				m_UnknownPoint[i].dx += X(2 * i, 0) / 1000;
				m_UnknownPoint[i].dy += X(2 * i + 1, 0) / 1000;
			}

			//求dXmax
			dXmax = fabs(X(0, 0));
			for (int i = 0; i < 2 * m_UnknownPointCount; i++)
			{
				dXmax = fabs(X(i, 0)) > dXmax ? fabs(X(i, 0)) : dXmax;
			}
		} while (dXmax > 0.1);

		//将零方向的角度改正归算到该测站其他照准边上
		double V0;//记录零方向的角度改正数
		for (int i = 0; i < m_ObsAngleCount; i++)
		{
			//如果是零方向
			if (m_ObsAngle[i].AngleValue(DEG) == 0)
			{
				V0 = V(i, 0);
				V(i, 0) = 0;//零方向改正数化为0

				//如果是该测站上其他照准方向
				for (int j = i + 1; j < m_ObsAngleCount; j++)
				{
					if (m_ObsAngle[j].AngleValue(DEG) != 0)
					{
						V(j, 0) = V(j, 0) - V0;
					}
					else
					{
						break;
					}
				}
			}
		}

		//计算边长和角度观测值的改正数以及改正后的值
		CalAdjustedAngleDist(V);

		//计算单位权中误差
		VT = ~V;
		CMatrix VV = VT * P * V;
		double m0 = sqrt(VV(0, 0) / double(V.Row() - X.Row())) / 1000;

		//精度评定，计算点位误差和误差椭圆参数
		PrecisionEvaluate(NBB1, m0);

		return TRUE;
	}
}



/*
//计算未知点的近似坐标
void ControlNetworkAdjust::CalUnknownPoint()
{
	do
	{
		//遍历未知点，计算未知点的近似坐标
		for (int i = 0; i < m_UnknownPointCount; i++)
		{
			//如果未知点的坐标未求出
			if (m_UnknownPoint[i].bState == FALSE)
			{
				//寻找未知点为照准点且测站点坐标已求出的角度观测值
				ObsAngle* LookedAngle1 = SearchObsAngle_ObjPtIsLookedPt(m_UnknownPoint[i].strID);

				//寻找测站点和LookedAngle1的相同且照准点坐标已求出的角度观测值
				ObsAngle* LookedAngle2 = SearchObsAngle_StationPtIsLookedPt(LookedAngle1->StationPoint->strID);

				//如果LookedAngle1和LookedAngle2存在
				if (LookedAngle1 != NULL && LookedAngle2 != NULL)
				{
					//寻找LookedAngle1照准边所在边边长观测值
					ObsDist* Dist = SearchObsDist(LookedAngle1->StationPoint->strID, LookedAngle1->ObjPoint->strID);

					//如果LookedAngle1照准边存在，证明该未知点的坐标可以求出
					if (Dist != NULL)
					{
						//计算LookedAngle1和LookedAngle2之间的夹角
						CAngle Angle = LookedAngle1->AngleValue - LookedAngle2->AngleValue;

						//计算LookedAngle2所在照准边的方位角
						CAngle LookedAngle2Azi = Azi(*LookedAngle2->StationPoint, *LookedAngle2->ObjPoint);

						//计算LookedAngle1所在照准边的方位角
						CAngle LookedAngle1Azi = LookedAngle2Azi + Angle;

						//如果求得LookedAngle1所在照准边的方位角大于等于360,减去360
						if (LookedAngle1Azi(DEG) > 360)
						{
							LookedAngle1Azi(DEG) = LookedAngle1Azi(DEG) - 360;
						}

						double dx, dy;//坐标增量

						//计算坐标增量
						dx = Dist->DistValue * cos(LookedAngle1Azi(RAD));
						dy = Dist->DistValue * sin(LookedAngle1Azi(RAD));

						//计算照准点的近似坐标
						LookedAngle1->ObjPoint->dx = LookedAngle1->StationPoint->dx + dx;
						LookedAngle1->ObjPoint->dy = LookedAngle1->StationPoint->dy + dy;
						LookedAngle1->ObjPoint->bState = TRUE;
					}
				}
			}
		}
	} while (!JudgeUnknownPoint());
}*/


/*
//计算未知点近似坐标
void ControlNetworkAdjust::CalUnknownPoint()
{
	CAngle StartAzi;//零方向边的方位角
	CAngle ObjAzi;//照准边的方位角
	double dx;//坐标x增量
	double dy;//坐标y增量
	ObsDist* Dist;//测站点和照准点之间的边长观测值，指针变量

	//计算未知点的近似坐标
	do
	{
		//遍历观测角度，求未知点近似坐标
		for (int i = 0; i < m_ObsAngleCount; i++)
		{
			//如果观测角度所在边是零方向边，且测站点和照准点的坐标均已求出
			//此时证明可根据该零方向边求出在该测站点上观测的其他未知点的坐标
			if (m_ObsAngle[i].AngleValue(DMS) == 0 && m_ObsAngle[i].StationPoint->bState == TRUE &&
				m_ObsAngle[i].ObjPoint->bState == TRUE)
			{
				//计算零方向边的方位角
				StartAzi = Azi(*m_ObsAngle[i].StationPoint, *m_ObsAngle[i].ObjPoint);

				//计算在该测站点上观测的其他未知点的坐标
				for (int j = i + 1; j < m_ObsAngleCount; j++)
				{
					//如果照准点是在该测站点上观测的其他未知点,且该未知点的坐标未求出
					if (m_ObsAngle[j].AngleValue(DMS) != 0 && m_ObsAngle[j].ObjPoint->bState == FALSE)
					{
						//计算照准边的方位角
						ObjAzi = StartAzi + m_ObsAngle[j].AngleValue;

						//如果求得照准边方位角大于等于360,减去360
						if (ObjAzi(DEG) >= 360)
						{
							ObjAzi(DEG) = ObjAzi(DEG) - 360;
						}

						//寻找测站点和照准点之间的边长观测值
						Dist = SearchObsDist(m_ObsAngle[j].StationPoint->strID, m_ObsAngle[j].ObjPoint->strID);

						//如果存在测站点和照准点之间的边长观测值
						if (Dist != NULL)
						{
							//计算坐标增量
							dx = Dist->DistValue * cos(ObjAzi(RAD));
							dy = Dist->DistValue * sin(ObjAzi(RAD));

							//计算照准点的近似坐标
							m_ObsAngle[j].ObjPoint->dx = m_ObsAngle[j].StationPoint->dx + dx;
							m_ObsAngle[j].ObjPoint->dy = m_ObsAngle[j].StationPoint->dy + dy;
							m_ObsAngle[j].ObjPoint->bState = TRUE;
						}
					}
					else if (m_ObsAngle[j].AngleValue(DMS) == 0)//如果是在其他测站点上观测的未知点
					{
						break;//跳出循环，说明在该测站点上观测的未知点的坐标已全部求出
					}
				}
			}
		}
	} while (!JudgeUnknownPoint());//如果未知点的坐标没有全部求出，继续循环
}
*/


