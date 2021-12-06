#include "pch.h"
#include "ControlNetworkAdjust.h"
#include "CommonSurveyFunctions.h"
#include <locale>

//���캯��
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

//�������캯��
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

//��������
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

//��֪�������Ƶ㣬��P1->P2�ķ�λ��
CAngle ControlNetworkAdjust::Azi(const CtrlPoint& P1, const CtrlPoint& P2)
{
	CAngle angAzi;

	angAzi(RAD) = Azimuth(P1.dx, P1.dy, P2.dx, P2.dy);

	return angAzi;
}

//��֪�������Ƶ㣬��P1��P2�ľ���
double ControlNetworkAdjust::DIST(const CtrlPoint& P1, const CtrlPoint& P2)
{
	return Dist(P1.dx, P1.dy, P2.dx, P2.dy);
}

//�ָ��ַ���
int ControlNetworkAdjust::SplitStringArray(CString str, char split, CStringArray& aStr)
{
	int startIdx = 0;
	int idx = str.Find(split, startIdx);
	aStr.RemoveAll();//�����
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



//���ݵ���Ѱ�ҷ��ϵ���֪���δ֪�㣬����CtrlPointָ�����
CtrlPoint* ControlNetworkAdjust::SearchPoint(CString strName)
{
	strName.Trim();
	CString strPointID;//������ĵ���

	//������֪��
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		strPointID = m_KnownPoint[i].strID.Trim();
		if (strName == strPointID)
		{
			return &m_KnownPoint[i];
		}
	}

	//����δ֪��
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		strPointID = m_UnknownPoint[i].strID.Trim();
		if (strPointID == strName)
		{
			return &m_UnknownPoint[i];
		}
	}
}


//�ж�δ֪���Ƿ�ȫ���������ȫ����������棬��֮���ؼ�
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

//Ѱ�Ҳ�վ��Ϊ���ҵ�����׼������������ĽǶȹ۲�ֵ
ObsAngle* ControlNetworkAdjust::SearchObsAngle_StationPtIsLookedPt(CString strName)
{
	strName.Trim();
	CString strStationPtID;//������Ƕȹ۲�ֵ��վ�����

	//�����Ƕȹ۲�ֵ
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

//Ѱ�Ҳ�վ�����׼��ĵ�ŷ��ϵĽǶȹ۲�ֵ
ObsAngle* ControlNetworkAdjust::SearchObsAngle(CString strStationName, CString strObjName)
{
	strStationName.Trim();
	strObjName.Trim();

	CString strStationID;//������Ĳ�վ��ĵ��
	CString strObjID;//���������׼��ĵ��

	//�����Ƕȹ۲�ֵ
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

//Ѱ���յ�Ϊ���ҵ����������������ı߳��۲�ֵ
//���룺strNameΪ���ҵ����
//��������������ı߳��۲�ֵ��ָ�����
ObsDist* ControlNetworkAdjust::SearchObsDist_EndPtIsLookedPt(CString strName)
{
	CString strEndID;//������߳��յ����

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

//Ѱ�����Ϊ���ҵ����յ�����������ı߳��۲�ֵ
//���룺strNameΪ���ҵ����
//��������������ı߳��۲�ֵ��ָ�����
ObsDist* ControlNetworkAdjust::SearchObsDist_StartPtIsLookedPt(CString strName)
{
	CString strStartID;//������߳��յ����

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

//Ѱ����һ�Ƕȹ۲�ֵ�в�վ���㷽������Ӧ�ĽǶȹ۲�ֵ
//���룺��һ�Ƕȹ۲�ֵ�۲�ֵ�������е�λ��iNum
//������ýǶȹ۲�ֵ�в�վ���㷽������Ӧ�ĽǶȹ۲�ֵ������ΪObsAngleָ�����
ObsAngle* ControlNetworkAdjust::SearchZeroDirection(int iNum)
{
	//�������Ƕȹ۲�ֵ
	for (int i = iNum; i >= 0; i--)
	{
		//����Ƕȹ۲�ֵΪ0��˵���ýǶȹ۲�ֵ����Ӧ�ķ���Ϊ�㷽��
		if (m_ObsAngle[i].AngleValue(DMS) == 0)
		{
			return &m_ObsAngle[i];
		}
	}

	return NULL;
}




//�����۲�Ƕ����̺�Ȩ����
//���룺��
//�����B��L��Ϊ����V=BX+L�нǶȺͱ߳��Ͼ���PΪ�ǶȺͱ߳���Ȩ��,����B��L����P�нǶ���ǰ���߳��ں�
void ControlNetworkAdjust::ObsAngleErrorMatrix(CMatrix& B, CMatrix& L, CMatrix& P)
{
	double p = 180 / PI * 3600;//ʹ������תΪ������Եĳ���

	//�����Ƕȹ۲�ֵ�������۲�Ƕ����̺�Ȩ����
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		//������׼����ķ�λ�ǽ���ֵ
		CAngle SightAzi = Azi(*m_ObsAngle[i].StationPoint, *m_ObsAngle[i].ObjPoint);

		//�����㷽��ķ�λ�ǽ���ֵ
		ObsAngle* ZeroDirection = SearchZeroDirection(i);
		CAngle ZeroAzi;
		if (ZeroDirection != NULL)
		{
			ZeroAzi = Azi(*ZeroDirection->StationPoint, *ZeroDirection->ObjPoint);
		}

		//������׼����ı߳�
		double dDist = DIST(*m_ObsAngle[i].StationPoint, *m_ObsAngle[i].ObjPoint);

		//��������B��ϵ��
		double dA = sin(SightAzi(RAD)) / dDist * p / 1000;
		double dB = cos(SightAzi(RAD)) / dDist * (-p) / 1000;

		//���㳣����
		CAngle Angle = SightAzi - ZeroAzi;//����н�
		//���������С��0
		if (Angle(DEG) < 0)
		{
			Angle(DEG) = Angle(DEG) + 360;
		}
		CAngle AngleError = Angle - m_ObsAngle[i].AngleValue;
		double dAngleErrorSec = AngleError(DEG) * 3600;//תΪ��

		//�������̺�Ȩ��
		{
			//����B
			//�����վ����δ֪��
            if (m_ObsAngle[i].StationPoint->ePointStyle == UnknownPoint)
			{
				//�����������������B�е�����
				int iNumx = 2 * m_ObsAngle[i].StationPoint->iNum;//����x��������B�е�����
				int iNumy = iNumx + 1;//����y��������B�е�����

				B(i, iNumx) = dA;
				B(i, iNumy) = dB;
			}
			//�����׼����δ֪��
			if (m_ObsAngle[i].ObjPoint->ePointStyle == UnknownPoint)
			{
				//�����������������B�е�����
				int iNumx = 2 * m_ObsAngle[i].ObjPoint->iNum;//����x��������B�е�����
				int iNumy = iNumx + 1;//����y��������B�е�����

				B(i, iNumx) = -dA;
				B(i, iNumy) = -dB;
			}

			//����L
			L(i, 0) = dAngleErrorSec;

			//����P
			P(i, i) = 1;
		}
	}
}

//�����۲�߳����̺�Ȩ����
//���룺��
//�����B��L��Ϊ����V=BX+L�нǶȺͱ߳��Ͼ���PΪ�ǶȺͱ߳���Ȩ��,����B��L����P�нǶ���ǰ���߳��ں�
void ControlNetworkAdjust::ObsDistErrorMatrix(CMatrix& B, CMatrix& L, CMatrix& P)
{
	//�����Ƕȹ۲�ֵ�������۲�Ƕ����̺�Ȩ����
	for (int i = 0; i < m_ObsDistCount; i++)
	{
		//����߳��ķ�λ�ǽ���ֵ
		CAngle DistAzi = Azi(*m_ObsDist[i].StartPoint, *m_ObsDist[i].EndPoint);

		//����߳�����ֵ
		double dDist = DIST(*m_ObsDist[i].StartPoint, *m_ObsDist[i].EndPoint);

		//���㳣����
		double dDistErrorMm = (dDist - m_ObsDist[i].DistValue) * 1000;//תΪ����

		//�����۲�߳����̺�Ȩ����
		{
			//����B
			//������Ϊδ֪��
			if (m_ObsDist[i].StartPoint->ePointStyle == UnknownPoint)
			{
				//�����������������B�е�����
				int iNumx = 2 * m_ObsDist[i].StartPoint->iNum;//����x��������B�е�����
				int iNumy = iNumx + 1;//����y��������B�е�����

				B(m_ObsAngleCount + i, iNumx) = -cos(DistAzi(RAD));
				B(m_ObsAngleCount + i, iNumy) = -sin(DistAzi(RAD));
			}

			//����յ�Ϊδ֪��
			if (m_ObsDist[i].EndPoint->ePointStyle == UnknownPoint)
			{
				//�����������������B�е�����
				int iNumx = 2 * m_ObsDist[i].EndPoint->iNum;//����x��������B�е�����
				int iNumy = iNumx + 1;//����y��������B�е�����

				B(m_ObsAngleCount + i, iNumx) = cos(DistAzi(RAD));
				B(m_ObsAngleCount + i, iNumy) = sin(DistAzi(RAD));
			}

			//����L
			L(m_ObsAngleCount + i, 0) = dDistErrorMm;

			//����P
			P(m_ObsAngleCount + i, m_ObsAngleCount + i) = 1000 / m_ObsDist[i].DistValue;
		}
	}
}

//����B��L��P
void ControlNetworkAdjust::CalBLP(CMatrix& B, CMatrix& L, CMatrix& P)
{
	//����B�нǶȹ۲�ֵ�㷽���ϵͳ���Z��ϵ��
	int iNumz = 2 * m_UnknownPointCount;//Z��B�е�����

	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		//������㷽��
		if (m_ObsAngle[i].AngleValue(DEG) == 0)
		{
			B(i, iNumz) = -1;

			//�Ѹò�վ�������۲�ֵ��Z��ϵ����ֵ
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

	//�������̺�Ȩ��
	ObsAngleErrorMatrix(B, L, P);
	ObsDistErrorMatrix(B, L, P);
}

//����ǶȺͱ߳��������Լ��������ֵ
void ControlNetworkAdjust::CalAdjustedAngleDist(CMatrix& V)
{
	//���������ĽǶȹ۲�ֵ
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		//����Ƕȸ�����,��λ�룬������λС��
		double dAngleError = V(i, 0);
		CString strAngleError;
		strAngleError.Format(_T("%.2f"), dAngleError);
		m_ObsAngle[i].AngleError = _tstof(strAngleError);

		//���������ĽǶ�ֵ
		//���ڽǶȸ�����һ�㲻����60�룬������ڣ���˵���۲�ֵ��������
		if (m_ObsAngle[i].AngleError < 60)
		{
			m_ObsAngle[i].AngleValue(DMS) = m_ObsAngle[i].AngleValue(DMS) + m_ObsAngle[i].AngleError / 10000;//���������Ƕ�ֵ
		}
		else
		{
			MessageBox(NULL, _T("�Ƕȸ���������60�룬�Ƕȹ۲�ֵ�������⣬���飡"), _T("����"), MB_OK);
		}
	}

	//���������ı߳��۲�ֵ
	for (int i = 0; i < m_ObsDistCount; i++)
	{
		//����߳�����������λm��������λС��
		double dDistError= V(m_ObsAngleCount + i, 0) / 1000;
		CString strDistError;
		strDistError.Format(_T("%.4f"), dDistError);
		m_ObsDist[i].DistError = _tstof(strDistError);

		//���������ı߳�
		m_ObsDist[i].DistValue = m_ObsDist[i].DistValue + m_ObsDist[i].DistError;
	}
}

//����������������λ���������Բ����
void ControlNetworkAdjust::PrecisionEvaluate(CMatrix& Qx, double m0)
{
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		//����δ֪������x,y�ھ���Qx�е�����
		int iNumx = 2 * m_UnknownPoint[i].iNum;
		int iNumy = iNumx + 1;

		double Qxx = Qx(iNumx, iNumx);
		double Qyy = Qx(iNumy, iNumy);
		double Qxy = Qx(iNumx, iNumy);

		//�����λ���
		m_UnknownPoint[i].dMx = m0 * sqrt(Qxx);
		m_UnknownPoint[i].dMy = m0 * sqrt(Qyy);
		m_UnknownPoint[i].dMk = sqrt(m_UnknownPoint[i].dMx * m_UnknownPoint[i].dMx + m_UnknownPoint[i].dMy * m_UnknownPoint[i].dMy);

		//���������Բ����
		double dAlfa;//���᷽λ��
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





//����dXmin,dYmin,dDetX,dDetY
void ControlNetworkAdjust::XYCal()
{
	double dXmax, dYmax;

	dXmin = m_KnownPoint[0].dx;
	dXmax = m_KnownPoint[0].dx;
	dYmin = m_KnownPoint[0].dy;
	dYmax = m_KnownPoint[0].dy;

	//������֪��
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		dXmin = m_KnownPoint[i].dx < dXmin ? m_KnownPoint[i].dx : dXmin;
		dXmax = m_KnownPoint[i].dx > dXmax ? m_KnownPoint[i].dx : dXmax;
		dYmin = m_KnownPoint[i].dy < dYmin ? m_KnownPoint[i].dy : dYmin;
		dYmax = m_KnownPoint[i].dy > dYmax ? m_KnownPoint[i].dy : dYmax;
	}

	//����δ֪��
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

//��������Ӧ���ڵĻ�ͼ����
void ControlNetworkAdjust::SetScale(CRect& rect)
{
	double a = rect.Width();
	double b = rect.Height();
	double ry = double(rect.Width()) * 3 / 4 / dDetY;//�������2/3���������ұ߻�������ͱ�����
	double rx = double(rect.Height()) * 3 / 4 / dDetX;
	dScale = rx < ry ? rx : ry;
}

//���߼�����ת�����ͻ�������
POINT ControlNetworkAdjust::LP2CP(const CtrlPoint& Point, CRect& rect)
{
	POINT P;
	P.y = rect.Height() - (Point.dx - dXmin) * dScale;
	P.x = (Point.dy - dYmin) * dScale;
	return P;
}

//��ָ����Ϊ���Ļ�С����������ʾ��֪��
void ControlNetworkAdjust::TrianglePointDrawing(CDC* pDC, const POINT& CentralPoint)
{
	POINT Vertice[3];//��������������

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

//���ƿ�����
void ControlNetworkAdjust::ControlNetworkDrawing(CDC* pDC, CRect& rect)
{
	XYCal();//����dXmin,dYmin,dDetX,dDetY
	SetScale(rect);//��������Ӧ���ڵĻ�ͼ����

	//��������ԭ�㣬���Կ���������������ƽ��
	double dOrgNetX = rect.Width() / 4;
	double dOrgNetY = -rect.Height() / 5;

	CPen pen(PS_SOLID, 1, RGB(0, 0, 0));
	CPen* pOldPen = pDC->SelectObject(&pen);

	POINT StartPoint;//�����ͼ�ϵ�����
	POINT EndPoint;//�յ���ͼ�ϵ�����

	//����֪���ϻ�������
	StartPoint = LP2CP(m_KnownPoint[0], rect);
	EndPoint = LP2CP(m_KnownPoint[1], rect);
	StartPoint.x = StartPoint.x + dOrgNetX;
	StartPoint.y = StartPoint.y + dOrgNetY;
	EndPoint.x = EndPoint.x + dOrgNetX;
	EndPoint.y = EndPoint.y + dOrgNetY;
	TrianglePointDrawing(pDC, StartPoint);
	TrianglePointDrawing(pDC, EndPoint);

	//������֪��֮���б��������һ��ƽ����
	int dx = StartPoint.x - EndPoint.x;
	int dy = StartPoint.y - EndPoint.y;
	double dAngle;//��ֱ����֪��֮��ֱ�ߵ�ֱ����X��ļн�
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

	//�����Ƕȹ۲�ֵ��������׼��
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
	

	//�������Բ
	
	//�����Բ������
	double dOrgX;
	double dOrgY;

	//����δ֪�㣬�������Բ
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		double dStartX, dStartY, dEndX, dEndY;

		//�����Բ����������
		POINT Org;
		Org = LP2CP(m_UnknownPoint[i], rect);//����ת��
		dOrgX = Org.x + dOrgNetX;
		dOrgY = Org.y + dOrgNetY;


		double dF = m_UnknownPoint[i].dF * 5000;
		double dE = m_UnknownPoint[i].dE * 5000;
		double dAlfa = m_UnknownPoint[i].dAlfa(RAD);

		//���ƶ̰���
		dStartX = (dF * sin(dAlfa)) * dScale + dOrgX;
		dStartY = (-dF * cos(dAlfa)) * dScale + dOrgY;
		dEndX = (-dF * sin(dAlfa)) * dScale + dOrgX;
		dEndY = (dF * cos(dAlfa)) * dScale + dOrgY;
		pDC->MoveTo(dStartX, dStartY);
		pDC->LineTo(dEndX, dEndY);

		//���Ƴ�����
		dStartX = (-dE * cos(dAlfa)) * dScale + dOrgX;
		dStartY = (-dE * sin(dAlfa)) * dScale + dOrgY;
		dEndX = (dE * cos(dAlfa)) * dScale + dOrgX;
		dEndY = (dE * sin(dAlfa)) * dScale + dOrgY;
		pDC->MoveTo(dStartX, dStartY);
		pDC->LineTo(dEndX, dEndY);

		double ex, fy;
		ex = dE;
		fy = 0;
		//ת���������᷽����
		dStartX = (ex * cos(dAlfa)
			- fy * sin(dAlfa)
			) * dScale + dOrgX;
		dStartY = (fy * cos(dAlfa)
			+ ex * sin(dAlfa)
			) * dScale + dOrgY;
		pDC->MoveTo(dStartX, dStartY);

		for (int i = 6; i <= 360; i += 6)
		{
			//�������᷽�������
			ex = dE * cos((i / 180.0) * PI);
			fy = dF * sin((i / 180.0) * PI);

			//ת���������᷽����
			dEndX = (ex * cos(dAlfa)
				- fy * sin(dAlfa)
				) * dScale + dOrgX;
			dEndY = (fy * cos(dAlfa)
				+ ex * sin(dAlfa)
				) * dScale + dOrgY;
			pDC->LineTo(dEndX, dEndY);
		}
	}
	



	//��������
	double dLenX = double(rect.Height()) / 10;//������X��ĳ���
	double dLenY = double(rect.Width())  / 10;//������Y��ĳ���

	//�������ԭ������
	double dOrgCoordX =  rect.left + double(rect.Width()) / 10;
	double dOrgCoordY = rect.bottom - double(rect.Height()) / 10;

	//��Y��
	pDC->MoveTo(dOrgCoordX, dOrgCoordY);
	pDC->LineTo(dOrgCoordX + dLenY, dOrgCoordY);
	pDC->LineTo(dOrgCoordX + dLenY - 5, dOrgCoordY - 5);
	pDC->MoveTo(dOrgCoordX + dLenY, dOrgCoordY);
	pDC->LineTo(dOrgCoordX + dLenY - 5, dOrgCoordY + 5);

	//��X��
	pDC->MoveTo(dOrgCoordX, dOrgCoordY);
	pDC->LineTo(dOrgCoordX, dOrgCoordY - dLenX);
	pDC->LineTo(dOrgCoordX - 5, dOrgCoordY - dLenX + 5);
	pDC->MoveTo(dOrgCoordX, dOrgCoordY - dLenX);
	pDC->LineTo(dOrgCoordX + 5, dOrgCoordY - dLenX + 5);



	//��������
	//�����ߵ����
	double dOrgScaleX = dOrgCoordX + dLenX + double(rect.Width()) / 10;
	double dOrgScaleY = dOrgCoordY;
	pDC->MoveTo(dOrgScaleX, dOrgScaleY);
	pDC->LineTo(dOrgScaleX + 200 * dScale, dOrgScaleY);
	pDC->LineTo(dOrgScaleX + 200 * dScale, dOrgScaleY - 2);
	pDC->MoveTo(dOrgScaleX, dOrgScaleY);
	pDC->LineTo(dOrgScaleX, dOrgScaleY - 2);


	pDC->SelectObject(pOldPen);
	pen.DeleteObject();




	CFont Font;//��������
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
		_T("����")));                 // lpszFacename

	//pDC->SetTextColor(RGB(0, 0, 0));///�ı�������ɫ

	CFont* pOldFont = pDC->SelectObject(&Font);

	//д����
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

	//д��������
	pDC->TextOutW(dOrgCoordX - 12, dOrgCoordY - dLenX - 2, _T("X"));
	pDC->TextOutW(dOrgCoordX + dLenY - 3, dOrgCoordY + 5, _T("Y"));

	pDC->SelectObject(pOldPen);
	Font.DeleteObject();


	//д��������
	pDC->TextOutW(dOrgScaleX + 70 * dScale, dOrgScaleY + 6, _T("200m"));
}




//���ļ��ж�ȡ����
void ControlNetworkAdjust::ReadFileData(CStdioFile& sf)
{
	CString strLine;
	CStringArray aStrTmp;

	//��ȡ��֪��
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

	//��ȡδ֪��
	sf.ReadString(strLine);
	m_UnknownPointCount = _ttoi(strLine);
	m_UnknownPoint = new CtrlPoint[m_UnknownPointCount];
	sf.ReadString(strLine);
	SplitStringArray(strLine, ',', aStrTmp);
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		m_UnknownPoint[i].strID = aStrTmp[i];
		m_UnknownPoint[i].dx = -1;
		m_UnknownPoint[i].dy = -1;//����x��y��ֵΪ-1����������ж��Ƿ��Ѿ����и���
		m_UnknownPoint[i].bState = FALSE;
		m_UnknownPoint[i].iNum = i;
		m_UnknownPoint[i].ePointStyle = UnknownPoint;
		m_UnknownPoint[i].dE = -1;//���ḳֵΪ-1����������ж��Ƿ��Ѿ�ƽ�����
	}

	//��ȡ�߳��۲�ֵ
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

	//��ȡ�Ƕȹ۲�ֵ
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

//����δ֪���������
bool ControlNetworkAdjust::CalUnknownPoint()
{
	if (m_UnknownPointCount == 0)
	{
		MessageBox(NULL, _T("δ��ȡ���ݣ�"), _T("����"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dx >= 0)
	{
		MessageBox(NULL, _T("�Ѿ����������㣡"), _T("����"), MB_OK);
		return FALSE;
	}
	else
	{
		do
		{
			//����δ֪�㣬����δ֪��Ľ�������
			for (int i = 0; i < m_UnknownPointCount; i++)
			{
				//�����δ֪�������δ���
				if (m_UnknownPoint[i].bState == FALSE)
				{
					ObsAngle* LookedAngle2 = NULL;//��վ�����׼����������ĽǶȹ۲�ֵ
					ObsAngle* LookedAngle1 = NULL;//��վ���������������׼��Ϊ����δ֪��ĽǶȹ۲�ֵ

					//Ѱ���յ�Ϊδ֪�����������������ı߳��۲�ֵ
					ObsDist* LookedDist = SearchObsDist_EndPtIsLookedPt(m_UnknownPoint[i].strID);
					
					//����߳��۲�ֵ����
					if (LookedDist != NULL)
					{
						//Ѱ�Ҳ�վ��ΪDist�������׼������������ĽǶȹ۲�ֵ
						LookedAngle2 = SearchObsAngle_StationPtIsLookedPt(LookedDist->StartPoint->strID);

						//Ѱ����׼��ΪDist���ڱߵĽǶȹ۲�ֵ
						LookedAngle1 = SearchObsAngle(LookedDist->StartPoint->strID, LookedDist->EndPoint->strID);
					}
					else
					{
						//Ѱ�����Ϊδ֪�����յ�����������ı߳��۲�ֵ
						LookedDist = SearchObsDist_StartPtIsLookedPt(m_UnknownPoint[i].strID);

						if (LookedDist != NULL)
						{
							//Ѱ�Ҳ�վ��ΪDist�������׼������������ĽǶȹ۲�ֵ
							LookedAngle2 = SearchObsAngle_StationPtIsLookedPt(LookedDist->EndPoint->strID);

							//Ѱ����׼��ΪDist���ڱߵĽǶȹ۲�ֵ
							LookedAngle1 = SearchObsAngle(LookedDist->EndPoint->strID, LookedDist->StartPoint->strID);
						}
					}
						
					//��������Ƕȹ۲�ֵ���ڣ�֤����δ֪����������
					if (LookedAngle1 != NULL && LookedAngle2 != NULL)
					{
						//����LookedAngle1��LookedAngle2֮��ļн�
						CAngle Angle = LookedAngle1->AngleValue - LookedAngle2->AngleValue;

						//����LookedAngle2������׼�ߵķ�λ��
						CAngle LookedAngle2Azi = Azi(*LookedAngle2->StationPoint, *LookedAngle2->ObjPoint);

						//����LookedAngle1������׼�ߵķ�λ��
						CAngle LookedAngle1Azi = LookedAngle2Azi + Angle;

						//������LookedAngle1������׼�ߵķ�λ�Ǵ��ڵ���360,��ȥ360
						if (LookedAngle1Azi(DEG) > 360)
						{
							LookedAngle1Azi(DEG) = LookedAngle1Azi(DEG) - 360;
						}

						double dx, dy;//��������

						//������������
						dx = LookedDist->DistValue * cos(LookedAngle1Azi(RAD));
						dy = LookedDist->DistValue * sin(LookedAngle1Azi(RAD));

						//������׼��Ľ�������
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

//��������ļ�
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

//���ı����������������
void ControlNetworkAdjust::PrintRoughCalResult(CString& strOut)
{
	strOut.Empty();

	CString strResult;

	strOut += (_T("*****************�������*****************\r\n\r\n"));

	strResult.Format(_T("%-20s%-22s%-22s\r\n\r\n"), _T("���"), _T("X����"), _T("Y����"));
	strOut += strResult;

	//�����֪������
	strOut += (_T("*****************��֪��*******************\r\n"));
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-20.6f%-20.6f\r\n"), m_KnownPoint[i].strID, m_KnownPoint[i].dx, m_KnownPoint[i].dy);
		strOut += strResult;
	}

	//���δ֪������
	strOut += (_T("\r\n*****************δ֪��*******************\r\n"));
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-20.6f%-20.6f\r\n"), m_UnknownPoint[i].strID, m_UnknownPoint[i].dx, m_UnknownPoint[i].dy);
		strOut += strResult;
	}

	return;
}

//���ı��������ƽ�����꣬�������ļ���
void ControlNetworkAdjust::PrintAdjustCalResult(CString& strOut)
{
	CString strFileData;
	CString strResult;

	strFileData += (_T("******************************ƽ�����************************************\n"));

	//�������۲�ɹ�
	strFileData += (_T("\n******************************����۲�ɹ�********************************\n"));
	strResult.Format(_T("%-13s%-13s%-13s%-10s%-17s%-10s\n"), _T("��վ"), _T("��׼"), _T("����ֵ(dms)"), _T("������(s)"),
		_T("ƽ���ֵ(dms)"), _T("��ע"));
	strFileData += strResult;
	for (int i = 0; i < m_ObsAngleCount; i++)
	{
		//������㷽��
		if (m_ObsAngle[i].AngleValue(DEG) == 0)
		{
			strResult.Format(_T("%-15s%-15s%-18.6f\n"), m_ObsAngle[i].StationPoint->strID,
				m_ObsAngle[i].ObjPoint->strID, m_ObsAngle[i].AngleValue(DMS));
			strFileData += strResult;
		}
		else//��������㷽��
		{
			strResult.Format(_T("%-15s%-15s%-18.6f%-15.2f%-18.6f\n"), m_ObsAngle[i].StationPoint->strID, m_ObsAngle[i].ObjPoint->strID,
				m_ObsAngle[i].AngleValue(DMS) - m_ObsAngle[i].AngleError / 10000, m_ObsAngle[i].AngleError, m_ObsAngle[i].AngleValue(DMS));
			strFileData += strResult;
		}
	}

	//�������۲�ɹ�
	strFileData += (_T("\n******************************����۲�ɹ�********************************\n"));
	strResult.Format(_T("%-13s%-15s%-14s%-10s%-13s%-20s\n"), _T("��վ"), _T("��׼"), _T("����(m)"), _T("������(m)"), _T("ƽ���ֵ(m)"),
		_T("��λ��(dms)"));
	strFileData += strResult;
	for (int i = 0; i < m_ObsDistCount; i++)
	{
		CAngle DistAzi = Azi(*m_ObsDist[i].StartPoint, *m_ObsDist[i].EndPoint);//���㷽λ��
		strResult.Format(_T("%-15s%-15s%-18.4f%-15.4f%-18.4f%-18.6f\n"), m_ObsDist[i].StartPoint->strID, m_ObsDist[i].EndPoint->strID,
			m_ObsDist[i].DistValue - m_ObsDist[i].DistError, m_ObsDist[i].DistError, m_ObsDist[i].DistValue, DistAzi(DMS));
		strFileData += strResult;
	}

	//���ƽ���λ���
	strFileData += (_T("\n******************************ƽ���λ���********************************\n"));
	strResult.Format(_T("%-13s%-12s%-10s%-15s%-16s%-20s\n"), _T("����"), _T("����(m)"), _T("����(m)"), _T("���᷽λ��(dms)"),
		_T("��λ�����(m)"), _T("��ע"));
	strFileData += strResult;
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-15.4f%-15.4f%-23.6f%-18.4f\n"), m_UnknownPoint[i].strID, m_UnknownPoint[i].dE,
			m_UnknownPoint[i].dF, m_UnknownPoint[i].dAlfa(DMS), m_UnknownPoint[i].dMk);
		strFileData += strResult;
	}

	//������Ƶ�ɹ�
	strFileData += (_T("\n******************************���Ƶ�ɹ�*********************************\n"));
	strResult.Format(_T("%-19s%-20s%-20s%-13s%-20s\n"), _T("����"), _T("X(m)"), _T("Y(m)"), _T("H(m)"), _T("��ע"));
	strFileData += strResult;
	for (int i = 0; i < m_KnownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-18.4f%-18.4f%-17s%-20s\n"), m_KnownPoint[i].strID, m_KnownPoint[i].dx, m_KnownPoint[i].dy, _T(" "), _T("δ֪��"));
		strFileData += strResult;
	}
	for (int i = 0; i < m_UnknownPointCount; i++)
	{
		strResult.Format(_T("%-15s%-18.4f%-18.4f%-17s%-20s\n"), m_UnknownPoint[i].strID, m_UnknownPoint[i].dx, m_UnknownPoint[i].dy, _T(" "), _T(" "));
		strFileData += strResult;
	}

	//������ļ���
	CString strFileName = _T("AdjustReult.dat");
	CStdioFile sf;
	if (!sf.Open(strFileName, CFile::modeReadWrite | CFile::modeCreate)) return;

	sf.WriteString(strFileData);
	MessageBox(NULL, _T("ƽ�����ѱ��浽��AdjustReult.dat���ļ��У�"), _T("��ʾ"), MB_OK);

	sf.SeekToBegin();//�ļ�ָ�붨λ����ͷ

	ReadFileContent(sf,strOut);//���ļ���ȡ����

	sf.Close();
}


//�ж��Ƿ�ﵽ��ͼ����
bool ControlNetworkAdjust::JudgeDrawConditon()
{
	if (m_UnknownPointCount == 0)
	{
		MessageBox(NULL, _T("δ��ȡ���ݣ�"), _T("����"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dx < 0)
	{
		MessageBox(NULL, _T("δ������㣡"), _T("����"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dE < 0)
	{
		MessageBox(NULL, _T("δƽ����㣡"), _T("����"), MB_OK);
		return FALSE;
	}
	else
	{
		return TRUE;
	}
}

//ƽ�����
bool ControlNetworkAdjust::AdjustCal()
{
	if (m_UnknownPointCount == 0)
	{
		MessageBox(NULL, _T("δ��ȡ���ݣ�"), _T("����"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dx < 0)
	{
		MessageBox(NULL, _T("δ������㣡"), _T("����"), MB_OK);
		return FALSE;
	}
	else if (m_UnknownPoint[0].dE >= 0)
	{
		MessageBox(NULL, _T("�Ѿ����ƽ����㣡"), _T("����"), MB_OK);
		return FALSE;
	}
	else
	{
		//�������̺�Ȩ��
		CMatrix X(m_UnknownPointCount * 3 + m_KnownPointCount, 1);
		CMatrix B(m_ObsAngleCount + m_ObsDistCount, m_UnknownPointCount * 3 + m_KnownPointCount);
		CMatrix L(m_ObsAngleCount + m_ObsDistCount, 1);
		CMatrix P(m_ObsAngleCount + m_ObsDistCount, m_ObsAngleCount + m_ObsDistCount);
		CMatrix V(m_ObsAngleCount + m_ObsDistCount, 1);


		CMatrix VT;
		CMatrix BT;//B��ת�þ���
		CMatrix NBB;
		CMatrix NBB1;//NBB�������

		//ƽ�����ֱ�����㾫��Ҫ��
		double dXmax;//X������x��y���������ֵ
		do
		{
			//����B��L��P
			CalBLP(B, L, P);

			BT = ~B;//��B��ת�þ���
			NBB = BT * P * B;//��NBB
			NBB1 = NBB.Inv();//��NBB�������
			X = -1 * NBB1 * BT * P * L;//��X
			V = B * X + L;//���������V

			//��δ֪������Ľ���ƽ��ֵ
			for (int i = 0; i < m_UnknownPointCount; i++)
			{
				m_UnknownPoint[i].dx += X(2 * i, 0) / 1000;
				m_UnknownPoint[i].dy += X(2 * i + 1, 0) / 1000;
			}

			//��dXmax
			dXmax = fabs(X(0, 0));
			for (int i = 0; i < 2 * m_UnknownPointCount; i++)
			{
				dXmax = fabs(X(i, 0)) > dXmax ? fabs(X(i, 0)) : dXmax;
			}
		} while (dXmax > 0.1);

		//���㷽��ĽǶȸ������㵽�ò�վ������׼����
		double V0;//��¼�㷽��ĽǶȸ�����
		for (int i = 0; i < m_ObsAngleCount; i++)
		{
			//������㷽��
			if (m_ObsAngle[i].AngleValue(DEG) == 0)
			{
				V0 = V(i, 0);
				V(i, 0) = 0;//�㷽���������Ϊ0

				//����Ǹò�վ��������׼����
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

		//����߳��ͽǶȹ۲�ֵ�ĸ������Լ��������ֵ
		CalAdjustedAngleDist(V);

		//���㵥λȨ�����
		VT = ~V;
		CMatrix VV = VT * P * V;
		double m0 = sqrt(VV(0, 0) / double(V.Row() - X.Row())) / 1000;

		//���������������λ���������Բ����
		PrecisionEvaluate(NBB1, m0);

		return TRUE;
	}
}



/*
//����δ֪��Ľ�������
void ControlNetworkAdjust::CalUnknownPoint()
{
	do
	{
		//����δ֪�㣬����δ֪��Ľ�������
		for (int i = 0; i < m_UnknownPointCount; i++)
		{
			//���δ֪�������δ���
			if (m_UnknownPoint[i].bState == FALSE)
			{
				//Ѱ��δ֪��Ϊ��׼���Ҳ�վ������������ĽǶȹ۲�ֵ
				ObsAngle* LookedAngle1 = SearchObsAngle_ObjPtIsLookedPt(m_UnknownPoint[i].strID);

				//Ѱ�Ҳ�վ���LookedAngle1����ͬ����׼������������ĽǶȹ۲�ֵ
				ObsAngle* LookedAngle2 = SearchObsAngle_StationPtIsLookedPt(LookedAngle1->StationPoint->strID);

				//���LookedAngle1��LookedAngle2����
				if (LookedAngle1 != NULL && LookedAngle2 != NULL)
				{
					//Ѱ��LookedAngle1��׼�����ڱ߱߳��۲�ֵ
					ObsDist* Dist = SearchObsDist(LookedAngle1->StationPoint->strID, LookedAngle1->ObjPoint->strID);

					//���LookedAngle1��׼�ߴ��ڣ�֤����δ֪�������������
					if (Dist != NULL)
					{
						//����LookedAngle1��LookedAngle2֮��ļн�
						CAngle Angle = LookedAngle1->AngleValue - LookedAngle2->AngleValue;

						//����LookedAngle2������׼�ߵķ�λ��
						CAngle LookedAngle2Azi = Azi(*LookedAngle2->StationPoint, *LookedAngle2->ObjPoint);

						//����LookedAngle1������׼�ߵķ�λ��
						CAngle LookedAngle1Azi = LookedAngle2Azi + Angle;

						//������LookedAngle1������׼�ߵķ�λ�Ǵ��ڵ���360,��ȥ360
						if (LookedAngle1Azi(DEG) > 360)
						{
							LookedAngle1Azi(DEG) = LookedAngle1Azi(DEG) - 360;
						}

						double dx, dy;//��������

						//������������
						dx = Dist->DistValue * cos(LookedAngle1Azi(RAD));
						dy = Dist->DistValue * sin(LookedAngle1Azi(RAD));

						//������׼��Ľ�������
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
//����δ֪���������
void ControlNetworkAdjust::CalUnknownPoint()
{
	CAngle StartAzi;//�㷽��ߵķ�λ��
	CAngle ObjAzi;//��׼�ߵķ�λ��
	double dx;//����x����
	double dy;//����y����
	ObsDist* Dist;//��վ�����׼��֮��ı߳��۲�ֵ��ָ�����

	//����δ֪��Ľ�������
	do
	{
		//�����۲�Ƕȣ���δ֪���������
		for (int i = 0; i < m_ObsAngleCount; i++)
		{
			//����۲�Ƕ����ڱ����㷽��ߣ��Ҳ�վ�����׼�������������
			//��ʱ֤���ɸ��ݸ��㷽�������ڸò�վ���Ϲ۲������δ֪�������
			if (m_ObsAngle[i].AngleValue(DMS) == 0 && m_ObsAngle[i].StationPoint->bState == TRUE &&
				m_ObsAngle[i].ObjPoint->bState == TRUE)
			{
				//�����㷽��ߵķ�λ��
				StartAzi = Azi(*m_ObsAngle[i].StationPoint, *m_ObsAngle[i].ObjPoint);

				//�����ڸò�վ���Ϲ۲������δ֪�������
				for (int j = i + 1; j < m_ObsAngleCount; j++)
				{
					//�����׼�����ڸò�վ���Ϲ۲������δ֪��,�Ҹ�δ֪�������δ���
					if (m_ObsAngle[j].AngleValue(DMS) != 0 && m_ObsAngle[j].ObjPoint->bState == FALSE)
					{
						//������׼�ߵķ�λ��
						ObjAzi = StartAzi + m_ObsAngle[j].AngleValue;

						//��������׼�߷�λ�Ǵ��ڵ���360,��ȥ360
						if (ObjAzi(DEG) >= 360)
						{
							ObjAzi(DEG) = ObjAzi(DEG) - 360;
						}

						//Ѱ�Ҳ�վ�����׼��֮��ı߳��۲�ֵ
						Dist = SearchObsDist(m_ObsAngle[j].StationPoint->strID, m_ObsAngle[j].ObjPoint->strID);

						//������ڲ�վ�����׼��֮��ı߳��۲�ֵ
						if (Dist != NULL)
						{
							//������������
							dx = Dist->DistValue * cos(ObjAzi(RAD));
							dy = Dist->DistValue * sin(ObjAzi(RAD));

							//������׼��Ľ�������
							m_ObsAngle[j].ObjPoint->dx = m_ObsAngle[j].StationPoint->dx + dx;
							m_ObsAngle[j].ObjPoint->dy = m_ObsAngle[j].StationPoint->dy + dy;
							m_ObsAngle[j].ObjPoint->bState = TRUE;
						}
					}
					else if (m_ObsAngle[j].AngleValue(DMS) == 0)//�������������վ���Ϲ۲��δ֪��
					{
						break;//����ѭ����˵���ڸò�վ���Ϲ۲��δ֪���������ȫ�����
					}
				}
			}
		}
	} while (!JudgeUnknownPoint());//���δ֪�������û��ȫ�����������ѭ��
}
*/


