#pragma once
#include "ValueName.h"
#include "Matrix.h"

const double EPSILON = 1.0E-12;
const double PI = 4.0 * atan(1.0);


/*******************************************************************
* ������ControlNetworkAdjust								       *
*																   *
* ������������ƽ����											   *
*																   *
* ʹ�÷�����													   *
*																   * 
* ��ʷ��**����** **����** **ǩ��**								   *
*       2021/7/6    ��     Է�ղ�								   *
*																   *
* �ⲿ�ࣺ��													   *
*******************************************************************/
class ControlNetworkAdjust
{
private:
	int m_KnownPointCount;//��֪����
	CtrlPoint* m_KnownPoint;//��֪��ָ��
	int m_UnknownPointCount;//δ֪����
	CtrlPoint* m_UnknownPoint;//δ֪��ָ��
	int m_ObsAngleCount;//�۲�Ƕ���
	ObsAngle* m_ObsAngle;//�۲�Ƕ�ָ��
	int m_ObsDistCount;//�۲�߳���
	ObsDist* m_ObsDist;//�۲�߳�ָ��

	//�����ǻ�ͼ�õ��ı���
	double dXmin;//����x����Сֵ
	double dYmin;//����y����Сֵ
	double dDetX;//����x�����ֵ����Сֵ�Ĳ�ֵ
	double dDetY;//����y�����ֵ����Сֵ�Ĳ�ֵ
	double dScale;//��ͼ����
public:
	ControlNetworkAdjust();//���캯��
	ControlNetworkAdjust(ControlNetworkAdjust& A);//�������캯��
	~ControlNetworkAdjust();//��������
private:
	friend void ReadFileContent(CStdioFile& sf, CString& strFileContent);//���ļ��ж�ȡ���ݵ�CString�ַ�����

	CAngle Azi(const CtrlPoint& P1, const CtrlPoint& P2);//��֪�������Ƶ㣬��P1->P2�ķ�λ��
	double DIST(const CtrlPoint& P1, const CtrlPoint& P2);//��֪�������Ƶ㣬��P1��P2�ľ���
	int SplitStringArray(CString str, char split, CStringArray& aStr);// �ָ��ַ������õ��ַ�������


	CtrlPoint* SearchPoint(CString strName);//���ݵ���Ѱ�ҷ��ϵ���֪���δ֪�㣬����CtrlPointָ�����
	bool JudgeUnknownPoint();//�ж�δ֪���Ƿ�ȫ���������ȫ����������棬��֮���ؼ�   
	ObsAngle* SearchObsAngle_StationPtIsLookedPt(CString strName);//Ѱ�Ҳ�վ��Ϊ���ҵ�����׼������������ĽǶȹ۲�ֵ
	ObsAngle* SearchObsAngle(CString strStationName, CString strObjName);//Ѱ�Ҳ�վ�����׼��ĵ�ŷ��ϵĽǶȹ۲�ֵ
	ObsDist* SearchObsDist_EndPtIsLookedPt(CString strName);//Ѱ���յ�Ϊ���ҵ����������������ı߳��۲�ֵ
	ObsDist* SearchObsDist_StartPtIsLookedPt(CString strName);//Ѱ�����Ϊ���ҵ����յ�����������ı߳��۲�ֵ
	ObsAngle* SearchZeroDirection(int iNum);//Ѱ����һ�Ƕȹ۲�ֵ�в�վ���㷽������Ӧ�ĽǶȹ۲�ֵ


	void ObsAngleErrorMatrix(CMatrix& B, CMatrix& L, CMatrix& P);//�����۲�Ƕ����̺�Ȩ����
	void ObsDistErrorMatrix(CMatrix& B, CMatrix& L, CMatrix& P);//�����۲�߳����̺�Ȩ����
	void CalBLP(CMatrix& B, CMatrix& L, CMatrix& P);//����B��L��P
	void CalAdjustedAngleDist(CMatrix& V);//����ǶȺͱ߳��������Լ�������߳��ͽǶ�
	void PrecisionEvaluate(CMatrix& Qx,double m0);//����������������λ���������Բ����


	void XYCal();//����dXmin,dYmin,dDetX,dDetY
	void SetScale(CRect& rect);//��������Ӧ���ڵĻ�ͼ����
	POINT LP2CP(const CtrlPoint& Point, CRect& rect);//���߼�����ת�����ͻ�������
	void TrianglePointDrawing(CDC* pDC, const POINT& CentralPoint);//��ָ����Ϊ���Ļ�С����������ʾ��֪��


	void PrintCMatrix(CMatrix& A, CString strFileName);//�������
public:
	void ReadFileData(CStdioFile& sf);//���ļ��ж�ȡ����
	bool CalUnknownPoint();//����δ֪���������
	void PrintRoughCalResult(CString& strOut);//���ı����������������
	bool AdjustCal();//ƽ�����
	void PrintAdjustCalResult(CString& strOut);//���ı��������ƽ�����꣬�������ļ���

	bool JudgeDrawConditon();//�ж��Ƿ�ﵽ��ͼ����
	void ControlNetworkDrawing(CDC* pDC, CRect& rect);//����������
};
