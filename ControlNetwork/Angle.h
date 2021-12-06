#pragma once

//ö���������ͣ����ڴ���Ƕ���ʽ
enum AngleStyle
{
	DEG,//ʮ���ƶ�
	DMS,//�ȷ���
	RAD//������
};

/*******************************************************************
* ������CAngle													   *
*																   *
* �������Ƕ���													   *
*																   *
* ʹ�÷�����1����ʼ����CAngle(double value=0,AngleStyle style=DMS) *
*			��Ĭ��ֵΪ0���Ƕ�Ϊ�ȷ��룬���valueΪa��AngleStyleΪ  *
*			Deg��RAD����������ֵ��Ϊa������Ϊ�Ȼ򻡶ȡ�          *
*			2����ȡָ�����͵ĽǶ�ֵ�ĺ���������					   *
*			��һ��Ϊdouble& operator() (AngleStyle style)��ʹ�ø�  *
*			�������������������������ֵ�����͸ı䣨ָ�����ͣ� *
*			�ڶ���Ϊdouble operator() (AngleStyle style) const��   *
*			const CAngle���ͱ������ã��ú������ص�����������ã� *
*			��������ֵ�����Ͳ��䣨ԭ�����ͣ�					   * 
*			3���ǶȼӼ�����CAngle operator + (const CAngle& m1,	   *
			const CAngle& m2)��CAngle operator - (const CAngle& m1,*
			const CAngle& m2)���ú������ص����������ΪRAD��       *
*																   *
* ��ʷ��**����** **����** **ǩ��**								   *
*       2021/7/6    ��     Է�ղ�								   *
*																   *
* �ⲿ�ࣺ��													   *
*******************************************************************/
class CAngle
{
public:
	CAngle(double value=0,AngleStyle style=DMS);
	~CAngle(void);
private:
	double dValue;//�Ƕ�ֵ
	AngleStyle  nCurStyle;//��ǰ�Ƕ�ֵ����
private:
	//���ó���Ա���������ã�1.���Ա���ᱻ�ı�
	//2.���Ա�����������
	double Deg(double dDms) const;
	double Dms(double dDeg) const;
	
public:
	//��ȡָ�������ͻ�ȡ�Ƕ�ֵ��
    //���ڷ��ص���dValue�����ã����Ը�ֵ��С���Ըı䣬�����Խ��и�ֵ
	double& operator() (AngleStyle style);   //CAanle ang1,   ang1.operator(DEG); angl1(DEG)

	//���أ���ȡָ�������ͻ�ȡ�Ƕ�ֵ����ֵ���ɸı䣬const CAngle���ͱ�������
	double operator() (AngleStyle style) const;
	//���������+/-
    friend CAngle operator + (const CAngle& m1,const CAngle& m2);
	friend CAngle operator - (const CAngle& m1,const CAngle& m2);
};
