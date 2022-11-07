#if !defined(AFX_CSPLINE_H__83C9C481_0CCE_4FFF_97B2_E87827E21301__INCLUDED_)
#define AFX_CSPLINE_H__83C9C481_0CCE_4FFF_97B2_E87827E21301__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "matrix.h"
#include "general_func.h"
#include <string>

class cSpline  // 3���X�v���C�����
{
	int N;
	
	enum { MAX_SIZE=250 };
	double _X[MAX_SIZE+1],_Y[MAX_SIZE+1];  // �����݌v�p�ɃA�h���X���Œ肷�邽�ߌŒ蒷�z��Ƃ���
	double *a,*b,*c,*d;

	void alloc();
	void erase();
	void calc();

	bool ref_public;  // (�����݌v�p��)�Q�Ƃ�Ԃ��֐� double& X(int i), double& Y(int i) �̂��߂ɕK�v�D
	                  // C#�̂悤�ȃv���p�e�B�������ł���΂悢���CC++�ł͕��G�ȕ��@�����Ȃ������D
public:
	cSpline();
	cSpline(const cSpline& x);
	virtual ~cSpline();
	cSpline& operator=(const cSpline &x);

	bool NotNull();
	cSpline reversed() const;
	void scale(double m);

	int  GetN() const;
	void SetN(int val);
	double GetX(int i);
	void   SetX(int i,double val);
	double GetY(int i);
	void   SetY(int i,double val);

	double& X(int i);
	double& Y(int i);

	int StartPointDerivativeZero, EndPointDerivativeZero;  // �[�_��1�����֐���0�ɂ��邩�ǂ����i�U�̂Ƃ��́C2�����֐���0�j

	double y(double x);
	double y1(double x);
	double y2(double x);

	std::string GetData();
	void        SetData(std::string data);
	void DoubleN();
	void SetXStep(double dx);
};

#endif // !defined(AFX_CSPLINE_H__83C9C481_0CCE_4FFF_97B2_E87827E21301__INCLUDED_)
