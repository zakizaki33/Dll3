#if !defined(AFX_CSPLINE_H__83C9C481_0CCE_4FFF_97B2_E87827E21301__INCLUDED_)
#define AFX_CSPLINE_H__83C9C481_0CCE_4FFF_97B2_E87827E21301__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "matrix.h"
#include "general_func.h"
#include <string>

class cSpline  // 3次スプライン補間
{
	int N;
	
	enum { MAX_SIZE=250 };
	double _X[MAX_SIZE+1],_Y[MAX_SIZE+1];  // 自動設計用にアドレスを固定するため固定長配列とする
	double *a,*b,*c,*d;

	void alloc();
	void erase();
	void calc();

	bool ref_public;  // (自動設計用の)参照を返す関数 double& X(int i), double& Y(int i) のために必要．
	                  // C#のようなプロパティが実装できればよいが，C++では複雑な方法しかなさそう．
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

	int StartPointDerivativeZero, EndPointDerivativeZero;  // 端点で1次導関数を0にするかどうか（偽のときは，2次導関数が0）

	double y(double x);
	double y1(double x);
	double y2(double x);

	std::string GetData();
	void        SetData(std::string data);
	void DoubleN();
	void SetXStep(double dx);
};

#endif // !defined(AFX_CSPLINE_H__83C9C481_0CCE_4FFF_97B2_E87827E21301__INCLUDED_)
