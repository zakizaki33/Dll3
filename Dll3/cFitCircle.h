#if !defined(AFX_CFITCIRCLE_H__3CAE4EBA_DB84_4480_BF72_F684D6151623__INCLUDED_)
#define AFX_CFITCIRCLE_H__3CAE4EBA_DB84_4480_BF72_F684D6151623__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#define _CRT_SECURE_NO_WARNINGS
#include "cLeastSquares.h" // Client.cpp Ç≈åƒÇ—èoÇµÇƒÇ¢ÇÈ
#include "general_func.h"  // Client.cpp Ç≈åƒÇ—èoÇµÇƒÇ¢ÇÈ

class cFitCircle : public cLeastSquares
{
	list<double> &x,&y,&z;
	double r,x0,y0;
	double a,b,c;
	double f(double x,double y,double dummy1,double dummy2);
	void SetA();
	void SetX();
	int NumberOfCoefficients();
public:
	cFitCircle();
	double GetR();  void SetR(double value);
	double GetX0(); void SetX0(double value);
	double GetY0(); void SetY0(double value);
	double GetX(int i); void SetX(int i,double value);
	double GetY(int i); void SetY(int i,double value);
	double yApproximate(double x,int sgn);

	int ModifyCoefficients();
};

#endif // !defined(AFX_CFITCIRCLE_H__3CAE4EBA_DB84_4480_BF72_F684D6151623__INCLUDED_)
