#if !defined(AFX_CMODLEGENDRE_H__FCE9369D_DD2C_4C92_B45B_646646555A9A__INCLUDED_)
#define AFX_CMODLEGENDRE_H__FCE9369D_DD2C_4C92_B45B_646646555A9A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include "string.h"
#include "general_func.h"

class cModLegendre
{
	enum {N=10};  // ‘½€®‚ÌÅ‘åŸ”
	double a[N+1][N+1];
public:
	static double P(int n,double x,int power=-1);    // Legendre‘½€® Pn(x)
	static double P1(int n,double x);   // Pn(x)‚Ì1ŠK”÷•ª
	static double P2(int n,double x);   // Pn(x)‚Ì2ŠK”÷•ª

	cModLegendre();
	cModLegendre(const cModLegendre& x);
	cModLegendre& operator=(const cModLegendre& x);
	
	double R0;	
	double& C(int i,int j);	

	void Clear();
	int NotNull();
	cModLegendre reversed();
	void scale(double m);

	double Z(double x,double y,int xpower=-1,int ypower=-1);
	double Zx(double x,double y);
	double Zy(double x,double y);
	double Zxx(double x,double y);
	double Zxy(double x,double y);
	double Zyy(double x,double y);

	std::string GetCoefficients();
	void        SetCoefficients(std::string s);
};

#endif // !defined(AFX_CMODLEGENDRE_H__FCE9369D_DD2C_4C92_B45B_646646555A9A__INCLUDED_)
