#if !defined(AFX_CSELFOC_H__7059589F_807B_4D63_A86C_D70B6139A85C__INCLUDED_)
#define AFX_CSELFOC_H__7059589F_807B_4D63_A86C_D70B6139A85C__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "vector.h"
#include "general_func.h"

class cSelfoc  
{
	// ã¸ê‹ó¶ N = no * ( 1 - rA*rA*(x*x+y*y)/2 )
	double no;
	double rA;
	double phi;
	double L;
	double Pn;

	double N(vector<double> r);
	vector<double> gradN(vector<double> r);
	vector<double> NgradN(vector<double> r);

public:
	static double N(double no,double rA,double x,double y);

	cSelfoc();
	cSelfoc(double no,double rA,double phi);
	virtual ~cSelfoc();

	double Get_no();
	void Set_no(double value);
	double Get_rA();
	void Set_rA(double value);
	double Get_phi();
	void Set_phi(double value);
	double Get_L();
	void Set_L(double value);
	double Get_Pn();
	void Set_Pn(double value);
	
	int trace(vector<double> &r,vector<double> &Q,double ds);
	int Trace(double &x,double &y,double &Qx,double &Qy,double &Qz,double N,double N1,double ds);
};

#endif // !defined(AFX_CSELFOC_H__7059589F_807B_4D63_A86C_D70B6139A85C__INCLUDED_)
