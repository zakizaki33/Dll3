#if !defined(AFX_CDCON_H__B0202B09_44C1_45DC_99C3_0330B15086FF__INCLUDED_)
#define AFX_CDCON_H__B0202B09_44C1_45DC_99C3_0330B15086FF__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#define _CRT_SECURE_NO_WARNINGS
#include "general_func.h"

class cDcon  
{
	static double QconCoefficient(int m,int i);
	static double Qcon(int m,double x,int me=-1);
	static double Qcon1(int m,double x);
	static double Qcon2(int m,double x);

	enum { M=8 };     //  ヘッダファイルでは int M=8; とはできないのでenumを使う．

public:
	cDcon();
	
	double A[M+1];

	double Rn;	
	void   SetRn(double val);
	
	double& A4();   // =A[0]
	double& A6();   // =A[1]
	double& A8();   // =A[2]
	double& A10();  // =A[3]
	double& A12();  // =A[4]
	double& A14();  // =A[5]
	double& A16();  // =A[6]
	double& A18();  // =A[7]
	double& A20();  // =A[8]

	void Clear();
	int NotNull();
	cDcon reversed();
	void scale(double m);

	double a4();
	double a6();
	double a8();
	double a10();
	double a12();
	double a14();
	double a16();
	double a18();
	double a20();

	void PSeriesToDcon(int n,double *a);

	double Z(double x,double y,int me=-1);
	double Zx(double x,double y);
	double Zy(double x,double y);
	double Zxx(double x,double y);
	double Zxy(double x,double y);
	double Zyy(double x,double y);

	std::string GetTerms();
	void        SetTerms(std::string terms);
};

#endif // !defined(AFX_CDCON_H__B0202B09_44C1_45DC_99C3_0330B15086FF__INCLUDED_)
