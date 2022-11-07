#if !defined(AFX_CLENEQ_H__441FB231_4BB5_11D6_A068_0001033BC35C__INCLUDED_)
#define AFX_CLENEQ_H__441FB231_4BB5_11D6_A068_0001033BC35C__INCLUDED_

#if _MSC_VER > 1000
#pragma once //インクルードを一度だけ行うという意味
#endif // _MSC_VER > 1000
#define _CRT_SECURE_NO_WARNINGS
#include "matrix.h"
#include "general_func.h"

class cLeneq  
{
	int NumberOfEq;
	int NumberOfVar;
	matrix<double> A;
	matrix<double> W;
	matrix<int> C; 
	matrix<double> F;
	matrix<double> X;
	bool calced;
	double calced_rho;

public:
	cLeneq();
	virtual ~cLeneq();
	void SetNumberOfEq(int m);
	int GetNumberOfEq();
	void SetNumberOfVar(int n);
	int GetNumberOfVar();
	void SetA(int i,int j,double value);
	double GetA(int i,int j) const;
	matrix<double> GetAMatrix() const;
	void SetB(int i,double value);
	double GetB(int i);
	void SetWeight(int i,double value);
	double GetWeight(int i);
	void SetConstrKind(int i,int value);
	int  GetConstrKind(int i);
	enum {DLS=0, EQ=1, LT=2, GT=3}; // コンストレインツの種類  DLS:減衰最小二乗法  EQ,LT,GT:一次式拘束(=,<,>)
	double GetDampedX(int j,double rho,int GradientDescent=0);
	double GetX(int j);
	double pvError(double rho);
	double pvError();
};

#endif // !defined(AFX_CLENEQ_H__441FB231_4BB5_11D6_A068_0001033BC35C__INCLUDED_)
