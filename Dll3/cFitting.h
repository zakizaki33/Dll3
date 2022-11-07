#if !defined(AFX_CFITTING_H__626D0810_9CB3_11D6_8AB6_E2716D88EA50__INCLUDED_)
#define AFX_CFITTING_H__626D0810_9CB3_11D6_8AB6_E2716D88EA50__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "matrix.h" // Client.cpp で呼び出している
#include "general_func.h"
#include "cLeastSquares.h" // Client.cpp で呼び出している
#include "cXYList.h"

class cFitting  : public cLeastSquares
{
	// (x1(i),x2(i),y(i))を
	//     y = Σ C(j)*(x1-x10)^x1dim(j)*(x2-x20)^x2dim(j) 
	// にフィッティングする．
	// Cとx10,x20を調整する．x10,x20は固定することも出来る．
	
	list<double> &y;
	int NumberOfTerms;
	int N();     // =NumberOfTerms+2(x10とx20)
	int *x1dim,*x2dim;
	matrix<double> C;
	double x10,x20;
	int x10Variable,x20Variable;   // x10,x20の固定，非固定
	bool CoefficientsIsNew;
	int *Digits;
	void alloc_terms();
	void free_terms();

	double f(double x1,double x2,double y,double dummy);
	void SetA();
	void SetX();
	int NumberOfCoefficients();

public:
	cFitting();
	virtual ~cFitting();

	void SetNumberOfTerms(int n);
	int GetNumberOfTerms() const;
	void SetOrder(int order);

	void dimensionSet(int j,int x1dimension,int x2dimension);
	int x1dimensionGet(int j) const;
	int x2dimensionGet(int j) const;
	void RemoveX1Odd();
	void RemoveX1Even();
	void RemoveX2Odd();
	void RemoveX2Even();

	void dataSet(int i,double x1_value,double x2_value,double y_value);
	void AddData(double x1_value,double x2_value,double y_value);
	void DataClear();
	double x1dataGet(int i) const;
	double x2dataGet(int i) const;
	double ydataGet(int i) const;
	
	void DigitsSet(int j,int value);
	int CalcCoefficients();
	double coefficientGet(int j,int CoefficientCalc);
	double coefficientGet(int x1_order,int x2_order,int CoefficientCalc);
	void   coefficientSet(int j,double value);
	double x10Get() const;
	void   x10Set(double value);
	double x20Get() const;
	void   x20Set(double value);
	int    x10VariableGet() const;
	void   x10VariableSet(int value);
	int    x20VariableGet() const;
	void   x20VariableSet(int value);
	double yApproximate(double x1value,double x2value,int EliminateConst,int CoefficientCalc);
	double yApproximate(double x1value,double x2value,int CoefficientCalc);
	double dyApproximate(double x1value,double x2value,int CoefficientCalc);
	double Error(int i,int CoefficientCalc);
	double rmsError(int CoefficientCalc);
	double pvError(int CoefficientCalc);
	double maxError(int CoefficientCalc);
	double R2(int CoeffecientCalc);
	void RemoveTerms();
	friend std::ostream& operator<<(std::ostream& to,const cFitting& x);
};

#endif // !defined(AFX_CFITTING_H__626D0810_9CB3_11D6_8AB6_E2716D88EA50__INCLUDED_)
