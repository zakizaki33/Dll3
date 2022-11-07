#if !defined(AFX_COMPLEX_H__DAF75BE6_6CB6_11D6_A068_0001033BC35C__INCLUDED_)
#define AFX_COMPLEX_H__DAF75BE6_6CB6_11D6_A068_0001033BC35C__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include <string>
#include <iostream>
#include <ostream>
#include <fstream>
#include "string.h"

class complex  
{
	double x,y;
public:
	complex();
	complex(double);
	complex(double,double);
	virtual ~complex();
	friend double& Re(complex&);
	// friend double Re(complex&);
	//friend double Re(const complex&);   //const で変更できないようにする？


	friend double& Im(complex&);
	friend double abs(const complex&);
// #if _MSC_VER < 1800 // Visual C++ 2013 より古いとき
// 	friend double abs(double);  // 【重要】新しいC++ではこの関数を削除すること（詳細は定義を参照）
// #endif // _MSC_VER < 1800
	friend double sqabs(const complex&);
	friend double arg(const complex&);
	friend complex conj(const complex&);       // conjugate number
	friend double  conj(const double&);
	friend complex operator*(const complex&);  // conjugate number
	friend bool operator==(const complex&,const complex&);
	friend bool operator!=(const complex&,const complex&);
	friend complex operator+(const complex&,const complex&);
	friend complex operator-(const complex&,const complex&);
	friend complex operator-(const complex&);
	friend complex operator*(const complex&,const complex&);
	friend complex operator/(const complex&,const complex&);
	complex& operator+=(const complex&);
	complex& operator-=(const complex&);
	complex& operator*=(const complex&);
	complex& operator/=(const complex&);
	friend complex exp(const complex&);
	friend complex sin(const complex&);
	friend complex cos(const complex&);
	friend complex sqrt(const complex&);
	friend std::ostream& operator<<(std::ostream& to,const complex&);
	friend std::istream& operator>>(std::istream& from,complex&);
};

#endif // !defined(AFX_COMPLEX_H__DAF75BE6_6CB6_11D6_A068_0001033BC35C__INCLUDED_)
