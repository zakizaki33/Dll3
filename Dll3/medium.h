#if !defined(AFX_MEDIUM_H__B9DC54E0_7E5A_11D7_BE64_AFA20D2D3150__INCLUDED_)
#define AFX_MEDIUM_H__B9DC54E0_7E5A_11D7_BE64_AFA20D2D3150__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include <fstream>
#include "string.h"
#include "cGlass.h"

class medium  
{
public:
	double d, delta_d; int dVariable;
	cGlass g; int gVariable;
public:
	static int file_ver;
	medium();
	medium(int cn);
	friend bool operator==(const medium& a,const medium& b);
	medium reversed();
	friend std::ostream& operator<<(std::ostream& to,medium x);
	friend std::istream& operator>>(std::istream& from,medium& x);
	void scale(double m);
};

#endif // !defined(AFX_MEDIUM_H__B9DC54E0_7E5A_11D7_BE64_AFA20D2D3150__INCLUDED_)
