#if !defined(AFX_CMATERIAL_H__2AC07BFA_7B81_4EFB_9859_CD5472A83582__INCLUDED_)
#define AFX_CMATERIAL_H__2AC07BFA_7B81_4EFB_9859_CD5472A83582__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include "complex.h"
#include "general_func.h"
#include "cXYList.h"

class cMaterial
{
public:
	static std::string GlassDataFile1;
	static std::string GlassDataFile2;
	static std::string ColorGlassDataFolder;

private:
	std::string name;    // “ª‚Ì"-"‚ðŠÜ‚Þ
	std::string name0;   // “ª‚Ì"-"‚ðœ‚¢‚½‚à‚Ì
	double sgn;
	double dN;
	int type; enum{NUL,AIR,CODE,FILE,N_NU,COMPLEX,REAL};
	int eq_type;
	int glass_code;
	double nd,nud;
	double a0,a1,a2,a3,a4,a5,a6,a7,a8;
	double b1,b2,b3,b4,b5,b6;
	cXYList n,k;
	cXYList lrA;
	int Open(std::string name0);

public:
	static std::string MaterialName(std::string gname);
	static double HerzEq(double Nd, double nud, double wl_nm,int digits=5);
	static double HerzEqDerivative(double Nd, double nud, double wl_nm);
	static double HerzEqDerivative2(double Nd, double nud, double wl_nm);
	static double HerzEqNdDerivative(double Nd, double nud, double wl_nm);
	static double HerzEqNUdDerivative(double Nd, double nud, double wl_nm);
	static double Nd(int glasscode);
	static double Nud(int glasscode);
	static double HerzEq(int glasscode, double wl_nm,int digits=5);
	static double HerzEqDerivative(int glasscode, double wl_nm);
	static double HerzEqDerivative2(int glasscode, double wl_nm);
	static complex Index(std::string name,double wl_nm,int digits=5);
	static double  GroupIndex(std::string name,double wl_nm);
	static complex IndexDerivative2(std::string name,double wl_nm);

	static double ColorGlassInnerTransmittance(std::string filename,double wl,double thickness);

	cMaterial();
	cMaterial(std::string name);
	virtual ~cMaterial();
	void SetName(std::string name);
	friend std::istream& operator>>(std::istream& from,cMaterial& x);
	complex Index(double wl_nm,int digits=5);
	complex IndexDerivative(double wl_nm);
	complex IndexDerivative2(double wl_nm);

	int IsGrin();

	double rA(double wl_nm);
	double GrinPhi;
};

#endif // !defined(AFX_CGLASS_H__C06EF99E_E23C_4930_A410_8401E0916897__INCLUDED_)
