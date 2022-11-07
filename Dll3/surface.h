#if !defined(AFX_SURFACE_H__7C06A9C0_7DD0_11D7_BE64_BA777A4A0750__INCLUDED_)
#define AFX_SURFACE_H__7C06A9C0_7DD0_11D7_BE64_BA777A4A0750__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include <fstream>
#include "string.h"
#include "cZernike.h"
#include "cDcon.h"
#include "cModLegendre.h"

struct freeform{
	enum { N=20 };     //  ヘッダファイルでは int N=20; などとできないのでenumを使う．
	double bb[N+1][N+1];
	double null;

	freeform();
	freeform(const freeform& x);
	freeform& operator=(const freeform& x);
	friend bool operator ==(const freeform &a,const freeform &b);

	double& b(int m,int n);
	int NotNull();
	freeform reversed();
	void scale(double m);

	std::string GetTerms();
	void        SetTerms(std::string terms);
};


class surface  
{
	double rr,cc, r0,c0;
	void match_r_c();
	double normh;
	double fi,pi, fi0,pi0;
	void match_fi_pi();
	double acoa,cacoa, acoa0,cacoa0;
	void match_acoa_cacoa();
	double bcoa,cbcoa, bcoa0,cbcoa0;
	void match_bcoa_cbcoa();
public:
	static int file_ver;
	double& r();
	double& c();
	double Newton,As0,As45,NewtonTol,AsTol; int rVariable;
	int asph_type; 
	double kp,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a18,a20;
	double& NormH();
	freeform b;
	cZernike zernike;
	cDcon Dcon;
	cModLegendre legendre;
	cSpline spline;
	double Apv,sfx,sfy;
	int UserDefSurf;
	int cylinder;
	int fresnel; double rbase;
	double ry,rx, kpy,kpx; int IsXToroid;
	double coneangle;
	double& fideal();
	double& pideal();
	double& aCOA(); double& caCOA();
	double& bCOA(); double& cbCOA();
	double tCOA;
	double SA0,CM0;   // 特性係数
	int grating,difforder; double gpitch,grx,gry,grz;
	double Diffusion;
	int Fish;
	int EAtype; double EAy,EAx,CHMy,CHMx,EAdy,EAdx;
	int decenter_type;
	double dx,dy,dz, rox,roy,roz; 
	int order;
	int ret;
	double dx1,dy1,dz1, rox1,roy1,roz1;
	int order1;
	std::string CoatName; int CoatReverse;
	std::string rem;

public:
	surface();
	surface(const surface& x);
	surface& operator=(const surface& x);
	friend std::ostream& operator<<(std::ostream& to,surface& x);
	friend std::istream& operator>>(std::istream& from,surface& x);
	friend bool operator==(surface& a,surface& b);
	surface reversed();
	void scale(double m,int with_EA,int with_decenter);
};

#endif // !defined(AFX_SURFACE_H__7C06A9C0_7DD0_11D7_BE64_BA777A4A0750__INCLUDED_)
