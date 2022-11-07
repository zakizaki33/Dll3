#if !defined(AFX_CZERNIKE_H__2E92D4E9_480D_4110_B685_846E96070E4D__INCLUDED_)
#define AFX_CZERNIKE_H__2E92D4E9_480D_4110_B685_846E96070E4D__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "matrix.h"
#include "list.h"
#include "general_func.h"
#include "vector.h"
#include "string.h"

class cZernike
{
	double r0;
	double r0fix;
	int normalize;
	int terms;
	int IsFringeOrder;
	int jBase;    // jBase==0のときピストンを第0項とする．それ以外では第1項とする．
	
	list<double> x,y,z;
	list<int> i,j;
	matrix<double> A,X,F;
	bool fitted;
	int Digits;
	void ConvertR0();

public:
	static int TotalTerms(int order,int IsFringeOrder);
	static int jNumber(int l,int n,int IsFringeOrder,int jBase);
	static int StandardNo(int FringeNo,int jBase);
	static int lNumber(int j,int IsFringeOrder,int jBase);
	static int mNumber(int j,int IsFringeOrder,int jBase);
	static int nNumber(int j,int IsFringeOrder,int jBase);
	static int Order(int j,int IsFringeOrder,int jBase);
	static void ConvertR0(double *a,int N,double r1,double r2,int Normalized,int IsFringeOrder);
	int TotalTerms(int order);
	int jNumber(int l,int n);
	int StandardNo(int FringeNo);
	int lNumber(int j) const;
	int mNumber(int j) const;
	int nNumber(int j) const;
	int Order(int j) const;
	
	cZernike();
	virtual ~cZernike();

	static double R(int l,int n,double rho);
	static double Rrho(int l,int n,double rho);
	static double Rrhorho(int l,int n,double rho);
	double U(int l,int n,double x,double y) const;
	double U(int j,double x,double y) const;
	double Ur(int l,int n,double x,double y) const;
	double Ur_pol(int l,int n,double r,double th) const;
	double Uth(int l,int n,double x,double y) const;
	double Urr(int l,int n,double x,double y) const;
	double Urr_pol(int l,int n,double r,double th) const;
	double Uthth(int l,int n,double x,double y) const;
	double Urth(int l,int n,double x,double y) const;
	double Urth_pol(int l,int n,double r,double th) const;
	double Ux(int l,int n,double x,double y) const;
	double Ux(int j,double x,double y) const;
	double Uy(int l,int n,double x,double y) const;
	double Uy(int j,double x,double y) const;
	double Uxx(int l,int n,double x,double y) const;
	double Uxx(int j,double x,double y) const;
	double Uyy(int l,int n,double x,double y) const;
	double Uyy(int j,double x,double y) const;
	double Uxy(int l,int n,double x,double y) const;
	double Uxy(int j,double x,double y) const;
	
	double GetR0() const;
	void   SetR0(double value);
	void   CalcR0();
	double GetR0Fix() const;
	void   SetR0Fix(double value);
	int  GetNormalize() const;
	void SetNormalize(int value);
	int  GetIsFringeOrder() const;
	void SetIsFringeOrder(int value);
	int  GetJBase() const;
	void SetJBase(int value);

	void DataClear();
	int NumberOfData();
	int  GetNumberOfTerms() const;
	void SetNumberOfTerms(int value);
	int  GetMaxOrder() const;
	void SetMaxOrder(int value);
	double GetXData(int i);
	double GetYData(int i);
	double GetZData(int i);
	int    GetIData(int i);
	int    GetJData(int i);
	void SetData(double x,double y,double z,int i,int j);
	void SetData(double x,double y,double z);
	void SetDigits(int value);

	int Fit();

	double GetC(int j);
	void   SetC(int j,double value);
	double& C(int j);
	int    NotNull();
	std::string GetCoefficients();
	void        SetCoefficients(std::string com);
	void   FitZToCoefficient();
	double Z(double x,double y,double Rref) const;
	double Z(double x,double y) const;
	double ZApproximate(double x,double y,double Rref);
	double ZApproximate(double x,double y);
	double Zr(double x,double y,double Rref,int RemoveTilt) const;
	double Zx(double x,double y) const;
	double Zy(double x,double y) const;
	double Zxx(double x,double y) const;
	double Zyy(double x,double y) const;
	double Zxy(double x,double y) const;
	double AxialR(double x,double y,double Rref,int RemoveTilt);
	void   ParaxialExpansion(double &cx,double &cxy,double &cy);
	void   ParaxialR(double &r1,double &r2,double &axis_deg,double Rref);
	double ParaxialR1(double Rref);
	double ParaxialR2(double Rref);
	double ParaxialRaxis();
	double Error(int i);
	double RMSError();
	double PVError();
	int RemoveTerms(int piston,int tilt,int sph);
	int AdjustTerms(int piston,int tilt,int sph);
	double RMS(int tilt,int sph,int cyl,int high_order);
	double Shear(double dx,double dy,int order);
};

#endif // !defined(AFX_CZERNIKE_H__2E92D4E9_480D_4110_B685_846E96070E4D__INCLUDED_)

