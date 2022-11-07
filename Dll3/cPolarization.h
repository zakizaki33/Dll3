#if !defined(AFX_CPOLARIZATION_H__E0D062C6_B364_474F_8757_F491FF9A8FDD__INCLUDED_)
#define AFX_CPOLARIZATION_H__E0D062C6_B364_474F_8757_F491FF9A8FDD__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "complex.h"
#include "matrix.h"
#include "list.h"
#include "vector.h"

class cPolarization  
{
	double S0,S1,S2,S3;    // Stokes parameter
	static void rotate(complex &Ex,complex &Ey,double xAngleDeg);

public:
	static void JonesToStokes(double& S0,double& S1,double& S2,double &S3,complex Ex,complex Ey);
	static void StokesToJones
		(complex &Ex,complex &Ey,double S1,double S2,double S3);
	static matrix<double> MullerMatrix
	    (double xAngleDeg, double Tx,double ArgExDeg, double Ty,double ArgEyDeg);

	static double DepolarizerEffect
		(double Angle1Deg, double dPhiRangeY1Deg,double dPhiRangeX1Deg,
		 double Angle2Deg, double dPhiRangeY2Deg,double dPhiRangeX2Deg,
		 double DetectorEfficiencyXYRatio, 
		 int IsGaussNotRect);
	static void EllipseShape
		(double& a,double& b,double& phi_deg,double a1,double a2,double delta_deg);
	static void Analyze
		(double& a1,double& a2,double& delta_deg,double a,double b,double phi_deg);
	
	cPolarization();
	cPolarization(const complex& Ex,const complex& Ey);
	cPolarization(const vector<complex>& E);
	void Initialize
	    (double xAngleDeg,double I,double IxRatio,double ArgExDeg,double ArgEyDeg,double PolDeg);
	void Apply(double xAngleDeg, double Tx,double ArgExDeg, double Ty,double ArgEyDeg);
	double Intensity() const;
	double PolarizedIntensity(double AngleDeg) const;
	double Amplitude(double AngleDeg) const;
	double EllipseLongAxisDirection() const;
	double EllipseLongShortRatio() const;
	double EllipseLongShortIntensityRatio() const;
	double Stokes(int i) const;
	double PolDeg() const;
	double PoincareLat() const;
	double PoincareLng() const;
};

#endif // !defined(AFX_CPOLARIZATION_H__E0D062C6_B364_474F_8757_F491FF9A8FDD__INCLUDED_)
