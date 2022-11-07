#if !defined(AFX_COPTICS_H__152EFA60_519B_11D7_BE64_8A53C09D0450__INCLUDED_)
#define AFX_COPTICS_H__152EFA60_519B_11D7_BE64_8A53C09D0450__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include <time.h>
#include "complex.h"
#include "geometry.h"
#include "cMaterial.h"
#include "general_func.h"
#include "string.h"
#include "cSpectrum.h"
#include "cGlass.h"

class cOptics  
{
	list<std::string> symbol, value;
public:
	cOptics();
	cOptics(const cOptics& x);
	cOptics& operator=(const cOptics& x);

	static std::string MaterialName(const std::string& gname);
	static double HerzEq(double Nd,double nud,const std::string& color);
	static complex Index(const std::string& s,double wl_nm);
	static double  GroupIndex(const std::string& s,double wl_nm);
	static complex IndexDerivative2(const std::string& s,double wl_nm);
	static double ReIndex(const std::string& s,double wl_nm);
	static double ImIndex(const std::string& s,double wl_nm);
	static double AbbeNumber(const std::string& glass,
		                     const std::string& color0,const std::string& color1,const std::string& color2);
	static double AbbeNumber(const std::string& glass);
	static std::string GlassCode(double nd,double nud);
	static std::string GlassCode(const std::string& s);
	static double Wavelength(const std::string& color);

	static double PolygonDutyFactor(double Facets,double ScanFullAngle);
	static double PolygonMinDiaNotIntensityVary
	              (double Facets,double BeamDia,double ScanFullAngle,double FeedAngle,double RollOffLength);
	static double BeamMaxDiaNotIntensityVary
	              (double Facets,double PolygonDia,double ScanFullAngle,double FeedAngle,double RollOffLength);
	static double BeamMaxDiaCenter(double Facets,double PolygonDia,double FeedAngle,double RollOffLength);
	static double PolygonEfficiency
	              (double Facets,double PolygonDia,double BeamDia,double FeedAngle,double RollOffLength,double ScanAngle);
	
	static point edge(double phi_x,double phi_y,int is_rect,double th_rad);
	static std::string Spherometer(double D,double H,double Rref=0,double Href=0);
	static std::string Plate(std::string glass,std::string color,double thickness,double incident_angle,
	                         double dndtemp,double alpha,std::string what);
	static double Deviation(std::string glass,std::string glass1,std::string color,double incident_angle);
	static std::string PrismDeviation(std::string glass,std::string color,double apex_angle,double incicdent_angle);
	static double MeasuredDistortion(double x1,double y1,double x2,double y2,double x3,double y3);
	static double BerekDOF(std::string &s,double wl_nm,double NAo,double TotalM,double VAmin);
	static double ScheimpflugTH1(double ms,double th_deg,double th_l_deg);
	static std::string Scheimpflug(double ms,double th_deg,double th_l_deg);
	static double BrewsterAngle(std::string glass_in,std::string glass_out,std::string color);
	static double ChirpedPulseWidth(double dt0_fs,double GDD_fs2);
	static double ODToT(double od);
	static double CoherenceLength(double wl_nm,double dwl_fwhm_nm);

	double Derivative(std::string command,double *x,double dx);
	std::string interpret_target(std::string sentence,
	                             std::string &com,std::string &sign,double &tol,double &target,double &weight,
	                             double set_weight=1);
	
	void symboltoval(std::string &s) const;
	static double arg(const complex& z);
	std::string arg(const std::string& sentence, int n);
	bool is_numeric(std::string& s);
	std::string cmd(const std::string& command,int val);
	virtual std::string scmd(std::string com,int val);
	virtual int push(){ return 0; };
	virtual int pop() { return 0; };
};

#endif // !defined(AFX_COPTICS_H__152EFA60_519B_11D7_BE64_8A53C09D0450__INCLUDED_)
