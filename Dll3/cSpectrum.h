#if !defined(AFX_CSPECTRUM_H__A83B28E4_4EBF_4F4F_8571_FA28844749CA__INCLUDED_)
#define AFX_CSPECTRUM_H__A83B28E4_4EBF_4F4F_8571_FA28844749CA__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <fstream>
#include <string>
#include <math.h>
#include "cXYList.h"
#include "matrix.h"
#include "cMaterial.h"

class cSpectrum
{
	cXYList spectrum;
	void make_constants();
	cXYList x_bar,y_bar,z_bar;
	cXYList A,C,D65,WHITE;
	double Tcp(double& duv);

public:
	static double WLpb1,WLpb2;
	static double StandardIlluminant(std::string name,double wl);
	static double BlackBodyRadiation(double temp,double wl_nm,double wl0_nm);
	static double LuminousEfficiency(double wl_nm);
	static void SpectrumLocus(double wl_nm,double& x,double &y);
	static double SpectrumLocusX(double wl_nm);
	static double SpectrumLocusY(double wl_nm);
	static void WhitePoint(double& xw,double& yw);
	static double XYToDominantWavelength(double x,double y);
	static double ExcitationPurity(double x,double y);
	static long XYToApproximateRGB(double x,double y);
	static long ApproximateColor(double wl_nm);
	static double LumenToWatt(double wl_nm,double flux_lumen);
	static double WattToLumen(double wl_nm,double flux_watt);

	cSpectrum();
	cSpectrum(const cSpectrum& x);
	cSpectrum& operator=(const cSpectrum& x);
	int Size();
	double WaveLength(int i);
	double Data(double wl);
	cSpectrum& AddData(double wl,double value);
	cSpectrum& AddConstantData(double wl_start,double wl_end,double wl_step,double value);
	cSpectrum& RemoveAllData();
	cSpectrum Copy() const;
	cSpectrum& Apply(cXYList &x);
	int Apply(const std::string& filename);
	cSpectrum& ApplyBlackBodyRadiation(double temp);
	cSpectrum& ApplyLuminousEfficiency();
	cSpectrum& ApplyColorGlass(std::string filename,double thickness);
	cSpectrum& Whiten();
	cSpectrum& Normalize(double new_max);
	cSpectrum& NormalizeByLocal(double wl_start,double wl_end,double new_max);
	cSpectrum& Shift(double multiplier);
	cSpectrum& ReverseTR();

	double x();
	double y();
	double u();
	double v();
	double Tcp();
	double duv();
	double LuminousFlux();
	double x(std::string IlluminantName,double temp);
	double y(std::string IlluminantName,double temp);
	double u(std::string IlluminantName,double temp);
	double v(std::string IlluminantName,double temp);
	double Tcp(std::string IlluminantName,double temp);
	double duv(std::string IlluminantName,double temp);
	double LuminousFluxRatio(std::string IlluminantName,double temp);
	double Integral(double wl_start,double wl_end);
	double Average(double wl_start,double wl_end);
		
	double Peak();
	double PeakLocal(double wl_start,double wl_end);
	double FW(double threshold);
	double FWStart(double threshold);
	double FWEnd(double threshold);
	double FWCenter(double threshold);
	double FWLocal(double wl_start,double wl_end,double threshold);
	double FWStartLocal(double wl_start,double wl_end,double threshold);
	double FWEndLocal(double wl_start,double wl_end,double threshold);
	double FWCenterLocal(double wl_start,double wl_end,double threshold);
	double FWHM();
	double FWHMLocal(double wl_start,double wl_end);
	double TransitionPoint(double val);
	double TransitionPointLocal(double wl_start,double wl_end,double val);
	double TransitionInterval(double val1,double val2);
	double TransitionIntervalLocal(double wl_start,double wl_end,double val1,double val2);

	std::string Str();
	friend std::istream& operator>>(std::istream& from, cSpectrum& x);
	friend std::ostream& operator<<(std::ostream& to, const cSpectrum& x);
	int Open(const std::string& filename);
	int Save(const std::string& filename);
};

#endif // !defined(AFX_CSPECTRUM_H__A83B28E4_4EBF_4F4F_8571_FA28844749CA__INCLUDED_)
