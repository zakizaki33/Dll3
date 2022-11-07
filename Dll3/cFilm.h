#if !defined(AFX_CFILM_H__57E20106_D1C5_11D6_8AB6_8CB66EF9E250__INCLUDED_)
#define AFX_CFILM_H__57E20106_D1C5_11D6_8AB6_8CB66EF9E250__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <string>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include "cOptics.h"
#include "matrix.h" // Client.cpp で呼び出している
#include "list.h" // Client.cpp で呼び出している
#include "cSpectrum.h"
#include "cLeneq.h" // Client.cpp で呼び出している
#include "string.h"
#include "cFitting.h"

class cFilm : public cOptics 
{
	std::string filename;
	int k;
	double wl0;    // nd測定波長．すなわちndのnは波長wl0での屈折率．
	               //   (膜厚がndで定義されたとき，nが波長に依存してしまうと
	               //   物理膜厚dが波長によって変化することになってしまう．)
	int tMode; enum {OPTICAL,PHYSICAL,OPTICALQW};
	double th_med_deg;
	std::string *mname;
	double *t; int *tVariable;
	double *ndTol, *dTol;
	double SubThickness;
	void create();
	void initialize();
	void erase();
	void assignment(const cFilm& x);
	stack<cFilm> Stack;
	list<cFilm> List;
	complex costh(int i,double th_med_deg,double wl_nm);
	complex Y(int i,double th_med_deg,double wl_nm,char s_or_p);
	matrix<complex> M;
	void calc_M(double th_med_deg,double wl_nm,char s_or_p);
	matrix<complex> Ms;
	void calc_Ms(double th_med_deg,double wl_nm,char s_or_p);
	double nd(int i,double wl);
	double nd(int i);
	double phase_thickness(double th_med_deg,double wl_nm);
	cXYList xylist;
	void make_xylist(std::string AboutWhat,double wl1,double wl2,double dwl);
	cSpectrum spectrum;
	void make_spectrum(std::string AboutWhat,double wl1,double wl2,double dwl);

public:
	static double NdOfSingleLayerARCoat(std::string material,double wl_nm,double th_deg);
	static double EquiFilmN(std::string SideMaterial1,std::string CenterMaterial2,
	                        double TotalPhi,double Phi2Phi1Ratio);
	static double EquiFilmN(std::string SideMaterial1,std::string CenterMaterial2,
	                        double TotalNd,double Nd2Nd1Ratio,double wl_nm);
	static double EquiFilmPhi(std::string SideMaterial1,std::string CenterMaterial2,
	                         double TotalPhi,double Phi2Phi1Ratio);
	static double EquiFilmPhi(std::string SideMaterial1,std::string CenterMaterial2,
	                         double TotalNd,double Nd2Nd1Ratio,double wl_nm);
	cFilm(int k=0);
	cFilm(const cFilm& x);
	virtual ~cFilm();
	void Reset();
	cFilm& operator=(const cFilm& x);
	int push(); int pop(); int undo(); int redo();
	void AddToList();
	void ClearList();
	int GetListData(int num);
	int  get_k();
	void set_k(int value);
	double get_wl0();
	void   set_wl0(double value_nm);
	std::string  get_tMode();
	void         set_tMode(std::string value);
	void tToOptical();
	void tToPhysical();
	void tToOpticalQW();
	double d(int i);
	void set_d(int i,double value);
	double get_th_med_deg();        
	void set_th_med_deg(double th_deg);
	std::string get_mname(int i);   
	void set_mname(int i,std::string name); 
	void ReplaceMaterial(std::string from,std::string to);
	double get_t(int i);           
	void   set_t(int i,double value);
	int  get_tVariable(int i);      
	void set_tVariable(int i,int value);
	double get_ndTol(int i);
	void   set_ndTol(int i,double value);
	double get_dTol(int i);
	void   set_dTol(int i,double value);
	double TotalNd();
	double TotalD();
	static int SubstrateAbsorptionEnable;
	double getSubThickness(); 
	void setSubThickness(double value);

	complex rs(double th_med_deg,double wl_nm);
	complex rs(double wl_nm);
	complex rp(double th_med_deg,double wl_nm);	
	complex rp(double wl_nm);	
	complex ts(double th_med_deg,double wl_nm);
	complex ts(double wl_nm);
	complex tp(double th_med_deg,double wl_nm);
	complex tp(double wl_nm);
	double Rs(double th_med_deg,double wl_nm);
	double Rs(double wl_nm);
	double Rp(double th_med_deg,double wl_nm);	
	double Rp(double wl_nm);	
	double Rave(double th_med_deg,double wl_nm);
	double Rave(double wl_nm);
	double Ts(double th_med_deg,double wl_nm);
	double Ts(double wl_nm);
	double Tp(double th_med_deg,double wl_nm);
	double Tp(double wl_nm);
	double Tave(double th_med_deg,double wl_nm);
	double Tave(double wl_nm);
	double As(double th_med_deg,double wl_nm);
	double As(double wl_nm);
	double Ap(double th_med_deg,double wl_nm);
	double Ap(double wl_nm);
	double Aave(double th_med_deg,double wl_nm);
	double Aave(double wl_nm);
	double RsBothSide(double th_med_deg,double wl_nm);
	double RsBothSide(double wl_nm);
	double RpBothSide(double th_med_deg,double wl_nm);
	double RpBothSide(double wl_nm);
	double RaveBothSide(double th_med_deg,double wl_nm);
	double RaveBothSide(double wl_nm);
	double TsBothSide(double th_med_deg,double wl_nm);
	double TsBothSide(double wl_nm);
	double TpBothSide(double th_med_deg,double wl_nm);
	double TpBothSide(double wl_nm);
	double TaveBothSide(double th_med_deg,double wl_nm);
	double TaveBothSide(double wl_nm);
	double ArgRs(double th_med_deg,double wl_nm);
	double ArgRs(double wl_nm);
	double ArgRp(double th_med_deg,double wl_nm);	
	double ArgRp(double wl_nm);
	double DArgR(double th_med_deg,double wl_nm);
	double DArgR(double wl_nm);
	double ArgTs(double th_med_deg,double wl_nm);
	double ArgTs(double wl_nm);
	double ArgTp(double th_med_deg,double wl_nm);
	double ArgTp(double wl_nm);
	double DArgT(double th_med_deg,double wl_nm);
	double DArgT(double wl_nm);
	double Tsub(double th_med_deg,double wl_nm);
	double Tsub(double wl_nm);
	
	double GD(double th_med_deg,std::string AboutWhat,double wl_nm,double dwl_nm);
	double GD(std::string AboutWhat,double wl_nm,double dwl_nm);
	double GDD(double th_med_deg,std::string AboutWhat,double wl_nm,double dwl_nm);
	double GDD(std::string AboutWhat,double wl_nm,double dwl_nm);
	
	double GDRs(double th_med_deg,double wl_nm,double dwl_nm);
	double GDRs(double wl_nm,double dwl_nm);
	double GDRp(double th_med_deg,double wl_nm,double dwl_nm);
	double GDRp(double wl_nm,double dwl_nm);
	double GDTs(double th_med_deg,double wl_nm,double dwl_nm);
	double GDTs(double wl_nm,double dwl_nm);
	double GDTp(double th_med_deg,double wl_nm,double dwl_nm);
	double GDTp(double wl_nm,double dwl_nm);

	double GDDRs(double th_med_deg,double wl_nm,double dwl_nm);
	double GDDRs(double wl_nm,double dwl_nm);
	double GDDRp(double th_med_deg,double wl_nm,double dwl_nm);
	double GDDRp(double wl_nm,double dwl_nm);
	double GDDTs(double th_med_deg,double wl_nm,double dwl_nm);
	double GDDTs(double wl_nm,double dwl_nm);
	double GDDTp(double th_med_deg,double wl_nm,double dwl_nm);
	double GDDTp(double wl_nm,double dwl_nm);

	std::string DispersionTable(double wl_start,double wl_end,int wl_points);
	
	std::string FileName();
	std::string FilmData(double monitor_wl_nm);
	std::string FilmDataZEMAX();
	friend std::ostream& operator<<(std::ostream& to,cFilm& x);
	friend std::istream& operator>>(std::istream& from,cFilm& x);
	int open(std::string filename);
	int save(std::string filename);
	void reverse(double wl_nm);
	void reverse();
	void AdjustThickness(double multiplier);
	void AdjustForIncidentAngle(double new_th_med_deg, double wl_nm);
	void ndPerturbe();
	cFilm ndPerturbed();
	void dPerturbe();
	cFilm dPerturbed();
	double xOfXYZ(std::string AboutWhat,std::string IlluminantName);
	double yOfXYZ(std::string AboutWhat,std::string IlluminantName);
	double uOfUCS(std::string AboutWhat,std::string IlluminantName);
	double vOfUCS(std::string AboutWhat,std::string IlluminantName);
	double Tcp(std::string AboutWhat,std::string IlluminantName);
	double duv(std::string AboutWhat,std::string IlluminantName);
	double DominantWavelength(std::string AboutWhat,std::string IlluminantName);
	double ExcitationPurity(std::string AboutWhat,std::string IlluminantName);
	double LuminousFluxRatio(std::string AboutWhat, std::string IlluminantName);
	double Average(std::string AboutWhat,double wl1,double wl2,double dwl);
	double Max(std::string AboutWhat,double wl1,double wl2,double dwl);
	double Min(std::string AboutWhat,double wl1,double wl2,double dwl);
	double FWHM(std::string AboutWhat,double wl1,double wl2,double dwl);
	double FWHMCenter(std::string AboutWhat,double wl1,double wl2,double dwl);
	double FW(std::string AboutWhat,double wl1,double wl2,double dwl,double val);
	double FWCenter(std::string AboutWhat,double wl1,double wl2,double dwl,double val);
	double TransitionPoint(std::string AboutWhat,double wl1,double wl2,double dwl,double val);
	double TransitionInterval(std::string AboutWhat,double wl1,double wl2,double dwl,double val1,double val2);
	double Integral(std::string AboutWhat,double wl1,double wl2,double dwl);
	double FiniteSlit(std::string AboutWhat,double wl,double dwl,int n);

	double tDerivative(std::string command,int i,double dt);
	
	double optimize(std::string command,double dt_ratio,double rho);
	void RemoveMinusLayers(double wl_nm);
	void RemoveMinusLayers();

	static list<cSpectrum> filters;
	static list<std::string> filternames;
	static std::string FilterDataFileLocation;
	static int AddFilter(std::string filename);
	static void RemoveFilters();
	static std::string FilterNames();
	static double Tfilters(double wl_nm);

	std::string scmd(std::string com,int val);
};

#endif // !defined(AFX_CFILM_H__57E20106_D1C5_11D6_8AB6_8CB66EF9E250__INCLUDED_)
