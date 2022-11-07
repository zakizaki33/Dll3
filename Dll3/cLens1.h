#if !defined(AFX_CLENS1_H__B050DAA9_7467_4AA8_8DFC_82CABBAF2484__INCLUDED_)
#define AFX_CLENS1_H__B050DAA9_7467_4AA8_8DFC_82CABBAF2484__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
// #define _AFXDLL

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include "cOptics.h"
#include "surface.h"
#include "medium.h"
#include "cLeneq.h"   // Client.cpp で呼び出している
#include "list.h"     // Client.cpp で呼び出している
#include "cFilm.h"
#include "cShapes.h"
#include "cSpectrum.h"
#include "geometry.h"
#include "vector.h"   // Client.cpp で呼び出している
#include "cPolarization.h"
#include "string.h"
#include "cZernike.h"
#include "cFitting.h"
#include "cSelfoc.h"
#include "cImage.h"
#include "AUserDef.h"
#include "cFitCircle.h"

struct group {
	int i1,i2;
	group(int i1=0,int i2=0);
	bool operator==(const group& x) const;
	bool operator>(const group& x) const;
	std::string str() const;
	std::string str(int k) const;
};

struct ray_data{
	int n;
	std::string color;
	double *x,*y,*z, *X,*Y,*Z, *X1,*Y1,*Z1;
	int *Ref;
	list<int> GrinRay_i;
	list< vector<double> > GrinRay;
	ray_data();
	ray_data(int n);
	ray_data(const ray_data& a);
	~ray_data();
	ray_data& operator=(const ray_data& a);
};

struct cPoint {
	point p;
	double wl;       // 波長(nm)
	double weight;
	double opl;      // 光路長
	cPoint();
	cPoint(double x,double y,double wl,double weight,double opl);
	cPoint(double x,double y,double wl,double weight);
	cPoint(double x,double y,double wl);
};

class cSpot {
	list<cPoint> spot;
	point origin;
public:
	cSpot();
	std::string Unit;
	int GetSize() const;
	void AddTail(const cPoint& p);
	int GetData(cPoint& p,int i);
	void RemoveAll();
	point GravityCenter();
	
	double InertiaPrincipalAxis();
	void PrincipalWidth(double &Wx,double &Wy);
	double PrincipalWx();
	double PrincipalWy();

	void OriginToGravityCenter();
	void OriginToNewCenter(double x,double y);
	complex OTF(double nu_y,double nu_x);
	double MTF(double nu_y,double nu_x);
	double TotalIntensity();
	double Intensity(double y,double x,double SensorPhi);
	double EncircledEnergy(double SensorPhi);
	double XIntensity(double x,double SensorFW);
	double YIntensity(double y,double SensorFW);
	double RmsPhi();
	double XYAbsMax();
	double XMax();
	double YMax();
	double XMin();
	double YMin();
	double XWidth();
	double YWidth();
	double XGravityCenter();
	double YGravityCenter();
	void RemoveTC();
};



class cLens1 : public cOptics
{
	friend class cLens;     // こうしないとcLensからcLens1のプライベートメンバにアクセスできない

	std::string filename;

	int k, mems;
	int cn, long_cn,short_cn,mid_cn;
	std::string *color; double* wl; double *colorweight;
	int IndexDigits;

	std::string Note;
	double FRW;             // ニュートンリング観察波長(nm)

	double s,t;

	surface *surf;
	double& r(int i); 
	double& rObj();
	double& rImage();
	double& c(int i);
	double& Newton(int i); double& As0(int i); double& As45(int i);
	double& NewtonTol(int i); double& AsTol(int i);
	int& rVariable(int i);
	int& asph_type(int i);
	double& kp(int i); double& NormH(int i); 
	double& a1(int i); double& a2(int i); double& a3(int i);
	double& a4(int i); double& a5(int i); double& a6(int i); double& a7(int i); double& a8(int i); double& a9(int i);
	double& a10(int i); double& a11(int i); double& a12(int i); double& a13(int i); double& a14(int i);
	double& a15(int i); double& a16(int i); double& a18(int i); double& a20(int i);
	double& b(int i,int m,int n);
	cZernike& zernike(int i);
	double& ZC(int i,int j);
	cDcon& Dcon(int i);
	double& DconRn(int i);
	double& DconA4(int i); double& DconA6(int i); double& DconA8(int i); double& DconA10(int i);
	double& DconA12(int i); double& DconA14(int i); double &DconA16(int i);
	double& DconA18(int i); double& DconA20(int i);
	cModLegendre& legendre(int i);
	double& LeC(int i,int m,int n);   // "LC" は使用済
	cSpline& spline(int i);
	double& SplineH(int i,int j);
	double& SplineZ(int i,int j);
	double& Apv(int i); double& sfx(int i); double& sfy(int i);
	int& UserDefSurf(int i);
	int& cylinder(int i);
	int& fresnel(int i);
	double& rbase(int i);
	double& ry(int i); double& rx(int i); double& kpy(int i); double& kpx(int i); int& IsXToroid(int i);
	double& coneangle(int i);
	double& fideal(int i); double& pideal(int i);
	double& aCOA(int i); double& caCOA(int i); double& bCOA(int i); double& cbCOA(int i); double& tCOA(int i);
	double& SA0(int i); double& CM0(int i);
	int& grating(int i); int& difforder(int i); double& gpitch(int i);
	double& grx(int i); double& gry(int i); double& grz(int i);
	double& Diffusion(int i);
	int& Fish(int i);
	int& EAtype(int i); double& EAy(int i); double& EAx(int i); double& CHMy(int i); double& CHMx(int i);
	double& EAdy(int i); double& EAdx(int i);
	int& decenter_type(int i);
	double& dx(int i); double& dy(int i); double& dz(int i);
	double& rox(int i); double& roy(int i); double& roz(int i);
	int& order(int i);
	int& ret(int i);
	double& dx1(int i); double& dy1(int i); double& dz1(int i);
	double& rox1(int i); double& roy1(int i); double& roz1(int i);
	int& order1(int i);
	std::string& CoatName(int i); int& CoatReverse(int i);
	std::string& rem(int i);

	medium *med;
	double& d(int i);
	double& delta_d(int i);
	int& dVariable(int i);
	std::string& gname(int i);
	int& gVariable(int i);
	double& Nd(int i);
	double& Nud(int i);
	double& dN(int i);
	double N(int i,int j);
	int IsGlass(int i);
	int IsGrin(int i);
	double rA(int i,int j);
	double GrinPhi(int i);
	double NNuActual(int i);

	double s1fix;
	double yObjectMax,xObjectMax;
	int stop;
	double EPD,EPDx, EPy,EPx;
	int nSpot;
	double SourcePhiY,SourcePhiX,SourceAxis; int SourceAreaType;

	double *phi, *u,*h,*up,*hp;
	double sum(const double *x, int i1, int i2);
	double rms(const double *x, int i1, int i2);
	double exponent(const double *x);
	double *hQ,*hQp;
	double *SA,*CM,*AS,*DS,*PT,*LC,*TC,*LC2,*TC2; 
	double SAt,CMt,ASt,DSt,PTt,LCt,TCt,LC2t,TC2t;
	double SAe,CMe,ASe,DSe,PTe,LCe,TCe,LC2e,TC2e;
	double *SA5,*CM41Z,*CM41,*CM23,*CM23P,*CM23Z,*CM41ALL,*CM23ALL; 
	double SA5t,CM41Zt,CM41t,CM23t,CM23Pt,CM23Zt,CM41ALLt,CM23ALLt;
	double SA5e,CM41Ze,CM41e,CM23e,CM23Pe,CM23Ze,CM41ALLe,CM23ALLe;
	double *SA32F,*SA32Z,*SA32,*SA32ALL;
	double SA32Ft,SA32Zt,SA32t,SA32ALLt;
	double SA32Fe,SA32Ze,SA32e,SA32ALLe;
	double *AS5,*SG5, *DS5;
	double AS5t,SG5t, DS5t;
	double AS5e,SG5e, DS5e;
	double *SAP,*CMP,*ASP,*DSP,*LCP;
	double SAPt,CMPt,ASPt,DSPt,LCPt;
	double SAPe,CMPe,ASPe,DSPe,LCPe;
	double *PRE,*DSE1,*DSE2,*ASE,*PTE,*CME;
	double PREt,DSE1t,DSE2t,ASEt,PTEt,CMEt;
	double PREe,DSE1e,DSE2e,ASEe,PTEe,CMEe;

	void create();
	void initialize();
	void erase();
	void assignment(const cLens1& x);

	void ParaxialExpansion(int i,double &cx,double &cxy,double &cy);
	double A[5][5], A2[5];
	void AMatrixCalc(int i1,int i2,int j=1);
	void InvertAMatrix();

	std::string NormalizeUnit;
	int NormalizeType;
	stack<cLens1> Stack;
	int is_plane(int i);
	int is_aspheric(int i);
	int is_ideallens(int i);
	int b_size(int i);
	int is_freeform_surface(int i);
	int is_zernike_surface(int i);
	int is_Dcon_surface(int i);
	int is_legendre_surface(int i);
	int is_spline_surface(int i);
	int is_periodic_surface(int i);
	int is_userdef_surface(int i);
	int ea_is_defined();
	int ea_all_regular();
	int surface_kind(int i);
	int is_refract_surf(int i);
	int is_reflect_surf(int i);
	int is_acting_surf(int i);
	int is_stop(int i);
	int is_mask(int i);
	int is_solid(int i);
	int is_dummy(int i);
	double phi_y(int i); double phi_x(int i); double ea_y(int i); double ea_x(int i); double ea_max(int i);
	
	double freeform_z(int i,double x,double y);
	double freeform_zx(int i,double x,double y);
	double freeform_zy(int i,double x,double y);
	double freeform_zxx(int i,double x,double y);
	double freeform_zyy(int i,double x,double y);
	double freeform_zxy(int i,double x,double y);

	double periodic_z(int i,double x,double y);
	double periodic_zx(int i,double x,double y);
	double periodic_zy(int i,double x,double y);
	double periodic_zxx(int i,double x,double y);
	double periodic_zyy(int i,double x,double y);
	double periodic_zxy(int i,double x,double y);

	double surface_sag(int i,double y,double x,int newton_enable,int& domain_error);
	double surface_sag(int i,double y,double x,int newton_enable);
	vector<double> surface_normal(int i,double y,double x,int newton_enable);
	void surface_curvature(int i,double& cx,double& cy,double& cxy,double x,double y,
	                       const matrix<double>& T,int newton_enable);

	void decenter_in(int i,double& x,double& y,double& z,int translate);
	void decenter_in(int i,complex& x,complex& y,complex& z,int translate);
	void decenter_in_rev(int i,double& x,double& y,double& z,int translate);
	void decenter_in_rev(int i,complex& x,complex& y,complex& z,int translate);
	void decenter_out(int i,double& x,double& y,double& z,int translate);
	void decenter_out(int i,complex& x,complex& y,complex& z,int translate);
	void decenter_out_rev(int i,double& x,double& y,double& z,int translate);
	void decenter_out_rev(int i,complex& x,complex& y,complex& z,int translate);
	vector<double>  ret_decenter(vector<double> v,int translae,int i,int inv);
	vector<complex> ret_decenter(vector<complex> v,int translate,int i,int inv);

	vector<double> *o,*o0, *ex,*ex0,*ey,*ey0,*ez,*ez0;
	void make_coordinate(double defocus);

	void transform(double& x,double& y,double& z,int i0,int i0_pre,int i,int i_pre,int translate) const;
	void transform(vector<double>& v,int i0,int i0_pre,int i,int i_pre,int translate) const;
	void transform(complex& x,complex& y,complex& z,int i0,int i0_pre,int i,int i_pre,int translate) const;
	
	void transform_plane(double &a,double &b,double &c,double &d,int i0,int i0_pre,int i,int i_pre);
	void image_plane(double &a,double &b,double &c,double &d,int i);

	void transform_line(double &x0,double &y0,double &z0,double &l,double &m,double &n,
	                    int i0,int i0_pre,int i,int i_pre);
	void image_line(double &x0,double &y0,double &z0,double &l,double &m,double &n,int i);
	
	int raytrace(int i,int j, double xp,double yp,double zp,double X,double Y,double Z,
		double& x,double& y,double& z,double & X1,double& Y1,double& Z1, 
		int EA_enable,int mask_enable,
		int as_trace,vector<double>& ha,vector<double>& dQa,vector<double>& hb,vector<double>& dQb,
		int E_trace, vector<complex>& E);

	double dk(double defocus);
	double x_pupil,y_pupil;
	double *x,*y,*z, *X,*Y,*Z, *X1,*Y1,*Z1, *xi,*xi1; //, hyxi1,uy1,hx,ux1;
	vector<double> *ha,*dQa,*hb,*dQb, *ha1,*dQa1,*hb1,*dQb1;
	vector<complex> *E,*E1;
	double apodization_amplitude;
	list<point> object_list;
	double *optical_path;
	double GrinDs;
	int GrinToListStep;
	list<int> GrinRay_i;
	list< vector<double> > GrinRay;

	list<point> gpoints;
	point gpoint_principal;
	int MakePupilGrid(double yObj,double xObj,int FindPupil,int n,int j,int Normalize,
	                  int Randomize,int makecoordinate=1,int ConsiderSymmetry=0);
	int MakePupilGrid(int n);

	double spot_area;

	matrix<double> rms_map;

	matrix<double> opd_map_x, opd_map_y, opd_map_w, opd_map_p, opd_map_a;
	matrix<int>    opd_map_e;
	cZernike opd_map_zernike;
	cImage psf_map;
	double psf_map_w_x,psf_map_w_y, psf_map_phi_x,psf_map_phi_y,psf_map_phi_threshold;
	cImage otf_map;
	double otf_map_w_x,otf_map_w_y;
	double otf_map_rayleigh_nu_x,otf_map_rayleigh_nu_y;

	double TcMinGlass,TcMinAir,TeMinGlass,TeMinAir;
	double CenterThickness(int i);
	list<group> groups;
	enum { SA_GRAPH, AS_GRAPH, DIST_GRAPH, DIST_CHART, DELTAH_GRAPH, LENSVIEW, SELECTEDVIEW, SPOT, SHAPES_SIZE };
	cShapes shapes[SHAPES_SIZE+1];
	cShapes* current_shapes;

	cLens1 *ref_sphere;
	int make_ref_sphere(double yObj,double xObj,
	                    double yPupil_principal,double xPupil_principal,double defocus,int j);
	int fno_for_unit_pupil(double& fno_y,double& fno_x,double yObj,double xObj,double defocus,int j,int IsExitpupil);
	
	void scale_condition(double m);

	cImage image,image0;

public:
	static double LargeNumber();
	static double NewtonToR(double r0,double EAphi,double FRW_nm,double rings);
	static double RToNewton(double r0,double EAphi,double FRW_nm,double r);
	static std::string RingsSpecRefSurf(double th_deg,double VirtualAS);
	static double lens_phi(double ea_phi);
	static double lens_ea(double phi);
	static double ZValue(double r1,double r2,double phi1,double phi2);
	
	static double ToroidZ(double x,double y,double rx,double ry,double kp,int IsXToroid,int &domain_error);
	static double ToroidZ(double x,double y,double rx,double ry,double kp,int IsXToroid);
	static void   ToroidDerivative(double &zx,double &zy,
	                               double x,double y,double rx,double ry,double kp,int IsXToroid);
	static void   Toroid2ndDerivative(double &zxx,double &zxy,double &zyy,double &zx,double &zy,
	                                  double x,double y,double rx,double ry,double kp,int IsXToroid);
	static double ToroidRa(double x,double y,double rx,double ry,double kp,int IsXToroid); 
		
	static double AASZ(double x,double y,double rx,double ry,double kpx,double kpy,int &domain_error);
	static double AASZ(double x,double y,double rx,double ry,double kpx,double kpy);
	static void   AASDerivative(double &zx,double &zy,
	                            double x,double y,double rx,double ry,double kpx,double kpy);
	static void   AAS2ndDerivative(double &zxx,double &zxy,double &zyy,double &zx,double &zy,
	                               double x,double y,double rx,double ry,double kpx,double kpy);

	static void   ConicOffAxis(double &z,double &zx,double &zy,double &zxx,double &zxy,double &zyy,
	                           double x,double y,double a,double b,double t,int &domain_error);
	static double ConicOffAxisZ(double x,double y,double a,double b,double t,int &domain_error);
	static double ConicOffAxisZ(double x,double y,double a,double b,double t);
	static void   ConicOffAxisDerivative(double &zx,double &zy,
		                                 double x,double y,double a,double b,double t);
	static void   ConicOffAxis2ndDerivative(double &zxx,double &zxy,double &zyy,double &zx,double &zy,
		                                    double x,double y,double a,double b,double t);

	static void rotate_angle(double& rox,double& roy,double Zx,double Zy,double Zz);
	static void rotate_angle2(double& rox,double& roy,double& roz,double Zx,double Zy,double Zz,double Xx,double Xy);
	static void rotate_angle3(double& rox,double& roy,double& roz,double Zx,double Zy,double Zz,double Yx,double Yy);
	static matrix<double> Tmatrix(const vector<double>& new_zdirection);
	static double ACoefficient(int cOrder,double phi,double N,double al,double h);
	static double BCoefficient(int cOrder,double phi,double N,double al,double h);
	static double SACoefficient(int cOrder,double phi,double N,double al,double h);
	static double CMCoefficient(int cOrder,double phi,double N,double al,double h,double hp);
	static double ASCoefficient(int cOrder,double phi,double N,double al,double h,double hp);
	static double DSCoefficient(int cOrder,double phi,double N,double al,double h,double hp);

	static cLens1 ThinDoubletOnNud2(double beta, 
	                                std::string color1,std::string color2,std::string color3, 
							        std::string glass1, double L,double B,double A, double nud2);
	static double ThinDoubletNd2OnNud2(double beta, 
	                                   std::string color1,std::string color2,std::string color3, 
								       std::string glass1, double L,double B,double A, double nud2);
	static double ThinDoubletSA5OnNud2(double beta, 
	                                   std::string color1,std::string color2,std::string color3, 
								       std::string glass1, double L,double B,double A, double nud2);
	static double ThinDoubletLC2OnNud2(double beta, 
	                                   std::string color1,std::string color2,std::string color3, 
								       std::string glass1, double L,double B,double A, double nud2);

	static cLens1 ThinDoublet(double beta, 
	                          std::string color1,std::string color2,std::string color3, 
							  std::string glass1,std::string glass2, double L,double B);
	static double ThinDoubletSA(double beta, 
	                            std::string color1,std::string color2,std::string color3, 
								std::string glass1,std::string glass2, double L,double B);
	static double ThinDoubletSA5(double beta, 
	                             std::string color1,std::string color2,std::string color3, 
								 std::string glass1,std::string glass2, double L,double B);
	static double ThinDoubletLC2(double beta, 
	                             std::string color1,std::string color2,std::string color3, 
								 std::string glass1,std::string glass2, double L,double B);

	static cLens1 ThinSeparateDoublet(double beta, 
	                                  std::string color1,std::string color2,std::string color3, 
							          std::string glass1,std::string glass2,
									  double L,double B,double A);
	static double ThinSeparateDoubletSA5(double beta, 
	                                     std::string color1,std::string color2,std::string color3, 
							             std::string glass1,std::string glass2,
									     double L,double B,double A);

	static cLens1 ThinTriplet(int FrontIsCemented, double beta,
	                          std::string color1,std::string color2,std::string color3, 
	                          std::string glass1,std::string glass2,std::string glass3,
	                          double PhiFront, double L,double B,double A);
	static double ThinTripletSA5(int FrontIsCemented, double beta,
	                             std::string color1,std::string color2,std::string color3, 
	                             std::string glass1,std::string glass2,std::string glass3,
	                             double PhiFront, double L,double B,double A);	
	static double ThinTripletLC2(int FrontIsCemented, double beta,
	                             std::string color1,std::string color2,std::string color3, 
	                             std::string glass1,std::string glass2,std::string glass3,
	                             double PhiFront, double L,double B,double A);

	static cLens1 DoubleThinDoublet(double beta,double t,double e,
	                                std::string color1,std::string color2,std::string color3,
	                                std::string glass1,std::string glass2,
						            std::string glass3,std::string glass4,
					                double Length, double L,double T,double CM,double SA,
								    int SolutionNo);
	static double DoubleThinDoubletSA5(double beta,double t,double e,
	                                   std::string color1,std::string color2,std::string color3,
	                                   std::string glass1,std::string glass2,
						               std::string glass3,std::string glass4,
					                   double Length, double L,double T,double CM,double SA,
								       int SolutionNo);
	static double DoubleThinDoubletAS(double beta,double t,double e,
	                                  std::string color1,std::string color2,std::string color3,
	                                  std::string glass1,std::string glass2,
						              std::string glass3,std::string glass4,
					                  double Length, double L,double T,double CM,double SA,
								      int SolutionNo);
	static double DoubleThinDoubletPT(double beta,double t,double e,
	                                  std::string color1,std::string color2,std::string color3,
	                                  std::string glass1,std::string glass2,
						              std::string glass3,std::string glass4,
					                  double Length, double L,double T,double CM,double SA,
								      int SolutionNo);
	static double DoubleThinDoubletDS(double beta,double t,double e,
	                                  std::string color1,std::string color2,std::string color3,
	                                  std::string glass1,std::string glass2,
						              std::string glass3,std::string glass4,
					                  double Length, double L,double T,double CM,double SA,
								      int SolutionNo);

	static void DoubleThinlens(double& A1,double& A2,double& B2, double& f1,double& f2, double& SigmaPhi,
	                           double& beta1,double& beta2, double& DS,
	                           double B1,
	                           double Length,double eLRatio,double beta,double t,double SA,double CM,double AS);
	static double DoubleThinlensA1(double B1,
	                           double Length,double eLRatio,double beta,double t,double SA,double CM,double AS);
	static double DoubleThinlensA2(double B1,
	                           double Length,double eLRatio,double beta,double t,double SA,double CM,double AS);
	static double DoubleThinlensB2(double B1,
	                           double Length,double eLRatio,double beta,double t,double SA,double CM,double AS);
	static double DoubleThinlensDS(double B1,
	                           double Length,double eLRatio,double beta,double t,double SA,double CM,double AS);
	static double DoubleThinlensF1(double Length,double eLRatio);
	static double DoubleThinlensF2(double Length,double eLRatio);
	static double DoubleThinlensSigmaPhi(double Length,double eLRatio);
	static double DoubleThinlensBeta1(double Length,double eLRatio,double beta);
	static double DoubleThinlensBeta2(double Length,double eLRatio,double beta);

	static cLens1 Triplet(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                     std::string color1,std::string color2,std::string color3,
						 std::string glass1,std::string glass2,std::string glass3,
						 double CM,double AS,double DS, double f,double d1,double d3,double d5);
	static double TripletSA(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                        std::string color1,std::string color2,std::string color3,
						    std::string glass1,std::string glass2,std::string glass3,
						    double CM,double AS,double DS, double f,double d1,double d3,double d5);
	static double TripletSA5(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                         std::string color1,std::string color2,std::string color3,
						     std::string glass1,std::string glass2,std::string glass3,
						     double CM,double AS,double DS, double f,double d1,double d3,double d5);
	static double TripletCM(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                        std::string color1,std::string color2,std::string color3,
						    std::string glass1,std::string glass2,std::string glass3,
						    double CM,double AS,double DS, double f,double d1,double d3,double d5);
	static double TripletAS(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                        std::string color1,std::string color2,std::string color3,
						    std::string glass1,std::string glass2,std::string glass3,
						    double CM,double AS,double DS, double f,double d1,double d3,double d5);
	static double TripletDS(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                        std::string color1,std::string color2,std::string color3,
						    std::string glass1,std::string glass2,std::string glass3,
						    double CM,double AS,double DS, double f,double d1,double d3,double d5);
	static double TripletLC(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                        std::string color1,std::string color2,std::string color3,
						    std::string glass1,std::string glass2,std::string glass3,
						    double CM,double AS,double DS, double f,double d1,double d3,double d5);
	static double TripletTC(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                        std::string color1,std::string color2,std::string color3,
						    std::string glass1,std::string glass2,std::string glass3,
						    double CM,double AS,double DS, double f,double d1,double d3,double d5);
	               
	static int ZernikeL(int term_no,int IsFringeOrder,int jBase);
	static int ZernikeN(int term_no,int IsFringeOrder,int jBase);

	static double GaussBeamTruncatedPower
	              (double BeamPhiX,double BeamPhiY,double AperturePhiX,double AperturePhiY,
                   double ApertureDX,double ApertureDY,int n);

	static std::string ZoomCamTable(double fv,double mo,
	                                double th1,double m_th1,double th2,double m_th2,
	                                double th_start,double th_end,double th_step,
	                                int variator_leading);

	static void Telephoto(double &p1,double &p2,double L,double e,double f);
	static double TelephotoF1(double L,double e,double f);
	static double TelephotoF2(double L,double e,double f);

	static void OffAxialConicAB(double &a,double &b,double t_deg,double LongR,double ShortR);

	static void MartinEq(double &Sph,double &Cyl,double Diopter,double th_deg,double N,int OriginalEq=0);
	static double MartinSph(double Diopter,double th_deg,double N);
	static double MartinCyl(double Diopter,double th_deg,double N);
	static void MartinEqGeneral(double &Sph,double &Cyl,double &Axis_deg, double Sph0,double Cyl0,double Axis0_deg,double N,
		                        double rox_deg,double roy_deg,int order,int inverse=0);

	//	cLens1(){ cLens1(2,3); }; とするとエクセルでエラー
	cLens1(int i=2,int j=3);
	cLens1(const cLens1& x);
	virtual ~cLens1();
	void Reset();
	cLens1& operator=(const cLens1& x);
	friend bool operator ==(const cLens1& a,const cLens1& b);
	int push(); int pop(); int undo(); int redo();

	double var,var1,var2,var3,var4,var5;
	// これら変数は他に影響しない．自動設計の変数として用いる．
	// 例：  (incidentangle pupil var 0 "" 0 0 1) を var を変数として θ(deg) に最適化すれば，
	//       画角θに対応する物体高がvarに代入される．

	std::string Get_Note();		   void Set_Note(std::string val);
	double Get_FRW();              void Set_FRW(double value);
	int Get_k();                   void Set_k(int val);
	int Get_cn();                  void Set_cn(int val);
	std::string Get_color(int j);  void Set_color(int j, std::string s);
	double Get_colorweight(int j); void Set_colorweight(int j,double value);
	double Get_s();                void Set_s(double value);
	double Get_t();                void Set_t(double value);
	double Get_r(int i);           void Set_r(int i,double value);
	double Get_rObj();             void Set_rObj(double value);
	double Get_rImage();           void Set_rImage(double value);
	double Get_Newton(int i);      void Set_Newton(int i,double value);
	double Get_As0(int i);         void Set_As0(int i,double value);
	double Get_As45(int i);        void Set_As45(int i,double value);
	int    Get_rVariable(int i);   void Set_rVariable(int i,int value);
	int    Get_asph_type(int i);   void Set_asph_type(int i,int value);
	enum{
		SPH=0,
		CONIC=1,
		TOROID=2,
		ANAMO=3,
		IDEAL=4,
		CONIC_OA=5
	};
	double Get_kp(int i);          void Set_kp(int i,double value);
	double Get_NormH(int i);       void Set_NormH(int i,double value);
	double Get_a1(int i);          void Set_a1(int i,double value);
	double Get_a2(int i);          void Set_a2(int i,double value);
	double Get_a3(int i);          void Set_a3(int i,double value);
	double Get_a4(int i);          void Set_a4(int i,double value);
	double Get_a5(int i);          void Set_a5(int i,double value);
	double Get_a6(int i);          void Set_a6(int i,double value);
	double Get_a7(int i);          void Set_a7(int i,double value);
	double Get_a8(int i);          void Set_a8(int i,double value);
	double Get_a9(int i);          void Set_a9(int i,double value);
	double Get_a10(int i);         void Set_a10(int i,double value);
	double Get_a11(int i);         void Set_a11(int i,double value);
	double Get_a12(int i);         void Set_a12(int i,double value);
	double Get_a13(int i);         void Set_a13(int i,double value);
	double Get_a14(int i);         void Set_a14(int i,double value);
	double Get_a15(int i);         void Set_a15(int i,double value);
	double Get_a16(int i);         void Set_a16(int i,double value);
	double Get_a18(int i);         void Set_a18(int i,double value);
	double Get_a20(int i);         void Set_a20(int i,double value);
	double Get_b20(int i);         void Set_b20(int i,double value);
	double Get_b11(int i);         void Set_b11(int i,double value);
	double Get_b02(int i);         void Set_b02(int i,double value);
	double Get_b30(int i);         void Set_b30(int i,double value);
	double Get_b21(int i);         void Set_b21(int i,double value);
	double Get_b12(int i);         void Set_b12(int i,double value);
	double Get_b03(int i);         void Set_b03(int i,double value);
	double Get_b40(int i);         void Set_b40(int i,double value);
	double Get_b31(int i);         void Set_b31(int i,double value);
	double Get_b22(int i);         void Set_b22(int i,double value);
	double Get_b13(int i);         void Set_b13(int i,double value);
	double Get_b04(int i);         void Set_b04(int i,double value);
	double Get_b(int i,int m,int n); void Set_b(int i,int m,int n,double value);
	std::string Get_bTerms(int i);
	void        Set_bTerms(int i,std::string s);
	cZernike    Get_ZSurface(int i);
	void        Set_ZSurface(int i, cZernike zernike);
	std::string Get_ZCoefficients(int i);
	void        Set_ZCoefficients(int i,std::string s);
	double      Get_ZernikeR0(int i);
	void        Set_ZernikeR0(int i,double value);
	int         Get_ZernikeMaxOrder(int i);
	void        Set_ZernikeMaxOrder(int i,int value);
	std::string Get_DconTerms(int i);
	void        Set_DconTerms(int i,std::string s);
	std::string Get_LCoefficients(int i);
	void        Set_LCoefficients(int i,std::string s);
	double      Get_LegendreR0(int i);
	void        Set_LegendreR0(int i,double value);
	std::string Get_SplineData(int i);
	void        Set_SplineData(int i,std::string s);
	double Get_Apv(int i);         void Set_Apv(int i,double value);
	double Get_sfx(int i);         void Set_sfx(int i,double value);
	double Get_sfy(int i);         void Set_sfy(int i,double value);
	int    Get_UserDefSurf(int i); void Set_UserDefSurf(int i,int value);
	int    Get_cylinder(int i);    void Set_cylinder(int i,int value);
	enum{
		GENY=1,
		GENX=2,
	};
	double Get_ry(int i);          void Set_ry(int i,double value);
	double Get_rx(int i);          void Set_rx(int i,double value);
	double Get_kpy(int i);         void Set_kpy(int i,double value);
	double Get_kpx(int i);         void Set_kpx(int i,double value);
	int    Get_IsXToroid(int i);   void Set_IsXToroid(int i,int value);
	double Get_fideal(int i);      void Set_fideal(int i,double value);
	double Get_aCOA(int i);        void Set_aCOA(int i,double value);
	double Get_bCOA(int i);        void Set_bCOA(int i,double value);
	double Get_tCOA(int i);        void Set_tCOA(int i,double value);
	double Get_SA0(int i);         void Set_SA0(int i,double value);
	double Get_CM0(int i);         void Set_CM0(int i,double value);
	int Get_grating(int i);        void Set_grating(int i,int value);
	int Get_difforder(int i);      void Set_difforder(int i,int value);
	double Get_gpitch(int i);      void Set_gpitch(int i,double value);
	double Get_grx(int i);         void Set_grx(int i,double value);
	double Get_gry(int i);         void Set_gry(int i,double value);
	double Get_grz(int i);         void Set_grz(int i,double value);
	double Get_Diffusion(int i);   void Set_Diffusion(int i,double value);
	int Get_Fish(int i);           void Set_Fish(int i,int value);
	double Get_d(int i);           void Set_d(int i,double value);
	double Get_delta_d(int i);     void Set_delta_d(int i,double value);
	int    Get_dVariable(int i);   void Set_dVariable(int i,int value);
	std::string Get_gname(int i);  void Set_gname(int i,std::string s);
	int    Get_gVariable(int i);   void Set_gVariable(int i,int value);
	double Get_N(int i,int j);
	int Get_EAtype(int i);         void Set_EAtype(int i,int value);
	double Get_EAy(int i);         void Set_EAy(int i,double value);
	double Get_EAx(int i);         void Set_EAx(int i,double value);
	double Get_CHMy(int i);        void Set_CHMy(int i,double value);
	double Get_CHMx(int i);        void Set_CHMx(int i,double value);
	double Get_EAdy(int i);        void Set_EAdy(int i,double value);
	double Get_EAdx(int i);        void Set_EAdx(int i,double value);
	int Get_decenter_type(int i);  void Set_decenter_type(int i,int value);
	double Get_dx(int i);          void Set_dx(int i,double value);
	double Get_dy(int i);          void Set_dy(int i,double value);
	double Get_dz(int i);          void Set_dz(int i,double value);
	double Get_rox(int i);         void Set_rox(int i,double value);
	double Get_roy(int i);         void Set_roy(int i,double value);
	double Get_roz(int i);         void Set_roz(int i,double value);
	int    Get_order(int i);       void Set_order(int i,int value);
	int    Get_ret(int i);         void Set_ret(int i,int value);
	double Get_dx1(int i);         void Set_dx1(int i,double value);
	double Get_dy1(int i);         void Set_dy1(int i,double value);
	double Get_dz1(int i);         void Set_dz1(int i,double value);
	double Get_rox1(int i);        void Set_rox1(int i,double value);
	double Get_roy1(int i);        void Set_roy1(int i,double value);
	double Get_roz1(int i);        void Set_roz1(int i,double value);
	int    Get_order1(int i);      void Set_order1(int i,int value);
	std::string Get_CoatName(int i);  void Set_CoatName(int i,std::string filename);
	int Get_CoatReverse(int i);    void Set_CoatReverse(int i,int value);
	std::string Get_rem(int i);    void Set_rem(int i,std::string value);
	double Get_s1fix();            void Set_s1fix(double value);
	std::string Get_AfocalMode();  void Set_AfocalMode(std::string value);
	double Get_yObjectMax();       void Set_yObjectMax(double value);
	double Get_xObjectMax();       void Set_xObjectMax(double value);
	double Get_yObjectMaxAng();    void Set_yObjectMaxAng(double value);
	double Get_xObjectMaxAng();    void Set_xObjectMaxAng(double value);
	int Get_stop();                void Set_stop(int i);
	double Get_EPD();              void Set_EPD(double value);
	double Get_EPDx();             void Set_EPDx(double value);
	double Get_EPy();              void Set_EPy(double value);
	double Get_EPx();              void Set_EPx(double value);
	int Get_nSpot();               void Set_nSpot(int value);
	double Get_SourcePhiY();       void Set_SourcePhiY(double value);
	double Get_SourcePhiX();       void Set_SourcePhiX(double value);
	double Get_SourceAxis();       void Set_SourceAxis(double value);
	int Get_SourceAreaType();      void Set_SourceAreaType(int value);
	double Get_GrinDs();           void Set_GrinDs(double value);
	int Get_GrinToListStep();      void Set_GrinToListStep(int value);

	int Get_ExcludeVirtualRay();    void Set_ExcludeVirtualRay(int value);
	int Get_ExcludeVirtualObject(); void Set_ExcludeVirtualObject(int value);

	int SurfNo(std::string rem);

	int Afocal;
	int AfocalMode; enum{ MIN=1, RAD=2, FUNDUS=3, VA=4 };
	double AfocalImageUnit();
	std::string AfocalImageUnitStr();
	std::string AfocalFreqUnitStr();
	int AfocalRotateEye;  // 真のとき，ディオプター値は主光線に沿った長さの逆数となる．
	                      // 偽のとき，ディオプター値はZ軸に沿った長さの逆数となる．
	                      // (松居“収差論” (3.2.30)式などは偽の場合に相当する．）

	int StopDominate;
	double tCalc();
	double EPDCalc();
	void EPCalc();
	double Mstop(); double Mstop1();

	int    ApodizationSurf;
	double ApodizationAmplitude();
	double ApodizationIntensity();
	double GaussianPhiX,GaussianPhiY;

	double Wl(int j);

	double power(int i1,int i2,int j=1);     double power(int j=1);
	double dpower(int i1,int i2,int j=1);    double dpower(int j=1);
	double f(int i1,int i2,int j=1);         double f(int j=1);
	double bf(int i1,int i2,int j=1);        double bf(int j=1);
	double ff(int i1,int i2,int j=1);        double ff(int j=1);
	double bfRatio(int i1,int i2);           double bfRatio();
	double ffRatio(int i1,int i2);           double ffRatio();
	double BFOverF2(int i1,int i2);          double BFOverF2();
	double FFOverF2(int i1,int i2);          double FFOverF2();
	double delta(int i1,int i2,int j=1);     double delta(int j=1);
	double delta1(int i1,int i2,int j=1);    double delta1(int j=1);
	double g_hat(int j=1);
	double g1_hat(int j=1);
	double g1(int j=1);
	double nodal(int i1,int i2,int j=1);     double nodal(int j=1);
	double nodal1(int i1,int i2,int j=1);    double nodal1(int j=1);
	double cc(int i);
	double cc1(int i);
	double vertex(int i,double si=0);
	double vertex1(int i,double si=0);
	double deadspace(int i1,int i2,int j=1); double deadspace(int j=1);
	double M(double s,int i1,int i2,int j);
	double M(int i1,int i2,int j=1);         double M(int j=1);
	void SetM(double value);
	double Mpupil(int i1,int i2,int j=1);    double Mpupil(int j=1);
	double gamma(int i1,int i2,int j=1);     double gamma(int j=1);
	double si(double s,int i,int j);
	double si(int i,int j=1);
	double s1i(double s,int i,int j);
	double s1i(int i,int j=1);
	double s1(double s,int j);
	double s1(int j=1);
	double vdiopter1i(double s,int i,int j);
	double vdiopter1i(int i,int j=1);
	double LCPar(int from,int to);
	double t1(int j=1);
	double ti(int i,int j=1);
	double t1i(int i,int j=1);
	double PupilToPupil(int i1=0,int i2=0);
	double ExitPupilZ(int findpupil);
	double ExitPupilDia(int findpupil);
	double ThinT(int i1,int i2);  double ThinT();
	double ImageInfinity();
	double TelecentricityObj();
	double TotalThickness(int i1,int i2); double TotalThickness();
	double TotalInAirThickness(int i1,int i2); double TotalInAirThickness();
	double ConjugateLength();
	double TotalOpticalThickness(int i1,int i2); double TotalOpticalThickness();
	double TotalOpticalThicknessDerivative2(int i1,int i2,int j);

	vector<double> VertexGlobal(int i);
	vector<double> exGlobal(int i);
	vector<double> eyGlobal(int i);
	vector<double> ezGlobal(int i);
	vector<double> IncidentPointGlobal(int i,double yObj,double xObj,std::string SetRay,int FindPupil);
	double Clearance(int iPoint,double yObjPoint,double xObjPoint,std::string SetRayPoint,
	                 int iLine,double yObjLine,double xObjLine,std::string SetRayLine,int FindPupil);

	double GDD(int i1,int i2,int j);
	double ChirpedPulseWidth(int i1,int i2,int j,double dt0_fs);
	
	std::string SurfaceSagTable(int i,double hStep);
	double SurfaceSagMax(int i,double hStep);
	double SurfaceSagMin(int i,double hStep);
	double SurfaceSlope(int i,double y,double x);
	double SurfaceSlopeMax(int i,double hStep);
	std::string SurfaceSlopeTable(int i,double hStep);
	double koba(int i);
	double ct(int i);
	double Steepness(int i);
	std::string DefaultBC(double min_koba=1,double min_ct=0.5,double max_steepness=0.8,double weight=1,int i1=0,int i2=0);
	double DAbsC(int i1,int i2);
	double ShapeFactor(int i);
	double ZValue(int i,int i1);
	double PreformKoba(int i);
	double Inflection(int i);
	double InflectionPoint(int i);
	double AxialRadiusY(int i, double y);
	double TangentialCurvature(int i,double y);
	double TangentialCurvatureMax(int i,double hStep);

	void DzToD(int i);
	void RotateSurface(int i,double Sx,double Sy,double Sz,double Rx,double Ry,double Rz,double th);
	void RotateBlock(int i1,int i2,double Sx,double Sy,double Sz,double Rx,double Ry,double Rz,double th);
	void RotateBlockX(int i1,int i2,double th);
	void RotateBlockY(int i1,int i2,double th);
	void RotateBlockAroundPupil(int i1,int i2,int stop,double rox);
	void RotateSurfaceXYZOrder
         (int i,double rox1,double roy1,double roz1,double rox2,double roy2,double roz2,int IsMirror);
	static void RotateSurfaceZYXOrder(double& roz,double& roy,double& rox,int IsMirror);
	void RotateSurfaceZYXOrder(int i,double roz,double roy,double rox,int IsMirror);
	
	double Scan(std::string ScanXY,std::string RotateAxisXY, double th_deg, int GetData);
	void   xScan(std::string RotateAxisXY, double th_deg);
	void   xScan(double th_deg);
	double xScan(std::string RotateAxisXY);
	double xScan();
	void   yScan(std::string RotateAxisXY, double th_deg);
	void   yScan(double th_deg);
	double yScan(std::string RotateAxisXY);
	double yScan();
	
	void   zScan(double dz);
	double zScan();

	static void InverseDecenter(double& dx,double& dy,double& dz,double& rox,double& roy,double& roz,int& order);
	void TransformRetDecenter(int i);
	void TransformRetDecenter();
	void DeleteDecenter();
	
	bool IsRotationallySymmetric();
	bool IsYAxisSymmetric();
	bool IsXAxisSymmetric();
	void Scale(double m,int i1,int i2,int WithD,int WithEA);
	void Scale(double m);
	void AdjustFocalLength(double fl,int i1,int i2,int KeepD=1,int ByAllSystem=0,int WithEA=0);
	void AdjustFocalLength(double fl,int KeepD=1);
	int cBend(int i1,int i2,double dc);
	double qValue(int i);
	int  qBend(int i,double q);
	void Add2ndOrderTerm(int i, double a2);
	void ToThinLens(int i1,int i2);
	void ToIdealLens(int i1,int i2);
	double ToAplanaticSurf(int i);
	double ToConcentricSurf(int i);
	double Get_e(int pre_i1,int pre_i2,int post_i1,int post_i2);
	void   Set_e(int pre_i1,int pre_i2,int post_i1,int post_i2,double value);
	double e(int i1,int i2);
	double e1(int i1,int i2);
	void TransformACoefficients(int i,double newNormH);
	void ReduceAsphTerms(int i,int max_order,double phi=0,double h_step=1);
	double SetSpheroid(int i,char majoraxis,int convex,double rshort,double ftof);
	void   SetSpheroid2(int i,int i1,int i2);
	void   SetSpheroid3(int i,double z1,double z2);
	void SetHyperboloid(int i,double f1,double f2);
	void SetHyperboloid2(int i,int i1,int i2);
	void SetParaboloid(int i,double f);
	void SetParaboloid2(int i,int i_f,int i_inf);
	void Chamfer(int i,double Wx,double Wy,double Cx,double Cy,double Offset);

	double DconSetRn(int i,double val);
	double DconSetRn(int i);
	void DconToPSeries(int i);
	void PSeriesToDcon(int i);
	void PSeriesToFreeForm(int i);
	void DconToFreeForm(int i);
	void LegendreToFreeForm(int i);
	void SetSpline(int i,int n,double step);
	int SplineN(int i);
	void SplineDoubleN(int i);
	void SplineSetHStep(int i,double dh);

	double dc(int i,double dz);
	double dkp(int i,double dz);
	double da(int i,int n,double dz);
	
	cLens1 SwapObjPupil();
	cLens1 Deleted(int i1,int i2);
	void   Delete(int i1,int i2);
	void   DeleteBetween(int i);
	cLens1 Trimed(int i1,int i2);
	void   Trim(int i,int j);
	cLens1 Added(int pre_surf,int n,double z);
	void   Add(int pre_surf,int n,double z);
	void   ToBlock(int i1,int i2);
	void   BreakCement(int i,std::string AdhesiveName="1");
	cLens1 Reversed(int z_reverse);
	void   Reverse();
	cLens1 zReversed(int ToMirrorImg);
	void   zReverse(int ToMirrorImg);
	void Replace(int i1,int i2,cLens1 X,int i1o=0,int i2o=0,int AdjustD=1);
	int  Replace(int i1,int i2,std::string filename,int i1o=0,int i2o=0,int AdjustD=1);
	void Merge(cLens1 X,double d,int Auto=1);
	int  Merge(std::string filename, double d);
	cLens1 Bended(int i);
	void   Bend(int i);
	cLens1 Turned(int i);
	void   Turn(int i);
	cLens1 DoubleTurned(int i1,int i2);
	void   DoubleTurn(int i1,int i2);
	cLens1 TripleTurned(int i1,int i2,int i3);
	void   TripleTurn(int i1,int i2,int i3);
	int Zoom(int i1,int i2,int i3,double MorF);

	void DispersionToZero(int i1,int i2);

	void rRound(int n);
	void dRound(int n);
	void RoundA(int n);  // aRound より RoundA の方が探しやすい 2016.03.28

	virtual std::string FileName();
	//  FileName()を仮想関数とすることにより，
	//  cLens1のメンバ，cLens1::f(){   FileName(); } を，
	//  cLens1から派生したクラスXから呼んだとき，すなわちX.f()として呼ぶとき，
	//  cLens1::FileName() ではなく， cLens::FileName()を呼ぶことができる．
	//  (このとき，cLens1::filename が設定されていなくとも，
	//   cLens::filenameが設定されていれば意図する結果が得られる．)
	std::string LensData(int colors=1);
	std::string LensDataZemax();
	std::string ParaxialValues();
	std::string Coefficients();
	std::string alpha_h_table();
	std::string Focallengths(int i1,int i2,int ForDrawing);

	friend std::ostream& operator<<(std::ostream& to,cLens1& x);
	friend std::istream& operator>>(std::istream& from,cLens1& x);
	int open(std::string filename);
	int save(std::string filename);
	int SaveAsRTPFile(std::string filename);

	int Get_NormalizeType();         void Set_NormalizeType(int value);
	std::string Get_NormalizeUnit(); void Set_NormalizeUnit(std::string value);
	void CalcCoefficients(int j2=0);
	double H(int i);
	double Hp(int i);
	double U(int i);
	double Up(int i);
	double Al(int i);
	double Alp(int i);

	double Coefficient(std::string name,int sum_rms,int exponent,int total,int i1,int i2);
	double Coefficient(std::string name);
	
	matrix<double> CMatrix(int i);
	matrix<double> CMatrix(int i1,int i2);
	double CMatrix(int i1,int i2,int i,int j);
	std::string CMatrixStr(int i1,int i2);

	double SA0(int i1,int i2); double SA0();
	double CM0(int i1,int i2); double CM0();

	double cDerivative(std::string command,int i,double dc);
	double dDerivative(std::string command,int i,double dd);

	void ToScheimpflugImagePlane();
	double ScheimpflugImagePlaneTiltX();
	double ScheimpflugImagePlaneTiltY();
	double ImagePlaneTiltX(int i,double yObj,double xObj,
	                       std::string SetRay,int FindPupil,double yPupil,double xPupil,double dyObj=0);
	double ImagePlaneTiltX(int i);

	double OffAxialMirrorComa(int i);

	int RayTrace(double yObj,double xObj,double yPupil,double xPupil,
		         double defocus,int j,int EA_enable,int mask_enable,int IsLambert,int as_trace,
				 int E_trace=0,const vector<complex>& E0=vector<complex>(), int MakeCoordinate=1);
	enum {
		NOT_THRU      =10000,
		NOT_INTERSECT =20000,
		CLIPPED       =30000,
		TOTAL_REF     =40000,
		INVALID       =50000,
		GRIN_OUT      =60000
	};
	int NotThruSurf(int errorcode);
	std::string ErrorMessage(int errorcode);
	int ExcludeVirtualRay;
	int ExcludeVirtualObject;
	
	int FindAThruRay(double& yPupil, double& xPupil,double yObj,double xObj,
	                 int j,int priority,int makecoordinate);
	double FindAThruRayTimeOut;
	int FindMarginalRay(double& yPupil,double& xPupil,double azimuth,
	                    double yObj,double xObj,double yPupilIni,double xPupilIni,int j,int makecoordinate);
	int FindMarginalRay2(double& yPupil,double& xPupil,int direction,
	                     double yObj,double xObj,double yPupilIni,double xPupilIni,int j,int makecoordinate);
	int FindPupilEdge(point& ymax,point& ymin, point& xmax,point& xmin,
	                  int& ymax_i,int& ymin_i, int& xmax_i,int& xmin_i,
	                  double yObj,double xObj, int j,int findpupil,int makecoordinate);
	int FindPupilEdge(point& ymax,point& ymin, point& xmax,point& xmin,
	                  double yObj,double xObj, int j,int findpupil,int makecoordinate);
	int FindPupilEdge(point& ymax,point& ymin, point& xmax,point& xmin,
     	              double yObj,double xObj, int j,int makecoordinate);
	int FindRay(point& p,double yObj,double xObj,std::string kind,int findpupil,int j);
	double PupilMargin;

	double OptimizedDefocus(double yObj,double xObj,int FindPupil,int j);
	double OptimizedDefocus(int FindPupil);
	double OptimizeS1fix(double yObj,double xObj,int FindPupil,int j);
	double OptimizeS1fix(int FindPupil);

	std::string SingleRayTrace(double yObj,double xObj,
	                           std::string SetRay,int findpupil,double yPupil,double xPupil,
		                       double defocus,int lastsurf,int j,
							   double pol_phi=90,double pol_ratio=0);
	std::string UsedRange(char y_or_x,double HObjMax,double HObjMin);
	std::string yUsedRange(double yObjMax,double yObjMin);
	std::string xUsedRange(double xObjMax,double xObjMin);

	vector<double> RayPos(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j);
	int RayPos(double* x,double* y,double* z,
		       int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j);
	double RayPosY(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j);
	double RayPosX(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j);
	double RayPosZ(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j);
	point RayHeight(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j);	
	
	void Footprint(double &ymax,double &ymin,double &xmax,double &xmin,
	              double &ydia,double &xdia,double &ycenter,double &xcenter,
	              double yObj,double xObj,int i,int j);
	double FootprintYdia(double yObj,double xObj,int i,int j);
	double FootprintXdia(double yObj,double xObj,int i,int j);

	double IncidentAngle(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                     int j,int InAir=0);
	double ExitAngle(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                 int j,int InAir=0);
	double IncidentAngleMax(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j);
	double IncidentAngleYFanMin(int i,double yObj,double xObj,int FindPupil,int j);
	double IncidentAnglePrinMin(int i,double yObj1,double yObj2,int FindPupil,int j);
	
	double IncidentAngleY(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                      int j,int InAir=0);
	double IncidentAngleX(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                      int j,int InAir=0);

	vector<double> IncidentDirectionCosine(int i,double yObj,double xObj,
				   std::string SetRay,int FindPupil,double yPupil,double xPupil,int j);
	vector<double> ExitDirectionCosine(int i,double yObj,double xObj,
				   std::string SetRay,int FindPupil,double yPupil,double xPupil,int j);
	double ExitDirectionCosineX(int i,double yObj,double xObj,
		                        std::string SetRay,int FindPupil,double yPupil,double xPupil,int j);
	double ExitDirectionCosineY(int i,double yObj,double xObj,
		                        std::string SetRay,int FindPupil,double yPupil,double xPupil,int j);
	double ExitDirectionCosineZ(int i,double yObj,double xObj,
	                            std::string SetRay,int FindPupil,double yPupil,double xPupil,int j);
	double ExitTanX(int i,double yObj,double xObj,
	                std::string SetRay,int FindPupil,double yPupil,double xPupil,int j);
	double ExitTanY(int i,double yObj,double xObj,
	                std::string SetRay,int FindPupil,double yPupil,double xPupil,int j);

	double DWDSag(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil);
	double DSagDFringe(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil);

	vector<complex> IncidentE(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j,
		                      double pol_phi=90,double pol_ratio=0);
	double IncidentAmplitude(std::string xyz,
	                         int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j,
		                     double pol_phi=90,double pol_ratio=0);
	double IncidentPhase(std::string xyz,
	                     int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j,
		                 double pol_phi=90,double pol_ratio=0);

	double IncidentIntensity(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j,
		                     double pol_phi=90,double pol_ratio=0);
	
	double FNumber(double yObj,double xObj,int findpupil);
	double FNumberXYAve(double yObj,double xObj,int findpupil);
	void   FNumberObj(double &FNumY,double &FNumX,double yObj,double xObj,int findpupil);
	double FNumberObj(double yObj,double xObj,int findpupil);
	double FNumberObjXYAve(double yObj,double xObj,int findpupil);
	double FNumberParaxial(int findpupil);
	double NA(double yObj,double xObj,int findpupil);
	double NAXYAve(double yObj,double xObj,int findpupil);
	void   NAObj(double &NAObjY,double &NAObjX,double yObj,double xObj,int findpupil);
	double NAObj(double yObj,double xObj,int findpupil);
	double NAObjXYAve(double yObj,double xObj,int findpupil);
	double NAParaxial(int findpupil);
	void   DiffractObj(double yObj,double xObj,int findpupil);
	void   DiffractObj(double yObj,double xObj,int findpupil,int nSpot,int DeletFrom);
	double FresnelNumber(int findpupil,double defocus);

	double AngleOfView(double h);

	double LSA(double yPupil,double xPupil,double yPupilPrincipal,double xPupilPrincipal,int j,int AbeTheory=0);
	double LSAAbe(double yPupil,double xPupil,double yPupilPrincipal,double xPupilPrincipal,int j,int order=3);
	double LSA(double yPupilNormalized,int findpupil,int j=1);
	double LSA();
	double LSA70();
	double LSA50();
	double LSA2nd(); 
	double LSA3rd();
	double DLSA2to3();
	double LSAp(double yObjNormalized,int findfield,int j=1);
	double LSAp();
	double LSAp70();
	double LSAp50();
	double LSAp2nd(double yObjNormalized,int findfield);
	double LSAp2nd();
	std::string LSAs(int findpupil,int colors,double weight=1);
	std::string LSAps(int findfield,int colors,double weight=1);

	double SC(double yPupil,double xPupil);
	double SC(double yPupilNormalized,int findpupil);
	double SC();
	double SC70();
	double SC50();
	double OSC(double yPupilNormalized,int findpupil);
	double OSC();
	double OSC70();
	double OSC50();
	double OSCp(double yObjNormalized,int findfield);
	double OSCp();

	double DeltaM(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j,int AbeTheory=0);
	double DeltaMAbe(double yObj,double xObj,int order=3);
	double DeltaM(double yObjNormalized,int findpupil,int j=1);
	double DeltaM();
	double DeltaM70();
	double DeltaM50();
	double DeltaS(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j,int AbeTheory=0);
	double DeltaSAbe(double yObj,double xObj,int order=3);
	double DeltaS(double yObjNormalized,int findpupil,int j=1);
	double DeltaS();
	double DeltaS70();
	double DeltaS50();
	double DSToDM(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j);
	double DSToDM(double yObjNormalized,int findpupil);
	double DSToDM();
	double DSToDM70();
	double DSToDM50();
	std::string DSDMs(int findpupil,double s_weight=1,double m_weight=1,int cols=1,int fine_pitch=0);

	int SCA(double& Sph,double& Cyl,double& Axis, 
	        double yObj,double xObj,double yPupil,double xPupil,double defocus,int i,int j);
	double Sph(double yObj,double xObj,double defocus,int findpupil=0,int i=0);
	double Sph();

	double Cyl(double yObj,double xObj,int findpupil=0,int i=0);
	double Cyl();
	double Cyl90();
	double Cyl70();
	double Cyl50();
	double CylSgn(double yObj,double xObj,int findpupil=0,int i=0);
	double Cyl0Deg(double yObj,double xObj,int findpupil=0,int i=0);
	double Cyl45Deg(double yObj,double xObj,int findpupil=0,int i=0);

	double Axis(double yObj,double xObj,int findpupil=0,int i=0);
	double Axis();
	double SphEq(double yObj,double xObj,double defocus,int findpupil=0,int i=0);
	double SphEq();
	double SphEq70();
	double SphEq50();
	double Mskew(char xy,double yObj,double xObj,double yPupil,double xPupil,int i1,int i2,int j,
	             int ObjImgTilt=1);
	double My(double yObj,double xObj,double yPupil,double xPupil,int i1,int i2,int j);
	double My(int i1,int i2);
	double My();
	double Mx(double yObj,double xObj,double yPupil,double xPupil,int i1,int i2,int j);
	double Mx(int i1,int i2);
	double Mx();
	double Mpskew(char xy,int i1,int i2,int j);
	double Mpy(int i1,int i2);
	double Mpy();
	double Mpx(int i1,int i2);
	double Mpx();
	double fskew(char xy,int j);
	double fx();
	double fy();
	double S1skew(char xy,double yObj,double xObj,double yPupil,double xPupil,int i,int j);
	double S1y(double yObj,double xObj,double yPupil,double xPupil,int i,int j);
	double S1y(int i);
	double S1y();
	double S1x(double yObj,double xObj,double yPupil,double xPupil,int i,int j);
	double S1x(int i);
	double S1x();
	double Dist(double yObj,double yPupil,double xPupil,int j);
	double Dist(double yObjNormalized,int findpupil);
	double Dist(); double Dist70();
	double DeltaH(int DYIs1DXIs2,
	              double yObj,double xObj,double yPupil,double xPupil,
	              double yPupilPrincipal,double xPupilPrincipal,double defocus,int j,int j0,int AbeTheory=0);
	double DeltaH(int DYIs1DXIs2,
	              double yObj,double xObj,double yPupil,double xPupil,
	              double yPupilPrincipal,double xPupilPrincipal,double defocus,int j,int j0,
				  int i1,int i2,int AbeTheory=0);
	double DeltaH(int DYIs1DXIs2,double yObj,double xObj,
				  int FindPupil,double ypNormalized,double xpNormalized,double defocus,int j,int j0);
	std::string SPO(double yObj,double xObj,int findpupil,int colors,int IgnoreTC=0,double weight=1,int yfan=0);
	std::string SPO2(double yObj,double xObj,int findpupil,int colors,int IgnoreTC=0,int n=0,double weight=1);
	double DYpr(double yObj,int j,int FindPupil);
	double DYpr(int j,int FindPupil);
	double DYpr75(int j,int FindPupil);
	double DYpr50(int j,int FindPupil);
	double DYymax(double yObj,int FindPupil=1); double DYymax(); 
	double DYymin(double yObj,int FindPupil=1); double DYymin();
	double DXxmax(double yObj); double DXxmax(); 
	double DYunsymmetric(double yObj,int FindPupil);
	
	int ImageHeight(double& yh,double& xh,double& zh,double yObj,double xObj,
	                std::string SetRay,double yPupil,double xPupil,
		            double defocus,int j,int EA_enable,int mask_enable,int IsLambert,int InAngle,
					int MakeCoordinate=1);
	int ImageHeightAbe(double& yh,double& xh,double yObj,double xObj,
	                   std::string SetRay,double yp, double xp,double defocus,int j,int order=3);
	int ImageHeight(point& height,const point& obj,const point& pupil,
	                double defocus,int j,int EA_enable,int mask_enable,int InAngle);
	double ImageHeightX(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                    double defocus,int j,int InAngle);
	double ImageHeightY(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                    double defocus,int j,int InAngle);
	double ImageHeightTh(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                     double defocus,int j);

	double DImageHeight(double yObj,double xObj,
	                std::string SetRay,double yPupil,double xPupil,
		            double defocus,int j,int EA_enable,int mask_enable,int IsLambert,int InAngle,
	                double DyObj,double DxObj);

	double Get_y(int i); double Get_x(int i); double Get_z(int i);
	
	int IgnoreTC;

	cSpot spot;
	point dObject();
	int MakeSpot(cSpot& spot,
		         double yObj,double xObj,int FindPupil,double defocus,int footprint,
				 int ColorStart,int ColorEnd,
	             int IsAreaSource,int IsLambert,int OriginIsGravityCenter,int Add,int Randomize,
				 int AutoInTangent=1);
	int MakeSpot(double yObj,double xObj,int FindPupil,double defocus,int footprint,
	             int ColorStart,int ColorEnd,
	             int IsAreaSource,int IsLambert,int OriginIsGravityCenter,int Add,
				 int AutoInTangent=1);
	double SpotTotalIntensity();
	double SpotIntensity(double y,double x,double SensorPhi);
	double SpotXIntensity(double x,double SensorFW);
	double SpotYIntensity(double y,double SensorFW);
	double SpotMTF(double nu_y,double nu_x);
	double SpotRmsPhi();
	double SpotXMax();
	double SpotYMax();
	double SpotXMin();
	double SpotYMin();
	double SpotXWidth();
	double SpotYWidth();
	double SpotXGravityCenter();
	double SpotYGravityCenter();
	double SpotInertiaPrincipalAxis();
	double SpotPrincipalWx();
	double SpotPrincipalWy();
	double SpotArea();
	int AddRayToSpot(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                 double defocus,int footprint,int color);
	void ClearSpot();
	
	double rmsphi(double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,int FindPupil,
	              int OptimizeDefocus);
	double rmsphiOn(int OptimizeDefocus);
	double rmsphiOff(int OptimizeDefocus);

	int MakeRmsMap(double yObjMax,double xObjMax,int IsRectObj,double defocus,
	               int ColorStart,int ColorEnd,int FindPupil);
	double RmsMap(int i,int j);

	double EncircledEnergy(double SensorPhi,double yObj,double xObj,double defocus,
		                   int ColorStart,int ColorEnd,int FindPupil,int OptimizeDefocus);
	
	complex OTF(double nu_y,double nu_x,double yObj,double xObj,double defocus,
				int IsGeometrical, int AberrationFree, int ColorStart,int ColorEnd, int FindPupil);
	complex OTFs(double nu,double yObj,double xObj,double defocus,
				 int IsGeometrical, int AberrationFree, int ColorStart,int ColorEnd, int FindPupil);
	complex OTFm(double nu,double yObj,double xObj,double defocus,
				 int IsGeometrical, int AberrationFree, int ColorStart,int ColorEnd, int FindPupil);
	
	double MTFy(double nu,double yObj,double xObj,double defocus,
	            int IsGeometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil=0);
	double MTFy(double nu,double yObj,double xObj);
	double MTFx(double nu,double yObj,double xObj,double defocus,
	            int IsGeometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil=0);
	double MTFx(double nu,double yObj,double xObj);
	double MTFxyave(double nu,double yObj,double xObj,double defocus,
	                int IsGeometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil=0);
	double MTFxyave(double nu,double yObj,double xObj);

	double MTFs(double nu,double yObj,double xObj,double defocus,
	            int IsGeometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil=0);
	double MTFs(double nu,double yObj,double xObj);
	double MTFm(double nu,double yObj,double xObj,double defocus,
	            int IsGeometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil=0);
	double MTFm(double nu,double yObj,double xObj);
	double MTFsmave(double nu,double yObj,double xObj,double defocus,
	                int IsGeometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil=0);
	double MTFsmave(double nu,double yObj,double xObj);

	double ResPower(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                int ColorStart,int ColorEnd,int FindPupil,double NuStepY,double NuStepX);
	double ResPowerY(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                 int ColorStart,int ColorEnd,int FindPupil,double NuStep);
	double ResPowerX(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                 int ColorStart,int ColorEnd,int FindPupil,double NuStep);
	double ResPowerS(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                 int ColorStart,int ColorEnd,int FindPupil,double NuStep);
	double ResPowerM(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                 int ColorStart,int ColorEnd,int FindPupil,double NuStep);

	double PTFy(double nu,double yObj,double xObj,double defocus,
	            int IsGeometrical, int ColorStart,int ColorEnd,int FindPupil=0);
	double PTFx(double nu,double yObj,double xObj,double defocus,
	            int IsGeometrical, int ColorStart,int ColorEnd,int FindPupil=0);

	list<double> nu_list;
	void ObjectListAdd(double y, double x);
	void ObjectListClear();
	int ObjectListSize();
	double ObjectListY(int i);
	double ObjectListX(int i);
	void SetDefaultObjects(double h1,double h2,double h3);
	double MinMTFAmongObjects(double nu,double defocus,int IsGeometrical=1);
	double MinMTF(double nu);
	double MTFArea(double defocus,int IsGeometrical);
	double MaximizeMTFArea(double defocus_step,int IsGeometrical);

	double OPL(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	           double defocus,int j,int i1=0,int i2=0,int makecoordinate=1);
	int OPD_calc(double& optical_path, int ref_sphere_make,double yObj,double xObj,double yPupil,double xPupil,
	             double yPupil_principal,double xPupil_principal,double defocus,int j0,int j);
	int OPD_calc(double& optical_path, double yObj,double xObj,double yPupil,double xPupil,
	             double yPupil_principal,double xPupil_principal,double defocus,int j0,int j);
	double OPD(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	           double yPupil_principal,double xPupil_principal,double defocus,int j0,int j);
	double OPD(double yObj,double xObj,
               int FindPupil,double ypNormalized,double xpNormalized,double defocus,int j0,int j);
	std::string OPDs(double yObj,double xObj,int findpupil,int colors,int IgnoreTC=0,int n=0,double weight=1);
	double OPD2(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	            double yPupil_principal,double xPupil_principal,double defocus,double wl_nm0,double wl_nm);
	std::string DispersionTable(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                            double yPupil_principal,double xPupil_principal,
								double wl_start,double wl_end,int wl_points);

	int MakeOPDMap(matrix<double>& X,matrix<double>& Y,
	               matrix<double>& W,matrix<double>& P,matrix<double>& A,matrix<int>& E,
				   cZernike& zernike,
	               double yObj,double xObj,double defocus,int j0,int j,int FindPupil,int InLambda,
				   int ZernikeOrder,int AdjustSph,int ZernikePupil,double ShearDx,double ShearDy,
				   std::string ZernikeFix);
	int MakeOPDMap(std::string filename,double yObj,double xObj,double defocus,int j0,int j,int FindPupil,int InLambda,
	               int AdjustSph,int ZernikePupil);
	int MakeOPDMap(double yObj,double xObj,double defocus,int j0,int j,int FindPupil,int InLambda,
	               int ZernikeOrder,int AdjustSph,int ZernikePupil,double ShearDx,double ShearDy,
				   std::string ZernikeFix);
	double OPDMap(int i,int j);
	double OPDMapRMS();
	double OPDMapPV();
	double OPDMapStrehlDef();
	double OPDMapZernikeR0();
	double GetOPDMapZernikeR0Fix();
	void   SetOPDMapZernikeR0Fix(double value);
	int  GetOPDMapZernikeNormalize();
	void SetOPDMapZernikeNormalize(int value);
	int  GetOPDMapZernikeIsFringeOrder();
	void SetOPDMapZernikeIsFringeOrder(int value);
	int  GetOPDMapZernikeJBase();
	void SetOPDMapZernikeJBase(int value);
	int OPDMapZernikeTotalTerms();
	double OPDMapZernikeCoefficient(int term_no);
	double OPDMapZernikeRMSError();
	double OPDMapZernikePVError();
	double OPDMapZernikeRMS(int tilt,int sph,int cyl,int high_order);
	int CopyOPDMapToImage(cImage& image,int PixelWidth);
	int SaveOPDMapAsBmp(std::string filename,int PixelWidth);
	
	double OPDRMS(double yObj,double xObj,double defocus,int j,int FindPupil,int InLambda,int AdjustSph,
		          int OptimizeDefocusOnAxis,std::string ZernikeFix);
	double OPDRMS0(double yObj,double xObj,int FindPupil);
	double OPDPV (double yObj,double xObj,double defocus,int j,int FindPupil,int InLambda,int AdjustSph,
	              int OptimizeDefocusOnAxis,std::string ZernikeFix);
	double StrehlDef(double yObj,double xObj,double defocus,int j,int FindPupil,int AdjustSph,
	                 int OptimizeDefocusOnAxis,std::string ZernikeFi);
	double StrehlDef(double yObj,double xObj);
	double ZernikeC(int term_no,double yObj,double xObj,double defocus,int j,int FindPupil,int AdjustSph,
	                int InLambda,int OptimizeDefocus);
	double ZernikeCHigh(double yObj,double xObj,double defocus,int j,int FindPupil,int AdjustSph,
	                    int InLambda,int OptimizeDefocus);

	int MakePSFMap(cImage& B,double& wx,double& wy,double& phix,double& phiy,double threshold,
	               double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,
				   int FindPupil,int StrehlNormalize,double Zoom,int AberrationFree,int AdjustSph,
				   std::string ZernikeFix);
	int MakePSFMap(double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,
	               int FindPupil,int StrehlNormalize,double Zoom,int AberrationFree,int AdjustSph,
				   std::string ZernikeFix);
	double PSFMap(int i,int j);
	double PSFMapWX();
	double PSFMapWY();
	double PSFMapPhiX();
	double PSFMapPhiY();
	double& PSFMapPhiThreshold();
	int CopyPSFMapToImage(cImage& image,int PixelWidth,double Width);
	int CopyPSFMapToImage(int PixelWidth,double Width);
	int SavePSFMapAsBmp(std::string filename,int PixelWidth,double Width,double gamma);
	double PSFMapFiberEfficiency(double MFD,int DoublePass,int Incoherent);
	int PSFMapToShadow(double MaskArea);
	int PSFMapAddAmplitude(double amp,double phase_deg);

	int MakeOTFMap(cImage& B,double& wx,double& wy,double& rayleigh_nu_x,double& rayleigh_nu_y,
	               double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,
				   int FindPupil,int AberrationFree,int AdjustSph,std::string ZernikeFix);
	int MakeOTFMap(double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,
	               int FindPupil,int AberrationFree,int AdjustSph,std::string ZernikeFix);
	double OTFMapMTF(int i,int j);
	double OTFMapWX();
	double OTFMapWY();
	double OTFMapRayleighNuX();
	double OTFMapRayleighNuY();

	double yVignetting(double yObj,int findpupil);
	double yVignetting(int findpupil=1);
	double AbyA0();
	double Transmittance(double yObject,double xObject,int IsAreaSource,int IsLambert,
	                     int MakeSurfaceList,list<int>& surfaces);
	double Transmittance(double yObject,double xObject,int IsAreaSource,int IsLambert);
	double Ton(); double Toff();
	std::string NotThruSurfaces(double yObject,double xObject,int IsAreaSource);
	double tLack();
	cLens1 perturbed(int IsEndNotUni=0);
	void perturbe_this(int IsEndNotUni=0);
	double PerturbDxDy(int i,double dr=0);
	double PerturbXYObjectMax(double dr=0);
	void AddGroup(int i1,int i2);
	void GroupsClear();
	int RowsNumber();
	int i1Row(int i);
	int i2Row(int i);
	std::string strRow(int i);

	void MakeThinDoubletOnNud2(double beta,
	                           std::string color1,std::string color2,std::string color3, 
	                           std::string glass1,
	                           double L,double B,double A, double nud2, double f);
	void MakeThinDoublet(double beta,
	                     std::string color1,std::string color2,std::string color3, 
						 std::string glass1,std::string glass2,
	                     double L,double B, double f);
	void MakeThinSeparateDoublet(double beta,
	                             std::string color1,std::string color2,std::string color3, 
						         std::string glass1,std::string glass2,
	                             double L,double B,double A, double f);
	void MakeThinTriplet(int FrontIsCemented, double beta,
	                     std::string color1,std::string color2,std::string color3, 
	                     std::string glass1,std::string glass2,std::string glass3,
	                     double Phi3Phi1Ratio, double L,double B,double A, double f);
	void MakeDoubleThinDoublet(double beta,double t,double e,
	                          std::string color1,std::string color2,std::string color3,
	                          std::string glass1,std::string glass2,
				              std::string glass3,std::string glass4,
			                  double Length, double L,double T,double CM,double SA,
						      int SolutionNo, double f,double d1,double d2,double d4,double d5);
	void MakeTriplet(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	             std::string color1,std::string color2,std::string color3,
	             std::string glass1,std::string glass2,std::string glass3,
                 double CM,double AS,double DS, double f,double d1,double d3,double d5);
	
	static double BeamFullDiv(double waist_dia,double wl_nm,double N,double M2);
	static double WaistDia(double beam_full_div,double wl_nm,double N,double M2);
	static double RaylaySemiDOF(double waist_dia,double wl_nm,double N,double M2);
	double GetM2(); void SetM2(double value); double M2;
	double WaistPosIn;
	double GetWaistDiaIn();    void SetWaistDiaIn(double value);    double waist_dia_in;
	double GetBeamFullDivIn(); void SetBeamFullDivIn(double value); double beam_full_div_in; 
	double RaylaySemiDOFIn();
	complex q0();
	complex q1();
	double WaistPosOut();
	double WaistDiaOut();
	double BeamDiaIn();
	double BeamDiaOut(double defocus);
	double BeamFullDivOut();
	double RaylaySemiDOFOut();

	int PolarizationTrace(double yObj,double xObj,
	                      std::string SetRay,int findpupil,double yPupil,double xPupil,
                          int j, double a,double b,double phi_deg);
	vector<complex> GetE(int i); double Ex(int i); double Ey(int i); double Edelta(int i);
	vector<complex> GetE1(int i);
	void Ellipse(int i,double& A,double& B,double& Phi_deg);
	int  Ellipse(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
	             int i,int j,double a,double b,double phi_deg,double& A,double& B,double& Phi_deg);
	double EllipseA(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
                    int i,int j, double a,double b,double phi_deg);
	double EllipseB(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
                    int i,int j, double a,double b,double phi_deg);
	double EllipseRatio(int i);
	double EllipseRatio(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
                        int i,int j, double a,double b,double phi_deg);
	double EllipseIntensityRatio(int i);
	double EllipsePhi(int i);
	double EllipsePhi(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
                      int i,int j, double a,double b,double phi_deg);

	cShapes* CurrentShapes(int n);
	void ShapesClear();
	void ShapesViewTopLeft(double x,double y);
	void ShapesViewWidth(double xw,double yw);
	int ShapesSize() const;
	int ShapeKind(int i);
	double ShapeX1(int i); double ShapeY1(int i);
	double ShapeX2(int i); double ShapeY2(int i);
	int ShapeLineStyle(int i);
	std::string ShapeText(int i);
	double ShapeFontSize(int i);
	long ShapeColor(int i);
	void ShapesZoom(double m);
	void ShapesTranslate(double dx,double dy);
	void ShapesAllView();
	int ShapesSaveAsBmp(std::string filename,int yPixels);

	int SAGraph();
	void MakeSAGraph(int FindPupil,double yPupilMax,double FullScale,int PupilSA,int AbeTheory=0);
	void MakeSAGraph(int FindPupil,double yPupilMax,double FullScale);
	int ASGraph();
	void MakeASGraph(double yObjMax,int FindPupil,double FullScale,double defocus,int colors=1,int AbeTheory=0,
		             int Scan=0,double yScanThMax=0,double yScanThMin=0);
	int DistGraph();
	void MakeDistGraph(double yObjMax,int FindPupil,double FullScale);
	int DistChart();
	std::string MakeDistChart(double yObjMax,double xObjMax,int DIVy,int DIVx,int FindPupil,double defocus,
	                          int Scan,
							  std::string yScanAxisXY,double yScanThMin,double yScanThMax,double yScanRatio,
							  std::string xScanAxisXY,double xScanThMin,double xScanThMax,double xScanRatio);
	std::string MakeDistChart(double yObjMax,double xObjMax,int DIVy,int DIVx,int FindPupil,double defocus);
	std::string MakeDistChart(double yObjMax,double xObjMax,int FindPupil,double defocus);
	std::string MakeScanDistChart(std::string yScanAxisXY,double yScanThMin,double yScanThMax,double yScanRatio,
		                          std::string xScanAxisXY,double xScanThMin,double xScanThMax,double xScanRatio,
							      int DIVy,int DIVx,int FindPupil,double defocus);
	std::string ExpandDist(double yObjMax,double xObjMax,int FindPupil,double defocus);
	int DeltaHGraph();
	void MakeDeltaHGraph(double yObjMax,double xObjMax,int FindPupil,double PupilMax,
	                     double defocus,int BothSide,double FullScale,std::string command,int i1,int i2,int AbeTheory=0);
	int LensView();
	void MakeLensView(std::string command,int FindPupil);
	int SelectedView();
	void MakeSelectedView(int WhichDraw,int FindPupil,double FullScale,std::string DHCommand,std::string LVCommand);
	int SpotDiagram();
	void MakeSpotDiagram(cShapes& shapes,cSpot& spot,double FullScale,int WithFrame);
	void MakeSpotDiagram(double yObj,double xObj,int FindPupil,double defocus,int footprint,
	                     int ColorStart,int ColorEnd,
	                     int IsAreaSource,int IsLambert,double OriginIsGravityCenter,int Add,
	                     double FullScale,
						 int defocusN,double defocusStep,std::string commands,
						 int FieldDiagram,int DIVy,int DIVx,int UseObjectList,int Scan=0);
	void MakeSpotDiagram(double yObj,double xObj,int FindPupil,double defocus,int footprint,
	                     int ColorStart,int ColorEnd,
	                     int IsAreaSource,int IsLambert,double OriginIsGravityCenter,int Add,
	                     double FullScale,
		                 int defocusN,double defocusStep,int UseObjectList);
	void MakeSpotDiagram2(double yObj,double xObj,int FindPupil,double defocus,int footprint,
	                      int ColorStart,int ColorEnd,
	                      int IsAreaSource,int IsLambert,double OriginIsGravityCenter,
	                      double FullScale,
						  int defocusN,double defocusStep,std::string commands,int UseObjectList);
	void MakeSpotDiagram3(int objects,double *yObj,double *xObj,int FindPupil,double defocus,int footprint,
	                      int ColorStart,int ColorEnd,
	                      int IsAreaSource,int IsLambert,double OriginIsGravityCenter,int Add,
	                      double FullScale,
		                  int defocusN,double defocusStep);
	void MakeSpotFieldDiagram(double yObj,double xObj,int DIVy,int DIVx,
	                          int FindPupil,double defocus,int footprint,int ColorStart,int ColorEnd,
							  int IsAreaSource,int IsLambert,double OriginIsGravityCenter,double FullScale,
							  std::string commands,int Scan);
	int AddSpotToImage(cImage& image,double weight);
	int AddSpotToImage(double weight);
	int CopySpotToImage(int xPixels,int yPixels,double xWidth);
	void Image0ToImage(int FindPupil,double defocus,int ColorStart,int ColorEnd,
	                   double ImageWidth=0,int ImageXPixels=0,int ImageYPixels=0,int RemoveMoire=0);
	void Image0ToImageScan(int FindPupil,double defocus,int ColorStart,int ColorEnd,
	                       double xScanStart,double xScanEnd,double yScanStart,double yScanEnd,
						   double xImageStart,double xImageEnd,double yImageStart,double yImageEnd,
	                       int xScanPoints,int yScanPoints,int RemoveMoire=0);

	double GetImage0Width(); void SetImage0Width(double value);
	double GetImageWidth();  void SetImageWidth(double value);
	int GetImageXPixels(); void SetImageXPixels(int value); 
	int GetImageYPixels(); void SetImageYPixels(int value);
	double GetImagePitch(); void SetImagePitch(double value);
	void ImageClear();
	int ImageAdd(double x,double y,double weight);
	void ImageNormalize(double new_max);
	void ImageReverseXY();
	double GetImageIntensity(int i,int j); void SetImageIntensity(int i,int j,double value);
	int LoadImage0(std::string filename);
	int LoadImage0FromBmp(std::string filename);
	int LoadImage(std::string filename);
	int LoadImageFromBmp(std::string filename);
	int SaveImage(std::string filename);
	int SaveImageAsBmp(std::string filename,double gamma);
	
	double EyeRefractivePower(double y,double VD,double SAm,int UseCorneaHeight,double PupilLocation);
	void EyeToAxial(double diopter,double VD);
	void EyeToRefractive(double diopter,double VD);
	double EyeZernikeCoefficient(double PupilZ,double RefSurfZ,double PupilPhi,int j);
	static double ModelEyeS1(double diopter,double VD);
	static double ModelEyeFundusR(double diopter,double VD);

	static cLens1 LMTestLens(double sph,double cyl);
	void ToLMTestLens(double sph,double cyl);
	double LensCenterPowerWorn(std::string SCA,double So,double Co,double Ao,
	                           double WrapAngle,double TiltAngle);

	bool is_numeric(std::string &s);
	int  IsCoefficient(const std::string& com);
	std::string scmd(std::string com,int val);
	double* property(std::string& s,int& args);
	double* property(std::string& s);
};

#endif // !defined(AFX_CLENS1_H__B050DAA9_7467_4AA8_8DFC_82CABBAF2484__INCLUDED_)
