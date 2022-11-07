#if !defined(AFX_CIMAGE_H__16186CDD_A7A6_4A2E_ADB2_70C8B7610251__INCLUDED_)
#define AFX_CIMAGE_H__16186CDD_A7A6_4A2E_ADB2_70C8B7610251__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#define _CRT_SECURE_NO_WARNINGS
#include "matrix.h"
#include "complex.h"
#include "cXYList.h"
#include "general_func.h"
#include "cBitmap.h"
#include "cFitting.h"

class cImage
{
	int ypixels,xpixels;
	double ypitch,xpitch;
	matrix<complex> A;
	void gravity_center(int& ig,int& jg);
	static cImage buf;

	enum{ MAX_ORDER=10 };  // const int MAX_ORDER=10; とするとコンパイルエラー．
	                       // "=10"を削除し，値はソースファイルで設定することになる．
	double Cx[MAX_ORDER+1][MAX_ORDER+1],Cy[MAX_ORDER+1][MAX_ORDER+1];
	// new, delete で動的にメモリを取る方法もあるが，
	// 代入演算子，コピーコンストラクタ等が必要となり面倒．
public:
	cImage();
	virtual ~cImage();
	cImage& ToSquare();
	cImage& ToSqrt();
	cImage& ToAbs();
	void SetData(const matrix<complex>& a);
	int  GetYpixels() const;  void SetYpixels(int value);
	int  GetXpixels() const;  void SetXpixels(int value);
	int  GetPixels() const;   void SetPixels(int value);
	double GetYpitch() const; void SetYpitch(double value);
	double GetXpitch() const; void SetXpitch(double value);
	double GetPitch() const;  void SetPitch(double value);
	double GetHeight() const; void SetHeight(double value);
	double GetWidth() const;  void SetWidth(double value);
	void Clear();
	void Resize(int ypixels,int xpixels);
	void ZeroFill(int ypixels,int xpixels);
	void Trim(int ypixels,int xpixels);
	void TrimCircle(int ypixels,int xpixels);

	complex GetComplexAmplitude(int i,int j) const;
	void    SetComplexAmplitude(int i,int j,complex value);
	double GetAmplitude(int i,int j) const;
	void   SetAmplitude(int i,int j,double value);
	double GetIntensity(int i,int j) const;
	void   SetIntensity(int i,int j,double value);
	void   AddIntensity(int i,int j,double value);
	double GetPhase(int i,int j) const;
	void   SetPhase(int i,int j,double value);
	int PixelLocation(int &i,int &j,double x,double y);
	int Add(double x,double y,double weight);
	int Add(double x,double y);
	int AddComplexAmplitude(double x,double y,complex value);
	void SetLine(double x1,double y1,double x2,double y2,double weight);
	void SetCircle(double xo,double yo,double Dia,double Intensity,
		           double GaussDia=0,double Tilt_deg=0,double TiltAzm_deg=0,double Wl_nm=632.8);

	void Mask(double InnerR,double OuterR);
	void MaskXCut();
	void MaskXPass();

	void NormalizeAmplitude(double new_max);
	void Normalize(double new_max);
	void NormalizeTotalIntensity(double total=1);

	void GCenterToCenter();
	void ToNeg(double new_max);
	void ReverseXY();
	void Zoom(double my,double mx,double io,double jo);
	void Zoom(double m,double io,double jo);
	void Zoom(double m);

	double Transform(int OriginIsCenter,int inv=0);
	double GetCx(int i,int j); void SetCx(int i,int j,double val);
	double GetCy(int i,int j); void SetCy(int i,int j,double val);
	int SetCxCy(std::string filename,std::string tagname);

	void DistCorrect(double dist);

	cImage DFTAmplitude(int optical,double zoom);
	cImage InvDFTAmplitude(int optical,double zoom);
	cImage DFTIntensity(int optical,double zoom);
	cImage InvDFTIntensity(int optical,double zoom);
	
	double Diffraction(double width,double width1,double wl,double z0,double z,double f,int IsFraunhofer);
	double DiffractionCircleAperture
	       (double diameter,int IsGauss,double width1,double wl,double z0,double z,double f,int IsFraunhofer);
	static double FWDiffractionCircleAperture
	       (int size,double diameter,int IsGauss,double wl,double z0,double z,double f,int IsFraunhofer,double threshold);
	double FWxIntensity(double threshold,double width);
	double FWyIntensity(double threshold,double width);
	complex Total();
	double TotalAmplitude();
	double TotalIntensity();
	double MaxIntensity();
	double GetBrightness();
	void   SetBrightness(double value);
	double GetContrast();
	void   SetContrast(double value);
	void   SetBrightnessContrast(double brightness,double contrast);
	void   Intensify(double multiplier);
	void   BiasIntensity(double offset);
	void   BiasAmplitude(double amp,double phase_deg);
	void   Merge(const cImage& X);
	int    MergeBmp(std::string filename);
	cImage MedianFilter();
	void   LowPassFilter(double RelativeR=1);
	void   SetGauss(double dia);
	cImage Edge();
	cImage Binarize(double threshold);

	int SaveIntensity(std::string filename) const;
	friend std::ostream& operator<<(std::ostream& to,const cImage& x);
	friend std::istream& operator>>(std::istream& from,cImage& x);
	int Open(std::string filename);
	int Save(std::string filename) const;
	
	cImage& FromBitmap(const cBitmap& x,int ByAmplitude=0);
	int OpenFromBmp(std::string filename,int ByAmplitude);
	int OpenFromBmp(std::string filename);
	int OpenFromBmpByAmplitude(std::string filename);
	
	cBitmap& ToBitmap(cBitmap& bitmap,int ByAmplitude=0) const;
	cBitmap ToBitmap(int ByAmplitude=0) const;
	static cBitmap ToBitmap(const cImage& r,const cImage& g,const cImage& b);
	int SaveAsBmp(std::string filename,double gamma,int ByAmplitude) const;
	int SaveAsBmp(std::string filename,double gamma=1) const;
	int SaveAsBmpByAmplitude(std::string filename,double gamma=1) const;
	
	friend cImage operator+(const cImage& A,const cImage& B);
	friend cImage SumOnIntensity(const cImage& A,const cImage& B);
	friend cImage operator*(const cImage& A,const cImage& B);
	friend cImage operator*(double x,const cImage& A);
	friend cImage ProductOnIntensity(const cImage& A,const cImage& B);
	friend cImage Convolution(const cImage& A,const cImage& B);
	friend cImage ConvolutionOnIntensity(const cImage& A,const cImage& B);
	static double NCC(const cImage& A,const cImage& B,int i1=0,int j1=0,int i2=0,int j2=0);
	static cImage POC(const cImage &A,const cImage &B);
	static int POC(std::string bmp_filename1,std::string bmp_filename2);
	
	// 以下は例えばVisualBasic上でbufを介して加算，乗算を行うための関数
	void ToBuf();
	void FromBuf();
	void ToBufPlus();
	void ToBufPlusIntensity();
	void ToBufMultiply();
	void ToBufConvolute();
	void ToBufConvoluteIntensity();
};

#endif // !defined(AFX_CIMAGE_H__16186CDD_A7A6_4A2E_ADB2_70C8B7610251__INCLUDED_)
