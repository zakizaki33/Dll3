#if !defined(AFX_CBITMAP_H__BC9C379C_6373_499F_A226_3C6629F24FB4__INCLUDED_)
#define AFX_CBITMAP_H__BC9C379C_6373_499F_A226_3C6629F24FB4__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <math.h>
#include "general_func.h"

typedef unsigned char  BYTE;
typedef unsigned short WORD;
typedef unsigned long  DWORD;
typedef long           LONG;

typedef struct tagBMPFILEHEADER {
        WORD  bfType;
        DWORD bfSize;
        WORD  bfReserved1;
        WORD  bfReserved2;
        DWORD bfOffBits;
} BMPFILEHEADER;

typedef struct tagBMPINFOHEADER{
        DWORD biSize;
        LONG  biWidth;
        LONG  biHeight;
        WORD  biPlanes;
        WORD  biBitCount;
        DWORD biCompression;
        DWORD biSizeImage;
        LONG  biXPelsPerMeter;
        LONG  biYPelsPerMeter;
        DWORD biClrUsed;
        DWORD biClrImportant;
} BMPINFOHEADER;

typedef struct tagPALLET{
	BYTE rgbBlue;
	BYTE rgbGreen;
	BYTE rgbRed;
	BYTE rgbReserved;
} PALLET;

class cBitmap  
{
	// âÊëf(i,j), ècï˚å¸ i=1Å`m
	//            â°ï˚å¸ j=1Å`n
	//            ç∂è„í[Ç™ (i,j)=(1,1)
	int m,n;
	BYTE **R,**G,**B;
	double gamma;
	BYTE a[256];    // gammaï‚ê≥ÉeÅ[ÉuÉã
	void alloc();
	void free();
	void make_gamma_table();
	static cBitmap buf,buf2;
public:
	cBitmap();
	cBitmap(const cBitmap& x);
	virtual ~cBitmap();
	cBitmap& operator=(const cBitmap& x);

	int  GetM() const; void SetM(int m);  // âÊëúçÇÇ≥
	int  GetN() const; void SetN(int n);  // âÊëúïù
	double GetGamma();
	void   SetGamma(double gamma);
	long GetRGB(int i,int j);
	void SetRGB(int i,int j,int R,int G,int B);
	void SetRGB(int i,int j,long col);
	void SetGrayScale(int i,int j,int val);
	
	void SetLine(double x1,double y1,double x2,double y2,int R,int G,int B);
	void SetCircle(double x0,double y0,double radius,int R,int G,int B,int fill=0);

	int  GetR(int i,int j) const;
	int  GetG(int i,int j) const;
	int  GetB(int i,int j) const;
	cBitmap& ToR();
	cBitmap& ToG();
	cBitmap& ToB();
	cBitmap& ToWhite();
	cBitmap& ToNeg();
	cBitmap& Expand();
	int Height();
	int Width();
	void Resize(int m,int n);
	void Resize(double ratio);
	void Trim(int i0,int j0,int height,int width);
	void TrimCircle();
	void GammaCompensation(double gamma);
	void Zoom(int i1,int j1,int i2,int j2);
	void SeeColorShift(int i,int j);
	void ColorShiftCorrect(std::string filename,int inv=0);
	void DistCorrect(double dist);
	void Edge();
	void Median();
	void Binarize(double threshold);
	void MarkSaturation();
	
	int Save(std::string filename,double gamma=1);
	int Open(std::string filename);
	int Merge(std::string filename_R,std::string filename_G,std::string filename_B);

	void ToBuf();
	void FromBuf();
	void ToBuf2();
	void FromBuf2();
};

#endif // !defined(AFX_CBITMAP_H__BC9C379C_6373_499F_A226_3C6629F24FB4__INCLUDED_)
