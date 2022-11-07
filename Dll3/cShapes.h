#if !defined(AFX_CSHAPES_H__A560E400_18B4_11D8_BE64_A01008E6B850__INCLUDED_)
#define AFX_CSHAPES_H__A560E400_18B4_11D8_BE64_A01008E6B850__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "list.h"
#include "geometry.h"
#include "general_func.h"
#include "cBitmap.h"
#include <string>

enum { LINE=1,TEXT=2, SOLID=1,DASH=2,DASHDOT=3 };

struct shape {
	int kind;
	double x1,y1, x2,y2;
	int linestyle;
	std::string text;
	double fontsize;;
	long color;

	shape() {};
	
	shape(double x1,double y1, double x2,double y2, long color, long linestyle) {
		kind=LINE;
		this->x1=x1; this->y1=y1; this->x2=x2; this->y2=y2;
		this->color=color;
		this->linestyle=linestyle;
	};
	
	shape(std::string text, double x1,double y1, double fontsize, long color) {
		kind=TEXT;
		this->text=text;
		this->x1=x1; this->y1=y1;
		this->fontsize=fontsize;
		this->color=color;
	};
	
	void offset(double dx,double dy){
		x1+=dx; y1+=dy;
		x2+=dx; y2+=dy;
	};
	
	int trim(const rect& r){
		line l; point p;
		switch(kind){
		case(LINE):
			l=line(x1,y1,x2,y2);
			if( l.trim(r) ){
				x1=l.p1.x; y1=l.p1.y; x2=l.p2.x; y2=l.p2.y;
				return 1;
			}
			else return 0;
		case(TEXT):
			p=point(x1,y1);
			if(p.IsIncluded(r,0)) return 1; else return 0;
		default:
			return 0;
		}
	};
};


class cShapes
{
	list<shape> shapes;	list<shape> shapes_in_window;
	rect window,view;
	matrix<double> T;
	void initialize();
	void Add(shape newdata);
	void ReSize(double m);

public:
	cShapes();
	void SetWindow(double x1,double y1,double x2,double y2);
	void SetWindow(double aspectratio);
	double WindowAspectRatio();
	void SetViewTopLeft(double x,double y);
	void SetViewWidth(double xw,double yw);
	
	void RemoveTail();
	void AddLine(double x1,double y1,double x2,double y2,long color=0,long linestyle=SOLID);
	void AddBox(double x1,double y1,double x2,double y2,long color=0,long linestyle=SOLID,
		        double ChamferX=0,double ChamferY=0);
	void AddEllipse(double x,double y,double rx,double ry,int div,long color=0,long linestyle=SOLID);
	void AddCircle(double x,double y,double r,int div,long color=0,long linestyle=SOLID);
	void AddText(std::string text,double x,double y,double fontsize,long color=0);
	void Add(cShapes newdata);
	cShapes AddInLine(cShapes newdata) const;

	shape Add3DLine(vector<double> v1,vector<double> v2,long color=0,long linestyle=SOLID);
	void ViewAngle(double rox,double roy,double roz);
	void TopView();
	void BottomView();
	void RightView();
	void LeftView();
	void BackView();
	void xView();
	void yView();
	void zView();
	matrix<double> Tmatrix();

	double OffsetX,OffsetY;
	void Reset();
	void Clear();
	int Size() const;
	int Kind(int i);
	double X1(int i);
	double Y1(int i);
	double X2(int i);
	double Y2(int i);
	int LineStyle(int i);
	std::string Text(int i);
	double FontSize(int i);
	long Color(int i);
	void Zoom(double m);
	void Translate(double dx,double dy);
	void AllView();
	int SaveAsBmp(std::string filename,int yPixels);
};

#endif // !defined(AFX_CSHAPES_H__A560E400_18B4_11D8_BE64_A01008E6B850__INCLUDED_)
