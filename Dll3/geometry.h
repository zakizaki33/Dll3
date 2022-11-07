#ifndef GEOMETRY_H_INCLUDED
#define GEOMETRY_H_INCLUDED
#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <math.h>
#include "general_func.h"   //S04フォルダーにある

struct point;
struct line;
struct rect;

struct point {
	double x,y;
	point();
	point(double x,double y);
	bool operator==(const point& p) const;
	bool operator!=(const point& p) const;
	friend double abs(const point& p);
	friend point operator+(const point& p1,const point& p2);
	friend point operator-(const point& p1,const point& p2);
	friend point operator-(const point& p);
	friend point operator*(const point& p,double a);
	friend point operator*(double a,const point& p);
	friend point operator/(const point& p,double a);
	point& operator+=(const point& p);
	point& operator-=(const point& p);
	point& operator*=(double a);
	point& operator/=(double a);
	point& rotate(double th_deg);
	friend double distance(const point& p1,const point& p2);
	friend double distance(const point& p,const line& l);
	int IsIncluded(const rect& r,int boundary_is_open) const;
	int IntersectPoint(const line& l1,const line& l2);
	std::string str() const;
	friend std::ostream& operator<<(std::ostream& to,const point& p);
};

struct line {
	point p1,p2;
	line();
	line(const point& p1,const point& p2);
	line(double x1,double y1,double x2,double y2);
	bool operator==(const line& l) const;
	int trim(const rect& r);
	std::string str() const;
};

struct rect {
	point p1,p2;
	rect();
	rect(const point& p1,const point& p2);
	rect(const line& l);
	rect(double x1,double y1,double x2,double y2);
	bool operator==(const rect& r) const;
	point TopLeft() const;
	point TopRight() const;
	point BottomLeft() const;
	point BottomRight() const;
	std::string str() const;
};

#endif // #ifndef GEOMETRY_H_INCLUDED