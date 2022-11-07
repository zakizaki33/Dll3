#include "stdafx.h"
#include "geometry.h"

////////  point  /////////////////////////////////////////////////
point::point() { x=y=0; }

point::point(double x,double y) { this->x=x; this->y=y; }

bool point::operator==(const point& p) const {
	return x==p.x && y==p.y;
}

bool point::operator!=(const point& p) const {
	return !(*this==p);
}

double abs(const point& p){
	return sqrt(p.x*p.x+p.y*p.y);
}

point operator+(const point& p1,const point& p2){
	return point( p1.x+p2.x, p1.y+p2.y );
}

point operator-(const point& p1,const point& p2){
	return point( p1.x-p2.x, p1.y-p2.y );
}

point operator-(const point& p){
	return point( -p.x, -p.y );
}

point operator*(const point& p,double a){
	return point( p.x*a, p.y*a );
}

point operator*(double a,const point& p){
	return point( p.x*a, p.y*a );
}

point operator/(const point& p,double a){
	return point( p.x/a, p.y/a );
}

point& point::operator+=(const point& p){
	x+=p.x; y+=p.y;
	return *this;
}

point& point::operator-=(const point& p){
	x-=p.x; y-=p.y;
	return *this;
}

point& point::operator*=(double a){
	x*=a; y*=a;
	return *this;
}

point& point::operator /=(double a){
	x/=a; y/=a;
	return *this;
}

point& point::rotate(double th_deg){
	double th, x0,y0;
	if(th_deg!=0){
		th=th_deg*PI/180;
		x0=x;
		y0=y;
		x=cos(th)*x0-sin(th)*y0;
		y=sin(th)*x0+cos(th)*y0;
	}
	return *this;
}

double distance(const point& p1,const point& p2){
	return sqrt( (p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y) );
}

double distance(const point& p,const line& l){
	// Ågê¸å`ë„êîÅhp.96
	double costh; 
	costh=((l.p2.x-l.p1.x)*(p.x-l.p1.x)+(l.p2.y-l.p1.y)*(p.y-l.p1.y))
	      /distance(l.p1,l.p2)/distance(p,l.p1);
	return distance(p,l.p1)*sqrt(1-costh*costh);
}

int point::IsIncluded(const rect& r,int boundary_is_open) const{
	if(boundary_is_open){
		return ( r.p1.x<r.p2.x ? (r.p1.x<x && x<r.p2.x):(r.p2.x<x && x<r.p1.x) ) &&
			   ( r.p1.y<r.p2.y ? (r.p1.y<y && y<r.p2.y):(r.p2.y<y && y<r.p1.y) );
	}
	else{
		return ( r.p1.x<r.p2.x ? (r.p1.x<=x && x<=r.p2.x):(r.p2.x<=x && x<=r.p1.x) ) &&
			   ( r.p1.y<r.p2.y ? (r.p1.y<=y && y<=r.p2.y):(r.p2.y<=y && y<=r.p1.y) );
	}
}

int point::IntersectPoint(const line& l1,const line& l2){
	double det;
	point p;
	det=(l1.p1.y-l1.p2.y)*(l2.p2.x-l2.p1.x)-(l2.p1.y-l2.p2.y)*(l1.p2.x-l1.p1.x);
	if(det==0) return 0;
	p.x=(1/det)*( (l2.p2.x-l2.p1.x)*(l1.p2.x*l1.p1.y-l1.p1.x*l1.p2.y)
	             +(l1.p1.x-l1.p2.x)*(l2.p2.x*l2.p1.y-l2.p1.x*l2.p2.y) );
	p.y=(1/det)*( (l2.p2.y-l2.p1.y)*(l1.p2.x*l1.p1.y-l1.p1.x*l1.p2.y)
	             +(l1.p1.y-l1.p2.y)*(l2.p2.x*l2.p1.y-l2.p1.x*l2.p2.y) );
	if( l1.p1.x!=l1.p2.x && (l1.p1.x-p.x)*(l1.p2.x-p.x)>0 ) return 0;
	if( l1.p1.y!=l1.p2.y && (l1.p1.y-p.y)*(l1.p2.y-p.y)>0 ) return 0;
	if( l2.p1.x!=l2.p2.x && (l2.p1.x-p.x)*(l2.p2.x-p.x)>0 ) return 0;
	if( l2.p1.y!=l2.p2.y && (l2.p1.y-p.y)*(l2.p2.y-p.y)>0 ) return 0;
	*this=p;
	return 1;
}

std::string point::str() const {
	char buf[100];
	sprintf(buf, "(%g,%g)", x,y);
	return buf;
}

std::ostream& operator<<(std::ostream& to,const point& p){
	to << '(' << p.x << ',' << p.y << ')';
	return to;
}


////////  line  /////////////////////////////////////////////////
line::line(){}

line::line(const point& p1,const point& p2){
	this->p1=p1;
	this->p2=p2;
}

line::line(double x1,double y1,double x2,double y2){
	p1=point(x1,y1);
	p2=point(x2,y2);
}

bool line::operator==(const line& l) const {
	return p1==l.p1 && p2==l.p2;
}

int line::trim(const rect& r){
	point q1,q2;
	bool q1_set=false, q2_set=false;
	line l1,l2,l3,l4;
	l1=line(r.p1,point(r.p1.x,r.p2.y));
	l2=line(point(r.p1.x,r.p2.y),r.p2);
	l3=line(r.p2,point(r.p2.x,r.p1.y));
	l4=line(point(r.p2.x,r.p1.y),r.p1);
	if( p1.IsIncluded(r,1) ){
		q1=p1;
		q1_set=true;
	}
	else{
		for(;;){
			if(q1.IntersectPoint(*this,l1)) { q1_set=true; break; }
			if(q1.IntersectPoint(*this,l2)) { q1_set=true; break; }
			if(q1.IntersectPoint(*this,l4)) { q1_set=true; break; }
			if(q1.IntersectPoint(*this,l3)) { q1_set=true; break; }
			break;
		}
	}
	if( p2.IsIncluded(r,1) ){
		q2=p2;
		q2_set=true;
	}
	else{
		for(;;){
			if(q2.IntersectPoint(*this,l3)) { q2_set=true; break; }
			if(q2.IntersectPoint(*this,l4)) { q2_set=true; break; }
			if(q2.IntersectPoint(*this,l2)) { q2_set=true; break; }
			if(q2.IntersectPoint(*this,l1)) { q2_set=true; break; }
			break;
		}
	}
	p1=q1; p2=q2;
	if( q1_set && q2_set ) return 1; else return 0;
}

std::string line::str() const {
	char buf[100];
	sprintf(buf, "(%g,%g)-(%g,%g)", p1.x,p1.y,p2.x,p2.y);
	return buf;
}

////////  rect  /////////////////////////////////////////////////
rect::rect(){}

rect::rect(const point& p1,const point& p2){
	this->p1=p1;
	this->p2=p2;
}

rect::rect(const line& l){
	this->p1=l.p1;
	this->p2=l.p2;
}

rect::rect(double x1,double y1,double x2,double y2){
	p1=point(x1,y1);
	p2=point(x2,y2);
}

bool rect::operator==(const rect& l) const {
	return p1==l.p1 && p2==l.p2;
}

point rect::TopLeft() const {
	return point( p1.x>p2.x ? p2.x:p1.x, p1.y>p2.y ? p1.y:p2.y );
}
point rect::TopRight() const {
	return point( p1.x>p2.x ? p1.x:p2.x, p1.y>p2.y ? p1.y:p2.y );
}
point rect::BottomLeft() const {
	return point( p1.x>p2.x ? p2.x:p1.x, p1.y>p2.y ? p2.y:p1.y );
}
point rect::BottomRight() const {
	return point( p1.x>p2.x ? p1.x:p2.x, p1.y>p2.y ? p2.y:p1.y );
}

std::string rect::str() const {
	char buf[100];
	sprintf(buf, "(%g,%g)-(%g,%g)", p1.x,p1.y,p2.x,p2.y);
	return buf;
}