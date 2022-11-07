#ifndef VECTOR_H_INCLUDED
#define VECTOR_H_INCLUDED
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1200   // コンパイラのバージョン VisualC++6.0 = 1200
#define TEMPLATECLASST template<class T>
#else
#define TEMPLATECLASST
#endif

#include "complex.h"
#include "matrix.h"
#include <math.h>

template<class T> class vector  
{
public:
	T x,y,z;
public:
	vector();
	vector(const T& x,const T& y,const T& z);

	// 【注意】
	//  Visual Studio 2013 でコンパイルするときは，
	//  friend 関数の先頭に "template<class T> " をつけないとエラーとなる．
	//  (Visual C++ 6.0 でコンパイルするときはつけるとエラーになる．)
	
	TEMPLATECLASST friend bool operator==(const vector<T>& a,const vector<T>& b);
	TEMPLATECLASST friend bool operator!=(const vector<T>& a,const vector<T>& b);
	TEMPLATECLASST friend double abs(const vector<T>& a);
	TEMPLATECLASST friend vector<T> DirectionCosine(const vector<T>& a);
	TEMPLATECLASST friend vector<T> operator+(const vector<T>& a,const vector<T>& b);
	TEMPLATECLASST friend vector<T> operator-(const vector<T>& a,const vector<T>& b);
	TEMPLATECLASST friend vector<T> operator-(const vector<T>& a);
	TEMPLATECLASST friend vector<T> operator*(const T& c,const vector<T>& a);
	TEMPLATECLASST friend vector<T> operator*(const vector<T>& a,const T& c);
	TEMPLATECLASST friend vector<T> operator*(const matrix<double>& A,const vector<T>& a);
	TEMPLATECLASST friend vector<T> operator/(const vector<T>& a,const T& c);
	vector<T>& operator+=(const vector<T>& a);
	vector<T>& operator-=(const vector<T>& a);
	vector<T>& operator*=(const T& c);
	TEMPLATECLASST friend T sProduct(const vector<T>& a,const vector<T>& b);
	TEMPLATECLASST friend T cosine(const vector<T>& a,const vector<T>& b);
	TEMPLATECLASST friend vector<T>  vProduct(const vector<T>& a,const vector<T>& b);
	vector<T> x_rotate(double th) const;
	vector<T> y_rotate(double th) const;
	vector<T> z_rotate(double th) const;
	TEMPLATECLASST friend void x_rotate(vector<T>& v,double th);
	TEMPLATECLASST friend void y_rotate(vector<T>& v,double th);
	TEMPLATECLASST friend void z_rotate(vector<T>& v,double th);
	TEMPLATECLASST friend std::ostream& operator<<(std::ostream& to,const vector<T>& a);
};



template<class T> vector<T>::vector() {
	// 高速化のため0に初期化しない（効果有り）．
	// double, int なども初期化されない．これに合わせる．
}

template<class T> vector<T>::vector(const T& x,const T& y,const T& z) {
	this->x=x; this->y=y; this->z=z;
}

template<class T> bool operator==(const vector<T>& a,const vector<T>& b){
	return ( a.x==b.x && a.y==b.y && a.z==b.z );
}

template<class T> bool operator!=(const vector<T>& a,const vector<T>& b){
	return !(a==b);
}

template<class T> double abs(const vector<T>& a) {
	// 2021-08-23 今回の目的には　conjは使わないのでいったん消す
	//return sqrt( abs(a.x*conj(a.x)+a.y*conj(a.y)+a.z*conj(a.z)) );
	return 0;
}

template<class T> vector<T> DirectionCosine(const vector<T>& a) {
	return a/abs(a);
}

template<class T> vector<T> operator+(const vector<T>& a,const vector<T>& b) {
	return vector<T>(a.x+b.x, a.y+b.y, a.z+b.z);
}

template<class T> vector<T> operator-(const vector<T>& a,const vector<T>& b) {
	return vector<T>(a.x-b.x, a.y-b.y, a.z-b.z);
}

template<class T> vector<T> operator-(const vector<T>& a) {
	return vector<T>(-a.x, -a.y, -a.z);
}

template<class T> vector<T> operator*(const T& c,const vector<T>& a) {
	return vector<T>(c*a.x, c*a.y, c*a.z);
}

template<class T> vector<T> operator*(const vector<T>& a,const T& c) {
	return vector<T>(c*a.x, c*a.y, c*a.z);
}

template<class T> vector<T> operator*(const matrix<double>& A,const vector<T>& a) {
	vector<T> v;
	if( A.rows()==3 && A.columns()==3 ){
		v.x=A[1][1]*a.x +A[1][2]*a.y +A[1][3]*a.z;
		v.y=A[2][1]*a.x +A[2][2]*a.y +A[2][3]*a.z;
		v.z=A[3][1]*a.x +A[3][2]*a.y +A[3][3]*a.z;
	}
	else{
		v=a;
	}
	return v;

//	以下のコードだと，コンパイラ最適化="実行速度"のときに，
//	cLens1 raytrace(...)のアスの追跡ができず，ASの収差図などがおかしくなる．
//	(戻り値がおかしな値になる．) 原因は分からない．        08.07.09
//
//	vector<T> v;
//	if( A.rows()==3 && A.columns()==3 ){
//		v.x=A[1][1]*a.x +A[1][2]*a.y +A[1][3]*a.z;
//		v.y=A[2][1]*a.x +A[2][2]*a.y +A[2][3]*a.z;
//		v.z=A[3][1]*a.x +A[3][2]*a.y +A[3][3]*a.z;
//		return v;
//	}
//	else{
//		v=a;
//		return v;
//	}
}

template<class T> vector<T> operator/(const vector<T>& a,const T& c) {
	return c==0 ? vector<T>(0,0,0) : vector<T>(a.x/c, a.y/c, a.z/c);
}

template<class T> vector<T>& vector<T>::operator+=(const vector<T>& a) {
	return *this=*this+a;
}

template<class T> vector<T>& vector<T>::operator-=(const vector<T>& a) {
	return *this=*this-a;
}

template<class T> vector<T>& vector<T>::operator*=(const T& c) {
	return *this=*this*c;
}

template<class T> T sProduct(const vector<T>& a,const vector<T>& b) {
	// 73行目と同じ理由
	// return a.x*conj(b.x)+a.y*conj(b.y)+a.z*conj(b.z);
	return 0;
}

template<class T> T cosine(const vector<T>& a,const vector<T>& b) {
	return sProduct(a,b)/abs(a)/abs(b);
}

template<class T> vector<T> vProduct(const vector<T>& a,const vector<T>& b) {
	return vector<T>(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
}

template<class T> vector<T> vector<T>::x_rotate(double th) const {
	double costh,sinth;
	if(th==0) return *this;
	costh=cos(th); sinth=sin(th);
	return vector<T>(x, costh*y-sinth*z, sinth*y+costh*z);	
}

template<class T> vector<T> vector<T>::y_rotate(double th) const {
	double costh,sinth;
	if(th==0) return *this;
	costh=cos(th); sinth=sin(th);
	return vector<T>(sinth*z+costh*x, y, costh*z-sinth*x);
}

template<class T> vector<T> vector<T>::z_rotate(double th) const {
	double costh,sinth;
	if(th==0) return *this;
	costh=cos(th); sinth=sin(th);
	return vector<T>(costh*x-sinth*y, sinth*x+costh*y, z);
}

template<class T> void x_rotate(vector<T>& v,double th) {
	static double cache=0;
	static double costh=1,sinth=0; 
	T y0,z0;
	if( th==0 ) return;
	if( th!=cache ) { costh=cos(th); sinth=sin(th); cache=th; }
	y0=v.y; z0=v.z;
	v.y=costh*y0-sinth*z0; v.z=sinth*y0+costh*z0;	
}

template<class T> void y_rotate(vector<T>& v,double th) {
	static double cache=0;
	static double costh=1,sinth=0; 
	T z0,x0;
	if( th==0 ) return;
	if( th!=cache ) { costh=cos(th); sinth=sin(th); cache=th; }
	z0=v.z; x0=v.x;
	v.z=costh*z0-sinth*x0; v.x=sinth*z0+costh*x0; 
}

template<class T> void z_rotate(vector<T>& v,double th) {
	static double cache=0;
	static double costh=1,sinth=0; 
	T x0,y0;
	if( th==0 ) return;
	if( th!=cache ) { costh=cos(th); sinth=sin(th); cache=th; }
	x0=v.x; y0=v.y;
	v.x=costh*x0-sinth*y0; v.y=sinth*x0+costh*y0;
}

template<class T> std::ostream& operator<<(std::ostream& to,const vector<T>& a){
	to << '(' << a.x << ',' << a.y << ',' << a.z << ')';
	return to;
}

#endif // #ifndef VECTOR_H_INCLUDED