#include "stdafx.h"
#include "complex.h"
const double PI=3.14159265358979;
const complex i=complex(0,1);

complex::complex() {
	// 高速化のため0に初期化しない．
	// double, int なども初期化されない．これに合わせる．
	// 追記 2021-05-12 Yamazaki
	this->x = 0; this->y = 0;
}

complex::complex(double x) {
	this->x=x; y=0;
}

complex::complex(double x,double y) {
	this->x=x; this->y=y;
}

complex::~complex() {
}

double& Re(complex& z) { return z.x; }
// double Re(const complex& z) { return z.x; }

double& Im(complex& z) { return z.y; }

double abs(const complex& z){
	return sqrt(z.x*z.x+z.y*z.y);
}

/*
#if _MSC_VER < 1800  // Visual C++ 2013 より古いとき
double abs(double x){
	return fabs(x);
	// テンプレートで T x; x=abs(x);
	// としたとき，本関数がないと，Tがcomplexであればよいが，doubleのとき，
	// int abs(int) が呼ばれてしまい，小数点以下が落ちてしまう．

	// 【重要】
	//  新しいC++(*)では，cmath.h で double abs(double)  が定義されており，
	//  コンパイルで再定義エラー(C2084)となるため，本関数を削除すること．
	//  (*) 例えば，Visual C++ 6.0 は古く，Visual Studio 2013 は新しい．
}
#endif // _MSC_VER < 1800
*/

double sqabs(const complex& z){
	return z.x*z.x+z.y*z.y;
}

double arg(const complex& z) {
	// -PI <= arg(z) < PI
	if(z.y==0){
		if(z.x>=0) return 0;
		else       return -PI;  // z.x<0
	}
	else if(z.y>0){
		if     (z.x==0) return PI/2;
		else if(z.x>0 ) return atan(z.y/z.x);
		else            return atan(z.y/z.x)+PI;  // z.x<0
	}
	else{  // z.y<0
		if     (z.x==0) return -PI/2;
		else if(z.x>0 ) return atan(z.y/z.x);
		else            return atan(z.y/z.x)-PI;  // z.x<0
	}
}

complex conj(const complex& z) { return complex(z.x,-z.y); }
double  conj(const double& x)  { return x; }

complex operator*(const complex& z) {
	return conj(z);
}

bool operator==(const complex& z1,const complex& z2) {
	if( z1.x==z2.x && z1.y==z2.y ) return true;
	else                           return false;
}

bool operator!=(const complex& z1,const complex& z2) {
	return !(z1==z2);
}

complex operator+(const complex& z1,const complex& z2){
	return complex(z1.x+z2.x, z1.y+z2.y);
}

complex operator-(const complex& z1,const complex& z2){
	return complex(z1.x-z2.x, z1.y-z2.y);
}

complex operator-(const complex& z){
	return complex(-z.x, -z.y);
}

complex operator*(const complex& z1,const complex& z2){
	return complex(z1.x*z2.x-z1.y*z2.y, z1.x*z2.y+z1.y*z2.x);
}

complex operator/(const complex& z1,const complex& z2){
	double a=z2.x*z2.x+z2.y*z2.y;
	if(a==0){
		return 0;
	}
	else {
		return complex( (z1.x*z2.x+z1.y*z2.y)/a, (-z1.x*z2.y+z1.y*z2.x)/a );
	}	
}

complex& complex::operator+=(const complex& z) { 
	return *this=*this+z;
}

complex& complex::operator -=(const complex& z) {
	return *this=*this-z;
}

complex& complex::operator *=(const complex& z) {
	return *this=*this*z;
}

complex& complex::operator /=(const complex& z) {
	return *this=*this/z;
}
	
complex exp(const complex& z){
	return exp(z.x)*complex(cos(z.y),sin(z.y));    
}

complex cos(const complex& z){
	return ( exp(i*z)+exp(-i*z) )/2;
}

complex sin(const complex& z){
	return ( exp(i*z)-exp(-i*z) )/2/i;
}

complex sqrt(const complex& z){
	return sqrt( abs(z) )*exp( i*arg(z)/2 );
}

std::ostream& operator<<(std::ostream& to,const complex& z){
	to << '(' << trim(str(z.x),0) << ',' << trim(str(z.y),0) << ')';
	// operator>> の都合により，空白が入り込まないようにする．
	return to;
}

std::istream& operator>>(std::istream& from,complex& z){
	std::string s;

	from>>s;                         // s="(0.12,5.43)"
	s=trim(s,1);                     // s="0.12,5.43"
	s=replace(s,","," ");            // s="0.12 5.43"
	z.x=atof(word(s,1,0).c_str());   // z.x=0.12
	z.y=atof(word(s,2,0).c_str());   // z.y=5.43
	return from;
}


