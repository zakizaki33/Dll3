#if !defined(AFX_MATRIX_H__441FB230_4BB5_11D6_A068_0001033BC35C__INCLUDED_)
#define AFX_MATRIX_H__441FB230_4BB5_11D6_A068_0001033BC35C__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if _MSC_VER > 1200   // コンパイラのバージョン VisualC++6.0 = 1200
#define TEMPLATECLASST template<class T>
#else
#define TEMPLATECLASST
#endif

#include <stdio.h>
#include <fstream>
#include <string>
#include "complex.h"  // abs()

template<class T> class matrix {
	int m,n;
public:
	T **a;
public:
	// 【注意】
	//  Visual Studio 2013 でコンパイルするときは，
	//  friend 関数の先頭に "template<class T> " をつけないとエラーとなる．
	//  (Visual C++ 6.0 でコンパイルするときはつけるとエラーになる．)
	matrix();
	matrix(int,int);
	matrix(const matrix<T>&);
	virtual ~matrix();
	matrix<T>& operator=(const matrix<T>&);
	matrix<T>& datacopy(const matrix<T>&);
	TEMPLATECLASST friend bool operator ==(const matrix<T>&, const matrix<T>&);
	int rows() const;
	int columns() const;
	void redim(int m,int n);
	matrix<T>& delrow(int);
	matrix<T>& delcol(int);
	TEMPLATECLASST friend matrix<T> inv(const matrix<T>&);             // inverse matrix
	TEMPLATECLASST friend int rank(const matrix<T>&);                  // rank of matrix
	TEMPLATECLASST friend int is_regular(const matrix<T>&);
	TEMPLATECLASST friend matrix<T> t(const matrix<T>&);               // transposed matrix
	TEMPLATECLASST friend matrix<T> operator*(const matrix<T>&,const matrix<T>&);
	TEMPLATECLASST friend matrix<T> operator*(T,const matrix<T>&);
	TEMPLATECLASST friend matrix<T> operator*(const matrix<T>&,T);
	TEMPLATECLASST friend matrix<T> operator/(const matrix<T>&,T);
	TEMPLATECLASST friend matrix<T> operator+(const matrix<T>&,const matrix<T>&);
	TEMPLATECLASST friend matrix<T> operator-(const matrix<T>&,const matrix<T>&);
	TEMPLATECLASST friend matrix<T> operator-(const matrix<T>&);
	TEMPLATECLASST friend matrix<T> zero(const matrix<T>&);            // zero matrix of the same size
	TEMPLATECLASST friend matrix<T> unit(const matrix<T>&);            // unit matrix of the same size
	T* operator[](int) const;
	TEMPLATECLASST friend void print(const matrix<T>&);
	TEMPLATECLASST friend std::string str(const matrix<T>&);
	TEMPLATECLASST friend std::ostream& operator<<(std::ostream& to,const matrix<T>&);
	TEMPLATECLASST friend std::istream& operator>>(std::istream& from,matrix<T>&);
	int save(std::string filename) const;
	int open(std::string filename);
};


template<class T> matrix<T>::matrix(){
	m=0; n=0;
	a=new T* [m+1];   // これをしないと，デストラクタでnewしてないのにdeleteすることになり，
	                  // エラーが発生する．
}

template<class T> matrix<T>::matrix(int m, int n){
	int i,j;
	this->m=m; this->n=n;
	a=new T* [m+1];
	for(i=1; i<=m; ++i) a[i]=new T [n+1];
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j) a[i][j]=0;
}

template<class T> matrix<T>::matrix(const matrix<T>& x) {
	int i,j;
	m=x.m; n=x.n;
	a=new T* [x.m+1];
	for(i=1; i<=x.m; ++i) a[i]=new T [x.n+1];
    for(i=1; i<=x.m; ++i) for(j=1; j<=x.n; ++j) a[i][j]=x.a[i][j];
}

template<class T> matrix<T>::~matrix() {
	int i;
	for(i=1; i<=m; ++i) delete [] a[i];
	delete [] a;
}

template<class T> matrix<T>& matrix<T>::operator=(const matrix<T>& x) { 
	int i,j;

	if(&x!=this){  // これがないと，自己代入でデータがなくなってしまう．
		redim(x.m,x.n);
		for(i=1; i<=x.m; ++i) for(j=1; j<=x.n; ++j) a[i][j]=x.a[i][j];
	}
	return *this;
}

template<class T> matrix<T>& matrix<T>::datacopy(const matrix<T>& from) {
	// サイズは変えずにデータのみコピーする
	int i,j;
    for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		if(i<=from.m && j<=from.n) a[i][j]=from.a[i][j];
	}
	return *this;
}

template<class T> bool operator ==(const matrix<T>& x,const matrix<T>& y) {
	int m,n, i,j;

	m=x.rows();    if(m!=y.rows()   ) return false;
	n=x.columns(); if(n!=y.columns()) return false;	
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		if(x[i][j]!=y[i][j]) return false;
	}
	return true;
}

template<class T> int matrix<T>::rows() const    { return m; }
template<class T> int matrix<T>::columns() const { return n; }

template<class T> void matrix<T>::redim(int m,int n) { 
	int i,j;

	if(m!=this->m || n!=this->n){
		for(i=1; i<=this->m; ++i) delete [] a[i];
		delete [] a;
		this->m=m; this->n=n;
		a=new T* [m+1];
		for(i=1; i<=m; ++i) a[i]=new T [n+1];
	}
    for(i=1; i<=m; ++i) for(j=1; j<=n; ++j) a[i][j]=0;
}

template<class T> matrix<T>& matrix<T>::delrow(int i){
//	int ii,jj;
//	matrix<T> x=*this;
//	
//	redim(m-1,n);
//	for(ii=1; ii<=m; ++ii) for(jj=1; jj<=n; ++jj){
//		if(ii>=i){
//			a[ii][jj]=x.a[ii+1][jj];
//		}
//		else{
//			a[ii][jj]=x.a[ii][jj];
//		}
//	}
//	return *this;

	// 以下の方が速い  2017.05.11
	int ii,jj;
	T** buf;
	
	for(ii=i; ii<=m-1; ++ii) for(jj=1; jj<=n; ++jj){
		a[ii][jj]=a[ii+1][jj];    // 削除する行より下の行を上へひとつずらして代入する
	}	
	delete [] a[m];               // 最下行を解放

	// ---- 各行先頭へのポインタ配列の最終要素を削除する -----------
	buf=new T* [m];
	for(ii=1; ii<=m-1; ++ii) buf[ii]=a[ii];
	delete [] a;
	a=new T* [m];
	for(ii=1; ii<=m-1; ++ii) a[ii]=buf[ii];
	m=m-1;
	delete [] buf;
	// -------------------------------------------------------------

	return *this;
}

template<class T> matrix<T>& matrix<T>::delcol(int j){
	int ii,jj;
	matrix<T> x=*this;
	
	redim(m,n-1);
	for(ii=1; ii<=m; ++ii) for(jj=1; jj<=n; ++jj){
		if(jj>=j){
			a[ii][jj]=x.a[ii][jj+1];
		}
		else{
			a[ii][jj]=x.a[ii][jj];
		}
	}
	return *this;
}

template<class T> matrix<T> inv(const matrix<T>& x) {
	if( x.m==x.n ) {
		T max,tmp; int ij,i,j, pivot;

		matrix<T> a=x, b(x.m,x.m);

		for(i=1; i<=x.m; ++i) b.a[i][i]=1; // 対角成分以外はコンストラクタにより既に0

		for(ij=1; ij<=x.m; ++ij){
			max=0;
			for(i=ij; i<=x.m; ++i){  // aのij列の要素が最大となる行をpivotとする
				if( abs(a.a[i][ij])>max ){ max=abs(a.a[i][ij]); pivot=i; }
			}
			
			if( max <= 1e-20 ) return zero(x);    // 非正則と判断
			
			for(j=1; j<=x.m; ++j){ // a,bのij行とpivot行を入れ替える (pivot行を行範囲ij～nの一番上に配置する)
				tmp=a.a[ij][j]; a.a[ij][j]=a.a[pivot][j]; a.a[pivot][j]=tmp; 
				tmp=b.a[ij][j]; b.a[ij][j]=b.a[pivot][j]; b.a[pivot][j]=tmp; 
			}
			
			for(i=ij+1; i<=x.m; ++i){ // ij行より下の行に行列基本操作を行い三角行列にしていく
				tmp=a.a[i][ij]/a.a[ij][ij];
				for(j=ij; j<=x.m; ++j) a.a[i][j]-=a.a[ij][j]*tmp;
				for(j=1;  j<=x.m; ++j) b.a[i][j]-=b.a[ij][j]*tmp;
			}
		}    
		
		for(ij=x.m; ij>=1; --ij){
			for(i=ij-1; i>=1; --i){
				tmp=a.a[i][ij]/a.a[ij][ij];
				for(j=1; j<=x.m; j++) b.a[i][j]-=b.a[ij][j]*tmp;
				// aに対する操作は逆行列を求める目的には不要なのでしない
			}
		}  
		
		for(i=1; i<=x.m; ++i) for(j=1; j<=x.m; ++j) b.a[i][j]=b.a[i][j]/a.a[i][i];  // aが単位行列になるときのbが解である
		
		return b;
	}
	else{
		return zero(x);
	}
}

/*
template<class T> matrix<T> inv(const matrix<T>& x) {
	if( x.m==x.n ) {
		T max,tmp; int ij,i,j, pivot;

		matrix<T> a=x, b(x.m,x.m);

		for(i=1; i<=x.m; ++i) b.a[i][i]=1; // 対角成分以外はコンストラクタにより既に0

		for(ij=1; ij<=x.m; ++ij){
			max=0; 
			for(i=ij; i<=x.m; ++i){
				if( abs(a.a[i][ij])>max ){ max=abs(a.a[i][ij]); pivot=i; }
			}
			if( max <= 1e-20 ) return zero(x);    // not regular
			for(j=1; j<=x.m; ++j){
				tmp=a.a[ij][j]; a.a[ij][j]=a.a[pivot][j]; a.a[pivot][j]=tmp; 
				tmp=b.a[ij][j]; b.a[ij][j]=b.a[pivot][j]; b.a[pivot][j]=tmp; 
			}
			for(i=ij+1; i<=x.m; ++i){
				for(j=1; j<=x.m; j++) b.a[i][j]-=b.a[ij][j]*a.a[i][ij]/a.a[ij][ij];
				for(j=x.m; j>=1; j--) a.a[i][j]-=a.a[ij][j]*a.a[i][ij]/a.a[ij][ij];
			}
		}    
		for(ij=x.m; ij>=1; --ij) for(i=ij-1; i>=1; --i){
			for(j=1; j<=x.m; j++) b.a[i][j]-=b.a[ij][j]*a.a[i][ij]/a.a[ij][ij];
			for(j=1; j<=x.m; j++) a.a[i][j]-=a.a[ij][j]*a.a[i][ij]/a.a[ij][ij];
		}  
		for(i=1; i<=x.m; ++i) for(j=1; j<=x.m; ++j) b.a[i][j]=b.a[i][j]/a.a[i][i];
		return b;
	}
	else{
		return zero(x);
	}
}
*/

template<class T> int rank(const matrix<T>& x) {
	if( x.m==x.n ) {
		T max,tmp; int ij,i,j, pivot;

		matrix<T> a=x; 

		for(ij=1; ij<=x.m; ++ij){
			
			max=0; 
			for(i=ij; i<=x.m; ++i){
				if( fabs(a.a[i][ij])>max ){ max=fabs(a.a[i][ij]); pivot=i; }
			}

			if( max <= 1e-20 ) return ij-1;
			
			for(j=1; j<=x.m; ++j){
				tmp=a.a[ij][j]; a.a[ij][j]=a.a[pivot][j]; a.a[pivot][j]=tmp; 
			}

			for(i=ij+1; i<=x.m; ++i){
				for(j=x.m; j>=1; j--) a.a[i][j]-=a.a[ij][j]*a.a[i][ij]/a.a[ij][ij];
			}
		}

		return x.m;
	}
	else{
		return -1;
	}
}

/*
template<class T> int rank(const matrix<T>& x) {
	if( x.m==x.n ) {
		T max,tmp; int ij,i,j, pivot;

		matrix<T> a=x; 

		for(ij=1; ij<=x.m; ++ij){
			max=0; 
			for(i=ij; i<=x.m; ++i){
				if( fabs(a.a[i][ij])>max ){ max=fabs(a.a[i][ij]); pivot=i; }
			}
			if( max <= 1e-20 ) return ij-1;    // return rank(x)
			for(j=1; j<=x.m; ++j){
				tmp=a.a[ij][j]; a.a[ij][j]=a.a[pivot][j]; a.a[pivot][j]=tmp; 
			}
			for(i=ij+1; i<=x.m; ++i){
				for(j=x.m; j>=1; j--) a.a[i][j]-=a.a[ij][j]*a.a[i][ij]/a.a[ij][ij];
			}
		}    
		for(ij=x.m; ij>=1; --ij) for(i=ij-1; i>=1; --i){
			for(j=1; j<=x.m; j++) a.a[i][j]-=a.a[ij][j]*a.a[i][ij]/a.a[ij][ij];
		}  
		return x.m;
	}
	else{
		return -1;
	}
}
*/

template<class T> int is_regular(const matrix<T>& x){
	// xが正則かどうか．
	if(rank(x)==x.m) return 1;
	else             return 0;
}

template<class T> matrix<T> t(const matrix<T>& x) {
	int i,j;
	matrix<T> p(x.n,x.m);
	for(i=1; i<=x.n; ++i) for(j=1; j<=x.m; ++j) p.a[i][j]=x.a[j][i];
	return p;
}

template<class T> matrix<T> operator*(const matrix<T>& x,const matrix<T>& y) {
	int i,j,k;
	matrix<T> c(x.m,y.n);  // コンストラクタにより零行列
	if(x.n==y.m){
		for(i=1; i<=x.m; ++i) for(j=1; j<=y.n; ++j){
			for(k=1; k<=x.n; ++k) c.a[i][j]+=x.a[i][k]*y.a[k][j];
		}
		return c;
	}
    else{
		return zero(c);
	}
}

template<class T> matrix<T> operator*(T x,const matrix<T>& y) {
	int i,j;
	matrix<T> c(y.m,y.n);
    for(i=1; i<=y.m; ++i) for(j=1; j<=y.n; ++j) c.a[i][j]=x*y.a[i][j];
	return c;
}

template<class T> matrix<T> operator*(const matrix<T>& x,T y) {
	int i,j;
	matrix<T> c(x.m,x.n);
    for(i=1; i<=x.m; ++i) for(j=1; j<=x.n; ++j) c.a[i][j]=x.a[i][j]*y;
	return c;
}

template<class T> matrix<T> operator/(const matrix<T>& x,T y) {
	int i,j;
	matrix<T> c(x.m,x.n);
	for(i=1; i<=x.m; ++i) for(j=1; j<=x.n; ++j) c.a[i][j]=x.a[i][j]/y;
	return c;
}

template<class T> matrix<T> operator+(const matrix<T>& x,const matrix<T>& y) {
	int i,j;
	matrix<T> c(x.m,x.n);
	if(x.m==y.m && x.n==y.n){
		for(i=1; i<=x.m; ++i) for(j=1; j<=x.n; ++j) {
			c.a[i][j]=x.a[i][j]+y.a[i][j];
		}
		return c;
	}
    else{
		return zero(c);
	}
}

template<class T> matrix<T> operator-(const matrix<T>& x,const matrix<T>& y) {
	int i,j;
	matrix<T> c(x.m,x.n);
	if(x.m==y.m && x.n==y.n){
		for(i=1; i<=x.m; ++i) for(j=1; j<=x.n; ++j) {
			c.a[i][j]=x.a[i][j]-y.a[i][j];
		}
		return c;
	}
    else{
		return zero(c);
	}
}

template<class T> matrix<T> operator-(const matrix<T>& x) {
	int i,j;
	matrix<T> c(x.m,x.n);
	for(i=1; i<=x.m; ++i) for(j=1; j<=x.n; ++j) c.a[i][j]=-x.a[i][j];
    return c;
}

template<class T> matrix<T> zero(const matrix<T>& x) {
	matrix<T> p(x.m,x.n);  // コンストラクタにより零行列になる．
	return p;
}

template<class T> matrix<T> unit(const matrix<T>& x) {
	int i,j;
	matrix<T> p(x.m,x.n);  // コンストラクタにより零行列

	if (x.m==x.n) {
		for(i=1; i<=x.m; ++i) for(j=1; j<=x.n; ++j)	{
			if(i==j) p.a[i][j]=1;
		}
		return p;
	}
	else {
		return zero(p);
	}
}

template<class T> matrix<T> unit(int m) {
	matrix<T> p(m,m);
	return unit(p);
}

template<class T> T* matrix<T>::operator[](int i) const{
	return a[i];

/*****************************************
	class X
	{
		double* a;

	public:
		X()
		{
			a=new double [10];	
		};

		~X()
		{
			delete [] a;
		};

		void f() const
		{
			a=0;      // error : lvalue is const
			a[1]=0;   // ok    : ( can't protect )
		};
	};
******************************************/
}

template<class T> void print(const matrix<T>& x) {
	int i,j;
	for(i=1; i<=x.m; ++i) {
		for(j=1; j<=x.n; ++j) printf("%7g  ",x.a[i][j]);
		printf("\n");
	}
}

template<class T> std::string str(const matrix<T>& x){
	int i,j;
	std::string s;
	for(i=1; i<=x.m; ++i) {
		for(j=1; j<=x.n; ++j){
			s+=str(x.a[i][j]);
			s+='\t';
		}
		s+='\n';
	}
	return s;
}

template<class T> std::ostream& operator<<(std::ostream& to,const matrix<T>& x){
	int i,j;
	to << x.m << ' ' << x.n << '\n';
	for(i=1; i<=x.m; ++i) {
		for(j=1; j<=x.n; ++j) to << x.a[i][j] << '\t';
		to << '\n';
	}
	return to;
}

template<class T> std::istream& operator>>(std::istream& from,matrix<T>& x){
	int i,j, m,n;
	from >> m >> n;
	x.redim(m,n);
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		from >> x.a[i][j];
	}
	return from;
}

template<class T> int matrix<T>::save(std::string filename) const{
	std::ofstream to(filename.c_str());
	if(to){
		to << *this;
		return 1;
	}
	else{
		return 0;
	}
}

template<class T> int matrix<T>::open(std::string filename){
	std::ifstream from(filename.c_str());
	if(from){
		from >> *this;
		return 1;
	}
	else{
		return 0;
	}
}

#endif // !defined(AFX_MATRIX_H__441FB230_4BB5_11D6_A068_0001033BC35C__INCLUDED_)
