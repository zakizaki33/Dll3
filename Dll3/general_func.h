#ifndef GENERAL_FUNC_H_INCLUDED
#define GENERAL_FUNC_H_INCLUDED
#define _CRT_SECURE_NO_WARNINGS
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "complex.h"
#include "matrix.h"
#include "vector.h"
#include "string.h"
#include "cBitmap.h"
#include "cSpline.h"

const double PI=3.14159265358979;

template <class T> void Dim(T** &p,int m,int n);
template <class T> void Erase(T** &p,int m);
template <class T> void Dim3(T*** &p,int n1,int n2,int n3);
template <class T> void Erase3(T*** &p,int n1,int n2);

int CreateConsole();
std::string RemoveExtension(std::string filename);

std::string Date();
std::string Time();
std::string LapTime(clock_t start);

void dwrite(std::ostream to,double x);
double dread(std::istream from);
void iwrite(std::ostream to,int x);
int iread(std::istream from);
void swrite(std::ostream to,std::string s);
std::string sread(std::istream from);

template <class T> void Swap(T &x,T &y){
	T buf;
	buf=y;
	y=x;
	x=buf;
}

bool is_even(const int& i);
bool is_odd(const int& i);
double factorial(const int& i);
double pw(const double& x,const int& n);
double fceil(const double& x);
double ffloor(const double& x);
double fceil_mantissa(const double& x);
double sgn(const double& x);
double Round(const double& x,const int& n);
int    ToInt(const double& x);
double ArcCos(const double& x);
void InverseMatrix(double a[4][4]);

void Srand();
double Random(const double& semi_width,const int& IsEndNotUni);
double Random(const double& limit1,const double& limit2,const int& IsEndNotUni);
double RandomGauss(const double& sigma);

long rgb(const long& r,const long& g,const long& b);
long Rrgb(const long& rgb);
long Grgb(const long& rgb);
long Brgb(const long& rgb);
long RGBComplement(const long& rgb);

double Max(const double& x1,const double& x2);
int    Max(const int& i1,const int& i2);
double Max(const double& x1,const double& x2,const double& x3);
double Max(const double& x1,const double& x2,const double& x3,const double& x4);
double Min(const double& x1,const double& x2);
int    Min(const int& i1,const int& i2);
double Min(const double& x1,const double& x2,const double& x3);
double Min(const double& x1,const double& x2,const double& x3,const double& x4);

template<class T> void Sort(T *buf,int n);
double Median(double *buf,int n);

int QuadraticEq(double& x1,double& x2,double a,double b,double c);
double QuadraticX1(double a,double b,double c);
double QuadraticX2(double a,double b,double c);

void DFT(complex *a,int n,int inv,int optical,double zoom);
void DFT(double *real,double *image,int n,int inv,int optical,double zoom);
void DFTRow(complex **a,int m,int n,int inv,int optical,double zoom);
void DFTColumn(complex **a,int m,int n,int inv,int optical,double zoom);
void DFT(complex **a,int m,int n,int inv,int optical,double zoom);

void HilbertT(complex *a,int n);
void HilbertT(double *real,double *image,int n);

void DispersionCor(double *a,int n,double A2,double A3);

void EnFace(int n,int m,std::string in_filename,std::string out_filename,double gamma=1);

double SphSag(double r,double h);
double OverlapAreaCircles(double r,double a);

double Distance(double &dx,double &dy,double &dz,double x1,double y1,double z1,double x2,double y2,double z2);
double DistanceX(double x1,double y1,double z1,double x2,double y2,double z2);
double DistanceY(double x1,double y1,double z1,double x2,double y2,double z2);
double DistanceZ(double x1,double y1,double z1,double x2,double y2,double z2);

void IntersectionsCircles(double &x1,double &y1,double &x2,double &y2,
                          double X1,double Y1,double R1,double X2,double Y2,double R2);
std::string IntersectionsCircles(double X1,double Y1,double R1,double X2,double Y2,double R2);
double Circle3Points(double& xc,double& yc,double& r,
	                 double x1,double y1,double x2,double y2,double x3,double y3);

double TriangleArea(double x1,double y1,double z1,double x2,double y2,double z2,
					double x3,double y3,double z3);
double QuadrangleArea(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4);
double PolygonArea(int n,double *x,double *y);
double PolygonPointDistance(int n,double *x,double *y,double X,double Y);

void ParaxialR(double &r1,double &r2,double &axis_deg,double axx,double ayy,double axy);

vector<double> NearestPointOnLine(vector<double> P0,vector<double> P1,vector<double> P2);
vector<double> NearestPointLineLine(vector<double> P,vector<double> V,vector<double> Po,vector<double> Vo);
double DistancePointLinesegment(vector<double> P0,vector<double> P1,vector<double> P2);
bool IsCrossLines(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4);

matrix<double> Tmatrix(double rox,double roy,double roz);

double Spline(double *x,double *y,int N,double xx,int SPDerivativeZero=0,int EPDerivativeZero=0);
double Lagrange(double *x,double *y,int N,double xx);

void Unwrap(double &ppre,double p0,double &psuc);
void Unwrap(double *th_deg,int N);

void AddSCA(double& S,double& C,double& A,
            double S1,double C1,double A1,char plus_minus,double S2,double C2,double A2);

double dBToT(double dB);

int CSVReadLine(std::ifstream& from,double* data,int datas,int& valid_line);

std::string scmd_general_func(const std::string& com,int val);


////  テンプレートの定義 //////////////////////////////////////////////////////

template <class T> void Dim(T** &p,int m,int n){
	// 2次元配列 p[0][0]～p[m][n] を生成し，各要素を0に初期化する．
	int i,j;

	p=new T*[m+1];           // p[0],p[1], ... p[m]
	for(i=0; i<=m; ++i){ 
		p[i]=new T[n+1];     // p[i][0],p[i][1], ... p[i][n]
		for(j=0 ;j<=n; ++j) p[i][j]=0;
	}
}

template <class T> void Erase(T** &p,int m){
	// 2次元配列 p[m][] を開放する．
	int i;

	for(i=0; i<=m; ++i) delete[] p[i];
	delete[] p;
}

template <class T> void Dim3(T*** &p,int n1,int n2,int n3){
	// 3次元配列 p[0][0][0]～p[n1][n2][n3] を生成し，各要素を0に初期化する．
	int i,j,k;

	p=new T**[n1+1];              // p[0],p[1], ... p[n1]
	for(i=0; i<=n1; ++i){ 
		p[i]=new T*[n2+1];        // p[i][0],p[i][1], ... p[i][n2]
		for(j=0; j<=n2; ++j){
			p[i][j]=new T[n3+1];  // p[i][j][0],p[i][j][1], ... p[i][j][n3]
			for(k=0; k<=n3; ++k) p[i][j][k]=0;
		}
	}
}

template <class T> void Erase3(T*** &p,int n1,int n2){
	// 3次元配列 p[n1][n2][] を開放する．
	int i,j;
	
	for(i=0; i<=n1; ++i){
		for(j=0; j<=n2; ++j) delete [] p[i][j];
		delete [] p[i];
	}
	delete[] p;
}


template<class T> void Sort(T *buf,int n){
	// buf[0]からbuf[n-1]までのデータを小さい順にソートし，
	// 結果を上書きする．
/*
	int i,j;

	for( i=0; i<=n-2; i++ ) {
		for( j=i+1; j<=n-1; j++ ) {
			if(buf[i]>buf[j]){
				Swap(buf[i],buf[j]);
			}
		}
	}
*/
	
	///// 高速化 ヒープソート /////////////////////////////
	
	int i,j,k,l;
	T *heap;

	heap=new T [n+1];  // heap[1]を根とするヒープ用配列

	for(l=1; l<=n; l++){    // ヒープに値を入れていく
		heap[l]=buf[l-1];
		i=l;
		j=i/2;                           // heap[i]の親はheap[i/2]
		while(j>=1 && heap[j]>heap[i]){  // 親の方が大きければ子と値を交換
			Swap(heap[i],heap[j]);
			i=j;
			j=i/2;
		}
	}

	for(l=1; l<=n; l++){     // 根の値をリストに加えていく
		buf[l-1]=heap[1];
		k=n-l;
		heap[1]=heap[k+1];   // ヒープの末尾の値を根に代入
		i=1; j=i*2;          // heap[i]の子はheap[i*2],heap[i*2+1]
		while(j<=k){
			if(j+1<=k && heap[j]>heap[j+1]) j++;   // 値が小さい方の子をheap[j]とする
			if(heap[i]>heap[j]){                   // 子の方が小さければ親と値を交換
				Swap(heap[i],heap[j]);
			}
			i=j;
			j=i*2;
		}
	}

	delete [] heap;
	
  //////////////////////////////////////////////////////
}

#endif // #ifndef GENERAL_FUNC_H_INCLUDED