#include "stdafx.h"
#include "cLeastSquares.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#include <math.h>

cLeastSquares::cLeastSquares(){}

cLeastSquares::~cLeastSquares(){}

int cLeastSquares::GetNumberOfData() const{
	return x1.GetSize();
}
void cLeastSquares::SetNumberOfData(int m){
	int i;

	x1.RemoveAll();
	x2.RemoveAll();
	x3.RemoveAll();
	x4.RemoveAll();
	w.RemoveAll();
	for(i=1; i<=m; i++){
		x1.AddTail(0);
		x2.AddTail(0);
		x3.AddTail(0);
		x4.AddTail(0);
		w.AddTail(1);
	}
}

double cLeastSquares::Get_x1(int i) const{ 
	if(1<=i && i<=GetNumberOfData()) return x1[i];
	else return 0; 
}
void cLeastSquares::Set_x1(int i,double value){
	if(1<=i && i<=GetNumberOfData()) x1[i]=value;
}

double cLeastSquares::Get_x2(int i) const{ 
	if(1<=i && i<=GetNumberOfData()) return x2[i];
	else return 0; 
}
void cLeastSquares::Set_x2(int i,double value){
	if(1<=i && i<=GetNumberOfData()) x2[i]=value;
}

double cLeastSquares::Get_x3(int i) const{ 
	if(1<=i && i<=GetNumberOfData()) return x3[i];
	else return 0; 
}
void cLeastSquares::Set_x3(int i,double value){
	if(1<=i && i<=GetNumberOfData()) x3[i]=value;
}

double cLeastSquares::Get_x4(int i) const{ 
	if(1<=i && i<=GetNumberOfData()) return x4[i];
	else return 0; 
}
void cLeastSquares::Set_x4(int i,double value){
	if(1<=i && i<=GetNumberOfData()) x4[i]=value;
}

double cLeastSquares::GetWeight(int i) const{ 
	if(1<=i && i<=GetNumberOfData()) return w[i];
	else return 0; 
}
void cLeastSquares::SetWeight(int i,double value){
	if(1<=i && i<=GetNumberOfData()) w[i]=value;
}
void cLeastSquares::SetWeight(double value){
	// 最終データのウエイトを設定する
	w[GetNumberOfData()]=value;
}


void cLeastSquares::AddData(double x1,double x2,double x3,double x4){
	this->x1.AddTail(x1);
	this->x2.AddTail(x2);
	this->x3.AddTail(x3);
	this->x4.AddTail(x4);
	this->w.AddTail(1);
}

void cLeastSquares::AddData(double x1,double x2,double x3){
	AddData(x1,x2,x3,0);
}

void cLeastSquares::AddData(double x1,double x2){
	AddData(x1,x2,0,0);
}

void cLeastSquares::DataClear(){
	SetNumberOfData(0);
}

double cLeastSquares::Error(int i){
	if(1<=i && i<=GetNumberOfData()) {
		return f(x1[i],x2[i],x3[i],x4[i]);
	}
	else return 0;
}

double cLeastSquares::rmsError(){
	int i;
	double err=0,a=0;
	for(i=1; i<=GetNumberOfData(); ++i){
		err+=f(x1[i],x2[i],x3[i],x4[i])*f(x1[i],x2[i],x3[i],x4[i]) *w[i]*w[i];
		a+=w[i];
	}
	err=sqrt(err/a);
	return err;
}

double cLeastSquares::pvError(){
	int i;
	double p=-1e30;
	double v= 1e30;
	for(i=1; i<=GetNumberOfData(); ++i){
		if( f(x1[i],x2[i],x3[i],x4[i])>p ) p=f(x1[i],x2[i],x3[i],x4[i]);
		if( f(x1[i],x2[i],x3[i],x4[i])<v ) v=f(x1[i],x2[i],x3[i],x4[i]);
	}
	return p-v;
}

double cLeastSquares::maxError(){
	// 全体のオフセットも評価するときは，pvエラーでは対応できない．
	// (オフセットが大きいとき，pvエラーも大きいとは限らない)
	int i;
	double a,max=-1e30;
	for(i=1; i<=GetNumberOfData(); ++i){
		a=fabs(f(x1[i],x2[i],x3[i],x4[i]));
		if( a>max ) max=a;
	}
	return max;
}

int cLeastSquares::ModifyCoefficients(double rho) {
	int i,j,count;
	int m=GetNumberOfData();
	int n=NumberOfCoefficients();
	double mf_start,mf_end, mf,mf1;

	A.redim(m,n);
	F.redim(m,1);
	X=new double* [n+1];
	SetX();

	// メリット関数の初期値
	//   メリット関数は Σ{w[i]^2 * f()^2} であるので，
	//   例えば，w[i]=2 の測定点1個は，w[i]=1 の測定点2個ではなく4個と等価であることに注意する．
	mf_start=0;
	for(i=1; i<=m; ++i){
		mf_start+=f(x1[i],x2[i],x3[i],x4[i])*f(x1[i],x2[i],x3[i],x4[i]) *w[i]*w[i];
	}

	for(count=1; count<=(rho==0 ? 1:20); ++count){  // とりあえず20回繰り返す
		                                            // ρ=0のときは線形とみなし1回とする
		SetA();
		
		for(i=1; i<=m; ++i) for(j=1; j<=n; ++j) A[i][j]*=w[i];
		// 注：
		//  matrix<double>W(m,m) とし，対角成分W[i][i]=w[i]とすれば，
		//  行列積A=W*Aが使えてコードがきれいに書けるが，
		//  このときWは対角成分以外が0なので，多くのメモリが無駄になる．
		//  特に，データ数が数万になると実行時にメモリ確保できずエラーが発生する．

		mf=0;
		for(i=1; i<=m; ++i){
			F[i][1]=-f(x1[i],x2[i],x3[i],x4[i]) *w[i];
			mf+=F[i][1]*F[i][1];
		}

		dX=inv( t(A)*(A)+rho*unit(t(A)*(A)) )*t(A)*F;

		for(i=1; i<=n; ++i){
			*X[i]+=dX[i][1];
		}
		
		mf1=0;
		for(i=1; i<=m; ++i){
			mf1+=f(x1[i],x2[i],x3[i],x4[i])*f(x1[i],x2[i],x3[i],x4[i]) *w[i]*w[i];
		}
		
		if(mf1<=mf){
			// メリット関数が増加しなければdamping factorを小さくする．
			// 【説明】
			//       mf1<mf とすると，ρが大きいとき，ρを小さくしてもmfが変化せず，
			//     十分最適化しないまま停止してしまうことがある．
			//       また，次のcountに移る時，*X[i]を元に戻さないため，
			//     一定の地点でmfが最小になるρを探すのではなく，移動しつつρも緩めていることになる．
			//       また，ρが大きくなる方向には探索しないため，引数ρは十分大きい
			//     ことが前提である．
			rho*=0.1;
		}
		else{
			// メリット関数が増加したらXを戻し，damping factorを大きくする．
			for(i=1; i<=n; ++i){
				*X[i]-=dX[i][1];
			}
			rho*=10;
		}
	}

	// メリット関数の終了値
	mf_end=0;
	for(i=1; i<=m; ++i){
		mf_end+=f(x1[i],x2[i],x3[i],x4[i])*f(x1[i],x2[i],x3[i],x4[i]) *w[i]*w[i];
	}

	delete [] X;

	// 改良されたかどうかを返す
	return mf_end<mf_start ? 1 : 0;
}
