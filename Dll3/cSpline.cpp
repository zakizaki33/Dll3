#include "stdafx.h"
#include "cSpline.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

void cSpline::alloc(){
	int i;

	for(i=1; i<=N; ++i) _X[i]=_Y[i]=0;
	
	a= new double[N];
	b= new double[N];
	c= new double[N];
	d= new double[N];

	for(i=1; i<=N-1; ++i) a[i]=b[i]=c[i]=d[i]=0;
}

void cSpline::erase(){
	delete [] a;
	delete [] b;
	delete [] c;
	delete [] d;
}

void cSpline::calc(){
	// 区分3次関数の係数 a[i],b[i],c[i],d[i] を計算する
	int i;
	double *h,*u, DX;

	if(N<3) return;

	u= new double[N+1];  // u[1]～u[N], 第1データから第Nデータ(データ数N)
	h= new double[N+1];  // h[1]～h[N-1], 第1区間から第N-1区間(区間数 N-1)

	matrix<double> A(N-2,N-2);   // 1からN-2, データ数Nより2少ない
	matrix<double> U(N-2,1);
	matrix<double> V(N-2,1);

	for(i=1; i<=N-1; i++) h[i]=_X[i+1]-_X[i];
	
	for(i=1; i<=N-2; i++) A.a[i][i]=2*(h[i]+h[i+1]);
	for(i=2; i<=N-2; i++) A.a[i][i-1]=h[i];
	for(i=1; i<=N-3; i++) A.a[i][i+1]=h[i+1];
	
	for(i=1; i<=N-2; i++) V.a[i][1]=6*( (_Y[i+2]-_Y[i+1])/h[i+1] - (_Y[i+1]-_Y[i])/h[i] );
	
	U=inv(A)*V;  // 両端の除く点での2次導関数が求まった

	for(i=1; i<=N; i++){ // 第1データから第Nデータ
		if(i==1 || i==N){
			u[i]=0;  // 両端の2次導関数は0
		}
		else{
			u[i]=U[i-1][1];
		}
	}

	for(i=1; i<=N-1; i++){ // 第1区間から第N-1区間
		a[i]=(u[i+1]-u[i])/6/(_X[i+1]-_X[i]);
		b[i]=u[i]/2;
		c[i]=(_Y[i+1]-_Y[i])/(_X[i+1]-_X[i])-(_X[i+1]-_X[i])*(2*u[i]+u[i+1])/6;
		d[i]=_Y[i];
	}

	if(StartPointDerivativeZero){
		// 始点で1次導関数が0のときはa[1],b[1],c[1],d[1]を再計算して上書きする
		DX=_X[2]-_X[1];
		d[1]=_Y[1];
		c[1]=0;
		b[1]=(6*(_Y[2]-_Y[1])-u[2]*DX*DX)/(4*DX*DX);
		a[1]=(u[2]-2*b[1])/(6*DX);
	}

	if(EndPointDerivativeZero){
		// 終点で1次導関数が0のときはa[N-1],b[N-1],c[N-1],d[N-1]を再計算して上書きする
		DX=_X[N]-_X[N-1];
		d[N-1]=_Y[N-1];
		b[N-1]=u[N-1]/2;
		c[N-1]=(3*_Y[N]-b[N-1]*DX*DX-3*d[N-1])/(2*DX);
		a[N-1]=-(2*b[N-1]*DX+c[N-1])/(3*DX*DX);
	}
	
	delete [] u;
	delete [] h;
}


// public menbers /////////////////////////////////////////////////////////

cSpline::cSpline(){
	N=1;             // NotNull() はfalseになる
	ref_public=false;
	alloc();
	StartPointDerivativeZero=EndPointDerivativeZero=0;  // 両端点では2次導関数が0（一般的な設定）
}

cSpline::cSpline(const cSpline& x){
	*this=x;
	ref_public=false;
}

cSpline::~cSpline(){
	erase();
}

cSpline& cSpline::operator =(const cSpline &x){
	int i;

	N=x.GetN(); alloc();
	// この行の代わりに SetN(x.GetN()); としてはいけない．new していないポインタに delete をすることになるため，
	// 実行時にエラーが出ることがある（出ないこともあるので厄介なバグとなる）．

	for(i=1; i<=N; ++i){
		_X[i]=x._X[i];
		_Y[i]=x._Y[i];
	}
	StartPointDerivativeZero=x.StartPointDerivativeZero;
	EndPointDerivativeZero=x.EndPointDerivativeZero;
	calc();

	return *this;
}

bool cSpline::NotNull(){
	return N>2 ? true : false;
}

cSpline cSpline::reversed() const{
	int i;
	cSpline x=*this;

	for(i=1; i<=x.N; ++i){
		x._Y[i]=-x._Y[i];
	}
	x.calc();
	return x;
}

void cSpline::scale(double m){
	int i;
	if(m==0) return;

	for(i=0; i<=N; ++i){
		_X[i]*=m;
		_Y[i]*=m;
	}
	calc();
}

int cSpline::GetN() const{
	return N;
}

void cSpline::SetN(int val){
	if(val<=MAX_SIZE){
		erase();
		N=val;
		alloc();
		calc();
	}
}

double cSpline::GetX(int i){
	if(1<=i && i<=N){
		return _X[i];
	}
	else{
		return 0;
	}
}

void cSpline::SetX(int i,double val){
	if(1<=i && i<=N){
		_X[i]=val;
		calc();
	}
}

double cSpline::GetY(int i){
	if(1<=i && i<=N){
		return _Y[i];
	}
	else{
		return 0;
	}
}

void cSpline::SetY(int i,double val){
	if(1<=i && i<=N){
		_Y[i]=val;
		calc();
	}
}

double& cSpline::X(int i){
	if(1<=i && i<=N){
		ref_public=true;  // 参照を一旦公開してしまうと，本関数を経ないで値を変更できてしまう
		return _X[i];
	}
	else{
		return _X[0];  // 使用していない
	}
}

double& cSpline::Y(int i){
	if(1<=i && i<=N){
		ref_public=true;  // 参照を一旦公開してしまうと，本関数を経ないで値を変更できてしまう
		return _Y[i];
	}
	else{
		return _Y[0];  // 使用していない
	}
}

double cSpline::y(double x){
	int i;

	if(ref_public){
		calc();
	}

	for(i=1; i<=N-1; i++){
		if( (_X[i]-x)*(_X[i+1]-x)<=0 ){
			return a[i]*(x-_X[i])*(x-_X[i])*(x-_X[i]) +b[i]*(x-_X[i])*(x-_X[i]) +c[i]*(x-_X[i]) +d[i];
		}
	}
	return 0;
}

double cSpline::y1(double x){
	// yの1次微分係数
	int i;

	if(ref_public){
		calc();
	}

	for(i=1; i<=N-1; i++){
		if( (_X[i]-x)*(_X[i+1]-x)<=0 ){
			return 3*a[i]*(x-_X[i])*(x-_X[i]) +2*b[i]*(x-_X[i]) +c[i];
		}
	}
	return 0;
}

double cSpline::y2(double x){
	// yの2次微分係数
	int i;

	if(ref_public){
		calc();
	}

	for(i=1; i<=N-1; i++){
		if( (_X[i]-x)*(_X[i+1]-x)<=0 ){
			return 6*a[i]*(x-_X[i]) +2*b[i];
		}
	}
	return 0;
}

std::string cSpline::GetData(){
	int i;
	char buf[100];
	std::string s;

	if(NotNull()){
		for(i=1; i<=N; ++i){
			sprintf(buf," %.8g %.8g;", _X[i],_Y[i]); s+=buf;
		}
	}
	return s;
}

void cSpline::SetData(std::string data){
	// 例： data = "0 0; 1 0.1; 2 0.4; 3 0.9; ... 10 10.0;"
	int i;
	std::string s,tmp;

	for(i=1; i<=sentences(data); ++i){
		tmp=sentence(data,i);
		if(words(tmp)==2) s+=(tmp+";");
	}

	SetN(sentences(s));

	for(i=1; i<=N; ++i){
		tmp=sentence(s,i);
		SetX(i,atof(word(tmp,1).c_str()));
		SetY(i,atof(word(tmp,2).c_str()));
	}

	calc();
}

void cSpline::DoubleN(){
	// データ数を2倍にする．隣接するデータの中間に補間データを挿入することによる．
	int i,new_n;
	double *_X1,*_Y1;

	new_n=N*2-1;

	if(new_n<=MAX_SIZE){
		_X1=new double [new_n+1];
		_Y1=new double [new_n+1];

		for(i=1; i<=new_n; i+=2){
			_X1[i]=_X[i/2+1];
			_Y1[i]=_Y[i/2+1];
		}

		for(i=2; i<=new_n-1; i+=2){
			_X1[i]=(_X1[i-1]+_X1[i+1])/2;
			_Y1[i]=y(_X1[i]);
		}

		SetN(new_n);
		
		for(i=1; i<=N; ++i){
			_X[i]=_X1[i];
			_Y[i]=_Y1[i];
		}

		delete [] _X1;
		delete [] _Y1;
	}

	calc();
}

void cSpline::SetXStep(double dx){
	// X[i]のステップをdxにする．ただし，両端(X[1]とX[N])は変えない．
	int i, new_n;
	double *_X1,*_Y1;
	
	for(i=1; i<=MAX_SIZE; ++i){
		if( Round(_X[1]+dx*(i-1),5) >= Round(_X[N],5) ) break;  // 計算誤差に備えてRound関数を使う
	}

	new_n=i;

	_X1=new double [new_n+1];
	_Y1=new double [new_n+1];

	_X1[1]=_X[1];     _Y1[1]=_Y[1];
	_X1[new_n]=_X[N]; _Y1[new_n]=_Y[N];

	for(i=2; i<=new_n-1; ++i){
		_X1[i]=_X1[1]+dx*(i-1);
		_Y1[i]=y(_X1[i]);
	}

	SetN(new_n);

	for(i=1; i<=N; ++i){
		_X[i]=_X1[i];
		_Y[i]=_Y1[i];
	}

	calc();

	delete [] _X1;
	delete [] _Y1;
}