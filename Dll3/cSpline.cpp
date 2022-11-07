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
	// �敪3���֐��̌W�� a[i],b[i],c[i],d[i] ���v�Z����
	int i;
	double *h,*u, DX;

	if(N<3) return;

	u= new double[N+1];  // u[1]�`u[N], ��1�f�[�^�����N�f�[�^(�f�[�^��N)
	h= new double[N+1];  // h[1]�`h[N-1], ��1��Ԃ����N-1���(��Ԑ� N-1)

	matrix<double> A(N-2,N-2);   // 1����N-2, �f�[�^��N���2���Ȃ�
	matrix<double> U(N-2,1);
	matrix<double> V(N-2,1);

	for(i=1; i<=N-1; i++) h[i]=_X[i+1]-_X[i];
	
	for(i=1; i<=N-2; i++) A.a[i][i]=2*(h[i]+h[i+1]);
	for(i=2; i<=N-2; i++) A.a[i][i-1]=h[i];
	for(i=1; i<=N-3; i++) A.a[i][i+1]=h[i+1];
	
	for(i=1; i<=N-2; i++) V.a[i][1]=6*( (_Y[i+2]-_Y[i+1])/h[i+1] - (_Y[i+1]-_Y[i])/h[i] );
	
	U=inv(A)*V;  // ���[�̏����_�ł�2�����֐������܂���

	for(i=1; i<=N; i++){ // ��1�f�[�^�����N�f�[�^
		if(i==1 || i==N){
			u[i]=0;  // ���[��2�����֐���0
		}
		else{
			u[i]=U[i-1][1];
		}
	}

	for(i=1; i<=N-1; i++){ // ��1��Ԃ����N-1���
		a[i]=(u[i+1]-u[i])/6/(_X[i+1]-_X[i]);
		b[i]=u[i]/2;
		c[i]=(_Y[i+1]-_Y[i])/(_X[i+1]-_X[i])-(_X[i+1]-_X[i])*(2*u[i]+u[i+1])/6;
		d[i]=_Y[i];
	}

	if(StartPointDerivativeZero){
		// �n�_��1�����֐���0�̂Ƃ���a[1],b[1],c[1],d[1]���Čv�Z���ď㏑������
		DX=_X[2]-_X[1];
		d[1]=_Y[1];
		c[1]=0;
		b[1]=(6*(_Y[2]-_Y[1])-u[2]*DX*DX)/(4*DX*DX);
		a[1]=(u[2]-2*b[1])/(6*DX);
	}

	if(EndPointDerivativeZero){
		// �I�_��1�����֐���0�̂Ƃ���a[N-1],b[N-1],c[N-1],d[N-1]���Čv�Z���ď㏑������
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
	N=1;             // NotNull() ��false�ɂȂ�
	ref_public=false;
	alloc();
	StartPointDerivativeZero=EndPointDerivativeZero=0;  // ���[�_�ł�2�����֐���0�i��ʓI�Ȑݒ�j
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
	// ���̍s�̑���� SetN(x.GetN()); �Ƃ��Ă͂����Ȃ��Dnew ���Ă��Ȃ��|�C���^�� delete �����邱�ƂɂȂ邽�߁C
	// ���s���ɃG���[���o�邱�Ƃ�����i�o�Ȃ����Ƃ�����̂Ŗ��ȃo�O�ƂȂ�j�D

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
		ref_public=true;  // �Q�Ƃ���U���J���Ă��܂��ƁC�{�֐����o�Ȃ��Œl��ύX�ł��Ă��܂�
		return _X[i];
	}
	else{
		return _X[0];  // �g�p���Ă��Ȃ�
	}
}

double& cSpline::Y(int i){
	if(1<=i && i<=N){
		ref_public=true;  // �Q�Ƃ���U���J���Ă��܂��ƁC�{�֐����o�Ȃ��Œl��ύX�ł��Ă��܂�
		return _Y[i];
	}
	else{
		return _Y[0];  // �g�p���Ă��Ȃ�
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
	// y��1�������W��
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
	// y��2�������W��
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
	// ��F data = "0 0; 1 0.1; 2 0.4; 3 0.9; ... 10 10.0;"
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
	// �f�[�^����2�{�ɂ���D�אڂ���f�[�^�̒��Ԃɕ�ԃf�[�^��}�����邱�Ƃɂ��D
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
	// X[i]�̃X�e�b�v��dx�ɂ���D�������C���[(X[1]��X[N])�͕ς��Ȃ��D
	int i, new_n;
	double *_X1,*_Y1;
	
	for(i=1; i<=MAX_SIZE; ++i){
		if( Round(_X[1]+dx*(i-1),5) >= Round(_X[N],5) ) break;  // �v�Z�덷�ɔ�����Round�֐����g��
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