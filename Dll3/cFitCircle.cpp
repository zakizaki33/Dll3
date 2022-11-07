#include "stdafx.h"
#include "MyDllOptics.h"
#include "cFitCircle.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

double cFitCircle::f(double x,double y,double dummy1,double dummy2){
	// f = (x-x0)^2 + (y-y0)^2 -r^2
	// �ŁCx0,y0,r��ϐ��ɂ���Ɛ��`�ł͂Ȃ��̂ŁC�œK���͒����ڋ߂ɂȂ�D
	// �������C f= x^2+2*x*x0+x0^2 +y^2+2*y*y0+y0^2 -r^2 �ɂ����āC
	//     a = -2*x0, b=-2*y0, c=x0^2+y0^2-r^2
	// �Ƃ����΁C f= x^2 +y^2 +a*x +b*y +c �ƂȂ�C�ϐ�a,b,c�ɂ��Đ��`�ƂȂ�D
	return x*x+y*y+a*x+b*y+c;
}

void cFitCircle::SetA(){
	for(int i=1; i<=GetNumberOfData(); ++i){
		A[i][1]=x[i];
		A[i][2]=y[i];
		A[i][3]=1;
	}
}

void cFitCircle::SetX(){
	X[1]=&a;
	X[2]=&b;
	X[3]=&c;
}

int cFitCircle::NumberOfCoefficients(){
	return 3;
}

cFitCircle::cFitCircle() : x(x1),y(x2),z(x3) {
	r=x0=y0=0;
	a=b=c=0;
}

double cFitCircle::GetR(){
	return r;
}
void cFitCircle::SetR(double value){
	r=value;
	c=x0*x0+y0*y0-r*r;
}

double cFitCircle::GetX0(){
	return x0;
}
void cFitCircle::SetX0(double value){
	x0=value;
	a=-2*x0;
	c=x0*x0+y0*y0-r*r;
}

double cFitCircle::GetY0(){
	return y0;
}
void cFitCircle::SetY0(double value){
	y0=value;
	b=-2*y0;
	c=x0*x0+y0*y0-r*r;
}

double cFitCircle::GetX(int i){ return x[i]; }
void cFitCircle::SetX(int i,double value){ x[i]=value; }

double cFitCircle::GetY(int i){ return y[i]; }
void cFitCircle::SetY(int i,double value){ y[i]=value; }

double cFitCircle::yApproximate(double x,int sgn){
	// sgn>=0 �̂Ƃ���y0��萳�����̉���Ԃ�
	// sgn<0  �̂Ƃ���y0��蕉�����̉���Ԃ�
	double v;
	
	v=r*r-(x-x0)*(x-x0);
	if(v>=0){
		return y0+(::sgn(sgn))*sqrt(v);
	}
	else{
		return 0;
	}
}

int cFitCircle::ModifyCoefficients(){
	int result;

	result=cLeastSquares::ModifyCoefficients(0);
	x0=-a/2;
	y0=-b/2;
	r=sqrt(a*a/4+b*b/4-c);
	return result;
}