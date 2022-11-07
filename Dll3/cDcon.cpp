#include "stdafx.h"
#include "MyDllOptics.h"
#include "cDcon.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

double cDcon::QconCoefficient(int m, int i){
	// Qcon(m,x) = Σ[i=0 to m] QconCoefficient(m,i) * x^i
	return (is_even(m-i) ? 1.0 : -1.0) *factorial(4+m+i)
	       /factorial(i) /factorial(4+i) /factorial(m-i);
}

double cDcon::Qcon(int m,double x,int me/*=-1*/){
	// Robust and fast computation for the polynomials of optics
	// G.W.Forbes 
	// 7 June 2010 / Vol. 18, No. 13 / OPTICS EXPRESS 13851
	// (1.7)式より
	// me>=0のときは m=meの項のみ抽出する．
	int i;
	double a;

	a=0;
	for(i=0; i<=m; ++i){
		if(me<0 || i==me){
			a+=QconCoefficient(m,i)*pw(x,i);
		}
	}
	return a;
}

double cDcon::Qcon1(int m,double x){
	// 1階微分
	int i;
	double a;

	a=0;
	for(i=1; i<=m; ++i){
		a+=QconCoefficient(m,i)*((double)i) *pw(x,i-1);
	}
	return a;
}

double cDcon::Qcon2(int m,double x){
	// 2階微分
	int i;
	double a;

	a=0;
	for(i=2; i<=m; ++i){
		a+=QconCoefficient(m,i)*((double)i)*((double)(i-1))*pw(x,i-2);
	}
	return a;
}

cDcon::cDcon(){
	int i;

	Rn=0;
	for(i=0; i<=M; ++i) A[i]=0; 
}

void cDcon::SetRn(double val){
	// Rnを新しい値valとして，それに応じてA[0],A[1],..A[M]を変換する．
	// 例えば，M=3として今のRnと新しいRnで，v=(R/Rn)^2を使って各次数(4,6,8,10..次)の式を書き出し，
	// 恒等式としてvのべきに関して係数比較を行なえば，
	// 少し面倒だが以下の式を帰納できる．
	int m,i;
	double *A1;
	double sum1,sum2;

	if(Rn==0){
		Rn=val;
	}
	else{
		A1=new double[M+1]; 
		for(m=M; m>=0; --m){
			sum1=0; for(i=m;   i<=M; ++i) sum1+=A[i]*QconCoefficient(i,m);
			sum2=0; for(i=m+1; i<=M; ++i) sum2+=A1[i]*QconCoefficient(i,m);
			A1[m]=pw(val,4+2*m)/QconCoefficient(m,m) * ( sum1/pw(Rn,4+2*m)-sum2/pw(val,4+2*m) );
		}
		for(m=0; m<=M; ++m) A[m]=A1[m];
		Rn=val;
		delete [] A1;
	}
}

double& cDcon::A4() { return A[0]; }
double& cDcon::A6() { return A[1]; }
double& cDcon::A8() { return A[2]; }
double& cDcon::A10(){ return A[3]; }
double& cDcon::A12(){ return A[4]; }
double& cDcon::A14(){ return A[5]; }
double& cDcon::A16(){ return A[6]; }
double& cDcon::A18(){ return A[7]; }
double& cDcon::A20(){ return A[8]; }

void cDcon::Clear(){
	int i;

	Rn=0;
	for(i=0; i<=M; ++i) A[i]=0;	
}

int cDcon::NotNull(){
	int i;
	
	if(Rn!=0) return true;
	for(i=0; i<=M; ++i){
		if(A[i]!=0) return true;
	}
	return false;
}

cDcon cDcon::reversed(){
	int i;
	cDcon x=*this;

	for(i=0; i<=M; ++i){
		x.A[i]=-x.A[i];
	}
	return x;
}

void cDcon::scale(double m){
	int i;
	if(m==0) return;

	Rn*=m;
	for(i=0; i<=M; ++i){
		A[i]*=m;
	}
}

double cDcon::a4() { return Z(1,0,0); }
double cDcon::a6() { return Z(1,0,1); }
double cDcon::a8() { return Z(1,0,2); }
double cDcon::a10(){ return Z(1,0,3); }
double cDcon::a12(){ return Z(1,0,4); }
double cDcon::a14(){ return Z(1,0,5); }
double cDcon::a16(){ return Z(1,0,6); }
double cDcon::a18(){ return Z(1,0,7); }
double cDcon::a20(){ return Z(1,0,8); }

void cDcon::PSeriesToDcon(int n,double *a){
	// z = a[0]*h^4 +a[1]*h^6 +a[2]*h^8 + .. +a[(n-4)/2]*h^n
	// であるとき，このzをcDconに加える．
	int M,m,i;
	double A1[this->M+1], z;

	M=(n-4)/2;

	for(m=M; m>=0; --m){
		z=a[m]*pw(Rn,4+m*2);
		for(i=M; i>m; --i){
			z-=A1[i]*Qcon(i,1,m);
		}
		A1[m]=z/Qcon(m,1,m);
	}
	for(m=0; m<=M; ++m){
		A[m]+=A1[m];
	}
}

double cDcon::Z(double x,double y,int me/*=-1*/){
	// Shape specification for axially symmetirc optical surfaces
	// G.W.Forbes
	// 16 Apr 2007 / Vol. 15, No. 8 / OPTICS EXPRESS 5218
	// (2)(3)式より（u^2=vとした）．
	// me>=0のときは i=meの項のみ抽出する．
	// (Δz=a4*h^4+a6*h^6 + ... としたときの a4,a6,...を求めるのに使用）
	int i;
	double v,a;
	
	if(Rn==0) return 0;
	v=(x*x+y*y)/Rn/Rn;
	a=0;
	for(i=0; i<=M; ++i){
		if(A[i]!=0){
			a+=v*v*A[i]*Qcon(i,v,me);  // この式より，最低次数はxやyの4乗となる．したがって近軸関係に影響しない
		}
	}
	return a;
}

double cDcon::Zx(double x,double y){
	int i;
	double v,vx,a;
	
	if(Rn==0) return 0;
	v=(x*x+y*y)/Rn/Rn;
	vx=2*x/Rn/Rn;
	a=0;
	for(i=0; i<=M; ++i){
		if(A[i]!=0){
			a+=2*v*vx*A[i]*Qcon(i,v)
			  +v*v*A[i]*Qcon1(i,v)*vx;
		}
	}
	return a;
}

double cDcon::Zy(double x,double y){
	int i;
	double v,vy,a;
	
	if(Rn==0) return 0;
	v=(x*x+y*y)/Rn/Rn;
	vy=2*y/Rn/Rn;
	a=0;
	for(i=0; i<=M; ++i){
		if(A[i]!=0){
			a+=2*v*vy*A[i]*Qcon(i,v)
			  +v*v*A[i]*Qcon1(i,v)*vy;
		}
	}
	return a;
}

double cDcon::Zxx(double x,double y){
	int i;
	double v,vx,vxx,a;
	
	if(Rn==0) return 0;
	v=(x*x+y*y)/Rn/Rn;
	vx=2*x/Rn/Rn;
	vxx=2/Rn/Rn;
	a=0;
	for(i=0; i<=M; ++i){
		if(A[i]!=0){
			a+=(2*vx*vx+2*v*vxx)*A[i]*Qcon(i,v)
			 +2*(2*v*vx*A[i]*Qcon1(i,v)*vx)   // (uv)" = u"v+2u'v'+uv"
			 +v*v*A[i]*(Qcon2(i,v)*vx*vx+Qcon1(i,v)*vxx);
		}
	}
	return a;
}

double cDcon::Zxy(double x,double y){
	int i;
	double v,vx,vy,a;
	
	if(Rn==0) return 0;
	v=(x*x+y*y)/Rn/Rn;
	vx=2*x/Rn/Rn;
	vy=2*y/Rn/Rn;
	a=0;
	for(i=0; i<=M; ++i){
		if(A[i]!=0){
			a+=2*vy*vx*A[i]*Qcon(i,v)
			  +2*v*vx*A[i]*Qcon1(i,v)*vy
			  +2*v*vy*A[i]*Qcon1(i,v)*vx
			  +v*v*A[i]*Qcon2(i,v)*vx;
		}
	}
	return a;
}

double cDcon::Zyy(double x,double y){
	int i;
	double v,vy,vyy,a;
	
	if(Rn==0) return 0;
	v=(x*x+y*y)/Rn/Rn;
	vy=2*y/Rn/Rn;
	vyy=2/Rn/Rn;
	a=0;
	for(i=0; i<=M; ++i){
		if(A[i]!=0){
			a+=(2*vy*vy+2*v*vyy)*A[i]*Qcon(i,v)
			  +2*(2*v*vy*A[i]*Qcon1(i,v)*vy)   // (uv)" = u"v+2u'v'+uv"
			  +v*v*A[i]*(Qcon2(i,v)*vy*vy+Qcon1(i,v)*vyy);
		}
	}
	return a;
}

std::string cDcon::GetTerms(){
	char buf[100];
	std::string s;

	if(NotNull()){
		sprintf(buf,"Rn %g; ",Rn); s+=buf;
		if(A4()!=0){ sprintf(buf,"A4 %g; ", A4()); s+=buf; }
		if(A6()!=0){ sprintf(buf,"A6 %g; ", A6()); s+=buf; }
		if(A8()!=0){ sprintf(buf,"A8 %g; ", A8()); s+=buf; }
		if(A10()!=0){ sprintf(buf,"A10 %g; ",A10()); s+=buf; }
		if(A12()!=0){ sprintf(buf,"A12 %g; ",A12()); s+=buf; }
		if(A14()!=0){ sprintf(buf,"A14 %g; ",A14()); s+=buf; }
		if(A16()!=0){ sprintf(buf,"A16 %g; ",A16()); s+=buf; }
		if(A18()!=0){ sprintf(buf,"A18 %g; ",A18()); s+=buf; }
		if(A20()!=0){ sprintf(buf,"A20 %g; ",A20()); s+=buf; }
	}
	return s;
}

void cDcon::SetTerms(std::string terms){
	std::string s,s0,s1;
	bool b1;
	int k;
	
	Clear();
	for(k=1; k<=sentences(terms); ++k){
		s=sentence(terms,k);
		s0=arg(s,0);
		s1=arg(s,1); b1=is_numeric(s1);

		if((s0=="Rn" || s0=="rn") && b1) Rn=atof(s1.c_str());
		if((s0=="A4" || s0=="a4") && b1) A4()=atof(s1.c_str());
		if((s0=="A6" || s0=="a6") && b1) A6()=atof(s1.c_str());
		if((s0=="A8" || s0=="a8") && b1) A8()=atof(s1.c_str());
		if((s0=="A10" || s0=="a10") && b1) A10()=atof(s1.c_str());
		if((s0=="A12" || s0=="a12") && b1) A12()=atof(s1.c_str());
		if((s0=="A14" || s0=="a14") && b1) A14()=atof(s1.c_str());
		if((s0=="A16" || s0=="a16") && b1) A16()=atof(s1.c_str());
		if((s0=="A18" || s0=="a18") && b1) A18()=atof(s1.c_str());
		if((s0=="A20" || s0=="a20") && b1) A20()=atof(s1.c_str());
	}
}

