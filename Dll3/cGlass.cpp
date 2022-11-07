#include "stdafx.h"
#include "MyDllOptics.h"
#include "cGlass.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

void cGlass::alloc(){
	wl=new double [cn+1];
	n=new double [cn+1];
	ra=new double [cn+1];
}

void cGlass::initialize(){
	int j;
	
	name=name_ini="1";
	nd=nd_ini=1;
	nud=nud_ini=0;
	dn=dn_ini=0;
	is_grin=0;
	for(j=1; j<=cn; ++j){
		wl[j]=546.07;    // e���̔g��
		n[j]=1;
	}

	coordinated=true;
}

void cGlass::free(){
	delete [] wl;
	delete [] n;
	delete [] ra;
}
                         
void cGlass::assignment(const cGlass &x){
	int j;

	name=x.name; name_ini=x.name_ini;
	nd=x.nd; nd_ini=x.nd_ini;
	nud=x.nud; nud_ini=x.nud_ini;
	dn=x.dn; dn_ini=x.dn_ini;
	is_grin=x.is_grin;
	grin_phi=x.grin_phi;
	for(j=1; j<=cn && j<=x.cn; ++j){
		wl[j]=x.wl[j];
		n[j]=x.n[j];
		ra[j]=x.ra[j];
	}
	coordinated=x.coordinated;
}

void cGlass::set_index(){
	// ���ܗ� n[j] �����X�V����Dname �͍X�V����Ă�����̂Ƃ���D
	int j;
	cMaterial m(name);

	for(j=1; j<=cn; ++j) n[j]=Re(m.Index(wl[j],Digits));  // N(i)�ŃA�N�Z�X����Ɩ������[�v
	if(is_grin=m.IsGrin()){
		for(j=1; j<=cn; ++j) ra[j]=m.rA(wl[j]);
		grin_phi=m.GrinPhi;
	}
}

void cGlass::set_index(int j){
	// j�Ԗڂ̔g���ɂ��Ă̂݋��ܗ� n[j] �����X�V����Dname �͍X�V����Ă�����̂Ƃ���D
	cMaterial m(name);

	n[j]=Re(m.Index(wl[j],Digits));  // N(i)�ŃA�N�Z�X����Ɩ������[�v
	if(is_grin=m.IsGrin()){
		ra[j]=m.rA(wl[j]);
		grin_phi=m.GrinPhi;
	}
}

std::string cGlass::values_to_name(double nd,double nud,double dn){
	return trim(str(Round(fabs(nd),-Digits)),0) 
	       + ":" + trim(str(Round(fabs(nud),-(Digits-3))),0) 
	       + " " + trim(str(Round(dn*100000,0)),0);
}

void cGlass::coordinate(){
	// Name, Nd, Nud, dN �𐮍������C���ܗ� n[j] �����X�V����D

	if(name!=name_ini){
		cMaterial m(name);
		double nC,nF;
		
		name_ini = name;
		nd_ini = nd = fabs(Re(m.Index(587.56,Digits)));
		nF=fabs(Re(m.Index(486.13,Digits)));
		nC=fabs(Re(m.Index(656.27,Digits)));
		nud_ini = nud = ( nC==nF ? 0 : (nd-1)/(nF-nC) );
		dn_ini  = dn  = atof(word(name,2).c_str())/100000;
		set_index();
	}
	else if(nd!=nd_ini){
		nd_ini   = nd;
		name_ini = name = values_to_name(nd,nud,dn);
		set_index();
	}
	else if(nud!=nud_ini){
		nud_ini  = nud;
		name_ini = name = values_to_name(nd,nud,dn);
		set_index();
	}
	else if(dn!=dn_ini){
		dn_ini   = dn;
		name_ini = name = word(name,1) + ( dn==0 ? "" : " " + str(Round(dn*100000,0)) );
		set_index();
	}

	coordinated=true;
}

// ---- public members --------------------------------------------
cGlass::cGlass(){
	cn=1; alloc();
	initialize();
}

cGlass::cGlass(int cn){
	this->cn=cn; alloc();
	initialize();
}

cGlass::cGlass(const cGlass& x){
	cn=x.cn; alloc();
	assignment(x);
}

cGlass::~cGlass(){
	free();
}

cGlass& cGlass::operator=(const cGlass& x){
	if(cn<x.cn){  // cLens1�̑���� cn>=x.cn �̏ꍇ�������̍Ċm�ۂ��s��Ȃ��D����ɍ��킹��D
		free();
		cn=x.cn; alloc();
	}
	assignment(x);
	return *this;
}

bool operator==(cGlass a,cGlass b){
	return ( a.Name()==b.Name() );
}

void cGlass::set_wl(int j,double val){
	if(1<=j && j<=cn){
		wl[j]=val;
		set_index(j);
	}
}

void cGlass::set_cn_wl(int cn,double *pwl){
	int j;
	
	free();
	this->cn=cn; alloc();
	for(j=1; j<=cn; ++j) wl[j]=pwl[j];
	set_index();
}

bool cGlass::UsePointer=false;
int  cGlass::Digits=5;

std::string& cGlass::Name() { if(!coordinated || UsePointer) coordinate(); coordinated=false; return name; }
double& cGlass::Nd() { if(!coordinated || UsePointer) coordinate(); coordinated=false; return nd;  }
double& cGlass::Nud(){ if(!coordinated || UsePointer) coordinate(); coordinated=false; return nud; }
double& cGlass::dN() { if(!coordinated || UsePointer) coordinate(); coordinated=false; return dn;  }

double cGlass::N(int j) { if(!coordinated || UsePointer) coordinate(); return n[j];     }

int cGlass::IsGlass(){ return (Nd()>1.1 && Nud()>1) ? 1 : 0; }

int cGlass::IsGrin()    { if(!coordinated || UsePointer) coordinate(); return is_grin;  }
double cGlass::rA(int j){ if(!coordinated || UsePointer) coordinate(); return ra[j];    }
double cGlass::GrinPhi(){ if(!coordinated || UsePointer) coordinate(); return grin_phi; }

double cGlass::NNuActual(double Nd,double Nud){
	// Nd-��d�}�ł̃K���X�̑��ݔ͈͂ɑ΂��āC
	// ���ɂ���Ƃ��͐��̒l���C�O�ɂ��鎞�͕��̒l���C���E��ɂ���Ƃ���0��Ԃ��D
	// �߂�l�̐�Βl�͋��E�܂ł̋����ɔ�Ⴗ��
	const int N=17;   // ���E���`����_�̐�
	int i;

	double x[N+1],y[N+1];

	x[ 1]=94.93; y[ 1]=1.43875;  // S-FPL53
	x[ 2]=81.54; y[ 2]=1.49700;  // S-FPL51
	x[ 3]=63.33; y[ 3]=1.61800;  // S-PHM52
	x[ 4]=54.68; y[ 4]=1.72916;  // S-LAL18
	x[ 5]=46.62; y[ 5]=1.81600;  // S-LAH59
	x[ 6]=40.76; y[ 6]=1.88300;  // S-LAH58
	x[ 7]=28.27; y[ 7]=2.00330;  // S-LAH79
	x[ 8]=18.90; y[ 8]=1.92286;  // S-NPH2
	x[ 9]=22.76; y[ 9]=1.80809;  // S-NPH1
	x[10]=26.52; y[10]=1.76182;  // S-TIH14
	x[11]=30.13; y[11]=1.69895;  // S-TIM35
	x[12]=34.46; y[12]=1.63980;  // S-TIM27
	x[13]=40.75; y[13]=1.58144;  // S-TIL25
	x[14]=45.79; y[14]=1.54814;  // S-TIL1
	x[15]=52.43; y[15]=1.51742;  // S-NSL36
	x[16]=64.14; y[16]=1.51633;  // S-BSL7
	x[17]=70.23; y[17]=1.48749;  // S-FSL5

	for(i=1; i<=N; ++i){
		y[i]*=100;  // �P�ʂ�Nd-Nud�`���[�g�̈�}�X��1�Ƃ���i��}�X��Nd��0.01, Nud��1�j
	}

	return PolygonPointDistance(N,x,y,Nud,Nd*100);
}

double cGlass::NNuActual(){
	return NNuActual(Nd(),Nud());
}

