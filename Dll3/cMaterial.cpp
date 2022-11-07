#include "stdafx.h"
#include "MyDllOptics.h"
#include "cMaterial.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

// �����ɋ��ܗ��f�[�^�t�@�C�����Ƃ��̏ꏊ���w�肷��D
//     �v���O������ GlassDataFile1,GlassDataFile2, �̏��ԂŃt�@�C����T���D	
std::string cMaterial::GlassDataFile1=(std::string)getenv("USERPROFILE")+"\\Documents\\GlassData\\glass.dat";
std::string cMaterial::GlassDataFile2="glass.dat";  // �J�����g�t�H���_�� glass.dat

// �����ɐF�K���X�t�B���^�̃f�[�^������t�H���_���w�肷��D
std::string cMaterial::ColorGlassDataFolder=(std::string)getenv("USERPROFILE")+"\\Documents\\GlassData\\ColorFilter";


int cMaterial::Open(std::string name0){
	int i;
	const int SIZE=100;
	static cMaterial buf[SIZE+1];
	static cMaterial buf2[SIZE+1];
	static int ip=1;
	static int ip2=1;
	cMaterial tmp;
	
	for(i=1; i<=SIZE; i++) {  //  �� SIZE�� ip-1,ip2-1 �ɒu�������Ă��������̌��ʂ͂Ȃ������������D 2017.02.10
		// �܂��C�o�b�t�@ buf1,buf2 �Ɉ�v������̂��Ȃ����T���D
		if(buf[i].name0==name0){
			tmp=buf[i];
			tmp.name=name; tmp.sgn=sgn; tmp.dN=dN; *this=tmp; // name,sgn,dN���ς�Ȃ��悤�ɂ���
			return 1;
		}
		if(buf2[i].name0==name0){
			return 0;
		}
	}

	// �o�b�t�@ buf1,buf2 �ɍ��v������̂��Ȃ������ꍇ
	std::ifstream from(GlassDataFile1.c_str());
	
	if(!from){
		from.clear();   // ���ꂪ�Ȃ��ƍēx��open���ł��Ȃ��l�q
		from.open(GlassDataFile2.c_str());
	}

	if(from){
		do{
			from>>tmp;
		} while(tmp.name0!=name0 && !from.eof());

		if(tmp.name0==name0 && ip<=SIZE){
			// �f�[�^�t�@�C���ɍ��v������̂��������ꍇ�� buf �ɓo�^����D
			buf[ip]=tmp;
			ip++;
			tmp.name=name; tmp.sgn=sgn; tmp.dN=dN; *this=tmp; // name,sgn,dN���ς�Ȃ��悤�ɂ���
			return 1;
		}
		if(tmp.name0!=name0 && ip2<=SIZE){
			// �f�[�^�t�@�C���ɍ��v������̂��Ȃ��ꍇ�� buf2 �ɓo�^����D
			buf2[ip2].name0=name0;
			ip2++;
			return 0;
		}
	}
	
	return 0;
}

// public members /////////////////////////////////////////////////////////////

double cMaterial::HerzEq(double Nd, double nud, double wl_nm,int digits/*=5*/){
	const double a=0.035;
	double A,B, N, x;
	x=wl_nm/1000;
	A=-1.294878 +0.088927*x*x +0.37349/(x*x-a) +0.005799/(x*x-a)/(x*x-a);
	B=0.00125 -0.007058*x*x +0.001071/(x*x-a) -0.000218/(x*x-a)/(x*x-a);
	N=1+(Nd-1)*(1+B+A/nud);
	N=Round(N,-digits);
	return N;
}

double cMaterial::HerzEqDerivative(double Nd, double nud, double wl_nm){
	// �g�������W��
	const double a=0.035;
	double A1,B1, N1, x;
	x=wl_nm/1000;
	A1=-1.294878+0.088927*(2*x) +0.37349*(-2*x/(x*x-a)/(x*x-a))
 	   +0.005799*(-4*x/(x*x-a)/(x*x-a)/(x*x-a));
	B1=0.00125-0.007058*(2*x) +0.001071*(-2*x/(x*x-a)/(x*x-a))
  	   -0.000218*(-4*x/(x*x-a)/(x*x-a)/(x*x-a));
	N1=(Nd-1)*(B1+A1/nud);
	N1/=1000;
	return N1;
}

double cMaterial::HerzEqDerivative2(double Nd, double nud, double wl_nm){
	// �g��2�K�����W��
	const double a=0.035;
	double A2,B2, N2, x;
	x=wl_nm/1000;
	A2=-1.294878+0.088927*2 +0.37349*((6*x*x+2*a)/(x*x-a)/(x*x-a)/(x*x-a))
 	   +0.005799*((20*x*x+4*a)/(x*x-a)/(x*x-a)/(x*x-a)/(x*x-a));
	B2=0.00125-0.007058*(2*x) +0.001071*((6*x*x+2*a)/(x*x-a)/(x*x-a)/(x*x-a))
  	   -0.000218*((20*x*x+4*a)/(x*x-a)/(x*x-a)/(x*x-a)/(x*x-a));
	N2=(Nd-1)*(B2+A2/nud);
	N2/=1000000;
	return N2;
}

double cMaterial::HerzEqNdDerivative(double Nd, double nud, double wl_nm){
	// Nd�����W��
	const double a=0.035;
	double A,B, x;
	x=wl_nm/1000;
	A=-1.294878 +0.088927*x*x +0.37349/(x*x-a) +0.005799/(x*x-a)/(x*x-a);
	B=0.00125 -0.007058*x*x +0.001071/(x*x-a) -0.000218/(x*x-a)/(x*x-a);
	return 	1+B+A/nud;
}

double cMaterial::HerzEqNUdDerivative(double Nd, double nud, double wl_nm){
	// nud�����W��
	const double a=0.035;
	double A, x;
	x=wl_nm/1000;
	A=-1.294878 +0.088927*x*x +0.37349/(x*x-a) +0.005799/(x*x-a)/(x*x-a);
	return -(Nd-1)*A/(nud*nud);
}

double cMaterial::Nd(int glasscode){
	double nd;
	nd=((double)(glasscode-glasscode%1000))/1000000+1;
	if(nd<1.3) nd+=1;   // Nd=1.3���Ⴂ�t��,�̂͂����炭�Ȃ��D
	                    // �C�̂��R�[�h�ŕ\�����Ƃ������炭�Ȃ��D
	                    // �R�[�h�̍ŏ���3����299�ȉ��̂Ƃ���
	                    // ���ܗ�2�ȏ�Ɋ��蓖�Ă�D
	return nd;
}

double cMaterial::Nud(int glasscode){
	return ((double)(glasscode%1000))/10;
}

std::string cMaterial::MaterialName(std::string gname){
	// �ޗ�����Ԃ�
	// ��F "-BK7" -> "BK7"
	//      "BK7 50" -> "BK7"
	if(gname=="") return "";
	if(gname[0]=='-'){
		gname=gname.erase(0,1);
	}
	if(gname=="1") return "";
	if(words(gname)>1){
		gname=word(gname,1);
	}
	return gname;
}

double cMaterial::HerzEq(int glasscode,double wl_nm,int digits/*=5*/){
	return HerzEq(Nd(glasscode),Nud(glasscode),wl_nm,digits);
}

double cMaterial::HerzEqDerivative(int glasscode,double wl_nm){
	return HerzEqDerivative(Nd(glasscode),Nud(glasscode),wl_nm);
}

double cMaterial::HerzEqDerivative2(int glasscode,double wl_nm){
	return HerzEqDerivative2(Nd(glasscode),Nud(glasscode),wl_nm);
}

complex cMaterial::Index(std::string name,double wl_nm,int digits/*=5*/){
	cMaterial g(name);
	return g.Index(wl_nm,digits);
}

double cMaterial::GroupIndex(std::string name,double wl_nm){
	// �g��wl_nm�ł̌Q���ܗ� ng = n^2/(n+��o(dn/d��o)),  ��o=�^�󒆋��ܗ�
	double n,n1;

	cMaterial g(name);
	n =Re(g.Index(wl_nm));
	n1=Re(g.IndexDerivative(wl_nm));
	return n*n/(n+wl_nm*n1);
}

complex cMaterial::IndexDerivative2(std::string name,double wl_nm){
	cMaterial g(name);
	return g.IndexDerivative2(wl_nm);
}

double cMaterial::ColorGlassInnerTransmittance(std::string filename,double wl,double thickness){
	// �t�@�C���̏���
	// 
	//   2.5              <- 1�s�� ����t0
	//   200  0           <- 2�s�ڈȉ��C�g��(nm) �� ����t0�ɂ����铧�ߗ� ���
	//   290  0
	//   300  .0038
    //   310  .087
	//   ..   ..

	int i;
	cXYData x;
	cXYList li;
	double t0;

	std::ifstream from( (ColorGlassDataFolder+'/'+filename).c_str() );
	
	if(from){
		from>>t0;
		do{
			from>>x;
			li.AddTail(x);
		} while( !from.eof() );
		for(i=1; i<=li.GetSize(); ++i) {
			li.GetData(x,i);
			x.y=pow( x.y, thickness/t0 );
			li.SetData(x,i);
		}
		return li.y(wl);
	}
	else return 1;
}


cMaterial::cMaterial(){
	type=NUL;
}

cMaterial::cMaterial(std::string name){	
	this->SetName(name);
}

cMaterial::~cMaterial(){}

void cMaterial::SetName(std::string name){
	// ��F
	//    "-BK7"   -> Nd=-1.51680
	//    "BK7 50" -> Nd=1.51680+0.00050=1.51730

	if(name=="") return;

	this->name=name;

	if(name[0]=='-'){
		name=name.erase(0,1);
		sgn=-1;
	}
	else{
		sgn=1;
	}

	if(words(name)>1){
		name0=word(name,1);
		dN=atof(word(name,2).c_str())/100000;
	}
	else{
		name0=name;
		dN=0;
	}

	type=NUL;

	if( name0=="1" ){
		type=AIR;
	}
	else if( name0.size()==6 && atoi(name0.c_str())>100 ){
		// >100 �̑���ɗႦ�� >100000 �͂悭�Ȃ��DS-LAH79�̃R�[�h�� 003283�ł���D
		type=CODE;
	}
	else if( is_numeric(arg(to_space(name0,':'),0)) && is_numeric(arg(to_space(name0,':'),1)) ){
		// "1.51680:64.17" -> N=1.51680,��=64.17
		// <����>
		//  �{�u���b�N�� else if(Open(name0)!=0) �̃u���b�N����Ɏ��s���邱�ƁD
		//  �������Ȃ��ƁCOpen()��buf2�ɃS�~�����܂邽�߂��C
		//  ("1.51633:64.14" �����o�b�t�@�ɓo�^����Ă��܂�)
		//  ���ܗ���ϐ��Ƃ��鎩���݌v�����ɒx���Ȃ�D
		type=N_NU;
	}
	else if( Open(name0)!=0 ){
		type=FILE;
	}
	else if( is_numeric(arg(to_space(name0,'_'),0)) && is_numeric(arg(to_space(name0,'_'),1)) ){
		// "0.07_-4.2" -> N=complex(0.07,-4.2)
		type=COMPLEX;
	}
	else if( is_numeric(arg(name0,0)) ){
		// "1.5" -> N=complex(1.5,0)
		type=REAL;
	}
}

std::istream& operator>>(std::istream& from,cMaterial& x){
	from>>x.name; 
	if(from.eof()) return from; // ���ꂪ�Ȃ��Ɩ������[�v�Ɋׂ�\��������D
	                            // �y�����z
	                            //  from�̍Ō�̃f�[�^�̌�ɋ󔒂���s������Ƃ��C
	                            //  from.eof()�̓f�[�^��ǂ݂��������_�ł�false�̂܂܂ł���C
	                            //  (�󔒂���s���Ȃ����true�ɂȂ�悤�ł���D)
	                            //  �f�[�^�̓ǂݏo���Ɏ��s�������_�ŏ��߂�true�ɂȂ�D
	                            //  ���������āCfrom��eq_type==6�̃f�[�^�ŏI����Ă��āC
	                            //  x.name��from�ɂȂ��Ƃ�("bk7"�Ȃ�)�C
	                            //       from�̃f�[�^��ǂ݂����Ă��{�֐���������x�Ă΂��
	                            //    -> from>>x.eq_type�Ɏ��s���Cx.eq_type�ɑO�̃f�[�^6���c��D
	                            //    -> case 6��do���[�v��wl==0�ƂȂ�Ȃ����ߖ������[�v�ƂȂ�D

	x.name0=x.name; x.sgn=1; x.dN=0;  // �t�@�C���� "-BK7" �̂悤�Ƀ}�C�i�X�̂����f�[�^��C
	                                  // "BK7 50" �̂悤�Ɍ덷�����܂ރf�[�^�͂Ȃ��Ƃ���D
	x.type=cMaterial::FILE;
	from>>x.eq_type;
	
	switch(x.eq_type){
	case 0:
		// �󔒂̂Ȃ��������0��������΁C�R�����g�Ƃ��Ďg�p�ł���D
		break;
	case 1:
		from>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6;
		break;
	case 2:
		from>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6;
		from>>x.b1>>x.b2>>x.b3>>x.b4>>x.b5>>x.b6;
		break;
	case 3:
		from>>x.a1>>x.a2>>x.a3;
		from>>x.b1>>x.b2>>x.b3;
		break;
	case 4:
		from>>x.a1>>x.b1;
		break;
	case 5:
		from>>x.glass_code;
		break;
	case 6:
		// ��
		// Ag 6
		// 147.6 1.032 -0.61    (�g�� ���� ����)
		// 151.2 0.993 -0.65
		//  :
		// 9537.0 12.210 -52.20
		// 9919.0 13.110 -53.70
		// 0                    (�g��0�͏I�[��\��)
		x.n.RemoveAll();
		x.k.RemoveAll();
		do{
			double wl,n,k;
			from>>wl;
			if(wl==0) break;
			from>>n>>k;
			x.n.AddData(wl,n);
			x.k.AddData(wl,k);
		}while(true);
		break;
	case 7:
		from>>x.a1>>x.a2>>x.a3;
		break;
	case 8:
		from>>x.nd>>x.nud;
		break;
	case 9:
		// ���K���X�Ђ̕��U��
		from>>x.a0>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6>>x.a7>>x.a8;
		break;

	case 11:
		// ��
		// selfoc300 11
		// 2                    (��)
		// 632.8 1.608 0.304    (�g�� no ��A)
		// 830   1.600 0.298
		// 1300  1.592 0.295
		// 1560  1.592 0.294
		// 0                    (�g��0�͏I�[��\��)
		from>>x.GrinPhi;
		x.n.RemoveAll();
		x.lrA.RemoveAll();
		do{
			double wl,no,rA;
			from>>wl;
			if(wl==0) break;
			from>>no>>rA;
			x.n.AddData(wl,no);
			x.lrA.AddData(wl,rA);
		}while(true);
		break;
	}
	return from;
}

complex cMaterial::Index(double wl_nm,int digits/*=5*/){
	// digits : �����_�ȉ�digits���Ŋۂ߂�D�قƂ�ǂ̏ꍇ��5���ł悢�D
	//          �������C�Ⴆ��SD-OCT�̕����v�̃Z���T��845�`913nm��1024��������D
	//          ���̂Ƃ������_�ȉ�5���ł͋��ܗ���1024�i�K�Ŋ��炩�ɕ\���Ȃ��D
	complex N;
	double x, x2,x4,x6,x8,x10,x12;

	switch(type){
	case AIR:
		N=1;
		break;
	case CODE:
		N=HerzEq(atoi(name0.c_str()),wl_nm,digits);
		break;
	case FILE:
		switch(eq_type){
		case 1:
			// ���[�����W�J(�������Ȃ�)
			x=wl_nm/1000;
			x2=x*x;
			x4=x2*x*x;
			x6=x4*x*x;
			x8=x6*x*x;
			Re(N)=sqrt( a1+a2*x2+a3/x2+a4/x4+a5/x6+a6/x8 );
			Im(N)=0;
			break;
		case 2:
			// ���[�����W�J(����������)
			x=wl_nm/1000;
			x2=x*x;
			x4=x2*x*x;
			x6=x4*x*x;
			x8=x6*x*x;
			Re(N)= sqrt( a1+a2*x2+a3/x2+a4/x4+a5/x6+a6/x8 );
			Im(N)=-sqrt( b1+b2*x2+b3/x2+b4/x4+b5/x6+b6/x8 );
			break;
		case 3:
			// �Z���}�C���[�̎�
			x=wl_nm/1000;
			Re(N)= sqrt( a1*x*x/(x*x-b1) +a2*x*x/(x*x-b2) +a3*x*x/(x*x-b3) +1 );
			Im(N)= 0;
			break;
		case 4:
			// ���ڎw��(���U�Ȃ�)
			Re(N)=a1;
			Im(N)=b1;
			break;
		case 5:
			// �K���X�R�[�h
			Re(N)=HerzEq(glass_code,wl_nm,digits);
			Im(N)=0;
			break;
		case 6:
			// �g���C���ܗ��̕\�̒������(�ߎ����ł��܂��\���ł��Ȃ������Ȃǂɗp����)
			Re(N)=n.y(wl_nm);
			Im(N)=k.y(wl_nm);
			break;
		case 7:
			// �R�[�V�[�̕��U��
			x=wl_nm/1000;
			Re(N)=a1+a2/x/x+a3/x/x/x/x;  //�y���Ӂza3 �̍��̎����� -4��i-3��ł͂Ȃ��j
			Im(N)=0;
			break;
		case 8:
			// �w���c�x���K�\���U��
			//�i�K���X�R�[�h�ł͂Ȃ��CNd�ƃ�d��^����i�K���X�R�[�h�͗L������3�������Ȃ��j�j�j
			Re(N)=HerzEq(nd,nud,wl_nm,digits);
			Im(N)=0;
			break;
		case 9:
			// ���K���X�Ђ̕��U��
			x=wl_nm/1000;
			x2=x*x;
			x4=x2*x*x;
			x6=x4*x*x;
			x8=x6*x*x;
			x10=x8*x*x;
			x12=x10*x*x;
			Re(N)=sqrt( a0+a1*x2+a2*x4+a3/x2+a4/x4+a5/x6+a6/x8+a7/x10+a8/x12 );
			Im(N)=0;
			break;

		case 11:
			// Selfoc�̒��S���ܗ�
			Re(N)=n.y(wl_nm);
			Im(N)=0;
			break;
		}
		
		Re(N)=Round( Re(N),-digits ); Im(N)=Round( Im(N),-digits );
		break;
	case N_NU:
		N=HerzEq(atof(arg(to_space(name0,':'),0).c_str()),
		         atof(arg(to_space(name0,':'),1).c_str()), wl_nm,digits);
		break;
	case COMPLEX:
		N=complex( Round(atof(arg(to_space(name0,'_'),0).c_str()), -digits),
		           atof(arg(to_space(name0,'_'),1).c_str()) );
		break;
	case REAL:
		N=Round( atof(name0.c_str()), -digits );
		break;

	default:
		return N=0;
	}

	return sgn*(N+dN);
}

complex cMaterial::IndexDerivative(double wl_nm){
	// ���ܗ��̔g�������W��
	// �P�ʂ� 1/nm
	complex N,N1;
	double x, x3,x5,x7,x9,x11,x13;

	switch(type){
	case AIR:
		N1=0;
		break;
	case CODE:
		N1=HerzEqDerivative(atoi(name0.c_str()),wl_nm);
		break;
	case FILE:
		switch(eq_type){
		case 1:
			// ���[�����W�J(�������Ȃ�)
			N=Index(wl_nm)*sgn;
			x=wl_nm/1000;
			x3=x*x*x;
			x5=x3*x*x;
			x7=x5*x*x;
			x9=x7*x*x;
			Re(N1)=(1.0/2/Re(N))*( 2*a2*x-2*a3/x3-4*a4/x5-6*a5/x7-8*a6/x9 );
			Re(N1)/=1000;  // �P�ʂ� 1/��m ���� 1/nm ��
			Im(N1)=0;
			break;
		case 2:
			// ���[�����W�J(����������)
			N=Index(wl_nm)*sgn;
			x=wl_nm/1000;
			x3=x*x*x;
			x5=x3*x*x;
			x7=x5*x*x;
			x9=x7*x*x;
			Re(N1)=(1.0/2/Re(N))*( 2*a2*x-2*a3/x3-4*a4/x5-6*a5/x7-8*a6/x9 );
			Re(N1)/=1000;
			Im(N1)=-(1.0/2/Im(N))*( 2*a2*x-2*a3/x3-4*a4/x5-6*a5/x7-8*a6/x9 );
			Im(N1)/=1000;
			break;
		case 3:
			// �Z���}�C���[�̎�
			N=Index(wl_nm)*sgn;
			x=wl_nm/1000;
			Re(N1)=(1.0/Re(N))
			      *(-a1*b1*x/(x*x-b1)/(x*x-b1)-a2*b2*x/(x*x-b2)/(x*x-b2)-a3*b3*x/(x*x-b3)/(x*x-b3));
			Re(N1)/=1000;
			Im(N1)=0;
			break;
		case 4:
			// ���ڎw��(���U�Ȃ�)
			Re(N1)=0;
			Im(N1)=0;
			break;
		case 5:
			// �K���X�R�[�h
			Re(N1)=HerzEqDerivative(glass_code,wl_nm);
			Im(N1)=0;
			break;
		case 6:
			// �g���C���ܗ��̕\�̒������(�ߎ����ł��܂��\���ł��Ȃ������Ȃǂɗp����)
			Re(N1)=n.dydx(wl_nm);
			Im(N1)=k.dydx(wl_nm);
			break;
		case 7:
			// �R�[�V�[�̕��U��
			x=wl_nm/1000;
			Re(N1)=-2*a2/x/x/x-4*a3/x/x/x/x/x;
			Re(N1)/=1000;
			Im(N1)=0;
			break;
		case 8:
			// �w���c�x���K�\���U��
			//�i�K���X�R�[�h�ł͂Ȃ��CNd�ƃ�d��^����i�K���X�R�[�h�͗L������3�������Ȃ��j�j�j
			Re(N1)=HerzEqDerivative(nd,nud,wl_nm);
			Im(N1)=0;
			break;
		case 9:
			// ���K���X�Ђ̕��U��
			N=Index(wl_nm)*sgn;
			x=wl_nm/1000;
			x3=x*x*x;
			x5=x3*x*x;
			x7=x5*x*x;
			x9=x7*x*x;
			x11=x9*x*x;
			x13=x11*x*x;
			Re(N1)=(1/Re(N))*( a1*x+2*a2*x3-a3/x3-2*a4/x5-3*a5/x7-4*a6/x9-5*a7/x11-6*a8/x13 );
			Re(N1)/=1000;  // �P�ʂ� 1/��m ���� 1/nm ��
			Im(N1)=0;
			break;

		case 11:
			// Selfoc�̒��S���ܗ�
			Re(N1)=n.dydx(wl_nm);
			Im(N1)=0;
			break;
		}
		
		break;
	case N_NU:
		N1=HerzEqDerivative(atof(arg(to_space(name0,':'),0).c_str()),
		                    atof(arg(to_space(name0,':'),1).c_str()), wl_nm);
		break;
	case COMPLEX:
		N1=0;
		break;
	case REAL:
		N1=0;
		break;

	default:
		return N1=0;
	}

	return N1*sgn;
}

complex cMaterial::IndexDerivative2(double wl_nm){
	// ���ܗ��̔g��2�K�����W�� (GDD(�Q���x�x�����U)�v�Z�̂��ߍ쐬 2011.10.28)
	// �P�ʂ� 1/nm^2
	complex N,N1,N2;
	double x, x2,x4,x6,x8,x10,x12,x14;

	switch(type){
	case AIR:
		N2=0;
		break;
	case CODE:
		N2=HerzEqDerivative2(atoi(name0.c_str()),wl_nm);
		break;
	case FILE:
		switch(eq_type){
		case 1:
			// ���[�����W�J(�������Ȃ�)
			N=Index(wl_nm)*sgn;
			N1=IndexDerivative(wl_nm)*sgn;
			x=wl_nm/1000;
			x4=x*x*x*x;
			x6=x4*x*x;
			x8=x6*x*x;
			x10=x8*x*x;
			Re(N2)=(1/Re(N))*( (a2+3*a3/x4+10*a4/x6+21*a5/x8+36*a6/x10)-Re(N1)*Re(N1) );
			Re(N2)/=1000000;  // �P�ʂ� 1/��m^2 ���� 1/nm^2 ��
			Im(N2)=0;
			break;
		case 2:
			// ���[�����W�J(����������)
			N=Index(wl_nm)*sgn;
			x=wl_nm/1000;
			x4=x*x*x*x;
			x6=x4*x*x;
			x8=x6*x*x;
			x10=x8*x*x;
			Re(N2)=(1/Re(N))*( (a2+3*a3/x4+10*a4/x6+21*a5/x8+36*a6/x10)-Re(N1)*Re(N1) );
			Re(N2)/=1000000;
			Im(N2)=(1/Im(N))*( (a2+3*a3/x4+10*a4/x6+21*a5/x8+36*a6/x10)-Im(N1)*Im(N1) );
			Im(N2)/=1000000;
			break;
		case 3:
			// �Z���}�C���[�̎�
			N=Index(wl_nm)*sgn;
			N1=IndexDerivative(wl_nm)*sgn;
			x=wl_nm/1000;
			Re(N2)=(1/Re(N))
			      *( (3*a1*b1*x*x-a1*b1*b1)/(x*x-b1)/(x*x-b1)/(x*x-b1)
			        +(3*a2*b2*x*x-a2*b2*b2)/(x*x-b2)/(x*x-b2)/(x*x-b2)
			        +(3*a3*b3*x*x-a3*b3*b3)/(x*x-b3)/(x*x-b3)/(x*x-b3) -Re(N1)*Re(N1) );
			Re(N2)/=1000000;
			Im(N2)=0;
			break;
		case 4:
			// ���ڎw��(���U�Ȃ�)
			Re(N2)=0;
			Im(N2)=0;
			break;
		case 5:
			// �K���X�R�[�h
			Re(N2)=HerzEqDerivative2(glass_code,wl_nm);
			Im(N2)=0;
			break;
		case 6:
			// �g���C���ܗ��̕\�̒������(�ߎ����ł��܂��\���ł��Ȃ������Ȃǂɗp����)
			Re(N2)=0;
			Im(N2)=0;
			break;
		case 7:
			// �R�[�V�[�̕��U��
			x=wl_nm/1000;
			Re(N2)=6*a2/x/x/x/x+20*a3/x/x/x/x/x/x;
			Re(N2)/=1000000;
			Im(N2)=0;
			break;
		case 8:
			// �w���c�x���K�\���U��
			//�i�K���X�R�[�h�ł͂Ȃ��CNd�ƃ�d��^����i�K���X�R�[�h�͗L������3�������Ȃ��j�j�j
			Re(N2)=HerzEqDerivative2(nd,nud,wl_nm);
			Im(N2)=0;
			break;
		case 9:
			// ���K���X�Ђ̕��U��
			N=Index(wl_nm)*sgn;
			N1=IndexDerivative(wl_nm)*sgn;
			x=wl_nm/1000;
			x2=x*x;
			x4=x2*x*x;
			x6=x4*x*x;
			x8=x6*x*x;
			x10=x8*x*x;
			x12=x10*x*x;
			x14=x12*x*x;
			Re(N2)=(1/Re(N))*( (a1+6*a2*x2+3*a3/x4+10*a4/x6+21*a5/x8+36*a6/x10+55*a7/x12+78*a8/x14)-Re(N1)*Re(N1) );
			Re(N2)/=1000000;  // �P�ʂ� 1/��m^2 ���� 1/nm^2 ��
			Im(N2)=0;
			break;

		case 11:
			// Selfoc�̒��S���ܗ�
			Re(N2)=0;
			Im(N2)=0;
			break;
		}
		
		break;
	case N_NU:
		N2=HerzEqDerivative2(atof(arg(to_space(name0,':'),0).c_str()),
		                     atof(arg(to_space(name0,':'),1).c_str()), wl_nm);
		break;
	case COMPLEX:
		N2=0;
		break;
	case REAL:
		N2=0;
		break;

	default:
		return N2=0;
	}

	return N2*sgn;
}

int cMaterial::IsGrin(){
	return (type==FILE && eq_type==11) ? 1 : 0;
}

double cMaterial::rA(double wl_nm){
	return lrA.y(wl_nm);
}

