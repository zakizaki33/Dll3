// #include "stdafx.h"
#include "cLens1.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
// TODO:�Ӗ��͕�����Ȃ�����U�R�����g�A�E�g���Ė������ʂ����F2022-10-06
// TODO:_AFXDLL����������`���Ă���F2022-10-06
// #define new DEBUG_NEW
#endif

const double LN=1e7;

///// protected members //////////////////////////////////////////////////////////////

// ---- �I�[�o�[�����h�~�̂��߂̈����`�F�b�N�̗�D�����炭�x���Ȃ�̂Ŏ����͕ۗ��D ----
// double err;
// double& cLens1::r(int i) { if(0<=i && i<=k+1) return surf[i].r(); else return err; }
// ------------------------------------------------------------------------------------

double& cLens1::r(int i) { return surf[i].r(); }
double& cLens1::rObj()   { return r(0); }
double& cLens1::rImage() { return r(k+1); }
double& cLens1::c(int i) { return surf[i].c(); }
double& cLens1::Newton(int i) { return surf[i].Newton; }
double& cLens1::As0(int i) { return surf[i].As0; }
double& cLens1::As45(int i) { return surf[i].As45; }
double& cLens1::NewtonTol(int i) { return surf[i].Newton; }
double& cLens1::AsTol(int i) { return surf[i].AsTol; }
int&    cLens1::rVariable(int i) { return surf[i].rVariable; }
int&    cLens1::asph_type(int i) { return surf[i].asph_type; }
double& cLens1::kp(int i) { return surf[i].kp; }
double& cLens1::NormH(int i) {return surf[i].NormH(); }
double& cLens1::a1(int i) { return surf[i].a1; }
double& cLens1::a2(int i) { return surf[i].a2; }
double& cLens1::a3(int i) { return surf[i].a3; }
double& cLens1::a4(int i) { return surf[i].a4; }
double& cLens1::a5(int i) { return surf[i].a5; }
double& cLens1::a6(int i) { return surf[i].a6; }
double& cLens1::a7(int i) { return surf[i].a7; }
double& cLens1::a8(int i) { return surf[i].a8; }
double& cLens1::a9(int i) { return surf[i].a9; }
double& cLens1::a10(int i) { return surf[i].a10; }
double& cLens1::a11(int i) { return surf[i].a11; }
double& cLens1::a12(int i) { return surf[i].a12; }
double& cLens1::a13(int i) { return surf[i].a13; }
double& cLens1::a14(int i) { return surf[i].a14; }
double& cLens1::a15(int i) { return surf[i].a15; }
double& cLens1::a16(int i) { return surf[i].a16; }
double& cLens1::a18(int i) { return surf[i].a18; }
double& cLens1::a20(int i) { return surf[i].a20; }
double& cLens1::b(int i,int m,int n) { return surf[i].b.b(m,n); }
cZernike& cLens1::zernike(int i){ return surf[i].zernike; }
double& cLens1::ZC(int i,int j){ return surf[i].zernike.C(j); }
cDcon&  cLens1::Dcon(int i){ return surf[i].Dcon; }
double& cLens1::DconRn(int i){ return surf[i].Dcon.Rn; }
double& cLens1::DconA4(int i){ return surf[i].Dcon.A4(); }
double& cLens1::DconA6(int i){ return surf[i].Dcon.A6(); }
double& cLens1::DconA8(int i){ return surf[i].Dcon.A8(); }
double& cLens1::DconA10(int i){ return surf[i].Dcon.A10(); }
double& cLens1::DconA12(int i){ return surf[i].Dcon.A12(); }
double& cLens1::DconA14(int i){ return surf[i].Dcon.A14(); }
double& cLens1::DconA16(int i){ return surf[i].Dcon.A16(); }
double& cLens1::DconA18(int i){ return surf[i].Dcon.A18(); }
double& cLens1::DconA20(int i){ return surf[i].Dcon.A20(); }
cModLegendre& cLens1::legendre(int i){ return surf[i].legendre; }
double& cLens1::LeC(int i,int m,int n){ return surf[i].legendre.C(m,n); }
cSpline& cLens1::spline(int i){ return surf[i].spline; }
double& cLens1::SplineH(int i,int j){ return surf[i].spline.X(j); }
double& cLens1::SplineZ(int i,int j){ return surf[i].spline.Y(j); }
double& cLens1::Apv(int i){ return surf[i].Apv; }
double& cLens1::sfx(int i){ return surf[i].sfx; }
double& cLens1::sfy(int i){ return surf[i].sfy; }
int&    cLens1::UserDefSurf(int i){ return surf[i].UserDefSurf; }
int&    cLens1::cylinder(int i) { return surf[i].cylinder; }
int&    cLens1::fresnel(int i) { return surf[i].fresnel; }
double& cLens1::rbase(int i) { return surf[i].rbase; }
double& cLens1::ry(int i) { return surf[i].ry; }
double& cLens1::rx(int i) { return surf[i].rx; }
double& cLens1::kpy(int i) { return surf[i].kpy; }
double& cLens1::kpx(int i) { return surf[i].kpx; }
int&    cLens1::IsXToroid(int i) { return surf[i].IsXToroid; }
double& cLens1::coneangle(int i) { return surf[i].coneangle; }
double& cLens1::fideal(int i) { return surf[i].fideal(); }
double& cLens1::pideal(int i) { return surf[i].pideal(); }
double& cLens1::aCOA(int i) { return surf[i].aCOA(); }
double& cLens1::caCOA(int i) { return surf[i].caCOA(); }
double& cLens1::bCOA(int i) { return surf[i].bCOA(); }
double& cLens1::cbCOA(int i) { return surf[i].cbCOA(); }
double& cLens1::tCOA(int i) { return surf[i].tCOA; }
double& cLens1::SA0(int i) { return surf[i].SA0; }
double& cLens1::CM0(int i) { return surf[i].CM0; }
int&    cLens1::grating(int i) { return surf[i].grating; }
int&    cLens1::difforder(int i) { return surf[i].difforder; }
double& cLens1::gpitch(int i) { return surf[i].gpitch; }
double& cLens1::grx(int i) { return surf[i].grx; }
double& cLens1::gry(int i) { return surf[i].gry; }
double& cLens1::grz(int i) { return surf[i].grz; }
double& cLens1::Diffusion(int i) { return surf[i].Diffusion; }
int&    cLens1::Fish(int i) { return surf[i].Fish; }
int&    cLens1::EAtype(int i) { return surf[i].EAtype; }
double& cLens1::EAy(int i) { return surf[i].EAy; }
double& cLens1::EAx(int i) { return surf[i].EAx; }
double& cLens1::CHMy(int i) { return surf[i].CHMy; }
double& cLens1::CHMx(int i) { return surf[i].CHMx; }
double& cLens1::EAdy(int i) { return surf[i].EAdy; }
double& cLens1::EAdx(int i) { return surf[i].EAdx; }
int&    cLens1::decenter_type(int i) { return surf[i].decenter_type; }
double& cLens1::dx(int i) { return surf[i].dx; }
double& cLens1::dy(int i) { return surf[i].dy; }
double& cLens1::dz(int i) { return surf[i].dz; }
double& cLens1::rox(int i) { return surf[i].rox; }
double& cLens1::roy(int i) { return surf[i].roy; }
double& cLens1::roz(int i) { return surf[i].roz; }
int&    cLens1::ret(int i) { return surf[i].ret; }
int&    cLens1::order(int i) {return surf[i].order; }
double& cLens1::dx1(int i) { return surf[i].dx1; }
double& cLens1::dy1(int i) { return surf[i].dy1; }
double& cLens1::dz1(int i) { return surf[i].dz1; }
double& cLens1::rox1(int i) { return surf[i].rox1; }
double& cLens1::roy1(int i) { return surf[i].roy1; }
double& cLens1::roz1(int i) { return surf[i].roz1; }
int&    cLens1::order1(int i) {return surf[i].order1; }
std::string& cLens1::CoatName(int i) { return surf[i].CoatName; }
int& cLens1::CoatReverse(int i) { return surf[i].CoatReverse; }
std::string& cLens1::rem(int i) { return surf[i].rem; }

double& cLens1::d(int i) { return med[i].d; }
double& cLens1::delta_d(int i) { return med[i].delta_d; }
int& cLens1::dVariable(int i) { return med[i].dVariable; }
std::string& cLens1::gname(int i) { return med[i].g.Name(); }
int& cLens1::gVariable(int i) { return med[i].gVariable; }
double& cLens1::Nd(int i) { return med[i].g.Nd(); }
double& cLens1::Nud(int i) { return med[i].g.Nud(); }
double& cLens1::dN(int i) { return med[i].g.dN(); }
double cLens1::N(int i,int j) { return med[i].g.N(j); }
int cLens1::IsGlass(int i) { return med[i].g.IsGlass(); }
int cLens1::IsGrin(int i) { return med[i].g.IsGrin(); }
double cLens1::rA(int i,int j) { return med[i].g.rA(j); }
double cLens1::GrinPhi(int i) { return med[i].g.GrinPhi(); }
double cLens1::NNuActual(int i) { return med[i].g.NNuActual(); }

double cLens1::sum(const double *x, int i1, int i2){
	int i; double sum;
	sum=0;
	if( 0<=i1 && i2<=k+1 ){
		sum+=x[i1];  // i1>=i2�̂Ƃ���x[i1]��Ԃ�
		for(i=i1+1; i<=i2; ++i) sum+=x[i];
		return sum;
	}
	else return 0;
}

double cLens1::rms(const double *x, int i1, int i2){
	int i,n; double a;
	n=0;
	a=0;
	if( 0<=i1 && i1<=i2 && i2<=k+1 ){
		for(i=i1; i<=i2; ++i){
			a+=x[i]*x[i];
			n+=1;
		}
		return sqrt(a/n);
	}
	else return 0;
}

double cLens1::exponent(const double *x){
	// x[1],x[2],...x[k], �y�� ��x[i] �𓯈�w���̎w���`���ŕ\���Ƃ��C
	// �������̏����_�ȏオ�ꌅ�ƂȂ�悤�Ȏw�������߂�D
	int i;
	double max,total;
	max=0; 
	for(i=0; i<=k+1; i++){ if( fabs(x[i])>max ) max=fabs(x[i]); }
	total=sum(x,0,k+1);
	if( fabs(total)>max ) max=fabs(total);
	return max==0 ? 0 : floor( log10(max) );
}

void cLens1::create(){
	int j;
	
	mems=k*2;   // �K�v�ʂ�葽�����������m�ۂ��C������Ƀ������̍Ċm��(�A�h���X���ς��)���N���ɂ�������D
	            // (������mems��傫����������Ǝ��s���x���x���Ȃ�̂Œ��ӁD�j
	            // �N���ɂ����Ȃ闝�R�͑�����Z�q�̃R�[�h���Q��(�����̃�����������Ă���� new ���s��Ȃ��j�D
	            // ����ɂ��COptimize�֐��ŕϐ��̃A�h���X x[] ���ς���Ă��܂���ʂ����Ȃ��Ȃ�̂ŗL�p�D
	            // �������C�������Ċm�ۂ��N�����Ă��܂��ꍇ�����蓾��̂ō��{�I�ȉ�����ł͂Ȃ��D             
	            // ������Ƃ��ėႦ�΁C
	            //     cLens *buf;
	            //     buf = this;      // buf�͓����������̈���Q�Ƃ���   
	            //     cLens(buf);      // this�̃��������Ċm��
	            //     ( create() �𔺂����� )
	            //     this = buf;      // ���̃C���X�^���X�ɖ߂�
	            // �Ƃ��āC������ɃA�h���X���܂߂Č��ʂ�ɂ��悤�Ƃ��Ă�
	            //  "this=buf" ���G���[�ƂȂ�ithis�͍��Ӓl�łȂ��j�D
	            // ���{�I�ȑ΍􂪂Ȃ��D 2016.05.13

	color= new std::string[cn+1]; wl=new double[cn+1]; colorweight= new double[cn+1];
	surf=new surface[mems+2];
	med=new medium[mems+2];  for(j=0; j<=mems+1; ++j) med[j]=medium(cn);
	phi=new double[mems+2];
	u=new double[mems+2]; h=new double[mems+2];
	up=new double[mems+2]; hp=new double[mems+2];
	hQ=new double[mems+2]; hQp=new double[mems+2];
	SA=new double[mems+2]; CM=new double[mems+2]; AS=new double[mems+2]; 
	DS=new double[mems+2]; PT=new double[mems+2];
	LC=new double[mems+2]; TC=new double[mems+2];	LC2=new double[mems+2]; TC2=new double[mems+2];
	SA5=new double[mems+2];
	CM41Z=new double[mems+2]; CM41=new double[mems+2]; 
	CM41ALL=new double[mems+2];
	CM23Z=new double[mems+2]; CM23P=new double[mems+2]; CM23=new double[mems+2]; 
	CM23ALL=new double[mems+2];
	SA32F=new double[mems+2]; SA32Z=new double[mems+2]; SA32=new double[mems+2];
	SA32ALL=new double[mems+2];
	AS5=new double[mems+2]; SG5=new double[mems+2];
	DS5=new double[mems+2];
	SAP=new double[mems+2]; CMP=new double[mems+2]; ASP=new double[mems+2]; DSP=new double[mems+2];
	LCP=new double[mems+2];
	PRE=new double[mems+2]; DSE1=new double[mems+2]; DSE2=new double[mems+2]; ASE=new double[mems+2]; 
	PTE=new double[mems+2]; CME=new double[mems+2];
	o=new vector<double>[mems+2];  o0=new vector<double>[mems+2];
	ex=new vector<double>[mems+2]; ex0=new vector<double>[mems+2];
	ey=new vector<double>[mems+2]; ey0=new vector<double>[mems+2];
	ez=new vector<double>[mems+2]; ez0=new vector<double>[mems+2];
	x=new double[mems+2]; y=new double[mems+2]; z=new double[mems+2];
	X=new double[mems+3]; Y=new double[mems+3]; Z=new double[mems+3]; 
	X1=new double[mems+3]; Y1=new double[mems+3]; Z1=new double[mems+3];
	xi=new double[mems+2]; xi1=new double[mems+2];
	ha =new vector<double>[mems+2]; dQa =new vector<double>[mems+2]; 
	hb =new vector<double>[mems+2]; dQb =new vector<double>[mems+2];
	ha1=new vector<double>[mems+2]; dQa1=new vector<double>[mems+2]; 
	hb1=new vector<double>[mems+2]; dQb1=new vector<double>[mems+2];
	E=new vector<complex>[mems+2]; E1=new vector<complex>[mems+2];
	optical_path=new double[mems+2];
}

void cLens1::initialize(){
	IndexDigits=5;
	Note="";
	FRW=550;
	int i,j; 
	for(j=1; j<=cn; j++){
		if     (j==1) color[j]="e";
		else if(j==2) color[j]="C'";
		else if(j==3) color[j]="F'";
		else if(j==4) color[j]="g";
		else if(j>=5) color[j]="e";
		wl[j]=Wl(j);
		colorweight[j]=1;
	}
	Set_color(0,"");
	for(i=0; i<=k+1; i++){
		E[i]=vector<complex>(0,0,0);
	}
	s=LN*10; t=0;
	ApodizationSurf=0;
	GaussianPhiX=GaussianPhiY=0;
	s1fix=0;
	Afocal=0; AfocalMode=MIN; AfocalRotateEye=1;
	StopDominate=0;
	yObjectMax=0; xObjectMax=0; EPD=EPDx=EPy=EPx=0;
	NormalizeUnit="f=1"; NormalizeType=1;
	Stack.SetMaxSize(100);
	stop=0;
	nSpot=20;
	SourcePhiY=SourcePhiX=SourceAxis=0; SourceAreaType=0;
	GrinDs=0.02;
	GrinToListStep=20;
	TcMinGlass=1; TcMinAir=0.5; TeMinGlass=1.2; TeMinAir=0;
	for(i=1; i<=SHAPES_SIZE-1; i++) shapes[i].Clear();
	ref_sphere=0;
	var=var1=var2=var3=var4=var5=0;
	M2=1;
	waist_dia_in=beam_full_div_in=WaistPosIn=0;
	current_shapes=&shapes[1];
	ExcludeVirtualRay=0;
	ExcludeVirtualObject=0;
	FindAThruRayTimeOut=0.1;
	PupilMargin=0.2;
	IgnoreTC=0;
}

void cLens1::erase(){
	delete[] color; delete[] wl; delete[] colorweight;
	delete[] surf;
	delete[] med;
	delete[] phi;
	delete[] u; delete[] h; delete[] up; delete[] hp;
	delete[] hQ; delete[] hQp;
	delete[] SA; delete[] CM; delete[] AS; delete[] DS; delete[] PT;  
	delete[] LC; delete[] TC; delete[] LC2; delete[] TC2;
	delete[] SA5; 
	delete[] CM41Z; delete[] CM41; delete[] CM41ALL;
	delete[] CM23Z; delete[] CM23P; delete[] CM23; delete[] CM23ALL;
	delete[] SA32F; delete[] SA32Z; delete[] SA32; delete[] SA32ALL;
	delete[] AS5; delete[] SG5;
	delete[] DS5;
	delete[] SAP; delete[] CMP; delete[] ASP; delete[] DSP; delete[] LCP;
	delete[] PRE; delete[] DSE1; delete[] DSE2; delete[] ASE;
	delete[] PTE; delete[] CME;
	delete[] o;  delete[] o0;
	delete[] ex; delete[] ex0;
	delete[] ey; delete[] ey0;
	delete[] ez; delete[] ez0;
	delete[] x; delete[] y; delete[] z;
	delete[] X; delete[] Y; delete[] Z;
	delete[] X1; delete[] Y1; delete[] Z1;
	delete[] xi; delete[] xi1;
	delete[] ha;  delete[] dQa;  delete[] hb;  delete[] dQb;
	delete[] ha1; delete[] dQa1; delete[] hb1; delete[] dQb1;
	delete[] E; delete[] E1;
	delete[] optical_path;
	if(ref_sphere!=0) { delete ref_sphere; ref_sphere=0; }
}

void cLens1::assignment(const cLens1& x){
	int i,j;
	filename=x.filename;
	IndexDigits=x.IndexDigits;
	Note=x.Note;
	FRW=x.FRW;
	s=x.s; t=x.t;
	for(j=1; j<=cn && j<=x.cn; j++) {
		color[j]=x.color[j];
		wl[j]=x.wl[j];
		colorweight[j]=x.colorweight[j];
	}
	Set_color(0,"");
	for(i=0; i<=k+1 && i<=x.k+1; i++) {
		surf[i]=x.surf[i];
		med[i]=x.med[i];
	}
	ApodizationSurf=x.ApodizationSurf;
	GaussianPhiX=x.GaussianPhiX; GaussianPhiY=x.GaussianPhiY;
	s1fix=x.s1fix;
	Afocal=x.Afocal; AfocalMode=x.AfocalMode; AfocalRotateEye=x.AfocalRotateEye;
	StopDominate=x.StopDominate;
	yObjectMax=x.yObjectMax; xObjectMax=x.xObjectMax; EPD=x.EPD; EPDx=x.EPDx; EPy=x.EPy; EPx=x.EPx;
	NormalizeUnit=x.NormalizeUnit; NormalizeType=x.NormalizeType;
	stop=x.stop;
	nSpot=x.nSpot;
	SourcePhiY=x.SourcePhiY; SourcePhiX=x.SourcePhiX;
	SourceAxis=x.SourceAxis; SourceAreaType=x.SourceAreaType;
	GrinDs=x.GrinDs;
	GrinToListStep=x.GrinToListStep;
	TcMinGlass=x.TcMinGlass; TcMinAir=x.TcMinAir;
	TeMinGlass=x.TeMinGlass; TeMinAir=x.TeMinAir;
	var=x.var; var1=x.var1; var2=x.var2; var3=x.var3; var4=x.var4; var5=x.var5;
	ExcludeVirtualRay=x.ExcludeVirtualRay;
	ExcludeVirtualObject=x.ExcludeVirtualObject;
	FindAThruRayTimeOut=x.FindAThruRayTimeOut;
	IgnoreTC=x.IgnoreTC;
	nu_list=x.nu_list;
	object_list=x.object_list;
}

void cLens1::ParaxialExpansion(int i,double &cx,double &cxy,double &cy){
	// i�ʌ`��̋ߎ��W�J�̌W�������߂�D z = (cx/2)x^2 +cxy*x*y +(cy/2)y^2
	// �y���Ӂz����ł�z������̉�]�ΐS�ɂ͑Ή����Ă��Ȃ��D
	double zcx,zcxy,zcy;

	cx=cy=c(i) + (asph_type(i)==CONIC ? a2(i)/NormH(i)/NormH(i) : 0);
	cxy=0;

	if(cylinder(i)==GENY) cy=0;
	if(cylinder(i)==GENX) cx=0;

	if(asph_type(i)==CONIC || asph_type(i)==TOROID || asph_type(i)==ANAMO || asph_type(i)==CONIC_OA){
		if(asph_type(i)==TOROID || asph_type(i)==ANAMO){
			cx= rx(i)==0 ? 0 : 1/rx(i);
			cy= ry(i)==0 ? 0 : 1/ry(i);
		}
		else if(asph_type(i)==CONIC_OA){
			matrix<double> T(3,3);
			surface_curvature(i,cx,cy,cxy,0,0,unit(T),0);
		}

		if(is_freeform_surface(i)){
			cx +=2*b(i,2,0);
			cxy+=b(i,1,1);
			cy +=2*b(i,0,2);
		}
		if(zernike(i).GetR0()>0){		
			zernike(i).ParaxialExpansion(zcx,zcxy,zcy);
			cx+=zcx;
			cxy+=zcxy;
			cy+=zcy;
		}
		if(is_spline_surface(i)){
			double zhh=surf[i].spline.y2(0); // ���_�ɂ�����T�O��2�������W�� = �ȗ�
			cx+=zhh;
			cy+=zhh;
		}
		if(UserDefSurf(i)){
			UserDefSurf2ndDerivative(zcx,zcxy,zcy,0,0);
			cx+=zcx;
			cxy+=zcxy;
			cy+=zcy;
		}
	}

	if(Newton(i) || As0(i) || As45(i)){
		double a,ea;

		ea=ea_max(i);
		if(ea>0){  // �����Ռ��̂Ƃ�ea<0
			a=(FRW/2/1000000)/(ea*ea/4);  // z=Newton*a*(x*x+y*y)
			                              // z=As0*a*(x*x-y*y)/2
			                              // z=As45*a*x*y
			cx +=Newton(i)*2*a+As0(i)*a;
			cxy+=As45(i)*a;
			cy +=Newton(i)*2*a-As0(i)*a;
		}
	}
}

void cLens1::AMatrixCalc(int i1,int i2,int j/*=1*/){
	int i,ii,jj;
	double ed, a1,a2,a3,a4, cx,cxy,cy, px,pxy,py, dN;
	double p[5][5],t[5][5],x[5][5];
	double p2[5],t2[5],x2[5];
	// �ߎ��s�� this->A,this->A2 ���v�Z����D
	//
	// | x' |       | x |
	// | y' | = A * | y | 
	// |��'x|       |��x|
	// |��'y|       |��y|
	//
	// 4x4�s��A�̓A�i�����t�B�b�N�n�ɑΉ��ł���D
	// 
	// |h' | = A2 * |h |
	// |��'| =      |��|
	//
    // 2x2�s��A2�͖{��������]�Ώ̂�O��Ƃ���D
	// �ߎ��p���[�Ƃ���x^2,y^2�̌W���̕��ς�p����D���̕��ς�x,y�ʓ��̉�]�ŕs��
	// �ł���C���Ȃ킿��o�������̃p���[�̕��ςł���D
	//
	// �ؖ��F
	//   �ʂ͋ߎ��ŁC
	//       z = axx*x^2 +axy*x*y +ayy*y^2
	//   �ƓW�J����C�K���ȍ��W�̉�]�ŁC
	//       z = axx'*x'^2 + ayy'*y'^2
	//   �ƂȂ�Caxx'��a'yy' ����ȗ��Ɗ֌W����D���̂Ƃ��C
	//       axx' = axx*(cos��)^2 +ayy*(sin��)^2 +axy*(sin��)(cos��)  ...(1)
	//       ayy' = axx*(sin��)^2 +ayy*(cos��)^2 -axy*(sin��)(cos��)  ...(2)
	//   �ł���D(1)(2)��������΁C
	//       axx' + ayy' = axx + ayy
	//   �ł��邩��D
	
	A[1][1]=1; A[1][2]=0; A[1][3]=0; A[1][4]=0;
	A[2][1]=0; A[2][2]=1; A[2][3]=0; A[2][4]=0;
	A[3][1]=0; A[3][2]=0; A[3][3]=1; A[3][4]=0;
	A[4][1]=0; A[4][2]=0; A[4][3]=0; A[4][4]=1;

	A2[1]=1; A2[2]=0;
	A2[3]=0; A2[4]=1;

	for(i=i1; i<=i2; i++){
		// ���܍s��p
		if(asph_type(i)==IDEAL){
			px=py= fideal(i)==0 ? 0 : 1/fideal(i);
			pxy=0;
		}
		else{
			ParaxialExpansion(i,cx,cxy,cy);
			dN=N(i,j)-N(i-1,j);
			px =dN*cx;
			pxy=dN*cxy;
			py =dN*cy;
		}

		p[1][1]=1;   p[1][2]=0;   p[1][3]=0; p[1][4]=0;
		p[2][1]=0;   p[2][2]=1;   p[2][3]=0; p[2][4]=0;
		p[3][1]=px;  p[3][2]=pxy; p[3][3]=1; p[3][4]=0;
		p[4][1]=pxy; p[4][2]=py;  p[4][3]=0; p[4][4]=1;

		p2[1]=1;         p2[2]=0;
		p2[3]=(px+py)/2; p2[4]=1;


		// �ڍs�s��t
		if(i<i2){
			if(IsGrin(i)){
				//  Selfoc�����Y�̏ꍇ�ߎ��ł� 
				//    y=a*sin((��A)z+��)  ...(1)
				//  �Ƃ���D����ƁC 
				//    y'=(��A)a*cos((��A)z+��)  ...(2)
				//  z=0��y=h,y'=-u�̏����������,
				//    sin��=h/a, cos��=-u/a/(��A)  ...(3)
				//  �ƂȂ�D
				//  (1)(2)�����@�藝�ŕό`���āC(3)��������Ύ��̓]�ڌ�����������D
				//      | h|=|      cos((��A)z)  -sin((��A)z)/(��A)/N || h|
				//      |Nu|=| (��A)sin((��A)z)N  cos((��A)z)         ||Nu|
				a1=cos(rA(i,j)*d(i));                a2=-sin(rA(i,j)*d(i))/rA(i,j)/N(i,j);
				a3=rA(i,j)*sin(rA(i,j)*d(i))*N(i,j); a4=cos(rA(i,j)*d(i));
				
				t[1][1]=a1; t[1][2]=0;  t[1][3]=a2; t[1][4]=0;
				t[2][1]=0;  t[2][2]=a1; t[2][3]=0;  t[2][4]=a2;
				t[3][1]=a3; t[3][2]=0;  t[3][3]=a4; t[3][4]=0;
				t[4][1]=0;  t[4][2]=a3; t[4][3]=0;  t[4][4]=a4;

				t2[1]=a1; t2[2]=a2;
				t2[3]=a3; t2[4]=a4;
			}
			else{
				ed=d(i)/N(i,j);
				t[1][1]=1; t[1][2]=0; t[1][3]=-ed; t[1][4]=0;
				t[2][1]=0; t[2][2]=1; t[2][3]=0;   t[2][4]=-ed;
				t[3][1]=0; t[3][2]=0; t[3][3]=1;   t[3][4]=0;
				t[4][1]=0; t[4][2]=0; t[4][3]=0;   t[4][4]=1;

				t2[1]=1; t2[2]=-ed;
				t2[3]=0; t2[4]=1;
			}
		}
		else{
			t[1][1]=1; t[1][2]=0; t[1][3]=0; t[1][4]=0;
			t[2][1]=0; t[2][2]=1; t[2][3]=0; t[2][4]=0;
			t[3][1]=0; t[3][2]=0; t[3][3]=1; t[3][4]=0;
			t[4][1]=0; t[4][2]=0; t[4][3]=0; t[4][4]=1;

			t2[1]=1; t2[2]=0;
			t2[3]=0; t2[4]=1;
		}
		
		for(ii=1; ii<=4; ii++) for(jj=1; jj<=4; jj++){
			// x = this->A
			x[ii][jj]=A[ii][jj];
		}
		for(ii=1; ii<=4; ii++) for(jj=1; jj<=4; jj++){
			// this->A = p*x
			A[ii][jj] =p[ii][1]*x[1][jj] +p[ii][2]*x[2][jj] +p[ii][3]*x[3][jj] +p[ii][4]*x[4][jj];
		}
		for(ii=1; ii<=4; ii++) for(jj=1; jj<=4; jj++){
			// x = this->A
			x[ii][jj]=A[ii][jj];
		}
		for(ii=1; ii<=4; ii++) for(jj=1; jj<=4; jj++){
			// this->A = t*x
			A[ii][jj] =t[ii][1]*x[1][jj] +t[ii][2]*x[2][jj] +t[ii][3]*x[3][jj] +t[ii][4]*x[4][jj];
		}

		for(ii=1; ii<=4; ii++) {
			// x2 = this->A2
			x2[ii]=A2[ii];
		}
		// this->A2 = p2*x2
		A2[1]=p2[1]*x2[1]+p2[2]*x2[3];
		A2[2]=p2[1]*x2[2]+p2[2]*x2[4];
		A2[3]=p2[3]*x2[1]+p2[4]*x2[3];
		A2[4]=p2[3]*x2[2]+p2[4]*x2[4];
		for(ii=1; ii<=4; ii++) {
			// x2 = this->A2
			x2[ii]=A2[ii];
		}
		// this->A2 = t2*x2
		A2[1]=t2[1]*x2[1]+t2[2]*x2[3];
		A2[2]=t2[1]*x2[2]+t2[2]*x2[4];
		A2[3]=t2[3]*x2[1]+t2[4]*x2[3];
		A2[4]=t2[3]*x2[2]+t2[4]*x2[4];
	}
}

void cLens1::InvertAMatrix(){
	int i,j;
	matrix<double> B;

	B.redim(4,4);
	for(i=1; i<=4; i++) for(j=1; j<=4; j++) B[i][j]=A[i][j];
	B=inv(B);
	for(i=1; i<=4; i++) for(j=1; j<=4; j++) A[i][j]=B[i][j];

	B.redim(2,2);
	B[1][1]=A2[1]; B[1][2]=A2[2];
	B[2][1]=A2[3]; B[2][2]=A2[4];
	B=inv(B);
	A2[1]=B[1][1]; A2[2]=B[1][2];
	A2[3]=B[2][1]; A2[4]=B[2][2];
}

int cLens1::is_plane(int i){
	return asph_type(i)==SPH && r(i)==0;
}

int cLens1::is_aspheric(int i){
	return asph_type(i)!=SPH || grating(i)!=0 || Diffusion(i)!=0;
}

int cLens1::is_ideallens(int i){
	if(asph_type(i)==IDEAL) return 1; else return 0;
}

int cLens1::b_size(int i) { return surf[i].b.N; }

int cLens1::is_freeform_surface(int i){
	switch(asph_type(i)){
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		return surf[i].b.NotNull();
	default:
		return 0;
	}
}

int cLens1::is_zernike_surface(int i){
	switch(asph_type(i)){
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		return Get_ZCoefficients(i)=="" ? 0 : 1;
	default:
		return 0;
	}
}

int cLens1::is_Dcon_surface(int i){
	switch(asph_type(i)){
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		return surf[i].Dcon.NotNull();
	default:
		return 0;
	}
}

int cLens1::is_legendre_surface(int i){
	switch(asph_type(i)){
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		return Get_LCoefficients(i)=="" ? 0 : 1;
	default:
		return 0;
	}
}

int cLens1::is_spline_surface(int i){
	switch(asph_type(i)){
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		return surf[i].spline.NotNull();
	default:
		return 0;
	}
}

int cLens1::is_periodic_surface(int i){
	int result=0;
	switch(asph_type(i)){
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		if(Apv(i)!=0) result=1;
		return result;
	default:
		return 0;
	}
}

int cLens1::is_userdef_surface(int i){
	int result=0;
	switch(asph_type(i)){
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		if(UserDefSurf(i)) result=1;
		return result;
	default:
		return 0;
	}
}

int cLens1::ea_is_defined(){
	int result=1;
	for(int i=1; i<=k; ++i) result*=( EAy(i)!=0 );
	return result;
}

int cLens1::ea_all_regular(){
	int result=1;
	for(int i=1; i<=k; ++i) result*=( EAx(i)==0 );
	return result;
}

int cLens1::surface_kind(int i){
	if( i<1 || k<i ) return 0;
	if( N(i-1,1)*N(i,1)>0 ) {
		if( EAy(i)<0 ) return 41;    // ���ܖʂł��邱�Ƃƃ}�X�N�ʂł��邱�Ƃ͔r���łȂ�
		if( asph_type(i)==IDEAL ) return 1;
		if( !is_solid(i-1) && !is_solid(i) && N(i-1,1)==N(i,1) ) return 3;
		else return 1;
	}
	else return 2;
}

int cLens1::is_refract_surf(int i){
	return surface_kind(i)%10==1 ? 1 : 0;
}

int cLens1::is_reflect_surf(int i){
	return surface_kind(i)==2 ? 1 : 0;
}

int cLens1::is_acting_surf(int i){
	return is_refract_surf(i) || is_reflect_surf(i) || grating(i) || i==0 || i==k+1;
	// ���̖ʂƑ��ʂ͖������ɐ^�Ƃ���D
	//   �������Ȃ��ƁCAllowVirtualRay���U�ł��Cs>0,���̋��N>0�̂Ƃ����_�����1��
	//   �Ɍ�������������r���ł��Ȃ��D
}

int cLens1::is_stop(int i){
	return surface_kind(i)==3 ? 1 : 0;
}

int cLens1::is_mask(int i){
	return surface_kind(i)==41 ? 1 : 0;
}

int cLens1::is_solid(int i){
	return fabs(N(i,1))>1.35 ? 1 : 0;
}

int cLens1::is_dummy(int i){
	return rem(i).substr(0,5)=="dummy" ? 1 : 0;  // rem(i)��"dummy"�Ŏn�܂邩
}

double cLens1::phi_y(int i){
	double phi;
	switch(EAtype(i)){
	case 0:
	case 1:
		phi= lens_phi(EAy(i));
		break;
	default:
		phi=0;
		break;
	}
	if(asph_type(i)==SPH){
		if( r(i)!=0 && fabs(phi)>=fabs(r(i))*2 ) phi=sgn(phi)*fabs(r(i))*2;
	}
	return phi;
}
double cLens1::phi_x(int i){
	double phi;
	switch(EAtype(i)){
	case 0:
	case 1:
		phi= EAx(i)==0 ? lens_phi(EAy(i)) : lens_phi(EAx(i));
		break;
	default:
		phi=0;
		break;
	}
	if(asph_type(i)==SPH){
		if( r(i)!=0 && fabs(phi)>=fabs(r(i))*2 ) phi=sgn(phi)*fabs(r(i))*2;
	}
	return phi;
}

double cLens1::ea_y(int i){
	double phi;
	switch(EAtype(i)){
	case 0:
	case 1:
		phi= EAy(i);
		break;
	default:
		phi=0;
		break;
	}
	if(asph_type(i)==SPH){
		if( r(i)!=0 && fabs(phi)>=fabs(r(i))*2 ) phi=sgn(phi)*fabs(r(i))*2;
	}
	return phi;
}
double cLens1::ea_x(int i){
	double phi;
	switch(EAtype(i)){
	case 0:
	case 1:
		phi= EAx(i)==0 ? EAy(i) : EAx(i);
		break;
	default:
		phi=0;
		break;
	}
	if(asph_type(i)==SPH){
		if( r(i)!=0 && fabs(phi)>=fabs(r(i))*2 ) phi=sgn(phi)*fabs(r(i))*2;
	}
	return phi;
}

double cLens1::ea_max(int i){
	// �j���[�g���Ȃ̌����͈͂̍ł��������a�͂���ɂȂ�D
	double phi;
	if( is_mask(i) ){
		phi=0;
	}
	else{
		switch(EAtype(i)){
			case 0:
				phi= EAy(i)>EAx(i) ? EAy(i) : EAx(i);
				break;
			case 1:
				if(EAx(i)==0){
					phi=sqrt(EAy(i)*EAy(i)*2);
				}
				else{
					phi=sqrt(EAy(i)*EAy(i)+EAx(i)*EAx(i));
				}
				break;
			default:
				phi=0;
				break;
		}
	}
	return phi;
}

double cLens1::freeform_z(int i,double x,double y){
	// ���R�Ȗʍ���z(x,y)
	// z(x,y)=��b(m,n)*x^m*y^n
	int m,n,N;
	double z=0;

	if(is_freeform_surface(i)){
		// �������̂���is_freeform_surface(i)���^�̏ꍇ�݈̂ȉ������s�D�傫�Ȍ��ʂ���D
		N=b_size(i);
		for(m=0; m<=N; m++) for(n=0; n<=N-m; n++){
			if(b(i,m,n)!=0) z+=b(i,m,n)*pw(x,m)*pw(y,n);
		}
	}
	return z;
}

double cLens1::freeform_zx(int i,double x,double y){
	// ���R�Ȗʍ���dz(x,y)/dx
	int m,n,N;
	double zx=0;
	
	if(is_freeform_surface(i)){
		// �������̂���is_freeform_surface(i)���^�̏ꍇ�݈̂ȉ������s
		N=b_size(i);
		for(m=1; m<=N; m++) for(n=0; n<=N-m; n++){
			if(b(i,m,n)!=0) zx+=b(i,m,n)*m*pw(x,m-1)*pw(y,n);
		}
	}
	return zx;
}

double cLens1::freeform_zy(int i,double x,double y){
	// ���R�Ȗʍ���dz(x,y)/dy
	int m,n,N;
	double zy=0;
	
	if(is_freeform_surface(i)){
		// �������̂���is_freeform_surface(i)���^�̏ꍇ�݈̂ȉ������s
		N=b_size(i);
		for(m=0; m<=N; m++) for(n=1; n<=N-m; n++){
			if(b(i,m,n)!=0) zy+=b(i,m,n)*pw(x,m)*n*pw(y,n-1);
		}
	}
	return zy;
}

double cLens1::freeform_zxx(int i,double x,double y){
	// ���R�Ȗʍ���dz(x,y)/dxdx
	int m,n,N;
	double zxx=0;

	if(is_freeform_surface(i)){
		// �������̂���is_freeform_surface(i)���^�̏ꍇ�݈̂ȉ������s
		N=b_size(i);
		for(m=2; m<=N; m++) for(n=0; n<=N-m; n++){
			if(b(i,m,n)!=0) zxx+=b(i,m,n)*m*(m-1)*pw(x,m-2)*pw(y,n);
		}
	}
	return zxx;
}

double cLens1::freeform_zyy(int i,double x,double y){
	// ���R�Ȗʍ���dz(x,y)/dydy
	int m,n,N;
	double zyy=0;

	if(is_freeform_surface(i)){
		// �������̂���is_freeform_surface(i)���^�̏ꍇ�݈̂ȉ������s
		N=b_size(i);
		for(m=0; m<=N; m++) for(n=2; n<=N-m; n++){
			if(b(i,m,n)!=0) zyy+=b(i,m,n)*pw(x,m)*n*(n-1)*pw(y,n-2);
		}
	}
	return zyy;
}

double cLens1::freeform_zxy(int i,double x,double y){
	// ���R�Ȗʍ���dz(x,y)/dxdy
	int m,n,N;
	double zxy=0;

	if(is_freeform_surface(i)){
		// �������̂���is_freeform_surface(i)���^�̏ꍇ�݈̂ȉ������s
		N=b_size(i);
		for(m=1; m<=N; m++) for(n=1; n<=N-m; n++){
			if(b(i,m,n)!=0) zxy+=b(i,m,n)*m*pw(x,m-1)*n*pw(y,n-1);
		}
	}
	return zxy;
}

double cLens1::periodic_z(int i,double x,double y){
	// periodic����z(x,y)
	// z(x,y)=-Apv[ (1/4)*{1+cos(2��*spx*x)}*{1+cos(2��*spy*y)} -1 ]  ...ZEMAX�Ɠ���(x=y=0 �� z=0)
	//    Apv = �U����peak to valley
	//    spx = x������Ԏ��g��
	//    spy = y������Ԏ��g��
	if(Apv(i)==0){
		return 0;
	}
	else{
		return -Apv(i)*( (1.0/4.0)*(1+cos(2*PI*sfx(i)*x))*(1+cos(2*PI*sfy(i)*y)) -1 );
	}
}

double cLens1::periodic_zx(int i,double x,double y){
	// periodic����dz(x,y)/dx
	if(Apv(i)==0){
		return 0;
	}
	else{
		return -Apv(i)*(1.0/4.0)*(-2*PI*sfx(i)*sin(2*PI*sfx(i)*x))*(1+cos(2*PI*sfy(i)*y));
	}
}

double cLens1::periodic_zy(int i,double x,double y){
	// periodic����dz(x,y)/dy
	if(Apv(i)==0){
		return 0;
	}
	else{
		return -Apv(i)*(1.0/4.0)*(1+cos(2*PI*sfx(i)*x))*(-2*PI*sfy(i)*sin(2*PI*sfy(i)*y));
	}
}

double cLens1::periodic_zxx(int i,double x,double y){
	// periodic����dz(x,y)/dxdx
	if(Apv(i)==0){
		return 0;
	}
	else{
		return -Apv(i)*(1.0/4.0)*(-4*PI*PI*sfx(i)*sfx(i)*cos(2*PI*sfx(i)*x))*(1+cos(2*PI*sfy(i)*y));
	}
}

double cLens1::periodic_zyy(int i,double x,double y){
	// periodic����dz(x,y)/dydy
	if(Apv(i)==0){
		return 0;
	}
	else{
		return -Apv(i)*(1.0/4.0)*(1+cos(2*PI*sfx(i)*x))*(-4*PI*PI*sfy(i)*sfy(i)*cos(2*PI*sfy(i)*y));
	}
}

double cLens1::periodic_zxy(int i,double x,double y){
	// periodic����dz(x,y)/dxdy
	if(Apv(i)==0){
		return 0;
	}
	else{
		return -Apv(i)*(1.0/4.0)*(-2*PI*sfx(i)*sin(2*PI*sfx(i)*x))*(-2*PI*sfy(i)*sin(2*PI*sfy(i)*y));
	}
}

inline double cLens1::surface_sag(int i,double y,double x,int newton_enable,int& domain_error) {
	// inline�w��͍������̌��ʂ��� (2013.01.04)
	double z,hh,H,HH,n,a,c;

	domain_error=0;
	c=this->c(i);
	
	switch(asph_type(i)){
	case SPH:
		if(c==0){
			z=0;
		}
		else{
			hh=y*y+x*x;
			a=1-hh*c*c;
			if(a>=0){  // a<0�ō����������ɂȂ���s����~����̂��������Da<0�ł̌��ʂ͕s��D
				z=hh*c/(1+sqrt(a));
			}
			else{
				domain_error=1; // (x,y)��`��G���[
			}
		}
		break;
	case CONIC: // ��]�Ώ̃R�[�j�b�N
		// z = chh/{1+sqrt(1-(kp+1)cchh)} +a4H^4 +a6H^6 +a8H^8 +a10H^10 + ... +a1H +a3H^3 + ... +a15H^15 +���R�Ȗʍ�   ...(1)
		//     �����ŁCh�͌������CH�͋K�i��������
		//
		// kp(i)�͂����� "Q�l", "�R�[�j�b�N�W��" �ł���D
		//   ���_�����_�ƈ�v���Cz���Ώ̂�2���Ȗʂ͈�ʂɁCyz���ʓ��ŁC
		//       y^2 + (kp+1)z^2 -2rz = 0         ...(2)
		//   �ƕ\�����Dr�͋ߎ��ȗ����a�ł���D�����(1)�̑�1���Ɠ����ł���D
		//   ����C�ȉ~�̕������́C
		//       (y/b)^2 + (z/a)^2 = 1     ...(3)
		//   �����C�����z�����ɕ��s�ړ����Ē��_�ƌ��_����v������΁C
		//      y^2 + ((b/a)^2)z^2 -2(b^2/a)z = 0    ...(4)
		//   �ƂȂ�D(2)(4)���C
		//      kp=(b/a)^2-1
		//   �ł���D
		//   ���l�ɁC�o�Ȗʂɂ��čl����D�������́C
		//       - (y/b)^2 + (z/a)^2 = 1   ...(5)
		//   ���_�ƌ��_����v������΁C
		//       y^2 - ((b/a)^2)z^2 �}2(b^2/a)z = 0   ...(6)
		//   (2)(6)���C
		//       kp=-(b/a)^2-1
		//   �ƂȂ�D
		//  
		//      kp<-1               �o�Ȗ� :  �œ_����̌��̔��˂͋��̗L�������ɏW��
		//      kp=-1    a=��       ������ :  �œ_����̌��̔��˂͕��s��
		//      -1<kp<0  0<b/a<1    �ȉ~�� :  �œ_����̌��̔��˂͎��̗L�������ɏW��
		//      kp=0     b=a        ����   :  �œ_����̌��̔��˂͏œ_�ɖ߂�
		//      0<kp     1<b/a      �΋��� :  �ȉ~��Z������ɉ񂵂��ʁD���ꂾ�����ʂ��steepening
		//      (������CodeV��Zemax�Ɠ����j
		//
		// ���Ȃ݂ɑȉ~�̗��S��   e=��(1-b^2/a^2) �Ƃ̊֌W�́C kp=-e^2
		//         �o�Ȗʂ̗��S�� e=��(1+b^2/a^2) �Ƃ̊֌W��,  kp=-e^2

		switch(cylinder(i)){
		case 0:    hh=y*y+x*x; break;
		case GENY: hh=x*x;     break;
		case GENX: hh=y*y;     break;
		}
		HH=hh/NormH(i)/NormH(i);  // �K�i�����ꂽ�������D
		                          // �����Y�a���傫���� �P�� an*h^n (h^n�����ł͂Ȃ��j���傫���l�ƂȂ�C
		                          // ��������̔g�����x�̐��x�Ƃ��邽�߂̌W��an�̕K�v�L���������傫���Ȃ�
		                          // �ꍇ�����邱�Ƃ��������� (2016/06/10)�D
		                          // �K�i�������� H=h/NormH �ɂ��, �ǂ̂��炢�̗L�������K�v����ڂŕ�����D
		H=sqrt(HH);  // ����񋅖ʗp
		a=1-(kp(i)+1)*c*c*hh;
		if(a>=0){  // a<0�ō����������ɂȂ���s����~����̂��������Da<0�ł̌��ʂ͕s��D
			n=sqrt(a);
			z=c*hh/(1+n)+a4(i)*HH*HH+a6(i)*HH*HH*HH+a8(i)*HH*HH*HH*HH+a10(i)*HH*HH*HH*HH*HH
			            +a12(i)*HH*HH*HH*HH*HH*HH+a14(i)*HH*HH*HH*HH*HH*HH*HH
						+a16(i)*HH*HH*HH*HH*HH*HH*HH*HH
						+a18(i)*HH*HH*HH*HH*HH*HH*HH*HH*HH
						+a20(i)*HH*HH*HH*HH*HH*HH*HH*HH*HH*HH
						+a1(i)*H +a2(i)*HH +a3(i)*HH*H +a5(i)*HH*HH*H +a7(i)*HH*HH*HH*H
						+a9(i)*HH*HH*HH*HH*H +a11(i)*HH*HH*HH*HH*HH*H +a13(i)*HH*HH*HH*HH*HH*HH*H +a15(i)*HH*HH*HH*HH*HH*HH*HH*H;
		}
		else{
			domain_error=1;
		}
		break;
	case TOROID: // �g���C�_��
		{	
			double kp;
			kp= IsXToroid(i) ? kpx(i) : kpy(i);
			z=ToroidZ(x,y,rx(i),ry(i),kp,IsXToroid(i),domain_error);
		}
		break;
	case ANAMO: // �A�i�����t�B�b�N�R�[�j�b�N
		z=AASZ(x,y,rx(i),ry(i),kpx(i),kpy(i),domain_error);
		break;
	case IDEAL:  // ���z�����Y
		z=0;
		break;
	case CONIC_OA:  // Off-Axial�񎟋Ȗ�
		z=this->ConicOffAxisZ(x,y,aCOA(i),bCOA(i),tCOA(i),domain_error);
		break;
	default: 
		break;
	}

	switch(asph_type(i)){
	case SPH:
		break;
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		z+=freeform_z(i,x,y);  // ���R�Ȗʐ���
		z+=zernike(i).Z(x,y)-zernike(i).Z(0,0);   // Zernike����(Z(0,0)�������Ȃ��ƌ��_��0�ɂȂ�Ȃ�) -> �s�X�g�����̍œK�����Ȃ���H
		z+=Dcon(i).Z(x,y);     // Frobes�ʐ���
		z+=legendre(i).Z(x,y)-legendre(i).Z(0,0); // Legendre������
		z+=spline(i).y(sqrt(x*x+y*y))-spline(i).y(0); // spline�ʐ���
		z+=periodic_z(i,x,y);  // periodic����
		if(UserDefSurf(i)) z+=UserDefSurfZ(x,y);  // ���[�U�[��`����
		break;
	}

	if(newton_enable) {
		switch(asph_type(i)){
		case IDEAL:
			break;
		default:
			if(Newton(i)!=0 || As0(i)!=0 || As45(i)!=0){
				double a,ea;

				ea=ea_max(i);
				if(ea>0){  // �����Ռ��̂Ƃ�ea<0
					a=(FRW/2/1000000)/(ea*ea/4);  // z=Newton*a*(x*x+y*y);
					                              // z=As0*a*(x*x-y*y)/2;
					                              // z=As45*a*x*y;
					z+= Newton(i)*a*(x*x+y*y) +As0(i)*a*(x*x-y*y)/2 +As45(i)*a*x*y;
				}
			}
			break;
		}
	}

	return z;
}

double cLens1::surface_sag(int i,double y,double x,int newton_enable){
	int dummy;
	return surface_sag(i,y,x,newton_enable,dummy);
}

inline vector<double> cLens1::surface_normal(int i,double y,double x,int newton_enable) {
	// |(l,m,n)|=1�͕ۏ؂��Ȃ��D �Ⴆ�΁Ccase1��c=0�̂Ƃ�n=1����, a4����0�Ȃ����m,n��0�łȂ��D
	// inline�w��͍������̌��ʂ��� (2013.01.04)
	double l,m,n, hh,H,HH,v;
	double c=this->c(i);

	switch(asph_type(i)){
	case SPH:
		if(c==0){
			n=1; l=m=0;
		}
		else{
			hh=y*y+x*x;
			n=sqrt( 1-hh*c*c );
			m=-y*c; l=-x*c;
		}
		break;
	case CONIC:
		// �ʂ̕�����z=f(x,y)
		// �A�֐��\���̓�(x,y,z)=z-f(x,y)=0
		// ��(x,y,z)=0�̖@��������(d��/dx,d��/dy,d��/dz)������C
		// �@������=(-df/dx,-df/dy,1)�ƂȂ�C(l,m,n)�͂����
		// n=sqrt( 1-(kp(i)+1)*hh*c*c ) ���悶�����́D
		switch(cylinder(i)){
		case 0:    hh=y*y+x*x; break;
		case GENY: hh=x*x;     break;
		case GENX: hh=y*y;     break;
		}
		HH=hh/NormH(i)/NormH(i);
		H=sqrt(HH);
		n=sqrt( 1-(kp(i)+1)*hh*c*c );
		v=c+n*(4*a4(i)*HH+6*a6(i)*HH*HH+8*a8(i)*HH*HH*HH+10*a10(i)*HH*HH*HH*HH
		       +12*a12(i)*HH*HH*HH*HH*HH+14*a14(i)*HH*HH*HH*HH*HH*HH
			   +16*a16(i)*HH*HH*HH*HH*HH*HH*HH
			   +18*a18(i)*HH*HH*HH*HH*HH*HH*HH*HH
			   +20*a20(i)*HH*HH*HH*HH*HH*HH*HH*HH*HH
			   +(H==0 ? 0 : a1(i)/H) +2*a2(i) +3*a3(i)*H +5*a5(i)*HH*H +7*a7(i)*HH*HH*H
			   +9*a9(i)*HH*HH*HH*H +11*a11(i)*HH*HH*HH*HH*H +13*a13(i)*HH*HH*HH*HH*HH*H +15*a15(i)*HH*HH*HH*HH*HH*HH*H)
			   /NormH(i)/NormH(i);
		switch(cylinder(i)){
		case 0:    m=-y*v; l=-x*v; break;
		case GENY: m=0;    l=-x*v; break;
		case GENX: m=-y*v; l=0;    break;
		}
		break;
	case TOROID:
		{
			double fx,fy,kp;
			kp= IsXToroid(i) ? kpx(i) : kpy(i);
			ToroidDerivative(fx,fy,x,y,rx(i),ry(i),kp,IsXToroid(i));
			l=-fx;
			m=-fy;
			n=1;
		}
		break;
	case ANAMO: 
		{
			double fx,fy;
			AASDerivative(fx,fy,x,y,rx(i),ry(i),kpx(i),kpy(i));
			l=-fx;
			m=-fy;
			n=1;
		}
		break;
	case IDEAL: 
		n=1; l=m=0;
		break;
	case CONIC_OA:
		{
			double fx,fy;
			ConicOffAxisDerivative(fx,fy,x,y,aCOA(i),bCOA(i),tCOA(i));
			l=-fx;
			m=-fy;
			n=1;
		}
	default: 
		break;
	}

	switch(asph_type(i)){
	case SPH:
		break;
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		l-= n*freeform_zx(i,x,y);  m-= n*freeform_zy(i,x,y);  // ���R�Ȗʐ���
		l-= n*zernike(i).Zx(x,y);  m-= n*zernike(i).Zy(x,y);  // zernike�ʐ���
		l-= n*Dcon(i).Zx(x,y);     m-= n*Dcon(i).Zy(x,y);     // Forbes�ʐ���
		l-= n*legendre(i).Zx(x,y); m-= n*legendre(i).Zy(x,y); // Legendre�ʐ���
		if(is_spline_surface(i)){  // spline����
			double h,zx,zy;
			h=sqrt(x*x+y*y);
			
			zx= h==0 ? 0 : spline(i).y1(h)*(x/h);
			zy= h==0 ? 0 : spline(i).y1(h)*(y/h);
			l-= n*zx; m-= n*zy;
		}
		l-= n*periodic_zx(i,x,y);  m-= n*periodic_zy(i,x,y);  // periodic����
		if(UserDefSurf(i)){
			double zx,zy;
			UserDefSurfDerivative(zx,zy,x,y);
			l-=n*zx; m-=n*zy;  // ���[�U�[��`����
		}
		break;
	}

	if(newton_enable) {
		switch(asph_type(i)){
		case IDEAL:
			break;
		default:
			if(Newton(i)!=0 || As0(i)!=0 || As45(i)!=0){
				double a,ea;

				ea=ea_max(i);
				if(ea>0){  // �����Ռ��̂Ƃ�ea<0
					a=(FRW/2/1000000)/(ea*ea/4); // z=Newton*a*(x*x+y*y);
					                             // z=As0*a*(x*x-y*y)/2;
					                             // z=As45*a*x*y;
					l-= n*(2*Newton(i)*a*x +As0(i)*a*x +As45(i)*a*y);
					m-= n*(2*Newton(i)*a*y -As0(i)*a*y +As45(i)*a*x);
				}
			}
			break;
		}
	}

	return vector<double>(l,m,n);
}

void cLens1::surface_curvature(int i,double& cx,double& cy,double& cxy,double x,double y,
                               const matrix<double>& T,int newton_enable){
	double fx,fy,fxx,fyy,fxy;
	double T11=T.a[1][1], T12=T.a[1][2], T13=T.a[1][3];
	double T21=T.a[2][1], T22=T.a[2][2], T23=T.a[2][3];
	double T31=T.a[3][1], T32=T.a[3][2], T33=T.a[3][3];

	switch(asph_type(i)){
	case SPH:
		if( newton_enable && (Newton(i)!=0 || As0(i)!=0 || As45(i)!=0) ){
			{
				double c,n, hh,hhx,hhy,hhxx,hhyy, zhh,zhhhh;
				c=this->c(i);
				hh=y*y+x*x;
				n=sqrt( 1-hh*c*c );
				zhh=c/2/n;            // z��hh�̂�����
				zhhhh=c*c*c/4/n/n/n;  // z��hh�ɂ��2�K����
				hhx=2*x; hhy=2*y; hhxx=2; hhyy=2;
				fx=zhh*hhx;
				fy=zhh*hhy;
				fxx=zhhhh*hhx*hhx+zhh*hhxx;
				fyy=zhhhh*hhy*hhy+zhh*hhyy;
				fxy=zhhhh*hhx*hhy;
				break;
			}			
		}
		else{
			cx=cy=this->c(i); cxy=0;
			return; // newton_enable�łȂ�SPH�ʂ͂�����return����D
		}
	case CONIC:
		{
			double c,n, hh,hhx,hhy,hhxx,hhyy, zhh,zhhhh, HH,H;
			c=this->c(i);
			switch(cylinder(i)){
			case 0:    hh=y*y+x*x; break;
			case GENY: hh=x*x;     break;
			case GENX: hh=y*y;     break;
			}
			HH=hh/NormH(i)/NormH(i);
			H=sqrt(HH);
			n=sqrt( 1-(kp(i)+1)*hh*c*c );
			zhh=c/2/n
				+(2*a4(i)*HH+3*a6(i)*HH*HH+4*a8(i)*HH*HH*HH+5*a10(i)*HH*HH*HH*HH
			    +6*a12(i)*HH*HH*HH*HH*HH
				+7*a14(i)*HH*HH*HH*HH*HH*HH
				+8*a16(i)*HH*HH*HH*HH*HH*HH*HH
				+9*a18(i)*HH*HH*HH*HH*HH*HH*HH*HH
				+10*a20(i)*HH*HH*HH*HH*HH*HH*HH*HH*HH
				+(H==0 ? 0 : a1(i)/2/H)
				+a2(i) +1.5*a3(i)*H +2.5*a5(i)*HH*H +3.5*a7(i)*HH*HH*H
				+4.5*a9(i)*HH*HH*HH*H +5.5*a11(i)*HH*HH*HH*HH*H +6.5*a13(i)*HH*HH*HH*HH*HH*H +7.5*a15(i)*HH*HH*HH*HH*HH*HH*H)
				/NormH(i)/NormH(i);   // z��hh�ɂ�����
			zhhhh=(kp(i)+1)*c*c*c/4/n/n/n
				 +(2*a4(i)+6*a6(i)*HH+12*a8(i)*HH*HH+20*a10(i)*HH*HH*HH
			     +30*a12(i)*HH*HH*HH*HH
				 +42*a14(i)*HH*HH*HH*HH*HH
				 +56*a16(i)*HH*HH*HH*HH*HH*HH
				 +72*a18(i)*HH*HH*HH*HH*HH*HH*HH
				 +90*a20(i)*HH*HH*HH*HH*HH*HH*HH*HH
				 +(H==0 ? 0 : -0.25*a1(i)/HH/H +0.75*a3(i)/H)
				 +3.75*a5(i)*H +8.75*a7(i)*HH*H +15.75*a9(i)*HH*HH*H +24.75*a11(i)*HH*HH*HH*H
				 +35.75*a13(i)*HH*HH*HH*HH*H +48.75*a15(i)*HH*HH*HH*HH*HH*H)
				 /NormH(i)/NormH(i)/NormH(i)/NormH(i); // z��hh�ɂ��2�K����
			switch(cylinder(i)){
			case 0:    hhx=2*x; hhy=2*y; hhxx=2; hhyy=2; break;
			case GENY: hhx=2*x; hhy=0;   hhxx=2; hhyy=0; break;
			case GENX: hhx=0;   hhy=2*y; hhxx=0; hhyy=2; break;
			}
			fx=zhh*hhx;
			fy=zhh*hhy;
			fxx=zhhhh*hhx*hhx+zhh*hhxx;
			fyy=zhhhh*hhy*hhy+zhh*hhyy;
			fxy=zhhhh*hhx*hhy;
			break;
		}
	case TOROID:
		{
			double kp;
			kp= IsXToroid(i) ? kpx(i) : kpy(i);
			Toroid2ndDerivative(fxx,fxy,fyy,fx,fy,x,y,rx(i),ry(i),kp,IsXToroid(i));
		}
		break;
	case ANAMO:
		AAS2ndDerivative(fxx,fxy,fyy,fx,fy,x,y,rx(i),ry(i),kpx(i),kpy(i));
		break;
	case CONIC_OA:
		ConicOffAxis2ndDerivative(fxx,fxy,fyy,fx,fy,x,y,aCOA(i),bCOA(i),tCOA(i));
		break;
	default:
		fxx=fyy=fxy=0;
		break;
	}

	switch(asph_type(i)){
	case SPH:
		break;
	case CONIC:
	case TOROID:
	case ANAMO:
	case CONIC_OA:
		// ���R�Ȗʍ��̕���������
		fx+=  freeform_zx(i,x,y);
		fy+=  freeform_zy(i,x,y);
		fxx+= freeform_zxx(i,x,y);
		fyy+= freeform_zyy(i,x,y);
		fxy+= freeform_zxy(i,x,y);
		// Forbes�ʐ�����������
		fx+=  Dcon(i).Zx(x,y);
		fy+=  Dcon(i).Zy(x,y);
		fxx+= Dcon(i).Zxx(x,y);
		fyy+= Dcon(i).Zyy(x,y);
		fxy+= Dcon(i).Zxy(x,y);
		// zernike�ʐ�����������
		fx+=  zernike(i).Zx(x,y);
		fy+=  zernike(i).Zy(x,y);
		fxx+= zernike(i).Zxx(x,y);
		fyy+= zernike(i).Zyy(x,y);
		fxy+= zernike(i).Zxy(x,y);
		// legendre�ʐ�����������
		fx+=  legendre(i).Zx(x,y);
		fy+=  legendre(i).Zy(x,y);
		fxx+= legendre(i).Zxx(x,y);
		fyy+= legendre(i).Zyy(x,y);
		fxy+= legendre(i).Zxy(x,y);
		if(is_spline_surface(i)){  // spline���̐�����������
			double zh,zhh, h,hx,hy,hxx,hxy,hyy, zx,zy,zxx,zxy,zyy;
			
			h=sqrt(x*x+y*y);
			zh=spline(i).y1(h);    // ��z/��h
			zhh=spline(i).y2(h);   // ��(zh)/��h

			if(h==0){
				zx=0;
				zy=0;
				zxx=zhh;
				zxy=0;        // ��]�Ώ̂�����
				zyy=zhh;
			}
			else{
				hx=x/h;
				hy=y/h;
				hxx=(h*h-x*x)/(h*h*h);
				hxy=-x*y/(h*h*h);
				hyy=(h*h-y*y)/(h*h*h);

				zx=zh*hx;
				zy=zh*hy;
				zxx=(zhh*hx)*hx+zh*hxx;  // zxx = ��(zh*hx)/��x = (��(zh)/��x)*hx + zh*hxx = (zhh*hx)*hx+zh*hxx
				zxy=(zhh*hy)*hx+zh*hxy;
				zyy=(zhh*hy)*hy+zh*hyy;
			}

			fx+=zx;
			fy+=zy;
			fxx+=zxx;
			fxy+=zxy;
			fyy+=zyy;
		}
		// periodic���̕���������
		fx+=  periodic_zx(i,x,y);
		fy+=  periodic_zy(i,x,y);
		fxx+= periodic_zxx(i,x,y);
		fyy+= periodic_zyy(i,x,y);
		fxy+= periodic_zxy(i,x,y);
		// ���[�U�[��`������������
		if(UserDefSurf(i)){
			double zx,zy,zxx,zxy,zyy;
			
			UserDefSurfDerivative(zx,zy,x,y);
			UserDefSurf2ndDerivative(zxx,zxy,zyy,x,y);
			fx+=zx;
			fy+=zy;
			fxx+=zxx;
			fxy+=zxy;
			fyy+=zyy;
		}
		break;
	}
	
	if(newton_enable) {
		switch(asph_type(i)){
		case IDEAL:
			break;
		default:
			if(Newton(i)!=0 || As0(i)!=0 || As45(i)!=0){
				double a,ea, A,B,C;

				ea=ea_max(i);
				if(ea>0){  // �����Ռ��̂Ƃ�ea<0
					a=(FRW/2/1000000)/(ea*ea/4); // z=Newton*a*(x*x+y*y);
					                             // z=As0*a*(x*x-y*y)/2;
					                             // z=As45*a*x*y;
					A=Newton(i)*a;
					B=As0(i)*a/2;
					C=As45(i)*a;
					// z�̑��� = A*(x*x+y*y) +B*(x*x-y*y) +C*x*y;
					fx +=2*A*x+2*B*x+C*y;
					fy +=2*A*y-2*B*y+C*x;
					fxx+=2*A+2*B;
					fyy+=2*A-2*B;
					fxy+=C;
				}
			}
			break;
		}
	}
	
	// �ȉ��́g�����Y���w(����)�h�� 3.4.4 ���Q�l�ɂ���Γ����o����D
	// ��̓I�ɂ� (3.4.54) �� x=... ���������D�������C�����ł�x���������ɂ��Ă���̂Ō��ʂ̎��͓����łȂ��D
	cx =( T12*T12*fyy +2*T12*T11*fxy +T11*T11*fxx )/( T33 -T32*fy -T31*fx );
	cy =( T22*T22*fyy +2*T21*T22*fxy +T21*T21*fxx )/( T33 -T32*fy -T31*fx );
	cxy=( T12*T22*fyy +(T22*T11+T21*T12)*fxy +T21*T11*fxx )/( T33 -T32*fy -T31*fx ); 

	return;
}


// �ΐS�ɂ��x�N�g���̕ϊ�
//   note 1 : �ȉ��Ƃ͋t�ɁC�Ⴆ�� decenter_in(int,double&, ..) ���� decenter_in(int,complex&, ..)
//            ���ĂԂ悤�ɍ��ƁC
//            (�������킴�킴���f���Ƃ��Čv�Z���邽��)���s���x�����ɒx���Ȃ�D
//   note 2 : ���x�N�g���ƕ��f�x�N�g�������ɑΉ����邽�߁C�e���v���[�g���g���l�������邪�C
//            ���̏ꍇ�C�֐���`���N���X�錾��(�܂� cLens1.h ��)�ɏ����Ȃ��ƂȂ���
//            �R���p�C���G���[�ƂȂ�D���ꂢ�łȂ��̂ł�߂��D
//
// order==0 �̂Ƃ��́C���s�ړ��Crox��]�Croy��]�Croz��]�̏��ɍs�Ȃ��D
// order!=0 �̂Ƃ��́Croz��]�Croy��]�Crox��]�C���s�ړ��̏��ɍs�Ȃ��D
//
// �� �ێ琫�ێ����߁C�edecenter_type(i)�̓��e�́C
//    decenter_in(_rev),decenter_out(_rev)�֐��̒������ɏ����悤�ɂ���D
//
void cLens1::decenter_in(int i,double& x,double& y,double& z,int translate){
	vector<double> v;
	switch( decenter_type(i) ){
	case 1:
	case 2:
	case 3:
		if(order(i)==0){
			if(translate){ x-=dx(i); y-=dy(i); z-=dz(i); }
			v.x=x; v.y=y; v.z=z;
			if(rox(i)!=0) x_rotate(v,-rox(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,-roy(i)*PI/180);
			if(roz(i)!=0) z_rotate(v,-roz(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
		}
		else{
			v.x=x; v.y=y; v.z=z;
			if(roz(i)!=0) z_rotate(v,-roz(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,-roy(i)*PI/180);
			if(rox(i)!=0) x_rotate(v,-rox(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
			if(translate){ x-=dx(i); y-=dy(i); z-=dz(i); }
		}
		break;
	}
}

void cLens1::decenter_in(int i,complex& x,complex& y,complex& z,int translate){
	vector<double> v1(Re(x),Re(y),Re(z)), v2(Im(x),Im(y),Im(z));
	decenter_in(i,v1.x,v1.y,v1.z,translate);
	decenter_in(i,v2.x,v2.y,v2.z,translate);
	x=complex(v1.x,v2.x);
	y=complex(v1.y,v2.y);
	z=complex(v1.z,v2.z);
}

void cLens1::decenter_in_rev(int i,double& x,double& y,double& z,int translate){
	vector<double> v;
	switch( decenter_type(i) ){
	case 1:
	case 2:
	case 3:
		if(order(i)==0){
			v.x=x; v.y=y; v.z=z;
			if(roz(i)!=0) z_rotate(v,roz(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,roy(i)*PI/180);
			if(rox(i)!=0) x_rotate(v,rox(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
			if(translate){ x+=dx(i); y+=dy(i); z+=dz(i); }
		}
		else{
			if(translate){ x+=dx(i); y+=dy(i); z+=dz(i); }
			v.x=x; v.y=y; v.z=z;
			if(rox(i)!=0) x_rotate(v,rox(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,roy(i)*PI/180);
			if(roz(i)!=0) z_rotate(v,roz(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
		}
		break;
	}
}

void cLens1::decenter_in_rev(int i,complex& x,complex& y,complex& z,int translate){
	vector<double> v1(Re(x),Re(y),Re(z)), v2(Im(x),Im(y),Im(z));
	decenter_in_rev(i,v1.x,v1.y,v1.z,translate);
	decenter_in_rev(i,v2.x,v2.y,v2.z,translate);
	x=complex(v1.x,v2.x);
	y=complex(v1.y,v2.y);
	z=complex(v1.z,v2.z);
}

void cLens1::decenter_out(int i,double& x,double& y,double& z,int translate){
	vector<double> v;
	switch( decenter_type(i) ){
	case 1:
		v.x=x; v.y=y; v.z=z;
		v=ret_decenter(v,translate,i,0);
		x=v.x; y=v.y; z=v.z;
		if(order1(i)==0){
			if(translate){ x-=dx1(i); y-=dy1(i); z-=dz1(i); }
			v.x=x; v.y=y; v.z=z;
			if(rox1(i)!=0) x_rotate(v,-rox1(i)*PI/180);
			if(roy1(i)!=0) y_rotate(v,-roy1(i)*PI/180);
			if(roz1(i)!=0) z_rotate(v,-roz1(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
		}
		else{
			v.x=x; v.y=y; v.z=z;
			if(roz1(i)!=0) z_rotate(v,-roz1(i)*PI/180);
			if(roy1(i)!=0) y_rotate(v,-roy1(i)*PI/180);
			if(rox1(i)!=0) x_rotate(v,-rox1(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
			if(translate){ x-=dx1(i); y-=dy1(i); z-=dz1(i); }
		}
		break;
	case 2:
		if(order(i)==0){
			v.x=x; v.y=y; v.z=z;
			if(roz(i)!=0) z_rotate(v,roz(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,roy(i)*PI/180);
			if(rox(i)!=0) x_rotate(v,rox(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
			if(translate){ x+=dx(i); y+=dy(i); z+=dz(i); }
		}
		else{
			if(translate){ x+=dx(i); y+=dy(i); z+=dz(i); }
			v.x=x; v.y=y; v.z=z;
			if(rox(i)!=0) x_rotate(v,rox(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,roy(i)*PI/180);
			if(roz(i)!=0) z_rotate(v,roz(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
		}
		break;
	case 3:
		if(order(i)==0){
			if(translate){ x-=dx(i); y-=dy(i); z-=dz(i); }
			v.x=x; v.y=y; v.z=z;
			if(rox(i)!=0) x_rotate(v,-rox(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,-roy(i)*PI/180);
			if(roz(i)!=0) z_rotate(v,-roz(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
		}
		else{
			v.x=x; v.y=y; v.z=z;
			if(roz(i)!=0) z_rotate(v,-roz(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,-roy(i)*PI/180);
			if(rox(i)!=0) x_rotate(v,-rox(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
			if(translate){ x-=dx(i); y-=dy(i); z-=dz(i); }
		}
		break;
	}
}

void cLens1::decenter_out(int i,complex& x,complex& y,complex& z,int translate){
	vector<double> v1(Re(x),Re(y),Re(z)), v2(Im(x),Im(y),Im(z));
	decenter_out(i,v1.x,v1.y,v1.z,translate);
	decenter_out(i,v2.x,v2.y,v2.z,translate);
	x=complex(v1.x,v2.x);
	y=complex(v1.y,v2.y);
	z=complex(v1.z,v2.z);
}

void cLens1::decenter_out_rev(int i,double& x,double& y,double& z,int translate){
	// ���̕��@�Ƃ��āC
	// ���ϊ��̕��s�ړ��ȊO������ v'=f(v)=Av �Ƃ���ƁCA�̗v�f�́C
	//  | a11 | 
	//  | a21 | = f(ex), ex��x�����P�ʃx�N�g��
	//  | a31 |
	// �ƁCy,z�ɂ��Ă̓��l�Ȍv�Z�ŋ��߂��C�t�ϊ��� v=inv(A)v' �ƂȂ�D
	// �������C���̕��@�̓R�[�h�͊Ȍ��Ȃ��̂́C�������Ԃ������邽�߂�߂��D 10.01.29

	vector<double> v;
	switch( decenter_type(i) ){
	case 1:
		if(order1(i)==0){
			v.x=x; v.y=y; v.z=z;
			if(roz1(i)!=0) z_rotate(v,roz1(i)*PI/180);
			if(roy1(i)!=0) y_rotate(v,roy1(i)*PI/180);
			if(rox1(i)!=0) x_rotate(v,rox1(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
			if(translate){ x+=dx1(i); y+=dy1(i); z+=dz1(i); }
		}
		else{
			if(translate){ x+=dx1(i); y+=dy1(i); z+=dz1(i); }
			v.x=x; v.y=y; v.z=z;
			if(rox1(i)!=0) x_rotate(v,rox1(i)*PI/180);
			if(roy1(i)!=0) y_rotate(v,roy1(i)*PI/180);
			if(roz1(i)!=0) z_rotate(v,roz1(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
		}
		v.x=x; v.y=y; v.z=z;
		v=ret_decenter(v,translate,i,1);
		x=v.x; y=v.y; z=v.z;
		break;
	case 2:
		if(order(i)==0){
			if(translate){ x-=dx(i); y-=dy(i); z-=dz(i); }
			v.x=x; v.y=y; v.z=z;
			if(rox(i)!=0) x_rotate(v,-rox(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,-roy(i)*PI/180);
			if(roz(i)!=0) z_rotate(v,-roz(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
		}
		else{
			v.x=x; v.y=y; v.z=z;
			if(roz(i)!=0) z_rotate(v,-roz(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,-roy(i)*PI/180);
			if(rox(i)!=0) x_rotate(v,-rox(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
			if(translate){ x-=dx(i); y-=dy(i); z-=dz(i); }
		}
		break;
	case 3:
		if(order(i)==0){
			v.x=x; v.y=y; v.z=z;
			if(roz(i)!=0) z_rotate(v,roz(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,roy(i)*PI/180);
			if(rox(i)!=0) x_rotate(v,rox(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
			if(translate){ x+=dx(i); y+=dy(i); z+=dz(i); }
		}
		else{
			if(translate){ x+=dx(i); y+=dy(i); z+=dz(i); }
			v.x=x; v.y=y; v.z=z;
			if(rox(i)!=0) x_rotate(v,rox(i)*PI/180);
			if(roy(i)!=0) y_rotate(v,roy(i)*PI/180);
			if(roz(i)!=0) z_rotate(v,roz(i)*PI/180);
			x=v.x; y=v.y; z=v.z;
		}
		break;
	}
}

void cLens1::decenter_out_rev(int i,complex& x,complex& y,complex& z,int translate){
	vector<double> v1(Re(x),Re(y),Re(z)), v2(Im(x),Im(y),Im(z));
	decenter_out_rev(i,v1.x,v1.y,v1.z,translate);
	decenter_out_rev(i,v2.x,v2.y,v2.z,translate);
	x=complex(v1.x,v2.x);
	y=complex(v1.y,v2.y);
	z=complex(v1.z,v2.z);
}

vector<double> cLens1::ret_decenter(vector<double> v,int translate,int i,int inv){
	// v��i�ʕΐS����W�����ret(i)�ʂ̕ΐS�O���W��ɖ߂�
	// translate =true  : v���ʒu�x�N�g���Ƃ��Ĉ���
	//                    (���_��i�ʌ��_���z������ret(i)�ʂ܂ł̖ʊԊu��ݐς����ʒu�ɂȂ�)
	//           =false : v������x�N�g���Ƃ��Ĉ���
	//
	// inv =true �̂Ƃ��͋t�ϊ����s���D

	int ii;

	if(decenter_type(i)!=1) return v;
	if(ret(i)<1 || i<ret(i) || i>k) return v;
		
	if(inv==0){
		decenter_in_rev(i,v.x,v.y,v.z,translate);

		for(ii=i-1; ii>=ret(i); ii--){
			if(translate) v.z+=d(ii);
			decenter_out_rev(ii,v.x,v.y,v.z,translate);
			decenter_in_rev(ii,v.x,v.y,v.z,translate);
		}

		if(translate) for(ii=ret(i); ii<=i-1; ii++) v.z-=d(ii);

		return v;
	}
	else{
		if(translate) for(ii=ret(i); ii<=i-1; ii++) v.z+=d(ii);
		
		for(ii=ret(i); ii<=i-1; ii++){
			decenter_in(ii,v.x,v.y,v.z,translate);
			decenter_out(ii,v.x,v.y,v.z,translate);
			if(translate) v.z-=d(ii);
		}

		decenter_in(i,v.x,v.y,v.z,translate);
		
		return v;
	}
}

vector<complex> cLens1::ret_decenter(vector<complex> v,int translate,int i,int inv){
	vector<double> x(Re(v.x),Re(v.y),Re(v.z));
	vector<double> y(Im(v.x),Im(v.y),Im(v.z));
	ret_decenter(x,translate,i,inv);
	ret_decenter(y,translate,i,inv);
	v.x=complex(x.x,y.x);
	v.y=complex(x.y,y.y);
	v.z=complex(x.z,y.z);
	return v;
}

void cLens1::make_coordinate(double defocus){
	// �e�ʂ̍��W����ݒ肷��D�O���[�o�����W��p����D
	//     o,o0 : �O���[�o�����W�ɂ�郍�[�J�����W���_�ʒu
	//     ex,ex0, ey,ey0, ez,ez0 : �O���[�o�����W�ɂ��x,y,z�e�������̒P�ʃx�N�g�� 
	//     (�Y��0�͕ΐS�O��\��) 
	
	// ���̊֐��������ǐՈ�{���Ɏ��s���邽�ߎ��s���x�̏�ő傫�ȏ�Q�ƂȂ��Ă���D
	// �������̂��߁C�x�N�g�����Z�𐬕����ɏ����Ă���D

	int i;
	vector<double> O,X,Y,Z, v, X0,Y0,Z0;

	const vector<double> Eo=vector<double>(0,0,0);
	const vector<double> Ex=vector<double>(1,0,0); 
	const vector<double> Ey=vector<double>(0,1,0); 
	const vector<double> Ez=vector<double>(0,0,1);
	                       
	d(0)= fabs(s)>=LN ? 0 : -s;
	d(k)= dk(defocus);

	// o0[0]=Eo+(0,0,-d(0)), ex0[0]=Ex ey0[0]=Ey ez0[0]=Ez
	//  "-d(0)"�̐���
	//     ����ő����̏ꍇ�C���_����1�ʂ̒��_�ɂȂ�D��O�͕ΐS�̂���ꍇ�ł���D
	//     ���������āC���̖ʂɌ��_��u���ꍇ�Ƃ͈قȂ�C
	//     ���̋�����ς��Ă������W�Ń����Y�̈ʒu�͕ς��Ȃ��D
	//     �Ⴆ�΁CMakeLensView()�ŕ��_��ς�������`�����Ƃ��ł���D
	o0[0].x=Eo.x;  o0[0].y=Eo.y;  o0[0].z=Eo.z-d(0);
	ex0[0].x=Ex.x; ex0[0].y=Ex.y; ex0[0].z=Ex.z;
	ey0[0].x=Ey.x; ey0[0].y=Ey.y; ey0[0].z=Ey.z;
	ez0[0].x=Ez.x; ez0[0].y=Ez.y; ez0[0].z=Ez.z;
	
	for(i=0; i<=k+1; ++i){	
		// ���ܑO�ΐS��̍��W
		if(decenter_type(i)!=0){
			// ���_�ʒu
			// v=Eo
			v.x=Eo.x; v.y=Eo.y; v.z=Eo.z;
			decenter_in_rev(i,v.x,v.y,v.z,1);
			// o[i]=o0[i]+v.x*ex0[i]+v.y*ey0[i]+v.z*ez0[i]
			o[i].x=o0[i].x+v.x*ex0[i].x+v.y*ey0[i].x+v.z*ez0[i].x;
			o[i].y=o0[i].y+v.x*ex0[i].y+v.y*ey0[i].y+v.z*ez0[i].y;
			o[i].z=o0[i].z+v.x*ex0[i].z+v.y*ey0[i].z+v.z*ez0[i].z;
			
			// ���W������
			// X0=ex0 Y0=ey0 Z0=ez0;
			X0.x=ex0[i].x; X0.y=ex0[i].y; X0.z=ex0[i].z;
			Y0.x=ey0[i].x; Y0.y=ey0[i].y; Y0.z=ey0[i].z;
			Z0.x=ez0[i].x; Z0.y=ez0[i].y; Z0.z=ez0[i].z;

			// v=Ex
			v.x=Ex.x; v.y=Ex.y; v.z=Ex.z;
			decenter_in_rev(i,v.x,v.y,v.z,0);
			// X=v.x*X0+v.y*Y0+v.z*Z0
			X.x=v.x*X0.x+v.y*Y0.x+v.z*Z0.x;
			X.y=v.x*X0.y+v.y*Y0.y+v.z*Z0.y;
			X.z=v.x*X0.z+v.y*Y0.z+v.z*Z0.z;
			
			// v=Ey
			v.x=Ey.x; v.y=Ey.y; v.z=Ey.z;
			decenter_in_rev(i,v.x,v.y,v.z,0);
			// Y=v.x*X0+v.y*Y0+v.z*Z0
			Y.x=v.x*X0.x+v.y*Y0.x+v.z*Z0.x;
			Y.y=v.x*X0.y+v.y*Y0.y+v.z*Z0.y;
			Y.z=v.x*X0.z+v.y*Y0.z+v.z*Z0.z;
			
			// Z=XxY
			Z.x=X.y*Y.z-X.z*Y.y;
			Z.y=X.z*Y.x-X.x*Y.z;
			Z.z=X.x*Y.y-X.y*Y.x;

			// ex[i]=X ey[i]=Y ez[i]=Z
			ex[i].x=X.x; ex[i].y=X.y; ex[i].z=X.z;
			ey[i].x=Y.x; ey[i].y=Y.y; ey[i].z=Y.z;
			ez[i].x=Z.x; ez[i].y=Z.y; ez[i].z=Z.z;
		}
		else{
			// o[i]=o0[i] ex[i]=ex0[i] ey[i]=ey0[i] ez[i]=ez0[i]
			o[i].x=o0[i].x; o[i].y=o0[i].y; o[i].z=o0[i].z;
			ex[i].x=ex0[i].x; ex[i].y=ex0[i].y; ex[i].z=ex0[i].z;
			ey[i].x=ey0[i].x; ey[i].y=ey0[i].y; ey[i].z=ey0[i].z;
			ez[i].x=ez0[i].x; ez[i].y=ez0[i].y; ez[i].z=ez0[i].z;
		}

		// ���܌�ΐS��̍��W
		if(i<=k){
			if(decenter_type(i)!=0){
				// ���_�ʒu
				// v=Eo
				v.x=Eo.x; v.y=Eo.y; v.z=Eo.z;
				decenter_out_rev(i,v.x,v.y,v.z,1);
				// o0[i+1]=o[i]+v.x*ex[i]+v.y*ey[i]+v.z*ez[i]
				o0[i+1].x=o[i].x+v.x*ex[i].x+v.y*ey[i].x+v.z*ez[i].x;
				o0[i+1].y=o[i].y+v.x*ex[i].y+v.y*ey[i].y+v.z*ez[i].y;
				o0[i+1].z=o[i].z+v.x*ex[i].z+v.y*ey[i].z+v.z*ez[i].z;

				// ���W������
				// X0=X Y0=Y Z0=Z
				X0.x=ex[i].x; X0.y=ex[i].y; X0.z=ex[i].z;
				Y0.x=ey[i].x; Y0.y=ey[i].y; Y0.z=ey[i].z;
				Z0.x=ez[i].x; Z0.y=ez[i].y; Z0.z=ez[i].z;

				// v=Ex
				v.x=Ex.x; v.y=Ex.y; v.z=Ex.z;
				decenter_out_rev(i,v.x,v.y,v.z,0);
				// X=v.x*X0+v.y*Y0+v.z*Z0
				X.x=v.x*X0.x+v.y*Y0.x+v.z*Z0.x;
				X.y=v.x*X0.y+v.y*Y0.y+v.z*Z0.y;
				X.z=v.x*X0.z+v.y*Y0.z+v.z*Z0.z;
				
				// v=Ey
				v.x=Ey.x; v.y=Ey.y; v.z=Ey.z;
				decenter_out_rev(i,v.x,v.y,v.z,0);
				// Y=v.x*X0+v.y*Y0+v.z*Z0
				Y.x=v.x*X0.x+v.y*Y0.x+v.z*Z0.x;
				Y.y=v.x*X0.y+v.y*Y0.y+v.z*Z0.y;
				Y.z=v.x*X0.z+v.y*Y0.z+v.z*Z0.z;
				
				// Z=XxY
				Z.x=X.y*Y.z-X.z*Y.y;
				Z.y=X.z*Y.x-X.x*Y.z;
				Z.z=X.x*Y.y-X.y*Y.x;

				// ���̖ʂ̕ΐS�O���W
				// ex0[i+1]=X ey0[i+1]=Y ez0[i+1]=Z
				ex0[i+1].x=X.x; ex0[i+1].y=X.y; ex0[i+1].z=X.z;
				ey0[i+1].x=Y.x; ey0[i+1].y=Y.y; ey0[i+1].z=Y.z;
				ez0[i+1].x=Z.x; ez0[i+1].y=Z.y; ez0[i+1].z=Z.z;
				// o0[i+1]+=d(i)*Z
				o0[i+1].x+=d(i)*Z.x; o0[i+1].y+=d(i)*Z.y; o0[i+1].z+=d(i)*Z.z;
			}
			else{
				// ���̖ʂ̕ΐS�O���W
				// o0[i+1]=o[i]+d(i)*ez[i]
				o0[i+1].x=o[i].x+d(i)*ez[i].x;
				o0[i+1].y=o[i].y+d(i)*ez[i].y;
				o0[i+1].z=o[i].z+d(i)*ez[i].z;
				// ex0[i+1]=ex[i] ey0[i+1]=ey[i] ez0[i+1]=ez[i]
				ex0[i+1].x=ex[i].x; ex0[i+1].y=ex[i].y; ex0[i+1].z=ex[i].z;
				ey0[i+1].x=ey[i].x; ey0[i+1].y=ey[i].y; ey0[i+1].z=ey[i].z;
				ez0[i+1].x=ez[i].x; ez0[i+1].y=ez[i].y; ez0[i+1].z=ez[i].z;
			}
		}
	}
}

inline void cLens1::transform(double& x,double& y,double& z,
					   int i0,int i0_pre,int i,int i_pre,int translate) const{
	// ��i0�ʂ̃��[�J�����W����C��i�ʂ̃��[�J�����W�ɕϊ�����D
	// i0_pre,i_pre : ��0�̂Ƃ��͕ΐS�O�̍��W�C0�̂Ƃ��͖ʍ��W

	// �������̂��߃x�N�g�����Z�͐������ɏ����Ă���D

	// inline�w��͍������̌��ʂ��� (2013.01.04)

	double px,py,pz;

	if(i0_pre){
		// p=x*ex0[i0]+y*ey0[i0+z*ez0[i0]
		px=x*ex0[i0].x+y*ey0[i0].x+z*ez0[i0].x;
		py=x*ex0[i0].y+y*ey0[i0].y+z*ez0[i0].y;
		pz=x*ex0[i0].z+y*ey0[i0].z+z*ez0[i0].z;
	}
	else{
		// p=x*ex[i0]+y*ey[i0]+z*ez[i0]
		px=x*ex[i0].x+y*ey[i0].x+z*ez[i0].x;
		py=x*ex[i0].y+y*ey[i0].y+z*ez[i0].y;
		pz=x*ex[i0].z+y*ey[i0].z+z*ez[i0].z;
	}

	if(translate!=0){
		// p-=o[i]-o[i0]
		if(i0_pre){
			px+=o0[i0].x;
			py+=o0[i0].y;
			pz+=o0[i0].z;
		}
		else{
			px+=o[i0].x;
			py+=o[i0].y;
			pz+=o[i0].z;
		}

		if(i_pre){
			px-=o0[i].x;
			py-=o0[i].y;
			pz-=o0[i].z;
		}
		else{
			px-=o[i].x;
			py-=o[i].y;
			pz-=o[i].z;
		}
	}

	if(i_pre){
		// x=sProduct(p,ex0[i]) y=sProcuct(p,ey0[i]) z=sProduct(p,ez0[i])
		x=px*ex0[i].x+py*ex0[i].y+pz*ex0[i].z;
		y=px*ey0[i].x+py*ey0[i].y+pz*ey0[i].z;
		z=px*ez0[i].x+py*ez0[i].y+pz*ez0[i].z;
	}
	else{
		// x=sProduct(p,ex[i]) y=sProduct(p,ey[i]) z=sProduct(p,ez[i])
		x=px*ex[i].x+py*ex[i].y+pz*ex[i].z;
		y=px*ey[i].x+py*ey[i].y+pz*ey[i].z;
		z=px*ez[i].x+py*ez[i].y+pz*ez[i].z;
	}
}

inline void cLens1::transform(vector<double>& v,int i0,int i0_pre,int i,int i_pre,int translate) const{
	transform(v.x,v.y,v.z,i0,i0_pre,i,i_pre,translate);
}

inline void cLens1::transform(complex& x,complex& y,complex& z,
					   int i0,int i0_pre,int i,int i_pre,int translate) const{
	vector<double> v1(Re(x),Re(y),Re(z)), v2(Im(x),Im(y),Im(z));
	transform(v1.x,v1.y,v1.z,i0,i0_pre,i,i_pre,translate);
	transform(v2.x,v2.y,v2.z,i0,i0_pre,i,i_pre,translate);
	x=complex(v1.x,v2.x);
	y=complex(v1.y,v2.y);
	z=complex(v1.z,v2.z);
}

void cLens1::transform_plane(double &a,double &b,double &c,double &d,int i0,int i0_pre,int i,int i_pre){
	// ����
	//   ax+by+cz+d=0
	// ��a,b,c,d���i0�ʂ̃��[�J�����W�����i�ʂ̃��[�J�����W�ł̕\���ɕϊ�����D
	// ���炩����make_coordinate()�����s���Ă������ƁD
	vector<double> n,r;

	// �@���x�N�g���̕ϊ�
	n=vector<double>(a,b,c);              // ax+by+cz+d=0�̖@���x�N�g����(a,b,c)
	transform(n,i0,i0_pre,i,i_pre,0);     // �V���W�֕ϊ�

	// �ʉߓ_�̕ϊ�
	if(a!=0){
		r=vector<double>(-d/a,0,0);
	}
	else if(b!=0){
		r=vector<double>(0,-d/b,0);
	}
	else if(c!=0){
		r=vector<double>(0,0,-d/c);
	}
	else{
		return;  // a=b=c=0 �͕��ʂ̕������ł͂Ȃ�
	}
	transform(r,i0,i0_pre,i,i_pre,1);     // �V���W�֕ϊ�
	
	// �V�������� nx(x-rx)+ny(y-ry)+nz(z-rz)=0
	a=n.x;
	b=n.y;
	c=n.z;
	d=-n.x*r.x-n.y*r.y-n.z*r.z;
}

void cLens1::image_plane(double& a,double &b,double &c,double &d,int i){
	// �����ʂ̑�i�ʂɂ�鑜�����߂�D
	// ���W�͑�i�ʂ̖ʍ��W�ŕ\������D
	// ���z����(���Ȃ킿�V���C���v���[�N����)�ɂ��D

	double s,h,u,u1;

	// ���̕��� ax+by+cz+d=0�C������ a'x+b'y+c'z+d=0 �� z=0(�啽��)�̌���͂��ꂼ��C
	//    ax+by+d=0
	//    a'x+b'y+d'=0
	// �ł���(a=b=0�͌����ɐ����Ȗʂ�\��)�C����瓯�꒼���ł��邩��C
	//    a'=a, b'=b, d'=d
	// �Ƃ���΂悢�D����������a,b,d�͕s�ςƂ���i�{���͈����ɓ����K�v�͂Ȃ��j�D
	// ���Ƃ͎啽�ʂƌ���(z��)�̌�_
	//    (0,0,s) = (0,0,-d/c)
	// ��肻�̑��_(0,0,s')�������֌W������Ƃ߁C
	//    c' = -d/s'
	// �Ƃ���΂悢�D
	// �������Cd=0�̂Ƃ��͕����ʂ����_��ʂ邩��
	// �����ʂ�z=0�̌����ɕ����ʂ�z���Ƃ̌�_������C
	// ���̕��@�ł͑����ʂ����߂邱�Ƃ��ł��Ȃ��D
	// ���̂Ƃ������������}���ł���Ε����ʂƑ����ʂ͓����ł���D
	// �}�����Ⴄ�Ƃ���
	//    c' = c*N(i-1,1)/N(i,1)
	// �ƂȂ�D����͗Ⴆ��y-z�ʂł� 
	//    y=-(c/b)z
	//    y=-(c'/b')z=-(c'/b)z
	// �̒����ƂȂ肱���͋��ܑO��̌����ƍl�����C
	// �ߎ����_��� N'c' = Nc ������D

	if(d!=0){
		s= c==0 ? LN : -d/c;
		if( fabs(s)>=LN ){
			u=0; h=1;
		}
		else{
			u=1/s; h=1;
		}
		u1=(N(i-1,1)*u+h*power(i,i,1))/N(i,1);
		c=-d*u1/h;
	}
	else{
		c=c*N(i-1,1)/N(i,1);
	}
}

void cLens1::transform_line(double &x0,double &y0,double &z0,double &l,double &m,double &n,
							int i0,int i0_pre,int i,int i_pre)
{
	// ����
	//   (x-x0)/l = (y-y0)/m = (z-z0)/n
	// ��x0,y0,z0,l,m,n���i0�ʂ̃��[�J�����W�����i�ʂ̃��[�J�����W�ł̕\���ɕϊ�����D
	// ���炩����make_coordinate()�����s���Ă������ƁD
	
	transform(x0,y0,z0,i0,i0_pre,i,i_pre,1);
	transform(l,m,n,i0,i0_pre,i,i_pre,0);
}

void cLens1::image_line(double &x0,double &y0,double &z0,double &l,double &m,double &n,int i){
	// ���� (x-x0)/l=(y-y0)/m=(z-z0)/n �̑�i�ʂɂ�鑜�����߂�D
	// ���W�͑�i�ʂ̖ʍ��W�ŕ\������D
	// ���z�����ɂ��D
	double x,y,xf,yf,o,ff;

	if(n==0) return;     // n==0�͎啽�ʂƕ��s�Ȓ���
	ff=this->ff(i,i,1);  // �p���[��0�̂Ƃ�ff()��0��Ԃ�

	// �啽��z=0�Ƃ̌�_	
	x=-z0*l/n+x0;
	y=-z0*m/n+y0;

	// �����œ_�ʂƂ̌�_
	xf=(ff-z0)*l/n+x0;
	yf=(ff-z0)*m/n+y0;
	
	// �V�����ʉߓ_(x0,y0,z0)�ƕ����x�N�g��(l,m,n)
	// ���ɕK�v�͂Ȃ����Cf��������n���傫�����ɂȂ�̂ŋK�i�����Ă����D
	x0=x;
	y0=y;
	z0=0;
	if(ff!=0){
		l=-xf;
		m=-yf;
		n=-ff*N(i,1)/N(i-1,1);
		o=sqrt(l*l+m*m+n*n);
		l/=o;
		m/=o;
		n/=o;
	}
	else{
		n*=N(i,1)/N(i-1,1);
	}
}


int cLens1::raytrace(int i,int j, double xp,double yp,double zp,double X,double Y,double Z,
                    double& x,double& y,double& z,double & X1,double& Y1,double& Z1, 
                    int EA_enable,int mask_enable,
                    int as_trace, 
					vector<double>& ha,vector<double>& dQa,vector<double>& hb,vector<double>& dQb,
                    int E_trace, vector<complex>& E) 
{
	// i-1�ʋ��܌�̏���(i-1�ʒ��_���W�)��i�ʋ��܌�(i�ʒ��_���W�)�ɕϊ�����D
	//  ���ˌ������� (X,Y,Z) �͒P�ʃx�N�g���ł��邱�ƁD�{�֐����ŋK�i���͂��Ă��Ȃ��D

	double val, c,p,Mz,MM,q,l,m,n,o,xi,xi1,G, N,N1,No,No1,H, EAx,EAy, EAdx,EAdy, D;
	int asph_type;
	
	//// �������̂��߂̏����D���w�n�f�[�^�����[�J���ϐ��ɕۑ�(�����o�֐��Ăяo���̍팸)�ȂǁD//////
	
	c  =this->c(i);
	N  =this->N(i-1,j);
	N1 =this->N(i,j);
	No =this->N(i-1,1);
	No1=this->N(i,1);
	EAx=this->EAx(i);
	EAy=this->EAy(i);
	EAdx=this->EAdx(i);
	EAdy=this->EAdy(i);
	asph_type=this->asph_type(i);

	H=fabs(N/N1);

	/////////////////////////////////////////////////////////////

	if(IsGrin(i-1)){
		cSelfoc s;
		double dz,dz0;
		int count;
		vector<double> r(xp,yp,zp),Q(X,Y,Z), r0,Q0;

		s=cSelfoc(N,rA(i-1,j),GrinPhi(i-1));

		// i-1�ʂ̖ʍ��W����i�ʂ̕ΐS�O���W(GRIN�}���̍��W)��
		transform(r,i-1,0,i,1,1);
		transform(Q,i-1,0,i,1,0);

		count=0;
		do{
			r0=r;
			Q0=Q;
			
			if( s.trace(r,Q,GrinDs)==0 ) return GRIN_OUT+i;

			count++;
			if(count%GrinToListStep==0){
				GrinRay_i.AddTail(i-1);
				GrinRay.AddTail(r);  // i��(�}���̏I�[)�̕ΐS�O���W
			}

			// i�ʂ̕ΐS�O���W����i�ʂ̖ʍ��W��
			transform(r,i,1,i,0,1);
			transform(r0,i,1,i,0,1);

			dz0=r0.z-surface_sag(i,r0.y,r0.x,1);
			dz=r.z-surface_sag(i,r.y,r.x,1);

			// i�ʂ̖ʍ��W����i�ʂ̕ΐS�O���W�ɖ߂�
			transform(r,i,0,i,1,1);
			transform(r0,i,0,i,1,1);

		} while( dz0*dz>0 );

		// i�ʂ̕ΐS�O���W����i-1�ʂ̖ʍ��W��
		transform(r0,i,1,i-1,0,1);
		transform(Q0,i,1,i-1,0,0);
		xp=r0.x; yp=r0.y; zp=r0.z;
		X=Q.x; Y=Q.y; Z=Q.z;
	}

	// i-1�ʂ̖ʍ��W����i�ʂ̖ʍ��W�֕ϊ�����
	transform(xp,yp,zp,i-1,0,i,0,1);
	transform(X,Y,Z,i-1,0,i,0,0);
	if(E_trace){
		transform(E.x,E.y,E.z,i-1,0,i,0,0);
	}
	if(as_trace){
		transform(ha.x,ha.y,ha.z,i-1,0,i,0,0);
		transform(dQa.x,dQa.y,dQa.z,i-1,0,i,0,0);
		transform(hb.x,hb.y,hb.z,i-1,0,i,0,0);
		transform(dQb.x,dQb.y,dQb.z,i-1,0,i,0,0);
	}

	{
		double cc;

		p=-(zp*Z+yp*Y+xp*X);
		Mz=zp+p*Z;
		MM=zp*zp+yp*yp+xp*xp-p*p;
		
		if(Z==0){
			// �Ⴆ�Ε��̖ʂɕΐS������ꍇ�Cs=t�łȂ��Ă��C
			// (s=t��Raytrace()�Ŕr�������)Z=0���N���肤��D
			// Z=0�ł͋��ʂƂ̌�_�͈�ӂɂ͌��܂�Ȃ�(�g�����Y�݌v�@�h�}3.2)
			// �̂ŁC�ǐՒ�~�Ƃ���D
			// (09.10.07)
			return NOT_INTERSECT+i;
		}
		else{
			if(asph_type==CONIC || asph_type==TOROID || asph_type==ANAMO || asph_type==CONIC_OA){
				int counter;

				switch(asph_type){
				case CONIC:
					cc=c*0.5;   // cc=c ����J�n����� ������(fish(i)=1) �Ō��ʂ� fish(i)=0 �Ɠ����ɂȂ�ꍇ�������i2018.02.16)
				                // ��1��ڂ̔񋅖ʂƂ̐ڐ��ƌ����Ƃ̌�_���C�����̔��Α��ɂȂ��Ă��܂����߂Ǝv����D
					
					if(EAy!=0){
						// cc��L���͈͂̒[��ʂ鋅�̋ȗ����a�Ƃ���
						// ���R�F
						//   �Ⴆ�΁C�����̍����܂ނȂǂɂ��L���͈͊O�ŃT�O���}���ɕω�����Ƃ��C
						//   �L���͈͊O���l������ƁC�ʂƌ����Ƃ̌�_��2�����ȏ゠��ꍇ������D
						//   ���̂Ƃ��Ccc�̋��ʂ��ʂ��痣��Ă����Ă���ƁC�Ӑ}���Ȃ����̌�_��I�����Ă��܂����Ƃ�����D
						//   (�X�}�z�����Y�̃J�����^�����Y�Ȃǁj
						double h=EAy/2;
						double z=surface_sag(i,h,0,0);

						if(z*z<h*h){  // �����𒴂��Ȃ����߂̏���
							cc=2*z/(z*z+h*h);							
						}
					}

					break;
				case TOROID:
				case ANAMO:
				case CONIC_OA:
					double cx,cxy,cy;
					ParaxialExpansion(i,cx,cxy,cy);
					cc=(cx+cy)*0.25;  // =(cx+cy)/2 *0.5;
					break;
				}
				
				for(counter=1; counter<=5; counter++){
					// ����ʂƌ����܂ŁC�ȗ����ɂ߂J��Ԃ��D
					// �ȑO(2018.02.15�ȑO)�� cc�����ߑł����Ă������C
					// �S�Ă̏ꍇ�ŋȗ��t�̋��ʂƌ����������Ƃ͌���Ȃ������D
					val=1-(cc/(Z*Z))*(MM*cc-2*Mz);
					if(val>0) break;
					cc*=0.5;
				}
			}
			else{
				cc=c;
				val=1-(cc/(Z*Z))*(MM*cc-2*Mz);
			}
		}

		if( val< 0 ){
			if( ExcludeVirtualObject!=0 && this->optical_path[i-1]==0 ){
				// ���̏ꍇ�͊���ʂƌ����Ȃ��Ă����s���ׂ��ł���D
				// �Ƃ肠����q=-1�Ƃ��đ��s����΁C��ɋ������Ɣ��肳���D
				q=-1;
				z=zp+q*Z; y=yp+q*Y; x=xp+q*X;
			}
			else{
				return NOT_INTERSECT+i;
			}
		}
		else{
			if(Fish(i)){
				if(cc==0){
					q=p-Mz/Z;
				}
				else{
					q=p+(MM*cc-2*Mz)/( Z*(1+sgn(N*Z)*sqrt(val)) );
					// �����͊���ʂ�2�����Ō����D
					// �����t�ɐi��(N*Z<0)�̂Ƃ���XY���ʂ��牓�������̗p����
					// (���჌���Y�̐擪�ʂȂǁj
				}
			}
			else{
				// �ʏ�̖� : �ǂ̂悤�ȏꍇ�ł�XY���ʂɋ߂������̗p����i"�����Y�݌v�@(����)" �̎�(3.9))
				q=p+(MM*cc-2*Mz)/( Z*(1+sqrt(val)) );
			}
			
			z=zp+q*Z; y=yp+q*Y; x=xp+q*X;
			n=1-z*cc; m=-y*cc; l=-x*cc;
		}
	}

	if( asph_type || (Newton(i)!=0 || As0(i)!=0 || As45(i)!=0) ) {
		// ���ˌ�����̓_(x,y,z)���ʂƈ�v����悤�ɁC�����ߎ����s�Ȃ��D
		double z01, we, x1,y1,z1, x0,y0,z0;
		int domain_error,count;
		const int MAXTIMES=10;
		vector<double> E;

		if(X*X+Y*Y>0){   // �����������ɕ��s�łȂ��Ƃ�
			// �قڊԈႢ�Ȃ�x,y���ʂ̒�`����ł��������̓_(x0,y0,z0)�����߂�D
			// ���ˌ�����XY�ʂւ̓��e�ւ̌��_����̐����̑���(x0,y0,0)�Ƃ��C
			// (x0,y0,0)�ɑΉ����������_��(x0,y0,z0)�Ƃ���D
			// (x,y)���ŏ������`����ł���Ζ��ʂȏ����ł���D
			x0=x-(x*X+y*Y)/(X*X+Y*Y)*X;
			y0=y-(x*X+y*Y)/(X*X+Y*Y)*Y;
			z0= Y!=0 ? z-(y-y0)*Z/Y : z-(x-x0)*Z/X;   // (x,y,z) : (x,y)����`����ɂ��������̓_
			surface_sag(i,y0,x0,1,domain_error);
			if(domain_error) return NOT_INTERSECT+i;
		}
		else{   // �����������ɕ��s�ȂƂ�
			surface_sag(i,y,x,1,domain_error);
			if(domain_error) return NOT_INTERSECT+i;
		}

		for(count=1; count<=MAXTIMES; ++count){ // �������[�v�Œ�~����̂�h��
			do{
				z01=surface_sag(i,y,x,1,domain_error);  // (x,y,z01) : �ʏ�̓_�i������ł͂Ȃ��j
				if(domain_error!=0){
					// (x,y)���ʂ̒�`��O�̂Ƃ��́C�O���`����ł������_�Ƃ̒��_�܂Ŗ߂��D
					y=y0+(y-y0)/2;
					x=x0+(x-x0)/2;
					z=z0+(z-z0)/2;
				}
			} while(domain_error!=0);

			E=surface_normal(i,y,x,1); l=E.x; m=E.y; n=E.z;
			we=n*(z01-z)/(X*l+Y*m+Z*n);
			x1=x+we*X; y1=y+we*Y; z1=z+we*Z;   // (x1,y1,z1) : (x,y,z01)�ł̐ڕ��ʂƌ����̌�_
			x0=x; y0=y; z0=z;
			x=x1; y=y1; z=z1;
			if(fabs(we)<1e-6) break;  // �g�����Y�݌v�@�i�����j�h�ł� |z1-z|<�� �Ŏ����𔻒肵�Ă��邪�C
			                          //  Z��0�i�������ق�XY�ʓ��j�̂Ƃ��덷��������D
			if(count==MAXTIMES) return NOT_INTERSECT+i;
		}
	}

	if( ( D=(x-xp)*X+(y-yp)*Y+(z-zp)*Z ) <0 ){  // �������̏ꍇ

		if ( ExcludeVirtualObject!=0 && this->optical_path[i-1]==0 ){
			if(i<k){
				// Schimpflug�J������팟�ᑤ����ǐՂ���Ƃ��C
				// ���_�̊�̉��s�������̈ʒu�ɂ���ĕ���Ԃ��قȂ�D
				// ����ɑΉ����邽�߁CExcludeVirtualObject=true �̂Ƃ��́C
				// ����Ԃ̌������������ɂȂ�܂ŁC
				// ���_�̋�Ԉʒu���Œ肵���܂܎��̋�Ԃֈڂ��D (09.10.01)
				//
				// �������Ci==k�ł����̏������s�Ȃ��ƁC�ǂ�Ȍ����ł��n��
				// �ʉ߂��邱�ƂɂȂ��Ă��܂��CFindMarginalRay()��
				// �������[�v����������D (09.10.06)
				x=xp; y=yp; z=zp;
				X1=X; Y1=Y; Z1=Z;
				xi=xi1=0;
				this->x[i]=x; this->y[i]=y; this->z[i]=z;
				this->X[i]=X; this->Y[i]=Y; this->Z[i]=Z;
				this->X1[i]=X1; this->Y1[i]=Y1; this->Z1[i]=Z1;
				this->xi[i]=xi; this->xi1[i]=xi1;
				return 0;
			}
			else{
				return NOT_INTERSECT;
			}
		}

		if ( is_acting_surf(i-1) && is_acting_surf(i) && ExcludeVirtualRay!=0 ){
			return NOT_INTERSECT+i;
		}
		
		// ������(Q�Ɣ��Ε����Ɍ���������)�����e����ƁC
		// �H�ɁC�����ʂ�Ȃ��͂��̏�ʂł��܂��܌������i�s�t�����Ŏ��̖ʂƌ�����
		// �ǐՂ��p������Ӑ}���Ȃ����ʂƂȂ�D 
		// �������S�Ĕr�����Ă��܂��ƁC�������̈ʒu�̉ˋ�ʂ��g���Ȃ��Ȃ�D
		// ���S�ɈӐ}�ǂ���ɂ���͓̂�����C���Ȃ��Ƃ���i,i-1�ʋ���
		// ���܁C���ˍ�p������Ƃ��͔r������D
		// �������ʒu�̍i��ł͏�L�s�s�����N����\�����c���Ă���D(09.08.27)

		// ��L�ł͐݌v�r���̌���0�̃����Y�ŃR�o�����ƂȂ�Ƃ������ǐՂ��o���Ȃ��D
		// �������̔r�����K�v�Ȃ̂͋H�ł��邩��C�V���� AllowVirtualRay ��݂��āC
		// �����l��1�Ƃ��C�s�x�ݒ肷�邱�Ƃɂ����D(09.09.16)
	
	}

	this->x[i]=x; this->y[i]=y; this->z[i]=z;
	this->X[i]=X; this->Y[i]=Y; this->Z[i]=Z;

	if(stop>0 && StopDominate!=0 && i!=stop){
		// StopDominate��true�̂Ƃ��Cstop�ʈȊO�ł�EA�𖳎�����D
		// ����́C�������ׂ��C���L���a�ɗ]�T�̂Ȃ��n�i�L�p�ڊ�Ȃǁj�̎����݌v�ŁC
		// ���O����������Ƃ���FindAThuruRay()�Ɏ��Ԃ��Ƃ��Ȃ����߂ɗL���ł���D
	}
	else if(i<k+1 && EA_enable){    // i==k+1(����)�ł�EA�ɂ�鐧���͂��Ȃ��D
		switch( EAtype(i) ){
		case 0:
			if( EAy>0 && EAx==0 ){
				if( (y-EAdy)*(y-EAdy)+(x-EAdx)*(x-EAdx)>EAy*EAy/4 ) return CLIPPED+i;
			}
			else if( EAy>0 && EAx>0 ){
				if( (y-EAdy)*(y-EAdy)/(EAy*EAy/4)+(x-EAdx)*(x-EAdx)/(EAx*EAx/4)>1 ) return CLIPPED+i;
			}
			else if( EAy==0 && EAx==0 ){
				return CLIPPED+i;
			}
			break;
		case 1:
			if( EAy>0 && EAx==0 ){
				if( (y-EAdy)*(y-EAdy)>EAy*EAy/4 || (x-EAdx)*(x-EAdx)>EAy*EAy/4 ) return CLIPPED+i;
			}
			else if( EAy>0 && EAx>0 ){
				if( (y-EAdy)*(y-EAdy)>EAy*EAy/4 || (x-EAdx)*(x-EAdx)>EAx*EAx/4 ) return CLIPPED+i;
			}
			else if( EAy==0 && EAx==0 ){
				return CLIPPED+i;
			}

			if(CHMx(i)>0 && CHMy(i)>0){ // �ʎ�肪����ꍇ
				double a,b;

				a=-CHMy(i)/CHMx(i);
				b=-a*( fabs((EAx==0 ? EAy:EAx)/2)-CHMx(i) ) +fabs(EAy/2); // ��1�ی��̖ʎ�΂ߐ��̕������� y=ax+b
				if((y-EAdy)> a*(x-EAdx)+b) return CLIPPED+i;  // ��1�ی�
				if((y-EAdy)>-a*(x-EAdx)+b) return CLIPPED+i;  // ��2�ی�
				if((y-EAdy)< a*(x-EAdx)-b) return CLIPPED+i;  // ��3�ی�
				if((y-EAdy)<-a*(x-EAdx)-b) return CLIPPED+i;  // ��4�ی�
			}
			break;
		}
	}
	
	if(EAy<0){  // ��������}��
		if(i<k+1 && mask_enable){
			switch( EAtype(i) ){
			case 0:
				if( EAy<0 && EAx==0 ){
					if( (y-EAdy)*(y-EAdy)+(x-EAdx)*(x-EAdx)<EAy*EAy/4 ) return CLIPPED+i;
				}
				else if( EAy<0 && EAx<0 ){
					if( (y-EAdy)*(y-EAdy)/(EAy*EAy/4)+(x-EAdx)*(x-EAdx)/(EAx*EAx/4)<1 ) return CLIPPED+i;
				}
				break;
			case 1:
				if( EAy<0 && EAx==0 ){
					if( (y-EAdy)*(y-EAdy)<EAy*EAy/4 && (x-EAdx)*(x-EAdx)<EAy*EAy/4 ) return CLIPPED+i;
				}
				else if( EAy<0 && EAx<0 ){
					if( (y-EAdy)*(y-EAdy)<EAy*EAy/4 && (x-EAdx)*(x-EAdx)<EAx*EAx/4 ) return CLIPPED+i;
				}

				if(CHMx(i)>0 && CHMy(i)>0){ // �ʎ�肪����ꍇ
					double a,b;

					a=-CHMy(i)/CHMx(i);
					b=-a*( fabs((EAx==0 ? EAy:EAx)/2)-CHMx(i) ) +fabs(EAy/2); // ��1�ی��̖ʎ�΂ߐ��̕������� y=ax+b
					if((x-EAdx)>=0 && (y-EAdy)>=0 && (y-EAdy)< a*(x-EAdx)+b) return CLIPPED+i;  // ��1�ی�
					if((x-EAdx)<=0 && (y-EAdy)>=0 && (y-EAdy)<-a*(x-EAdx)+b) return CLIPPED+i;  // ��2�ی�
					if((x-EAdx)<=0 && (y-EAdy)<=0 && (y-EAdy)> a*(x-EAdx)-b) return CLIPPED+i;  // ��3�ی�
					if((x-EAdx)>=0 && (y-EAdy)<=0 && (y-EAdy)>-a*(x-EAdx)-b) return CLIPPED+i;  // ��4�ی�
				}
				break;
			}
		}
	}

	// i-1�ʌ�_����i�ʌ�_�܂ł̌��H����������D
	optical_path[i]=optical_path[i-1]+fabs(N)*D;  // D = X*(x-xp)+Y*(y-yp)+Z*(z-zp)
	//	optical_path[i]=optical_path[i-1]+fabs(N[i-1][j])*sqrt( (x-xp)*(x-xp)+(y-yp)*(y-yp)+(z-zp)*(z-zp) ) ���ƁC
	//	�ˋ�ʂȂǂɂ��(x-xp,y-yp,z-zp)��(X,Y,Z)�Ɣ��Ό����ɂȂ镔��(���łȂ�����)������ƕs����N����D
	
	// ���ܓ_�ɂ�����GRIN�}���̋��ܗ���ݒ肷��
	if(IsGrin(i-1)) N =cSelfoc::N(N,rA(i-1,j),x,y);
	if(IsGrin(i)  ) N1=cSelfoc::N(N1,rA(i,j),x,y);

	o=1/sqrt(l*l+m*m+n*n); l*=o; m*=o; n*=o;
	
	switch(asph_type){
	case IDEAL: 
		{
			double _x,_y;
			
			xi=X*l+Y*m+Z*n;
			if( ( val=1-H*H*(1-xi*xi) ) < 0 ) return TOTAL_REF+i;
			xi1=sgn(N*N1*xi)*sqrt(val);
			G=xi1-H*xi;
			X1=H*X+G*l;
			Y1=H*Y+G*m;
			Z1=H*Z+G*n;
			switch( cylinder(i) ) {
				// �g�􉽌��w(�O��)�h���C
				//      y' = -f(y/z) (=-f*tan��) (3.6)
				//      n'/f' = -n/f (=fideal)   (3.17)
				// ���������āC
				//      f' = fideal*n'
				//      y' = fideal*n*tan��
			case GENY:
				_x= fideal(i)==0 ? x+X1/Z1 : fideal(i)*No*X/Z;
				_y= fideal(i)==0 ? x+X1/Z1 : y+fideal(i)*No1*Y1/Z1;
				break;
			case GENX:
				_x= fideal(i)==0 ? y+Y1/Z1 : x+fideal(i)*No1*X1/Z1;
				_y= fideal(i)==0 ? y+Y1/Z1 : fideal(i)*No*Y/Z;
				break;
			default:
				_x= fideal(i)==0 ? x+X1/Z1 : fideal(i)*No*X/Z;
				_y= fideal(i)==0 ? y+Y1/Z1 : fideal(i)*No*Y/Z;
				break;
			}
			X1= _x-x;
			Y1= _y-y;
			Z1= fideal(i)==0 ? 1 : No1*fideal(i);
			o=1/sqrt(X1*X1+Y1*Y1+Z1*Z1);
			X1*=o; Y1*=o; Z1*=o;
		}
		break;
	default:
		xi=X*l+Y*m+Z*n;
		if(grating(i)){
			// ���܂̖@������܊i�q�֊g�����g�����Y�݌v�@�h�̏������C��
			// (3.5)  -> |N'|(ExQ')=|N|(ExQ)+(m��/d)(Exg)
			// (3.21) -> ��'=( (NN'��)/|NN'��| ) * sqrt(1-|ExQ'|^2)
			// (3.23) -> X' += (1/|N'|)(m��/d)(gx-(E,g)X)
			//           Y' += (1/|N'|)(m��/d)(gy-(E,g)Y)
			//           Z' += (1/|N'|)(m��/d)(gz-(E,g)Z)
			//
			//  �ɁF�^�󒆔g��
			//  m : ��܎���
			//  d : �i�q�Ԋu
			//  g=(gx,gy,gz) : �i�q�����̒P�ʃx�N�g���D�i�q��g��@���Ƃ���Ԋud�̕��ʌQ��
			//                 ���ܖʂƂ̌���ł���D
			vector<double> C(l,m,n);
			vector<double> Q(X,Y,Z);
			vector<double> g(grx(i),gry(i),grz(i)); g=g/abs(g);
			vector<double> v;
			double a;

			a=(double)difforder(i)*(wl[j]*1e-6)/gpitch(i);
			v=(1/fabs(N1))*( fabs(N)*vProduct(C,Q)+a*vProduct(C,g) ); 
			// (��)�㎮�͔g�����܂ނ��ߒP�ʌn�Ɉˑ�
			if( ( val=1-abs(v)*abs(v) ) < 0 ) return TOTAL_REF+i;
			xi1=sgn(N*N1*xi)*sqrt(val);
			G=xi1-H*xi;
			val=sProduct(C,g);
			X1=H*X+G*l +a/fabs(N1)*(g.x-val*Q.x);
			Y1=H*Y+G*m +a/fabs(N1)*(g.y-val*Q.y);
			Z1=H*Z+G*n +a/fabs(N1)*(g.z-val*Q.z);
		}
		else if(Diffusion(i)!=0){     // �g�U��
			// Diffusion(i)>0 �̏ꍇ�͊��S�g�U�ł���C�g����p�̍ő��180���ł���D
			// Diffusion(i)<0 �̏ꍇ�́C
			// ���S�g�U���ł��邪�C|Diffusion(i)| ���L����p(�S�p)�̌��E�Ƃ���i�p�x0�͖ʖ@���j�D			
			// ���w�n�ɉ����� |Diffusion(i)|<0 ���w�肷�邱�Ƃɂ��v�Z���ԒZ�k��}�邱�Ƃ��ł���D
			vector<double> Ex,Ey,Ez;
			double th,phi, x,y,z, semi_angle;

			semi_angle= Diffusion(i)<0 ? -Diffusion(i)/2 : 90;
			
			xi=sgn(X*l+Y*m+Z*n);  // ���ˌ����x�N�g���Ɩʖ@���x�N�g�������������D�����Y�݌v�@(3.20)�ɑΉ��D
			xi1=sgn(N*N1*xi);     // �o�ˌ������ʖ@���x�N�g���Ɠ��������D�����Y�݌v�@(3.21)�ɑΉ��D
			X1=xi1*l;
			Y1=xi1*m;
			Z1=xi1*n;
			Ez=vector<double>(X1,Y1,Z1);  // ���̐i�s�������������ʖ@���P�ʃx�N�g�������߂�ꂽ

			// Ez�ɐ����ȂP�̒P�ʃx�N�g��Ex�𐶐�����
			// (�O��QxEz���g����Q//Ez�̂Ƃ��ɑΉ��ł��Ȃ�)
			if( Ez.x*Ez.x<=Ez.y*Ez.y && Ez.x*Ez.x<=Ez.z*Ez.z ){ // x������Βl���ŏ��̂Ƃ�
				Ex=vector<double>(0, Ez.z, -Ez.y);
				Ex=Ex/abs(Ex);
			}
			else if( Ez.y*Ez.y<=Ez.x*Ez.x && Ez.y*Ez.y<=Ez.z*Ez.z ){
				Ex=vector<double>(Ez.z,0,-Ez.x);
				Ex=Ex/abs(Ex);
			}
			else{
				Ex=vector<double>(Ez.y,-Ez.x,0);
				Ex=Ex/abs(Ex);
			}
			// Ex,Ez�ɐ����ȃx�N�g��Ey�D�����z�����ʖ@�����������������[�J�����W�n(Ex,Ey,Ez)���ł����D
			Ey=vProduct(Ex,Ez);
			
			// �����ŏo�ˌ��������̋ɍ��Wth,phi�𐶐�����D
			//     x=sin(th)cos(phi)
			//     y=sin(th)sin(phi)
			//     z=cos(th)
			// th�̕��z��cos(th)�ɔ�Ⴗ��(�����x���g��)�Dphi�Ɋւ��Ă͈�l�D
			// (3�����̏ꍇ�C��l�ȕ��z��sin(th)�ɔ�Ⴗ�邱�Ƃɒ���)
			do{
				th=Random(0,90,0);
			} while( Random(0,1,0)>sin(th*PI/180)*cos(th*PI/180) || th>semi_angle );
			phi=Random(0,360,0);

			x=sin(th*PI/180)*cos(phi*PI/180);
			y=sin(th*PI/180)*sin(phi*PI/180);
			z=cos(th*PI/180);

			X1=x*Ex.x+y*Ey.x+z*Ez.x;
			Y1=x*Ex.y+y*Ey.y+z*Ez.y;
			Z1=x*Ex.z+y*Ey.z+z*Ez.z;
			o=1/sqrt(X1*X1+Y1*Y1+Z1*Z1);
			X1*=o; Y1*=o; Z1*=o;
		}
		else{
			if( ( val=1-H*H*(1-xi*xi) ) < 0 ) return TOTAL_REF+i;
			xi1=sgn(N*N1*xi)*sqrt(val);
			G=xi1-H*xi;
			X1=H*X+G*l;
			Y1=H*Y+G*m;
			Z1=H*Z+G*n;	
		}
		break;
	}

	if(as_trace) {
		switch(asph_type) {
		case IDEAL:
			{
				// h��Q
				// h'y=hy, h'x=hx �Ƃ��邱�Ƃ��ł��Ch'��Q' ��� h'z�����܂�D
				// t,t'�͌����ɉ������אڌ�����y�����������ɂȂ�ʒu�܂ł̋����Ƃ���ƁC
				// [h + t (Q +dQ )]y = [t Q ]y
				// [h'+ t'(Q'+dQ')]y = [t'Q']y
				// �܂��Cz���ɉ�����������s�Ƃ���ƁC
				// s = tQz
				// s'= t'Qz'
				// N'/s' = N/s + 1/f
				// �ȏ���C
				// dQ'y�����܂�DdQ'x�����l�D
				// dQ'��Q'���dQ'z�����܂�D
				double tau, py,px;
				vector<double> Q(X,Y,Z), Q1(X1,Y1,Z1);

				switch(cylinder(i)) {
				case GENY:
					py=0; 
					px= fideal(i)==0 ? 0 : 1/fideal(i);
					break;
				case GENX:
					py= fideal(i)==0 ? 0 : 1/fideal(i); 
					px=0;
					break;
				default:
					py=px= fideal(i)==0 ? 0 : 1/fideal(i);
					break;
				}
				tau=(z-zp)/Z;

				ha=ha+tau*dQa;
				hb=hb+tau*dQb;
				this->ha[i]=ha; this->dQa[i]=dQa; this->hb[i]=hb; this->dQb[i]=dQb;
				
				ha.z=-(ha.x*Q1.x+ha.y*Q1.y)/Q1.z;
				dQa.y=(N/N1)*(Q1.z/Q.z)*dQa.y -ha.y*Q1.z/N1*py;
				dQa.x=(N/N1)*(Q1.z/Q.z)*dQa.x -ha.x*Q1.z/N1*px;
				dQa.z=-(Q1.x*dQa.x+Q1.y*dQa.y)/Q1.z;

				hb.z=-(hb.x*Q1.x+hb.y*Q1.y)/Q1.z;
				dQb.y=(N/N1)*(Q1.z/Q.z)*dQb.y -hb.y*Q1.z/N1*py;
				dQb.x=(N/N1)*(Q1.z/Q.z)*dQb.x -hb.x*Q1.z/N1*px;
				dQb.z=-(Q1.x*dQb.x+Q1.y*dQb.y)/Q1.z;

				this->ha1[i]=ha; this->dQa1[i]=dQa; this->hb1[i]=hb; this->dQb1[i]=dQb;
			}
			break;
		default:
			{
				// �yha,hb,dQa,dQb �̒�`�z �ɂ��Ắg�����Y���w(����O)�h��(3.4.11)���̐������Q��
				double tau;
				matrix<double> T(3,3);
				vector<double> Q(X,Y,Z), Q1(X1,Y1,Z1), dQa1,dQb1;
				double cx,cy,cxy, G;

				// h,dQ�͒��_���W�ɂ��
				tau=(z-zp)/Z;  // ��_�����ǐՂ͏�Ɍ����i�s�����𕄍��̊�Ƃ��Ă��邱�Ƃɒ���
				               // (Z�Ŋ����Ă���)
				ha=ha+tau*dQa; hb=hb+tau*dQb;      // �ʂƎ�����̌�_��ʂ���ˎ�����ɐ����ȖʂƋߖT�����̌�_
				this->ha[i]=ha; this->dQa[i]=dQa; this->hb[i]=hb; this->dQb[i]=dQb;

				// �@�����W�n�ֈڍs
				T=Tmatrix(vector<double>(l,m,n));  // ���_���W��������ʉߓ_�̖@�����W�֕ϊ�����s��
				Q=T*Q; Q1=T*Q1;
				dQa=T*dQa; dQb=T*dQb;
				ha=T*ha; hb=T*hb;
				ha=ha+( -ha.z/Q.z )*Q; hb=hb+( -hb.z/Q.z )*Q;  // (3.4.44) =�ʂƎ�����̌�_��ʂ�@���ɐ����ȖʂƋߖT�����̌�_
				G=Q1.z-H*Q.z;    // (3.4.36)

				surface_curvature(i, cx,cy,cxy, x,y, T, 1);  // �@�����W�ł̃e�[���[�W�J��2���W��

				dQa1.x=H*dQa.x-G*(cx*ha.x+cxy*ha.y);   // (3.4.35)
				dQb1.x=H*dQb.x-G*(cx*hb.x+cxy*hb.y);   // (3.4.35)
				dQa1.y=H*dQa.y-G*(cy*ha.y+cxy*ha.x);   // (3.4.35)
				dQb1.y=H*dQb.y-G*(cy*hb.y+cxy*hb.x);   // (3.4.35)
				dQa1.z=-( Q1.x*dQa1.x+Q1.y*dQa1.y )/Q1.z;   // (3.4.47) ��� dQ1��Q1
				dQb1.z=-( Q1.x*dQb1.x+Q1.y*dQb1.y )/Q1.z;   // (3.4.47) ��� dQ1��Q1

				// ���_���W�֖߂�
				T=inv(T);
				dQa1=T*dQa1; dQb1=T*dQb1;
				ha=T*ha; hb=T*hb;
				Q1=vector<double>(X1,Y1,Z1);
				ha=ha-sProduct(ha,Q1)*Q1; hb=hb-sProduct(hb,Q1)*Q1;  // �ʂƎ�����̌�_��ʂ�o�ˎ�����ɐ����ȖʂƋߖT�����̌�_
				
				dQa=dQa1; dQb=dQb1;

				this->ha1[i]=ha; this->dQa1[i]=dQa; this->hb1[i]=hb; this->dQb1[i]=dQb;
			}
			break;
		}
	}

	if( E_trace ){
		vector<complex> Q(X,Y,Z), Q1(X1,Y1,Z1), C(l,m,n), S,P,P1;
		cFilm film(0);
		complex tp,ts; double Tp,Ts; 
		double th_deg;
		std::string mname,mname1;

		E=E-sProduct(E,Q)*Q;  // E����E��Q������������������(�d��E�͐i�s����Q�ɐ����j�D
		                      // �璷��������Ȃ����C���ꂪ�Ȃ��Ə����l��^����̂��ʓ|�ɂȂ�D

		this->E[i]=E;  // i�ʒ��_���W�ŕ\�������ˌ��d��
		
		if(CoatName(i)!="") film.open(CoatName(i).c_str());
		if(CoatReverse(i)) film.reverse();
		
		// �����f�[�^�̓��ˑ��}����ݒ�
		mname=gname(i-1);   // i�ʒ��O�̔}��
		if( mname[0]=='-' ) mname.erase(0,1);
		film.set_mname(film.get_k()+1,mname);
		
		// �����f�[�^�̏o�ˑ��}����ݒ�
		if(No*No1>0){              // ���߂̂Ƃ�
			mname1=gname(i);       // i�ʒ���̔}��
			if( mname1[0]=='-' ) mname1.erase(0,1);
		}
		else if(CoatName(i)==""){  // ���ˁC���R�[�g�̐ݒ肪�Ȃ��Ƃ��D
			// ���˂̊�́C�t�@�C���w��̂Ȃ��Ƃ��ȉ��Ƃ���D
			// (�Ⴆ�΁CcFilm�̃f�t�H���g�R���X�g���N�^�����"1"�Ƃ���Ƃ��C
			//  �}����"1"�ł���Ɣ��˂�0�ɂȂ�ǐՂ��~�܂��Ă��܂�)
			mname1= mname=="1" ? "1.5" : "1";
		}
		film.set_mname(0,mname1);

		th_deg=asin( abs(vProduct(Q,C)) )*180/PI;
		// th_deg=acos( Re(sProduct(Q,C)) ) ...  �Ƃ���ƁC���ˊp��0�̂Ƃ�
		// �v�Z�덷�ɂ�� acos �̈�����1���킸���ɒ�����ƒ�`��G���[�ƂȂ��Ă��܂��D

		if( No*No1>0 ){
			tp=film.tp(th_deg,wl[j]); ts=film.ts(th_deg,wl[j]);
			Tp=film.Tp(th_deg,wl[j]); Ts=film.Ts(th_deg,wl[j]);
			tp=tp*sqrt( Tp/(sqabs(tp)) );
			ts=ts*sqrt( Ts/(sqabs(ts)) );
		}
		else{
			tp=film.rp(th_deg,wl[j]); ts=film.rs(th_deg,wl[j]);
		}

		if( abs(vProduct(C,Q))!=0 ){
			// s,p�̐������́C
			//   S=(CxQ)/|CxQ|
			//   P=SxQ
			//   P1= SxQ1 (����)
			//      -SxQ1 (����)
			// �Ƃ���D
			// ����́C�Ⴆ�΁C
			//    (1)�g���w�����Ɛ����Z�p ������ �A���o�b�N�h
			//    (2)�g���w�����t�B���^�[�f�U�C�� ���w�R���M �I�v�g���j�N�X�h
			// �Ɠ��������ł���D
			//    (3)�g�����p�Z�p1988 �g�����w �c���r�� JOEM�h
			//    (4)�g���w���� �����j�Y�� �����o�Łh
			// �Ƃ͈قȂ�D
			// (3)��P1�̎��������ߔ��ˋ� P1=SxQ1 �ƂȂ�C���̓_�ȒP�����C
			// ���˂ɂ����Ĉʑ����^�����Ȃ���Ԃ��Carg(rs)=0, arg(rp)=180��
			// �ƂȂ�C�����I�ɔc�����ɂ����D�i2010.4.8�܂ł͂��̕��@�ł������D�j

			S=vProduct(C,Q)/complex(abs(vProduct(C,Q)));  // S�����P�ʃx�N�g��	
			P=vProduct(S,Q);                              // P�����P�ʃx�N�g��(����)
			if(No*No1>0){
				P1= vProduct(S,Q1);  // P�����P�ʃx�N�g��(���ߌ�)
			}
			else{
				P1=-vProduct(S,Q1);  // P�����P�ʃx�N�g��(���ˌ�)
			}
			E=tp*sProduct(E,P)*P1+ts*sProduct(E,S)*S;
		}
		else{
			E=tp*E;
		}
		
		this->E1[i]=E;  // i�ʒ��_���W�ŕ\�����o�ˌ��d��
	}

	this->X1[i]=X1; this->Y1[i]=Y1; this->Z1[i]=Z1;
	this->xi[i]=xi; this->xi1[i]=xi1;

	return  0;
}



double cLens1::dk(double defocus) {
	if(!Afocal){
		if(s1fix==0) return s1()+defocus;
		else         return s1fix+defocus;
	}
	else{
		return 0;
	}
}

int cLens1::MakePupilGrid(double yObj,double xObj,int FindPupil,int n,int j,int Normalize,
						  int Randomize,int makecoordinate/*=1*/,int ConsiderSymmetry/*=0*/){
	///////  ���i�q�̐ݒ�����郋�[�`�����������̊֐��Ɋ܂܂�邽�߁C���o���s�Ȃ����D /////////
	///////  �{�֐��ł́C�����o list<point> gpoint �Ɋi�q�_��ݒ肷��D
	///////  ����C�e�֐��ɂ����鋤�ʕ�����{�֐��Œu�������čs���\��D     2015.05.08  /////////

	// FindPupil % 10 ==0 : �����̋��E��EPD,EPDx�Ŏw�肷����˓��Ƃ���D
	// FindPupil % 10 !=0 : �����̋��E�͑S�Ă̗L���a��ʂ邱�Ƃ̂ł���ő�̌����̂���ł���D
	// FindPupil / 10 ==0 : �����S��
	// FindPupil / 10 !=0 : ���������̂�
	//    ��F FindPupil== 0 : EPD,EPDx�Ō��܂�����D�����S�́D
	//         FindPupil== 1 : �e�ʗL���a�Ō��܂�����D�����S�́D
	//         FindPupil==10 : EPD,EPDx�Ō��܂�����D���������̂݁D
	//         FindPupil==11 : �e�ʗL���a�Ō��܂�����D���������̂݁D
	//
	// n : �i�q�_�̍s���Ɨ񐔁Dn=0�̂Ƃ���this->nSpot�Ƃ���D
	// j : ��j�g����ΏۂƂ���D
	// Normalize : �^�̂Ƃ��C�����W�� |ymax-ymin|/2, |xmax-xmin|/2 ��1�Ƃ��ċK�i������D
	//             ���̂Ƃ��C�����W�̌��_�͎�����Ƃ���D
	// Randomize !=0 �Ƃ��邱�Ƃɂ��C�œ_���牓���f�t�H�[�J�X�ʒu�ł�OTF�v�Z�Ŋi�q�s�b�`�̉e����h���D
	//
	// �߂�l : �쐬���ꂽ�i�q�_�̑���
	
	int iy,ix;
	double yPupil,xPupil, th, rx,ry, a,b;
	point dr;
	cPoint p;
	point ymax,ymin,xmax,xmin, yw,xw, p0;
	const double dEP=0.000001;
	bool rotationally_symmetric, y_axis_symmetric, valid;

	if(ConsiderSymmetry){
		rotationally_symmetric=IsRotationallySymmetric();
		y_axis_symmetric=IsYAxisSymmetric();
	}

	this->gpoints.RemoveAll();
	this->gpoint_principal=point(0,0);

	if( !FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,j,FindPupil,makecoordinate) ) return 0;
	
	p0=(ymax+ymin)/2;
	yw=ymax-ymin;
	xw=xmax-xmin;
	ry=abs(yw)/2;
	rx=abs(xw)/2;

	if((FindPupil%10 !=0) && (!Normalize)){ // �����̋��E�͊e�ʗL���a�Ō��肵�C���� Normalize �����Ȃ�
		ymax+=yw*PupilMargin/2;  // PupilMargin ��������
		ymin-=yw*PupilMargin/2;
		xmax+=xw*PupilMargin/2;
		xmin-=xw*PupilMargin/2;
	}

	if(n==0) n=this->nSpot;

	for(iy=1; iy<=n; iy++) for(ix=1; ix<=n; ix++){
		if(FindPupil/10 !=0){ // ���������̂�
			if(ix>1) break;   // ���Ȃ킿�Ciy=1�`n �ł̂݌J��Ԃ��D
			th=360.0/n*iy;

			if(FindPupil%10 !=0){
				FindMarginalRay(yPupil,xPupil,th,yObj,xObj,p0.y,p0.x, j,makecoordinate);
			}
			else{
				if(EPDx==0){
					yPupil=EPD/2*sin(th*PI/180);
					xPupil=EPD/2*cos(th*PI/180);
				}
				else{
					// �ȉ~�̕�����  x^2/a^2 + y^2/b^2 = (rcos��)^2 / a^2 + (rsin��)^2 / b^2 = 1
					// ��� r �����܂�̂ŁC x=rcos��, y=rsin�� �ɂ�� x,y �����܂�D
					double r, cs,sn;

					a=EPDx/2;
					b=EPD/2;
					cs=cos(th*PI/180);
					sn=sin(th*PI/180);
					r=sqrt( 1/(cs*cs/a/a+sn*sn/b/b) );
					xPupil=r*cs;
					yPupil=r*sn;
				}
			}
		}
		else{ // ���������݂̂łȂ�
			yPupil=ymin.y+(ymax-ymin).y/n*(iy-0.5+(Randomize ? Random(0.5,0):0));
			xPupil=xmin.x+(xmax-xmin).x/n*(ix-0.5+(Randomize ? Random(0.5,0):0));
		}

		a= rx+dEP/2;
		b= ry+dEP/2;
		valid=false;

		if(FindPupil){
			if(Normalize){
				if((yPupil-p0.y)*(yPupil-p0.y)/b/b+(xPupil-p0.x)*(xPupil-p0.x)/a/a<=1){
					xPupil = rx==0 ? 0 : (xPupil-p0.x)/rx;
					yPupil = ry==0 ? 0 : (yPupil-p0.y)/ry;
					valid=true;
					// �y���z p0��(xPupil,yPupil)�̌��_�ɂȂ�
				}
			}
			else{
				valid=true;
				// �y���z PupilMargin�ɂ�����(xPupil,yPupil)�����w�n��ʉ߂ł���Ƃ͌���Ȃ��D
			}
		}
		else{
			if(yPupil*yPupil/b/b+xPupil*xPupil/a/a<=1){
				if(Normalize){
					xPupil = rx==0 ? 0 : xPupil/rx;
					yPupil = ry==0 ? 0 : yPupil/ry;
				}
				valid=true;
			}
		}
		
		if(valid){
			if(ConsiderSymmetry){
				if(rotationally_symmetric && yObj==0 && xObj==0){
					if(yPupil>=0 && xPupil>=0){
						this->gpoints.AddTail(point(xPupil,yPupil));
					}
				}
				else if(y_axis_symmetric && xObj==0){
					if(xPupil>=0){
						this->gpoints.AddTail(point(xPupil,yPupil));
					}
				}
				else{
					this->gpoints.AddTail(point(xPupil,yPupil));
				}
			}
			else{
				this->gpoints.AddTail(point(xPupil,yPupil));
			}
		}
	}
	
	return this->gpoints.GetSize();
}

int cLens1::MakePupilGrid(int n){
	return MakePupilGrid(0,0,0,n,0,1,0,0);  // normalize=1
}

double cLens1::CenterThickness(int i){
	if(1<=i && i<=k-1) return N(0,1)*N(i,1)>0 ? d(i) : -d(i);
	else               return 0;
}

int cLens1::make_ref_sphere(double yObj,double xObj,
                           double yPupil_principal,double xPupil_principal,double defocus,int j){
	// ����     : this->k+1�� == ref_sphere->k��
	// �Q�Ƌ��� : this->k+2�� == ref_sphere->k+1��
	// 
	// �Q�Ɩʂ̈ʒu��OPD�v�Z���x�ɂ���   06.10.18
	//     �ȑO�͎Q�ƖʂƂ��đ��ʂ��猩�ĎQ�Ƌ��̌��w�n�ŏI�ʂƔ��Α��̖ʂ��g���Ă����D
	//     �������C���̏ꍇ���̗��OPD�ɑ傫�Ȍ덷�𐶂����D
	//
	//             R     d    glass    EA
	//         1   ��          1        1  (��1mm�s���z�[��)
	//                         1
	//         s=667  s1fix=2000  (1D�̃f�t�H�[�J�X)
	//         �ɂ��Ď������������OPD���v�Z����D
	//
	//     �Q�Ƌ��̌��w�n�ŏI�ʑ��̖ʂ��Q�ƖʂƂ��Ď��݂��Ƃ���덷�͉��P���ꂽ���߁C
	//     �{�֐����C�������D
	//
	// 2010.12.09
	//   �S�Ă̌������Q�Ƌ��̒��S�Ɍ������Ă���i�Q�ƕ��ʂł͕��s�j�̂łȂ���΁C
	//   �قȂ��̎Q�Ƌ��ʁi���ʁj�̊Ԃ̌����̒����͈Ⴄ����C
	//   �g�ʎ����͎Q�Ƌ���R�i���ʂ̈ʒu�j�Ɉˑ�����D
	//   ���������āC�œ_�͂��ꂪ�傫���Ƃ��͎Q�Ƌ���R�i�܂��͎Q�ƕ��ʂ̈ʒu�j�ɒ��ӂ���D
	//
	// 2017.05.17
	// ���̗�ł́Cs1fix(�t���l����)�ɂ���ĉ�ܑ��͕ω����Ȃ���΂Ȃ�Ȃ��D
	// ���̂��߂ɂ́C�Q�Ƌ��ʂ͑�1�ʂƈ�v���Ȃ��Ă͂Ȃ�Ȃ��D
	//
	//             R     d    glass    EA
	//         1   ��          1        1  (��1mm�s���z�[��)
	//                         1
	//         s=��  s1fix=100
	// ��ʂɁC�Q�Ƌ��ʂ͍i��ʒu�i�܂��͍i�葜���ˏo��)�Ɉ�v���Ȃ���΂Ȃ�Ȃ��D


	double R;
	if(ref_sphere!=0) delete ref_sphere;
	ref_sphere=new cLens1(k+1,cn);
	ref_sphere->assignment(*this);
	if( !Afocal ) {
		ref_sphere->d(this->k)=dk(defocus);  // ����
		//R= fabs(dk(defocus))<100 ? 100*sgn(dk(defocus)) : dk(defocus);
		// �� �̂悤��R��100�ȏ�ɐ�������ƁC�Ⴆ�΁g1mm�J������50�o�̈ʒu�̋��x���z�h���v�Z�ł��Ȃ��Ȃ�D
		R= dk(defocus)-t1();  // �Q�Ƌ����a�� ������-�ˏo���ʒu
		ref_sphere->EAtype(this->k+1)=0; ref_sphere->EAy(this->k+1)=1e+30;  // ���ʂ̗L���a�͑傫��
		// ���ʋ��܌�ΐS������΍폜 [�J�n]
		if(ref_sphere->decenter_type(this->k+1)==2) ref_sphere->decenter_type(this->k+1)=1;
		ref_sphere->ret(this->k+1)=0;
		ref_sphere->dx1(this->k+1)=0;
		ref_sphere->dy1(this->k+1)=0;
		ref_sphere->dz1(this->k+1)=0;
		ref_sphere->rox1(this->k+1)=0;
		ref_sphere->roy1(this->k+1)=0;
		ref_sphere->roz1(this->k+1)=0;
		// ���ʋ��܌�ΐS������΍폜 [�I��]
		ref_sphere->s1fix=-R;  // �Q�Ƌ��ʂ�z�ʒu
		ref_sphere->Set_gname(this->k+1,gname(k)); // ���ʌ�̔}���͑���Ԕ}���Ɠ���
		// �Q�Ƌ����S���j�g��������Ƒ��ʂ̌�_�Ɉ�v������
		if( ref_sphere->RayTrace(yObj,xObj,yPupil_principal,xPupil_principal,defocus,j,1,0,0,0)!=0 ){
			return 0;
		}
		ref_sphere->decenter_type(this->k+2)=1;
		ref_sphere->dy(this->k+2)=ref_sphere->y[this->k+1];
		ref_sphere->dx(this->k+2)=ref_sphere->x[this->k+1];
		ref_sphere->r(this->k+2)=R+ref_sphere->z[this->k+1];
	}
	else {
		// Afocal�̏ꍇ�Q�Ɩ�(this->k+2��)�͍ŏI��(this->k��)�̈ʒu�ɂ���D
		// �Q�Ɩʂ̈ʒu���ŏI�ʂ��炸�炵�����Ƃ��͌��f�[�^�� *this�̍Ō�ɉˋ�ʂ�ǉ�����K�v������D
		// this->k+1�ʂ͓��ɈӖ��������Ȃ��D
		ref_sphere->d(this->k)=0;  // this->k+1�ʂ͍ŏI�ʂ̈ʒu�ɂ���
		ref_sphere->EAtype(this->k+1)=0; ref_sphere->EAy(this->k+1)=1e+30;  // this->k+1�ʂ̗L���a�͑傫��
		ref_sphere->decenter_type(this->k+1)=0;  // this->k+1�ʂ͕ΐS���Ȃ�
		//ref_sphere->s1fix=0.000001;  // �Q�Ƌ��� this->k+1�ʂ̈ʒu (=�ŏI�ʂ̈ʒu) �ɂ���
		ref_sphere->s1fix=t1();  // �Q�Ƌ��͎ˏo���ʂɂ���i2017.05.17)
		ref_sphere->Set_gname(this->k+1,gname(this->k));  // �Q�Ƌ��̂����Ԃ̔}���͑���Ԕ}��
		// ��j�g���̎�������g���D
		if( ref_sphere->RayTrace(yObj,xObj,yPupil_principal,xPupil_principal,defocus,j,1,0,0,0)!=0 ){
			return 0;
		}
		// �Q�Ɩʂ̌��_������������Ɉړ�����D�������Ȃ��ƃf�t�H�[�J�X�ɂ���ĎQ�Ɩʂ�R�������Ƃ��C
		// �����������0�łȂ�(�ˏo���ƍŏI�ʂ�����Ă���Ƃ���)�Ɣg�ʎ����Ƀ`���g��������D
		// ����ɂ��f�t�H�[�J�X�����Ƃ���OTF�v�Z���ʂɖ����ł��Ȃ��덷�������邱�Ƃ��������D
		// (�����I�ɂ͍��͐����Ȃ���?)
		// �Ⴆ�΃X���b�g�����v�ōŏI�ʂ��ˏo���łȂ��ڊ�ŏI�ʂɂ����Ƃ��̃X���[�t�H�[�J�XOTF�v�Z�D
		// 08.02.14
		ref_sphere->decenter_type(this->k+2)=1;
		ref_sphere->dy(this->k+2)=ref_sphere->y[this->k+1];
		ref_sphere->dx(this->k+2)=ref_sphere->x[this->k+1];
		// �Q�Ɩʂ�������Ɛ����ɂ���D
		// ("*sgn(ref_sphere->Z[k+1]"���Ȃ������̂Ō�����-z�����֐i�s���Ă���ꍇrox��180������Ă��܂����D080403)
		rotate_angle(ref_sphere->rox(this->k+2),ref_sphere->roy(this->k+2),
		             ref_sphere->X[this->k+1]*sgn(ref_sphere->Z[this->k+1]),
					 ref_sphere->Y[this->k+1]*sgn(ref_sphere->Z[this->k+1]),
					 ref_sphere->Z[this->k+1]*sgn(ref_sphere->Z[this->k+1]));
		// defocus>0�̂Ƃ��f�t�H�[�J�X�����Q�Ɩʂ͌����i�s�����ɉ�(������)�Ƃ���D
		// (���x�␳�@�\���t�������w�@��͎��x�␳�l>0�̂Ƃ���ɑ΂��Ď�����������D)
		ref_sphere->r(this->k+2)= defocus==0 ? 0 : 1000/defocus*sgn(ref_sphere->Z[this->k+1]);
	}
	return 1;
}

int cLens1::fno_for_unit_pupil(double& fno_y,double& fno_x,double yObj,double xObj,double defocus,int j,int IsExitPupil){
	// ���㒷��1�ɑ΂��鑜��F�i���o
	//     ����F�i���o = �Q�Ƌ����a / �Q�Ƌ���̌�������(���_�ߖT��)���ʂւ̓��e
	//       ���˓�(IsExitPupil=false)�̂Ƃ��C �Q�Ƌ����a / ���ܗ� / ���˓��㒷��1�ɑΉ�����Q�Ƌ��ɉ����������̑��ʓ��e
	//       �ˏo��(IsExitPupil=true )�̂Ƃ��C �Q�Ƌ����a / ���ܗ�
	//      <�Q�l>  ��������F�i���o�́C
	//                  �Q�Ƌ����a / �Q�Ƌ���̌�����
	//              �ł���D����𑜖ʂ֓��e����̂́C
	//              ���O�ȂǂŎ�����Ƒ��ʂ������łȂ��̂�␳����̂��{���ƍl������D
	//              (�����̏ꍇ���ʂ͕��ʂł��邽��xy�ʂƕ��s�ŁC���e��xy�ʂւ̂��̂ƂȂ�)
	// ���v�Z����D
	// afocal�̂Ƃ���
	//     1 / ���ܗ� / ���˓��㒷��1�ɑΉ�����Q�Ɩʂɉ���������
	//        �� ��Ŕ`�����w�n�̏ꍇ�C��̉�����l����΂���������I�ł���D
	//           (afocal�̂Ƃ��͎Q�Ɩʂ�����������ɕΐS����̂ł���͑���xy�ʓ��e�ł���D)
	// ���v�Z����D
	// <����> ���˓��㒷��x�ɑ΂���F�i���o�� fno_x*x �ł͂Ȃ� fno_x/x �ƂȂ�D
	double r;                   // �Q�Ƌ����a
	point py1,py2,px1,px2;      // ���˓�������ʒu
	double dydpy,dxdpx;
	double N;
	double r0;
	vector<double> vy1,vy2,vx1,vx2, vdy,vdx, v,n;

	if(this->FindPupilEdge(py1,py2,px1,px2,yObj,xObj,j,1)){
		make_ref_sphere(yObj,xObj,((py1+py2)/2).y,((py1+py2)/2).x,defocus,j);
		r=ref_sphere->r(ref_sphere->k+1);
		N=ref_sphere->N(ref_sphere->k,j);

		if(IsExitPupil){
			// �ˏo��(�a�͂��悻�̒l)�ł̎�������4�{�ƎQ�Ƌ��ʏ�̈ʒu vy1,vy2,vx1,vx2
			// x,y�����͎Q�Ƌ����[�J�����W�̂���ł���D
			r0=(distance(py1,py2)+distance(px1,px2))/4;
			vy1=vector<double>(   0, r0,ref_sphere->surface_sag(ref_sphere->k+1, r0,  0,0) );
			vy2=vector<double>(   0,-r0,ref_sphere->surface_sag(ref_sphere->k+1,-r0,  0,0) );
			vx1=vector<double>(  r0,  0,ref_sphere->surface_sag(ref_sphere->k+1,  0, r0,0) );
			vx2=vector<double>( -r0,  0,ref_sphere->surface_sag(ref_sphere->k+1,  0,-r0,0) );
			vdy=vy1-vy2;
			vdx=vx1-vx2;
		}
		else{
			// ��������4�{�ƎQ�Ƌ��ʏ�̈ʒu vy1,vy2,vx1,vx2
			// x,y�����͓��˓����W�̂���ł���D
			// ���������āCz��]�ΐS�ɂ�葜�ʂł͂˂���Ă��邩������Ȃ��D
			vy1=ref_sphere->RayPos(ref_sphere->k+1,yObj,xObj,"",py1.y,py1.x,j);
			vy2=ref_sphere->RayPos(ref_sphere->k+1,yObj,xObj,"",py2.y,py2.x,j);
			vx1=ref_sphere->RayPos(ref_sphere->k+1,yObj,xObj,"",px1.y,px1.x,j);
			vx2=ref_sphere->RayPos(ref_sphere->k+1,yObj,xObj,"",px2.y,px2.x,j);
			vdy=vy1-vy2;
			vdx=vx1-vx2;
		}

		if(Afocal){
			// Afocal�̂Ƃ��͑��ʂ͎�����ɐ����ƍl����D(��Ŕ`�����w�n�̏ꍇ�C��̉�����l����΍����I)
			// dydpy = ����̒P�ʋ����ɑΉ�����Q�Ɩʏ�̒���
			if(IsExitPupil){
				dydpy=abs(vdy)/(vy1.y-vy2.y);
				dxdpx=abs(vdx)/(vx1.x-vx2.x);
			}
			else{
				dydpy=abs(vdy)/(py1.y-py2.y);
				dxdpx=abs(vdx)/(px1.x-px2.x);
			}

			fno_y= fabs(1/N/dydpy);
			fno_x= fabs(1/N/dxdpx);
		}
		else{
			// ���ʂƎ�����̌�_�ɂ����鑜�ʂ̖@���P�ʃx�N�g�� n ���v�Z����
			if(asph_type(this->k+1)==SPH && c(this->k+1)==0){
				// ���ʂ����ʂł���Ƃ�
				n.x=0; n.y=0; n.z=1;
			}
			else{
				// ���ʂ����ʂłȂ��Ƃ�
				v=ref_sphere->RayPos(this->k+1,yObj,xObj,"",((py1+py2)/2).y,((py1+py2)/2).x,j);
				n=ref_sphere->surface_normal(this->k+1,v.y,v.x,0);
				n=n/abs(n);
			}

			// n �𑜖ʃ��[�J�����W����Q�Ƌ����[�J�����W�֕ϊ�����D
			// (afocal�łȂ���� make_ref_sphere() ��藼���W�͕��s�Ȃ̂ŏ璷�����C�֐��̓Ɨ��������߂邽�߁j
			ref_sphere->transform(n,this->k+1,0,ref_sphere->k+1,0,0);

			// vdy,vdx �� n �ɐ����Ȗʂɓ��e����C���Ȃ킿���_�ߖT�̑��ʂɓ��e���C
			// ����̒����Ŋ���D
			// dydpy = ����̒P�ʋ����ɑΉ�����Q�Ɩʏ�̒����̑��ʂւ̓��e
			if(IsExitPupil){
				dydpy=abs( vdy-sProduct(vdy,n)*n )/(vy1.y-vy2.y);
				dxdpx=abs( vdx-sProduct(vdx,n)*n )/(vx1.x-vx2.x);
			}
			else{
				dydpy=abs( vdy-sProduct(vdy,n)*n )/(py1.y-py2.y);
				dxdpx=abs( vdx-sProduct(vdx,n)*n )/(px1.x-px2.x);
			}

			fno_y=fabs(r/N/dydpy);
			fno_x=fabs(r/N/dxdpx);
		}

		return 1;
	}
	else{
		fno_y=fno_x=0;
		return 0;
	}
}
/*   **** ���R�[�h�i�O�̂��߁C���΂炭�ۑ�����D 20140616) ****
int cLens1::fno_for_unit_pupil(double& fno_y,double& fno_x,double yObj,double xObj,double defocus,int j,int IsExitPupil){
	// ���㒷��1�ɑ΂��鑜��F�i���o
	//     ����F�i���o = �Q�Ƌ����a / �Q�Ƌ���̒�����xy�ʂւ̓��e (�g�����Y�݌v�̂��߂̔g�ʌ��w(����)�h (2-13) )
	//     ���˓�(IsExitPupil=false) : �Q�Ƌ����a / ���ܗ� / ���˓��㒷��1�ɑΉ�����Q�Ƌ��ɉ�����������xy�ʓ��e
	//     �ˏo��(IsExitPupil=true ) : �Q�Ƌ����a / ���ܗ�
	//      <�Q�l> �Q�Ƌ����a = ���ɉ������������^cos��
	//             �Q�Ƌ��㒷����xy�ʓ��e = ����������݂����˓��a * (cos��)^2
	//             �ƂȂ�C�q�ߕ������OF�i���o�͈�ʂɌ�����悤�Ɏ���� 1/(cos��)^3 �ƂȂ�D
	//             ����́g�����Y�݌v�̂��߂̔g�ʌ��w(����)�h�� (4.7)���i�Q�Ƌ�R�ƎQ�Ƌ��Ƃ̌�_
	//             �̃ăŖʓ��e�ŉ�ܑ��̑傫�������܂�j�Ƃ���v����D
	//             �i�ʂ̌������������ (cos��)^3 ��3��̈Ӗ��́C
	//                 �E���_���猩��Ɠ��a�� cos�� �{�ɂȂ�
	//                 �E�����瑜�_�̋���������� 1/cos�� �{�ɂȂ�
	//                 �E������Ƒ��ʂ̂Ȃ��p�ɂ��PSF�� 1/cos�� �{�ɐL�т�
	//              �ł���j
	// ���v�Z����D
	// afocal�̂Ƃ���
	//     1 / ���ܗ� / ���˓��㒷��1�ɑΉ�����Q�Ɩʂɉ���������
	//        �� ��Ŕ`�����w�n�̏ꍇ�C��̉�����l����΂���������I�ł���D
	//           (afocal�̂Ƃ��͎Q�Ɩʂ�����������ɕΐS����̂ł���͑���xy�ʓ��e�ł���D)
	// ���v�Z����D
	// <����> ���˓��㒷��x�ɑ΂���F�i���o�� fno_x*x �ł͂Ȃ� fno_x/x �ƂȂ�D
	double r;                   // �Q�Ƌ����a
	point py1,py2,px1,px2;      // ���˓�������ʒu
	point y1,y2,x1,x2;          // �Q�Ƌ��ʏ�����ʒu
	double dydpy,dxdpx;
	double N;
	if(this->FindPupilEdge(py1,py2,px1,px2,yObj,xObj,j,1)){
		make_ref_sphere(yObj,xObj,((py1+py2)/2).y,((py1+py2)/2).x,defocus,j);
		r=ref_sphere->r(ref_sphere->k+1);
		N=ref_sphere->N[ref_sphere->k][j];
		if(IsExitPupil){
			fno_y= !Afocal ? fabs(r/N) : fabs(1/N);
			fno_x= fno_y;
		}
		else{
			y1=ref_sphere->RayHeight(ref_sphere->k+1,yObj,xObj,"",py1.y,py1.x,j);
			y2=ref_sphere->RayHeight(ref_sphere->k+1,yObj,xObj,"",py2.y,py2.x,j);
			x1=ref_sphere->RayHeight(ref_sphere->k+1,yObj,xObj,"",px1.y,px1.x,j);
			x2=ref_sphere->RayHeight(ref_sphere->k+1,yObj,xObj,"",px2.y,px2.x,j);
			// x1��x2�̋���dx��|x1.x-x2.x|�ł͂Ȃ��ꍇ������(�r����z������̕ΐS�������ꍇ�Ȃ�)
			// ����distance�֐����g���D
			dydpy=distance(y1,y2)/(py1.y-py2.y);
			dxdpx=distance(x1,x2)/(px1.x-px2.x);
			fno_y= !Afocal ? fabs(r/N/dydpy) : fabs(1/N/dydpy);
			fno_x= !Afocal ? fabs(r/N/dxdpx) : fabs(1/N/dxdpx);
		}
		return 1;
	}
	else{
		fno_y=fno_x=0;
		return 0;
	}
}
*/

void cLens1::scale_condition(double m){
	if(m==0) return;
	if(fabs(s)<LN) s*=m; 
	if(fabs(t)<LN) t*=m;
	s1fix*=m;
	if(fabs(s)<LN) { yObjectMax*=m; xObjectMax*=m; }
	if(fabs(t)<LN) { EPD*=m; EPDx*=m; EPy*=m; EPx*=m; }
}



///// public members ////////////////////////////////////////////////////////////////////////

double cLens1::LargeNumber(){ return LN;}

double cLens1::NewtonToR(double r0,double EAphi,double FRW_nm,double rings){
	double dz,z;
	if(rings==0) return r0;
	dz=(FRW_nm/2/1000000)*rings;
	if (r0==0) { z=dz; }
	if (r0> 0) { z=r0-sqrt(r0*r0-(EAphi/2)*(EAphi/2)); z+=dz; }
	if (r0< 0) { z=r0+sqrt(r0*r0-(EAphi/2)*(EAphi/2)); z+=dz; }
	return ( z*z+(EAphi/2)*(EAphi/2) )/2/z; 
}

double cLens1::RToNewton(double r0,double EAphi,double FRW_nm,double r){
	double z0,z;
	if (r0==0) z0=0;
	if (r0> 0) z0=r0-sqrt(r0*r0-(EAphi/2)*(EAphi/2));
	if (r0< 0) z0=r0+sqrt(r0*r0-(EAphi/2)*(EAphi/2));
	if (r==0) z=0;
	if (r> 0) z=r-sqrt(r*r-(EAphi/2)*(EAphi/2));
	if (r< 0) z=r+sqrt(r*r-(EAphi/2)*(EAphi/2));
	return (z-z0)/(FRW_nm/2/1000000);
}

std::string cLens1::RingsSpecRefSurf(double th_deg,double VirtualAS){
	// ���ˊpth�̔��˖ʂ̌����ɉ������ő�A�X�����ʂ�,
	// ���ˊp0�x�C���ܗ�1.5�C�}����C�C�A�X�{���� VirtualAS �̓��ߖ�
	//  (�g�ʎ����� VirtualAS/4 �g���ƂȂ�C��ʓI�ȃA�X1�{�̂Ƃ��g�ʎ����̓�/4�ɂȂ�D)
	// �̃A�X�����ʈȉ��Ƃ��������̂Ƃ��C
	// ���̔��˖ʂ� �j���[�g���{���K�iN, �A�X�{���K�iAS����肤��͈͂́C
	// N��x���CAS��y���ɂƂ�ƁC
	//    N/cos(th) - (N - AS)cos(th) = VirtualAS/4  (�Ō�� /4 �͔��ˁC���܂̈Ⴂ���j
	//    AS=2N  (AS�K�i��N�K�i��2�{�ȏ�ɂȂ肦�Ȃ�)
	// ��2�{�̒�����x�����͂ގO�p�`�ƂȂ�D
	// �����ł́C
	//     (1) N�̎�肤��ő�lN1�i���̂Ƃ�AS��0)
	//     (2) AS�̎�肤��ő�lAS2�Ƃ��̂Ƃ���N2
	//     (3) AS=N�̂Ƃ��̍ő�l N3
	//     (4) AS=N/2�̂Ƃ��̍ő�l N4
	//     (5) AS=N/3�̂Ƃ��̍ő�l N5
	//     (6) AS=N/4�̂Ƃ��̍ő�l N6
	//     (7) AS=N/5�̂Ƃ��̍ő�l N7
	// ��3�_���v�Z���\������D(1)(2)�Ɍ��_��������3�_���O�p�`�̒��_�ƂȂ�D
	//
	// tips > VirtualAS �ɁC(�ʐ��x�����a/�����a)^2 / �}�����ܗ� ���悶���l��^�����
	//        ���Z�̎�Ԃ��Ȃ��Ȃ�D

	std::string s;
	double th, N1,N2,N3,N4,N5,N6,N7,AS2;
	char buf[100];
	
	th=th_deg*PI/180;
	N1=VirtualAS/4/(1/cos(th)-cos(th));
	N2=(VirtualAS/4/cos(th))/(1/cos(th)/cos(th)+1);
	AS2=N2*2;
	N3=VirtualAS/4/(1/cos(th));
	N4=VirtualAS/4/(1/cos(th)-(1.0-1.0/2.0)*cos(th));
	N5=VirtualAS/4/(1/cos(th)-(1.0-1.0/3.0)*cos(th));
	N6=VirtualAS/4/(1/cos(th)-(1.0-1.0/4.0)*cos(th));
	N7=VirtualAS/4/(1/cos(th)-(1.0-1.0/5.0)*cos(th));

	if(th!=0){ 
		sprintf(buf,"N<=%g (N=%g then AS=0)\n",N1,N1); s+=buf;
		sprintf(buf,"AS<=%g (AS=%g then N=%g)\n",AS2,AS2,N2); s+=buf;
		sprintf(buf,"N<=%g (AS=%g)\n",N3,N3); s+=buf;
		sprintf(buf,"N<=%g (AS=%g)\n",N4,N4/2); s+=buf;
		sprintf(buf,"N<=%g (AS=%g)\n",N5,N5/3); s+=buf;
		sprintf(buf,"N<=%g (AS=%g)\n",N6,N6/4); s+=buf;
		sprintf(buf,"N<=%g (AS=%g)\n",N7,N7/5); s+=buf;
	}
	else{
		sprintf(buf,"AS<=%g (N>=2AS)\n", VirtualAS/4); s+=buf;
	}
	return s;
}

double cLens1::lens_phi(double ea_phi){
	// REG401D �\4�̒����R�����̒l�i����1�ۂ�荶�R�����̒l�ł͎���̌������������Ƃ̏��L)
	// �i�����s������ea_phi�ɉ�����Ă���萔�𗼕ӂ�������ΘR��͂Ȃ����Ƃ��m�F�ł���D�j
	if( ea_phi<=0                       ) return ea_phi;      // �����Ռ��̏ꍇ
	if( ea_phi+1.5<=30                  ) return ea_phi+1.5;
	if( 30<ea_phi+2 && ea_phi+2<=50     ) return ea_phi+2;
	if( 50<ea_phi+2.5 && ea_phi+2.5<=80 ) return ea_phi+2.5;
	if( 80<ea_phi+3                     ) return ea_phi+3;
	return 0;
}

double cLens1::lens_ea(double phi){
	if( phi<=30           ) return phi-1.5;
	if( 30<phi && phi<=50 ) return phi-2;
	if( 50<phi && phi<=80 ) return phi-2.5;
	if( 80<phi            ) return phi-3;
	return 0;
}

double cLens1::ZValue(double r1,double r2,double phi1,double phi2){
	// �x���N�����v���S�o���\���]��
	// Z��0.1        ����
	// 0.1��Z��0.15  �\���L��
	// 0.15��Z       �\
	double c1,c2;
	c1= r1==0 ? 0 : 1/r1;
	c2= r2==0 ? 0 : 1/r2;
	return fabs( phi1*c1-phi2*c2 )/4;
}

double cLens1::ToroidZ(double x,double y,double rx,double ry,double kp,int IsXToroid,int &domain_error){
	// �g���C�_���ʂ�sag��z���v�Z����
	//   x,y   : x,y���W
	//   rx,ry : ��ȗ����a(��o����x,y�����Ƃ���)
	//   kp    : ��]������~���Ȑ��̃R�[�j�b�N�W��
	//   IsXToroid : true  X�g���C�h�C���arx�̉~���Ȑ��𔼌ary�ŉ�
	//               false Y�g���C�h�C���ary�̉~���Ȑ��𔼌arx�ŉ�
	double z,cx,cy,al,n;
	
	cx= rx==0 ? 0 : 1/rx;
	cy= ry==0 ? 0 : 1/ry;	
	domain_error=0;

	if(IsXToroid){
		n=1-(kp+1)*x*x*cx*cx;
		if(n<0){ domain_error=1; return 0; }
		al=1+sqrt(n);
		n=(1-x*x*cy*cx/al)*(1-x*x*cy*cx/al)-y*y*cy*cy;
		if(n<0){ domain_error=1; return 0; }
		z=ry-ry*sqrt(n);
	}
	else{
		// x,y������
		n=1-(kp+1)*y*y*cy*cy;
		if(n<0){ domain_error=1; return 0; }
		al=1+sqrt(n);
		n=(1-y*y*cx*cy/al)*(1-y*y*cx*cy/al)-x*x*cx*cx;
		if(n<0){ domain_error=1; return 0; }
		z=rx-rx*sqrt(n);
	}
	
	return z;
}

double cLens1::ToroidZ(double x,double y,double rx,double ry,double kp,int IsXToroid){
	int dummy;
	return ToroidZ(x,y,rx,ry,kp,IsXToroid,dummy);
}

void cLens1::ToroidDerivative(double &zx,double &zy,double x,double y,double rx,double ry,double kp,int IsXToroid){
	// �g���C�_����z(x,y)�� zx=dz/dx, zy=dz/dy���v�Z����D
	//   x,y   : x,y���W
	//   rx,ry : ��ȗ����a(��o����x,y�����Ƃ���)
	//   kp    : ��]������~���Ȑ��̃R�[�j�b�N�W��
	//   IsXToroid : true  X�g���C�h�C���arx�̉~���Ȑ��𔼌ary�ŉ�
	//               false Y�g���C�h�C���ary�̉~���Ȑ��𔼌arx�ŉ�

	//  �Ⴆ��X�g���C�h�̂Ƃ��C
	//      (z-ry)^2+y^2 = {cx*x^2/(1+sqrt(1-(kp+1)*cx^2*x^2)) -ry}^2 = {A-ry}^2
	//  �̗��ӂ�x,y�ŕΔ������邱�Ƃ�zx,zy���v�Z����D
	double cx,cy,A,Ay,Ax,root,z;
	
	cx= rx==0 ? 0 : 1/rx;
	cy= ry==0 ? 0 : 1/ry;	

	if(IsXToroid){
		root=sqrt(1-(kp+1)*cx*cx*x*x);
		A=cx*x*x/(1+root);
		Ax=cx*x/root;
		z=ry-ry*sqrt((cy*A-1)*(cy*A-1)-cy*cy*y*y);
		zy=y*cy/(1-z*cy);
		zx=(A*cy-1)*Ax/(z*cy-1);
	}
	else{
		// x,y������
		root=sqrt(1-(kp+1)*cy*cy*y*y);
		A=cy*y*y/(1+root);
		Ay=cy*y/root;
		z=rx-rx*sqrt((cx*A-1)*(cx*A-1)-cx*cx*x*x);
		zx=x*cx/(1-z*cx);
		zy=(A*cx-1)*Ay/(z*cx-1);
	}
}

void cLens1::Toroid2ndDerivative(double &zxx,double &zxy,double &zyy,double &zx,double &zy,
                                 double x,double y,double rx,double ry,double kp,int IsXToroid){
	// �g�[���b�N��z(x,y)��2�������W�� zxx,zxy,zyy ���v�Z����D
	//   x,y   : x,y���W
	//   rx,ry : ��ȗ����a(��o����x,y�����Ƃ���)
	//   kp    : ��]������~���Ȑ��̃R�[�j�b�N�W��
	//   IsXToroid : true  X�g���C�h�C���arx�̉~���Ȑ��𔼌ary�ŉ�
	//               false Y�g���C�h�C���ary�̉~���Ȑ��𔼌arx�ŉ�

	//  �Ⴆ��X�g���C�h�̂Ƃ��C
	//      (z-ry)^2+y^2 = {cx*x^2/(1+sqrt(1-(kp+1)*cx^2*x^2)) -ry}^2 = {A-ry}^2
	//  �̗��ӂ�x,y�ŕΔ������邱�Ƃ�zx,zy���v�Z����D
	double cx,cy,A,Ay,Ayy,Ax,Axx,root,z;
	
	cx= rx==0 ? 0 : 1/rx;
	cy= ry==0 ? 0 : 1/ry;	

	if(IsXToroid){
		root=sqrt(1-(kp+1)*cx*cx*x*x);
		A=cx*x*x/(1+root);
		Ax=cx*x/root;
		Axx=(cx*root+cx*cx*cx*x*x*(kp+1)/root)/root/root;
		z=ry-ry*sqrt((cy*A-1)*(cy*A-1)-cy*cy*y*y);
		zy=y*cy/(1-z*cy);
		zx=(A*cy-1)*Ax/(z*cy-1);
		zyy=(cy+cy*cy*(y*zy-z))/(1-z*cy)/(1-z*cy);
		zxy=y*zx*cy*cy/(1-z*cy)/(1-z*cy);
		zxx=((cy*Ax*Ax+(A*cy-1)*Axx)*(z*cy-1)-(A*cy-1)*Ax*cy*zx)/(z*cy-1)/(z*cy-1);
	}
	else{
		// x,y������
		root=sqrt(1-(kp+1)*cy*cy*y*y);
		A=cy*y*y/(1+root);
		Ay=cy*y/root;
		Ayy=(cy*root+cy*cy*cy*y*y*(kp+1)/root)/root/root;
		z=rx-rx*sqrt((cx*A-1)*(cx*A-1)-cx*cx*x*x);
		zx=x*cx/(1-z*cx);
		zy=(A*cx-1)*Ay/(z*cx-1);
		zxx=(cx+cx*cx*(x*zx-z))/(1-z*cx)/(1-z*cx);
		zxy=x*zy*cx*cx/(1-z*cx)/(1-z*cx);
		zyy=((cx*Ay*Ay+(A*cx-1)*Ayy)*(z*cx-1)-(A*cx-1)*Ay*cx*zy)/(z*cx-1)/(z*cx-1);
	}	
}

double cLens1::ToroidRa(double x,double y,double rx,double ry,double kp,int IsXToroid){
	// �g�[���b�N�ʂ�axial�ȗ����aRa���v�Z����
	//   x,y   : x,y���W
	//   rx,ry : ��ȗ����a(��o����x,y�����Ƃ���)
	//   IsXToroid : true  X�g���C�h�C���arx�̉~���Ȑ��𔼌ary�ŉ�
	//               false Y�g���C�h�C���ary�̉~���Ȑ��𔼌arx�ŉ�
	// Ra�Ƃ́Cz�����܂ޒf�ʂɂ����āC�ʏ�̓_A�ł̖@����z���ƌ����_B�̋���A-B�ł���D
	// (3�����ōl���Ă��܂��ƁC��ʂɖʂ̖@����z���ƌ����Ȃ�����Ra���`�ł��Ȃ��D�j
	double zx,zy,r,zr,ra;

	r=sqrt(x*x+y*y);
	ToroidDerivative(zx,zy,x,y,rx,ry,kp,IsXToroid);
	zr= r==0 ? 0 : zx*x/r+zy*y/r;   // z��r�Δ���=(dz/dx)(dx/dr)+(dz/dy)(dy/dr)
	ra= zr==0 ? 0 : r/sin(atan(zr));
	return ra;
}

double cLens1::AASZ(double x,double y,double rx,double ry,double kpx,double kpy,int &domain_error){
	// �A�i�����t�B�b�N�񋅖�(AAS)�̃T�O���v�Z����D
	// z = (cx*x*x+cy*y*y)/[ 1+sqrt{1-(kpx+1)*cx*cx*x*x-(kpy+1)*cy*cy*y*y} ]
	double cx,cy,n;

	cx= rx==0 ? 0 : 1/rx;
	cy= ry==0 ? 0 : 1/ry;
	n=1-(kpx+1)*cx*cx*x*x-(kpy+1)*cy*cy*y*y;
	if(n>=0){
		domain_error=0;
		return (cx*x*x+cy*y*y)/(1+sqrt(n));
	}
	else{
		domain_error=1;
		return 0;
	}
}

double cLens1::AASZ(double x,double y,double rx,double ry,double kpx,double kpy){
	int dummy;
	return AASZ(x,y,rx,ry,kpx,kpy,dummy);
}

void cLens1::AASDerivative(double &zx,double &zy,
                           double x,double y,double rx,double ry,double kpx,double kpy){
	double cx,cy, a,b,rt, ax,ay,bx,by;

	cx= rx==0 ? 0 : 1/rx;
	cy= ry==0 ? 0 : 1/ry;
	// z=a/b
	a=cx*x*x+cy*y*y;
	rt=sqrt(1-(kpx+1)*cx*cx*x*x-(kpy+1)*cy*cy*y*y);
	b=1+rt;
	
	ax=2*cx*x;
	ay=2*cy*y;
	bx=-(kpx+1)*cx*cx*x/rt;
	by=-(kpy+1)*cy*cy*y/rt;

	zx=(ax*b-a*bx)/b/b;
	zy=(ay*b-a*by)/b/b;
}

void cLens1::AAS2ndDerivative(double &zxx,double &zxy,double &zyy,double &zx,double &zy,
                              double x,double y,double rx,double ry,double kpx,double kpy){
	double cx,cy, a,b,rt, ax,ay,bx,by, axx,axy,ayy,bxx,bxy,byy;

	cx= rx==0 ? 0 : 1/rx;
	cy= ry==0 ? 0 : 1/ry;
	a=cx*x*x+cy*y*y;
	rt=sqrt(1-(kpx+1)*cx*cx*x*x-(kpy+1)*cy*cy*y*y);
	b=1+rt;

	ax=2*cx*x;
	ay=2*cy*y;
	axx=2*cx;
	axy=0;
	ayy=2*cy;

	bx=-(kpx+1)*cx*cx*x/rt;
	by=-(kpy+1)*cy*cy*y/rt;
	bxx=(-(kpx+1)*cx*cx*rt-(kpx+1)*(kpx+1)*cx*cx*cx*cx*x*x/rt)/(rt*rt);
	bxy=-(kpx+1)*cx*cx*x*(kpy+1)*cy*cy*y/(rt*rt*rt);
	byy=(-(kpy+1)*cy*cy*rt-(kpy+1)*(kpy+1)*cy*cy*cy*cy*y*y/rt)/(rt*rt);
	
	zx=(ax*b-a*bx)/b/b;
	zy=(ay*b-a*by)/b/b;
	zxx=((axx*b-a*bxx)*b*b-(ax*b-a*bx)*2*b*bx)/(b*b*b*b);
	zxy=((axy*b+ax*by-ay*bx-a*bxy)*b*b-(ax*b-a*bx)*2*b*by)/(b*b*b*b);
	zyy=((ayy*b-a*byy)*b*b-(ay*b-a*by)*2*b*by)/(b*b*b*b);
}

void cLens1::ConicOffAxis(double &z,double &zx,double &zy,double &zxx,double &zxy,double &zyy,
						  double x,double y,double a,double b,double t,int &domain_error){
	// �gOff-Axial���w�n�̋ߎ��E�����_�I���(�r�،h�� ����w�ʘ_��)�h��(2.2)���ɂ��C
	// �gOff-Axial�񎟋Ȗʁh�̃T�O�ʂƂ��̔����W�����v�Z����D
	//  �������C�ȉ��̂悤��a,b���`�������߁C�_���� a-b �� b-a �ɒu��������D

	//  a = yz>0 �̗̈�i��1�܂��͑�3�ی��j�ɂ���œ_�܂ł̌��_����̋����Dz�̐��̑��ɂ���Ƃ����Ƃ���D
	//  b = yz<0 �̗̈�i��2�܂��͑�4�ی��j�ɂ���œ_�܂ł̌��_����̋����Dz�̐��̑��ɂ���Ƃ����Ƃ���D
    //  t = ���_�Əœ_�����Ԓ�����z���̂Ȃ��p�x(deg)�D���̋��p�ŕ\���΂悢�D
	double p, cost,sint,tant, ca,cb;
	double u, ux,uy, uxx,uxy,uyy, v, vx,vy, vxx,vxy,vyy;
	double z1,z2,z3, z1x,z1y,z2x,z2y,z3x,z3y, z1xx,z1xy,z1yy,z2xx,z2xy,z2yy,z3xx,z3xy,z3yy, sqrtz3;

	p=1;   // �_���̎��ɂ���p��1�ɌŒ肷��D���Ȃ킿�C2�œ_�����Ԓ����Ɋւ��ĉ�]�Ώ̂Ƃ���D

	ca= a==0 ? 0 : 1/a;
	cb= b==0 ? 0 : 1/b;
	t=t*PI/180;
	cost=cos(t);
	sint=sin(t);
	tant=tan(t);

	z3=1.0 +(cb-ca)*sint*y -ca*cb*y*y -(ca*cb+0.25*tant*tant*(ca+cb)*(ca+cb))*p*x*x;  // ������

	if(z3>=0){
		domain_error=0;

		// �T�O�ʂ̌v�Z z = u/v, u=z1, v=(cost*2)*(z2+��z3)
		z1=(ca+cb)*(cost*cost*y*y+p*x*x);  // �_���̎��̕��q
		z2=1.0+0.5*(cb-ca)*sint*y;         // �_���̎��̕��ꊇ�ʓ���1��
		sqrtz3=sqrt(z3);                   // �_���̎��̕��ꊇ�ʓ���2��
		u=z1;
		v=(cost*2)*(z2+sqrtz3);
		z=u/v;                             // �T�O��
		
		// ��K�����̌v�Z
		z1x=(ca+cb)*(p*x*2);
		z1y=(ca+cb)*(cost*cost*y*2);
		z2x=0;
		z2y=0.5*(cb-ca)*sint;
		z3x=-(ca*cb+0.25*tant*tant*(ca+cb)*(ca+cb))*p*x*2;
		z3y=(cb-ca)*sint -ca*cb*y*2;
		
		ux=z1x;                    
		uy=z1y; 
		vx=(cost*2)*( z2x+z3x/(sqrtz3*2) );
		vy=(cost*2)*( z2y+z3y/(sqrtz3*2) );

		zx=(ux*v-u*vx)/(v*v);
		zy=(uy*v-u*vy)/(v*v);

		// ��K�����̌v�Z
		z1xx=(ca+cb)*(p*2);
		z1xy=0;
		z1yy=(ca+cb)*(cost*cost*2);
		z2xx=0;
		z2xy=0;
		z2yy=0;
		z3xx=-(ca*cb+0.25*tant*tant*(ca+cb)*(ca+cb))*p*2;
		z3xy=0;
		z3yy=-ca*cb*2;

		uxx=z1xx;
		uxy=z1xy;
		uyy=z1yy;
		vxx=(cost*2)*(z2xx+(z3xx*sqrtz3*2 -z3x*z3x/sqrtz3))/(z3*4);
		vxy=(cost*2)*(z2xy+(z3xy*sqrtz3*2 -z3x*z3y/sqrtz3))/(z3*4);
		vyy=(cost*2)*(z2yy+(z3yy*sqrtz3*2 -z3y*z3y/sqrtz3))/(z3*4);

		zxx=( (uxx*v-u*vxx)*v*v-(ux*v-u*vx)*2*v*vx )/(v*v*v*v);
		zxy=( (uxy*v+ux*vy-uy*vx-u*vxy)*v*v-(ux*v-u*vx)*2*v*vy )/(v*v*v*v);
		zyy=( (uyy*v-u*vyy)*v*v-(uy*v-u*vy)*2*v*vy )/(v*v*v*v);
	}
	else{
		domain_error=1;
		z=zx=zy=zxx=zxy=zyy=0;
	}
}

double cLens1::ConicOffAxisZ(double x,double y,double a,double b,double t,int &domain_error){
	double z,zx,zy,zxx,zxy,zyy;

	ConicOffAxis(z,zx,zy,zxx,zxy,zyy,x,y,a,b,t,domain_error);
	return z;
}

double cLens1::ConicOffAxisZ(double x,double y,double a,double b,double t){
	int domain_error;

	return ConicOffAxisZ(x,y,a,b,t,domain_error);
}

void cLens1::ConicOffAxisDerivative(double &zx,double &zy,
                                    double x,double y,double a,double b,double t){
	double z,zxx,zxy,zyy;
	int domain_error;

	ConicOffAxis(z,zx,zy,zxx,zxy,zyy,x,y,a,b,t,domain_error);
}

void cLens1::ConicOffAxis2ndDerivative(double &zxx,double &zxy,double &zyy,double &zx,double &zy,
                                       double x,double y,double a,double b,double t){
	double z;
	int domain_error;

	ConicOffAxis(z,zx,zy,zxx,zxy,zyy,x,y,a,b,t,domain_error);
}

void cLens1::rotate_angle(double& rox,double& roy,double Zx,double Zy,double Zz){
	// (Zx,Zy,Zz) ��]��Z������
	double o=sqrt(Zx*Zx+Zy*Zy+Zz*Zz);
	Zx/=o; Zy/=o; Zz/=o;
	rox=-arg(complex(Zz,Zy));   // ZX�ʂ�VZ���Ɉ�v������D
	roy=asin(Zx);               // Z����VZ���ƈ�v������irox��]�ɂ��V��Z���͉s�p���Ȃ��Ă���̂ł���ł悢�j
	rox*=180/PI; roy*=180/PI;
}
void cLens1::rotate_angle2(double& rox,double& roy,double& roz,double Zx,double Zy,double Zz,double Xx,double Xy){
	// (Zx,Zy,Zz) ��]��Z������
	// (Xx,Xy,* ) ��]��X������
	// ���_�F
	//   ��]��Z����XY�ʓ��̂Ƃ�(roz�ɂ��X���̋O�Ղ�Z���������猩�ĖʂɂȂ�Ȃ��Ƃ��j
	//   roz�̉���(Xz��^���Ȃ���)��ӂłȂ��D
	double o=sqrt(Zx*Zx+Zy*Zy+Zz*Zz);
	Zx/=o; Zy/=o; Zz/=o;
	rox=-arg(complex(Zz,Zy));
	roy=asin(Zx);
	// roz��]�O��Y���P�ʃx�N�g���CZ���P�ʃx�N�g����ey,ez�Ƃ���ƁC
	//     ey=( 0, cos(rox), sin(rox) )
	//     ez=( sin(roy), -cos(roy)sin(rox), cos(roy)cos(rox) )
	// X���P�ʃx�N�g��ex=ey x ez
	//     ex=( cos(roy), sin(rox)sin(roy), -cos(rox)sin(roy) ) 
	// roz��]���X���P�ʃx�N�g���́Ccos(roz)ex+sin(roz)ey �ł��邱�Ƃ���C
	//     Xx=cos(roz)cos(roy)
	//     Xy=cos(roz)sin(rox)sin(roy)+sin(roz)cos(rox)
	// �v�Z�����
	//     cos(roz)cos(rox)cos(roy)= cos(rox)Xx
	//     sin(roz)cos(rox)cos(roy)=-sin(rox)sin(roy)Xx+cos(roy)Xy
	// ����������
	//     roz=arg( complex( cos(roz),sin(roz) ) )
	// ���C
	roz=arg( sgn(cos(rox)*cos(roy)) * complex( cos(rox)*Xx, -sin(rox)*sin(roy)*Xx+cos(roy)*Xy ) );
	// ( Xx,Xy��0(Z���ɋ߂�)�̂Ƃ��㎮��arg�̈������ق�complex(0,0)�ɂȂ�roz�̐��x�������Ȃ�D)
	// ( sgn(...) ���Ȃ��Ɓ}180���̕s�萫�������� �j
	rox*=180/PI; roy*=180/PI; roz*=180/PI;
}
void cLens1::rotate_angle3(double& rox,double& roy,double& roz,double Zx,double Zy,double Zz,double Yx,double Yy){
	// (Zx,Zy,Zz) ��]��Z������
	// (Yx,Yy,* ) ��]��Y������
	double o=sqrt(Zx*Zx+Zy*Zy+Zz*Zz), sin_roz,cos_roz;
	Zx/=o; Zy/=o; Zz/=o;
	roy=asin(Zx);
	rox=-arg(complex(Zz,Zy));
	// roz��]�O��Z���P�ʃx�N�g���CY���P�ʃx�N�g����ez,ey�Ƃ���ƁC
	// X���P�ʃx�N�g��ex=ey x ez
	// roz��]���Y���P�ʃx�N�g���́C-sin(roz)ex+cos(roz)ey �ł��邱�Ƃ���C
	sin_roz=-Yx/cos(roy);
	cos_roz=(Yy+sin(rox)*sin(roy)*sin_roz)/cos(rox);
	roz=arg( complex(cos_roz,sin_roz) );
	rox*=180/PI; roy*=180/PI; roz*=180/PI;
}

matrix<double> cLens1::Tmatrix(const vector<double>& new_zdirection){
	// ��ԃx�N�g��������z���� new_zdirection ���������W�n�ɕϊ�����s��i�g�����Y���w�h3.3.1 )
	// �V�C�����W�n��y���ƐV����z���͓��ꕽ�ʏ�ɂ���Ƃ���D
	double rox,roy;
	matrix<double> Rx(3,3),Ry(3,3);
	rotate_angle(rox,roy,new_zdirection.x,new_zdirection.y,new_zdirection.z);
	return ::Tmatrix(rox,roy,0);
}


double cLens1::ACoefficient(int cOrder,double phi,double N,double al,double h) {
	// "�����Y�݌v�@" eq.(4.29). Return the coefficient of front curvature c^cOrder.
	if(N==1) return 0;
	switch(cOrder) {
	case 0 :
		return (N/(N-1))*(N/(N-1))*phi*phi*phi +((3*N+1)/(N-1))*phi*phi*al/h 
		       +((3*N+2)/N)*phi*al*al/h/h;
	case 1 :
		return -((2*N+1)/(N-1))*phi*phi - (4*(N+1)/N)*phi*al/h;
    case 2 :
		return  ( (N+2)/N )*phi;
	default :
		return 0;
	}
}
double cLens1::BCoefficient(int cOrder,double phi,double N,double al,double h) {
	// "�����Y�݌v�@" eq.(4.29). Return the coefficient of front curvature c^cOrder.
	if(N==1) return 0;
	switch(cOrder) {
    case 0:
		return -(N/(N-1))*phi*phi -((2*N+1)/N)*phi*al/h;
	case 1:
		return ((N+1)/N)*phi;
	default:
		return 0;
	}
}
double cLens1::SACoefficient(int cOrder,double phi,double N,double al,double h) { 
	// "�����Y�݌v�@" eq.(4.30)
	double A2,A1,A0;
	A2=ACoefficient(2,phi,N,al,h); A1=ACoefficient(1,phi,N,al,h); A0=ACoefficient(0,phi,N,al,h);
	switch(cOrder) {
	case 2:  return h*h*h*h*A2;
	case 1:  return h*h*h*h*A1;
	case 0:  return h*h*h*h*A0;
	default: return 0;
	}
}
double cLens1::CMCoefficient(int cOrder,double phi,double N,double al,double h,double hp) { 
	// "�����Y�݌v�@" eq.(4.30)
	double A2,A1,A0,B1,B0;
	A2=ACoefficient(2,phi,N,al,h); A1=ACoefficient(1,phi,N,al,h); A0=ACoefficient(0,phi,N,al,h);
	B1=BCoefficient(1,phi,N,al,h); B0=BCoefficient(0,phi,N,al,h);
	switch(cOrder) {
	case 2:  return h*h*h*hp*A2;
	case 1:  return h*h*h*hp*A1 +h*h*B1;
	case 0:  return h*h*h*hp*A0 +h*h*B0;
	default: return 0;
	}
}
double cLens1::ASCoefficient(int cOrder,double phi,double N,double al,double h,double hp) { 
	// "�����Y�݌v�@" eq.(4.30)
	double A2,A1,A0,B1,B0;
	A2=ACoefficient(2,phi,N,al,h); A1=ACoefficient(1,phi,N,al,h); A0=ACoefficient(0,phi,N,al,h);
	B1=BCoefficient(1,phi,N,al,h); B0=BCoefficient(0,phi,N,al,h);
	switch(cOrder) {
	case 2:  return h*h*hp*hp*A2;
	case 1:  return h*h*hp*hp*A1 +2*h*hp*B1;
	case 0:  return h*h*hp*hp*A0 +2*h*hp*B0 +phi;
	default: return 0;
	}
}
double cLens1::DSCoefficient(int cOrder,double phi,double N,double al,double h,double hp) { 
	// "�����Y�݌v�@" eq.(4.30)
	double A2,A1,A0,B1,B0;
	A2=ACoefficient(2,phi,N,al,h); A1=ACoefficient(1,phi,N,al,h); A0=ACoefficient(0,phi,N,al,h);
	B1=BCoefficient(1,phi,N,al,h); B0=BCoefficient(0,phi,N,al,h);
	switch(cOrder) {
	case 2:  return h*hp*hp*hp*A2;
	case 1:  return h*hp*hp*hp*A1 +3*hp*hp*B1;
	case 0:  return h*hp*hp*hp*A0 +3*hp*hp*B0 +(hp/h)*(3*phi+phi/N);
	default: return 0;
	}
}

cLens1 cLens1::ThinDoubletOnNud2(double beta, 
	                             std::string color1,std::string color2,std::string color3, 
						         std::string glass1, double L,double B,double A, double nud2)
{
	//   ��2�����Y�ގ��̕��U��nud2
	//   beta = �{��
	//   L    = �c�F����(f=1)
	//   B   = �R�}�����W��(f=1 t=0)
	//   A   = ���ʎ����W��(f=1)
	// �̏����𖞂��������ڍ��_�u���b�g(f=1)���쐬����
	cLens1 x(3,3),x_pre(3,3);
	double wl_ref,wl_long,wl_short;
	double n1,nu1, nd2,n2,nu2;
	const double N1=1.3, N2=2;
	double p1,p2, al1,al2, c1,c2,c3, sa1,sa2,sa,sa_pre;

	if(beta==0) x.Set_s(LN*10); else x.Set_s(-1+1/beta);
	x.Set_color(1,color1); x.Set_color(2,color2); x.Set_color(3,color3);
	wl_ref=x.wl[1]; wl_long=x.wl[x.long_cn]; wl_short=x.wl[x.short_cn];
	n1=Re(Index(glass1,wl_ref));
	nu1=(n1-1)/( Re(Index(glass1,wl_short))-Re(Index(glass1,wl_long)) );
	for(nd2=N1; nd2<=N2; nd2+=0.001) {
		n2=cMaterial::HerzEq(nd2,nud2,wl_ref);
		nu2=(n2-1)/( cMaterial::HerzEq(nd2,nud2,wl_short)-cMaterial::HerzEq(nd2,nud2,wl_long) );

		p1=( 1/nu2-L)/(1/nu2-1/nu1);
		p2=(-1/nu1+L)/(1/nu2-1/nu1);
		al1=beta/(1-beta); al2=al1+p1;
		// �ȉ��g�����Y�݌v�@ �h(4.29) ���
		c1=( B -BCoefficient(0,p1,n1,al1,1) -BCoefficient(0,p2,n2,al2,1)
				+BCoefficient(1,p2,n2,al2,1)*p1/(n1-1) )
		   / ( BCoefficient(1,p1,n1,al1,1) +BCoefficient(1,p2,n2,al2,1) );
		c2=c1-p1/(n1-1);
		sa1= ACoefficient(0,p1,n1,al1,1)
			+ACoefficient(1,p1,n1,al1,1)*c1
			+ACoefficient(2,p1,n1,al1,1)*c1*c1;
		sa2= ACoefficient(0,p2,n2,al2,1)
			+ACoefficient(1,p2,n2,al2,1)*c2
			+ACoefficient(2,p2,n2,al2,1)*c2*c2;
		sa=sa1+sa2;

		if(c1==0) x.r(1)=0; else x.r(1)=1/c1;
		if(c2==0) x.r(2)=0; else x.r(2)=1/c2;
		c3=c2-p2/(n2-1);
		if(c3==0) x.r(3)=0; else x.r(3)=1/c3;
		x.Set_gname(1,glass1); x.Set_gname(2,GlassCode(nd2,nud2)); 
		if( nd2!=N1 && (sa-A)*(sa_pre-A)<0 ) {
			return fabs(sa_pre-A)<fabs(sa-A) ? x_pre : x;
		}
		sa_pre=sa; x_pre=x;
	}
	return cLens1(3,3);
}
double cLens1::ThinDoubletNd2OnNud2(double beta, 
	                                std::string color1,std::string color2,std::string color3, 
							        std::string glass1, double L,double B,double A, double nud2)
{
	cLens1 x=ThinDoubletOnNud2(beta, color1,color2,color3, glass1, L,B,A, nud2);
	return Re( Index(x.gname(2),Wavelength("d")) );
}
double cLens1::ThinDoubletSA5OnNud2(double beta, 
	                                std::string color1,std::string color2,std::string color3, 
							        std::string glass1, double L,double B,double A, double nud2)
{
	cLens1 x=ThinDoubletOnNud2(beta, color1,color2,color3, glass1, L,B,A, nud2);
	return x.Coefficient("SA5");
}
double cLens1::ThinDoubletLC2OnNud2(double beta, 
	                                std::string color1,std::string color2,std::string color3, 
							        std::string glass1, double L,double B,double A, double nud2)
{
	cLens1 x=ThinDoubletOnNud2(beta, color1,color2,color3, glass1, L,B,A, nud2);
	return x.Coefficient("LC2");
}

cLens1 cLens1::ThinDoublet(double beta, 
                           std::string color1,std::string color2,std::string color3, 
                           std::string glass1,std::string glass2, double L,double B)
{
	//   beta = �{��
	//   L    = �c�F����(f=1)
	//   CM   = �R�}�����W��(f=1) (���Ȃ킿�����W��B)
	// �̏����𖞂��������ڍ��_�u���b�g(f=1)���쐬����
	// �Е��̃����Y�̍ގ�����C�ȂǕ��U�̂Ȃ��Ƃ��C���̃����Y�̃p���[��0�Ƃ��C
	// ���Ȃ킿�P�����Y�ɂ��Čv�Z���邱�Ƃ��ł���D
	cLens1 x(3,3);
	double n1,dn1,nu1, n2,dn2,nu2;
	double p1,p2, al1,al2, c1,c2,c3;

	if(beta==0) x.Set_s(LN*10); else x.Set_s(-1+1/beta);
	x.Set_color(1,color1); x.Set_color(2,color2); x.Set_color(3,color3);
	x.Set_gname(1,glass1); x.Set_gname(2,glass2);

	n1=x.N(1,1);
	n2=x.N(2,1);
	dn1=x.N(1,x.short_cn)-x.N(1,x.long_cn);
	dn2=x.N(2,x.short_cn)-x.N(2,x.long_cn);
	
	if(n1==0 || n2==0){
		return x;
	}
	else if(dn1==0){
		p1=0; p2=1;
	}
	else if(dn2==0){
		p1=1; p2=0;
	}
	else if(dn1==0 && dn2==0){
		return x;
	}
	else{
		nu1=(n1-1)/dn1;
		nu2=(n2-1)/dn2;
		if(nu1==nu2){
			return x;
		}
		else{
			p1=( 1/nu2-L)/(1/nu2-1/nu1);
			p2=(-1/nu1+L)/(1/nu2-1/nu1);
		}
	}

	al1=beta/(1-beta); al2=al1+p1;
	
	// �ȉ��g�����Y�݌v�@ �h(4.29) ���
	c1=( B -BCoefficient(0,p1,n1,al1,1) -BCoefficient(0,p2,n2,al2,1)
		+BCoefficient(1,p2,n2,al2,1)*(p1/(n1-1)) )
	   / ( BCoefficient(1,p1,n1,al1,1) +BCoefficient(1,p2,n2,al2,1) );
	if(c1==0) x.r(1)=0; else x.r(1)=1/c1;
	c2=c1-(p1==0 ? 0 : p1/(n1-1));
	if(c2==0) x.r(2)=0; else x.r(2)=1/c2;
	c3=c2-(p2==0 ? 0 : p2/(n2-1));
	if(c3==0) x.r(3)=0; else x.r(3)=1/c3;
	return x;
}
double cLens1::ThinDoubletSA(double beta, 
                            std::string color1,std::string color2,std::string color3, 
                            std::string glass1,std::string glass2, double L,double B)
{
	cLens1 x=ThinDoublet(beta, color1,color2,color3, glass1,glass2, L,B);
	return x.Coefficient("SA");
}
double cLens1::ThinDoubletSA5(double beta, 
                             std::string color1,std::string color2,std::string color3, 
                             std::string glass1,std::string glass2, double L,double B)
{
	cLens1 x=ThinDoublet(beta, color1,color2,color3, glass1,glass2, L,B);
	return x.Coefficient("SA5");
}
double cLens1::ThinDoubletLC2(double beta, 
                             std::string color1,std::string color2,std::string color3, 
                             std::string glass1,std::string glass2, double L,double B)
{
	cLens1 x=ThinDoublet(beta, color1,color2,color3, glass1,glass2, L,B);
	return x.Coefficient("LC2");
}


cLens1 cLens1::ThinSeparateDoublet(double beta, 
                                   std::string color1,std::string color2,std::string color3, 
                                   std::string glass1,std::string glass2,
								   double L,double B,double A)
{
	//   beta = �{��
	//   L    = �c�F����(f=1)
	//   B    = �R�}�����W��(f=1 t=0)
	//   A    = ���ʎ����W��(f=1)
	// �̏����𖞂���������ڍ��_�u���b�g(f=1)���쐬����
	cLens1 x(4,3), x1,x2;
	double n1,dn1,nu1, n2,dn2,nu2;
	double p1,p2, al1,al2, c1,c2,c3,c4;
	double b11,b01,b12,b02, a21,a11,a01,a22,a12,a02;
	double E,F,G,H,I;

	if(beta==0) x.Set_s(LN*10); else x.Set_s(-1+1/beta);
	x.Set_color(1,color1); x.Set_color(2,color2); x.Set_color(3,color3);
	x.Set_gname(1,glass1); x.Set_gname(3,glass2); 

	n1=x.N(1,1); dn1=x.N(1,x.short_cn)-x.N(1,x.long_cn);
	n2=x.N(3,1); dn2=x.N(3,x.short_cn)-x.N(3,x.long_cn);
	
	if(n1==0 || n2==0 || dn1==0 || dn2==0){
		return x;
	}
	else{
		nu1=(n1-1)/dn1;
		nu2=(n2-1)/dn2;
		if(nu1==nu2){
			return x;
		}
		else{
			p1=( 1/nu2-L)/(1/nu2-1/nu1);
			p2=(-1/nu1+L)/(1/nu2-1/nu1);
		}
	}

	al1=beta/(1-beta); al2=al1+p1;
	
	// �ȉ��g�����Y�݌v�@ �h(4.29) ���

	// B1=b11*c1+b01
	// B2=b12*c3+b02
	// CM=B1+B2
	// -> c3=E*c1+F E=-b11/b12 F=(CM-b01-b02)/b12
	b11=BCoefficient(1,p1,n1,al1,1);
	b01=BCoefficient(0,p1,n1,al1,1);
	b12=BCoefficient(1,p2,n2,al2,1);
	b02=BCoefficient(0,p2,n2,al2,1);
	E=-b11/b12;
	F=(B-b01-b02)/b12;
	
	// A1=a21*c1^2+a11*c1+a01
	// A2=a22*c3^2+a12*c3+a02
	// SA=A1+A2
	// -> G*c1^2+H*c1+I=0
	//      G=a21+a22*E^2
	//      H=a11+2*a22*E*F+a12*E
	//      I=a01+a22*F^2+a12*F+a02-SA
	a21=ACoefficient(2,p1,n1,al1,1);
	a11=ACoefficient(1,p1,n1,al1,1);
	a01=ACoefficient(0,p1,n1,al1,1);
	a22=ACoefficient(2,p2,n2,al2,1);
	a12=ACoefficient(1,p2,n2,al2,1);
	a02=ACoefficient(0,p2,n2,al2,1);
	G=a21+a22*E*E;
	H=a11+2*a22*E*F+a12*E;
	I=a01+a22*F*F+a12*F+a02-A;

	if(H*H-4*G*I<0) return x;

	x1=x;
	c1=(-H+sqrt(H*H-4*G*I))/2/G;   if(c1==0) x1.r(1)=0; else x1.r(1)=1/c1;
	c2=c1-(p1==0 ? 0 : p1/(n1-1)); if(c2==0) x1.r(2)=0; else x1.r(2)=1/c2;
	c3=E*c1+F;                     if(c3==0) x1.r(3)=0; else x1.r(3)=1/c3;
	c4=c3-(p2==0 ? 0 : p2/(n2-1)); if(c4==0) x1.r(4)=0; else x1.r(4)=1/c4;

	x2=x;
	c1=(-H-sqrt(H*H-4*G*I))/2/G;   if(c1==0) x2.r(1)=0; else x2.r(1)=1/c1;
	c2=c1-(p1==0 ? 0 : p1/(n1-1)); if(c2==0) x2.r(2)=0; else x2.r(2)=1/c2;
	c3=E*c1+F;                     if(c3==0) x2.r(3)=0; else x2.r(3)=1/c3;
	c4=c3-(p2==0 ? 0 : p2/(n2-1)); if(c4==0) x2.r(4)=0; else x2.r(4)=1/c4;

	// SA5�̐�Βl�̏�������(�����炭���p��)��Ԃ�
	x1.CalcCoefficients();
	x2.CalcCoefficients();
	return fabs(x1.SA5t) > fabs(x2.SA5t) ? x2 : x1;
}
double cLens1::ThinSeparateDoubletSA5(double beta, 
                                      std::string color1,std::string color2,std::string color3, 
                                      std::string glass1,std::string glass2,
									  double L,double B,double A)
{
	cLens1 x=ThinSeparateDoublet(beta, color1,color2,color3, glass1,glass2, L,B,A);
	return x.Coefficient("SA5");
}


cLens1 cLens1::ThinTriplet(int FrontIsCemented, double beta,
	                     std::string color1,std::string color2,std::string color3, 
	                     std::string glass1,std::string glass2,std::string glass3,
                         double PhiFront, double L,double B,double A)
{
	double wl1,wl2,wl3, n1,n2,n3, nu1,nu2,nu3;
	double det, p1,p2,p3, al1,al2,al3, c1,c2,c3;
	cLens1 x(5,3), x1,x2;
	double A21,A11,A01,A22,A12,A02,A23,A13,A03, B11,B01,B12,B02,B13,B03;
	double E,F,G,H,I;
	
	if(beta==0) x.Set_s(LN*10); else x.Set_s(-1+1/beta);	
	x.Set_color(1,color1); x.Set_color(2,color2); x.Set_color(3,color3);
	if(FrontIsCemented){
		x.Set_gname(1,glass1); x.Set_gname(2,glass2); x.Set_gname(4,glass3);
	}
	else{
		x.Set_gname(1,glass1); x.Set_gname(3,glass2); x.Set_gname(4,glass3);
	}

	wl1=x.wl[1],wl2=x.wl[x.long_cn],wl3=x.wl[x.short_cn];
	n1=Re(Index(glass1,wl1)); n2=Re(Index(glass2,wl1)); n3=Re(Index(glass3,wl1));
	nu1=(n1-1)/( Re(Index(glass1,wl3))-Re(Index(glass1,wl2)) );
	nu2=(n2-1)/( Re(Index(glass2,wl3))-Re(Index(glass2,wl2)) );
	nu3=(n3-1)/( Re(Index(glass3,wl3))-Re(Index(glass3,wl2)) );

	if(FrontIsCemented){
		det=1/nu2-1/nu1;
		if( det==0 ) return cLens1(5,3);
		p1=( PhiFront/nu2-L+(1-PhiFront)/nu3 )/det;
		p2=PhiFront-p1;
		p3=1-p1-p2;	
	}
	else{
		det=1/nu3-1/nu2;
		if( det==0 ) return cLens1(5,3);
		p1=PhiFront;
		p2=( (1-PhiFront)/nu3-L+PhiFront/nu1 )/det;
		p3=1-p1-p2;
	}
	al1=beta/(1-beta); al2=al1+p1; al3=al2+p2;

	A21=ACoefficient(2,p1,n1,al1,1); A11=ACoefficient(1,p1,n1,al1,1); A01=ACoefficient(0,p1,n1,al1,1);
	A22=ACoefficient(2,p2,n2,al2,1); A12=ACoefficient(1,p2,n2,al2,1); A02=ACoefficient(0,p2,n2,al2,1);
	A23=ACoefficient(2,p3,n3,al3,1); A13=ACoefficient(1,p3,n3,al3,1); A03=ACoefficient(0,p3,n3,al3,1);
	B11=BCoefficient(1,p1,n1,al1,1); B01=BCoefficient(0,p1,n1,al1,1);
	B12=BCoefficient(1,p2,n2,al2,1); B02=BCoefficient(0,p2,n2,al2,1);
	B13=BCoefficient(1,p3,n3,al3,1); B03=BCoefficient(0,p3,n3,al3,1);

	if(FrontIsCemented){
		// c3=E*c1+F;
		E=(-B11-B12)/B13;
		F=(B-B01+B12*p1/(n1-1)-B02-B03)/B13;
		// A=G*c1*c1+H*c1+I
		G= A21 +A22 +A23*E*E;
		H= A11 -2*A22*p1/(n1-1) +A12 +2*A23*E*F +A13*E;
		I= A01 +A22*p1*p1/(n1-1)/(n1-1) -A12*p1/(n1-1) +A02 +A23*F*F +A13*F +A03;
		if( H*H-4*G*(I-A)<0 ) return cLens1(5,3);

		x1=x;
		c1=( -H+sqrt(H*H-4*G*(I-A)) )/2/G;
		c2=c1-p1/(n1-1);
		c3=E*c1+F;
		x1.r(1)= c1==0 ? 0 : 1/c1;
		x1.r(2)= c2==0 ? 0 : 1/c2;
		x1.r(3)= c2-p2/(n2-1)==0 ? 0 : 1/( c2-p2/(n2-1) );
		x1.r(4)= c3==0 ? 0 : 1/c3;
		x1.r(5)= c3-p3/(n3-1)==0 ? 0 : 1/( c3-p3/(n3-1) );

		x2=x;
		c1=( -H-sqrt(H*H-4*G*(I-A)) )/2/G;
		c2=c1-p1/(n1-1);
		c3=E*c1+F;
		x2.r(1)= c1==0 ? 0 : 1/c1;
		x2.r(2)= c2==0 ? 0 : 1/c2;
		x2.r(3)= c2-p2/(n2-1)==0 ? 0 : 1/( c2-p2/(n2-1) );
		x2.r(4)= c3==0 ? 0 : 1/c3;
		x2.r(5)= c3-p3/(n3-1)==0 ? 0 : 1/( c3-p3/(n3-1) );
	}
	else{
		// c1=E*c2+F
		E=(-B12-B13)/B11;
		F=(B-B01-B02+B13*p2/(n2-1)-B03)/B11;
		// A=G*c2*c2+H*c2+I
		G= A21*E*E +A22 +A23;
		H= 2*A21*E*F +A11*E +A12 -2*A23*p2/(n2-1) +A13;
		I= A21*F*F +A11*F +A01 +A02 +A23*p2*p2/(n2-1)/(n2-1) -A13*p2/(n2-1) +A03;
		if( H*H-4*G*(I-A)<0 ) return cLens1(5,3);

		x1=x;
		c2=( -H+sqrt(H*H-4*G*(I-A)) )/2/G;
		c1=E*c2+F;
		c3=c2-p2/(n2-1);
		x1.r(1)= c1==0 ? 0 : 1/c1;
		x1.r(2)= c1-p1/(n1-1)==0 ? 0 : 1/( c1-p1/(n1-1) );
		x1.r(3)= c2==0 ? 0 : 1/c2;
		x1.r(4)= c3==0 ? 0 : 1/c3;
		x1.r(5)= c3-p3/(n3-1)==0 ? 0 : 1/( c3-p3/(n3-1) );

		x2=x;
		c2=( -H-sqrt(H*H-4*G*(I-A)) )/2/G;
		c1=E*c2+F;
		c3=c2-p2/(n2-1);
		x2.r(1)= c1==0 ? 0 : 1/c1;
		x2.r(2)= c1-p1/(n1-1)==0 ? 0 : 1/( c1-p1/(n1-1) );
		x2.r(3)= c2==0 ? 0 : 1/c2;
		x2.r(4)= c3==0 ? 0 : 1/c3;
		x2.r(5)= c3-p3/(n3-1)==0 ? 0 : 1/( c3-p3/(n3-1) );
	}

	return fabs(x1.Coefficient("SA5")) < fabs(x2.Coefficient("SA5")) ? x1 : x2;
}
double cLens1::ThinTripletSA5(int FrontIsCemented, double beta,
	                         std::string color1,std::string color2,std::string color3, 
	                         std::string glass1,std::string glass2,std::string glass3,
	                         double PhiFront, double L,double B,double A)
{
	cLens1 x=ThinTriplet(FrontIsCemented, beta, 
	                    color1,color2,color3, glass1,glass2,glass3, PhiFront,L,B,A);
	return x.Coefficient("SA5");
}
double cLens1::ThinTripletLC2(int FrontIsCemented, double beta,
	                         std::string color1,std::string color2,std::string color3, 
	                         std::string glass1,std::string glass2,std::string glass3,
	                         double PhiFront, double L,double B,double A)
{
	cLens1 x=ThinTriplet(FrontIsCemented, beta, 
	                    color1,color2,color3, glass1,glass2,glass3, PhiFront,L,B,A);
	return x.Coefficient("LC2");
}

cLens1 cLens1::DoubleThinDoublet(double beta,double t,double e,
	                    std::string color1,std::string color2,std::string color3,
	                    std::string glass1,std::string glass2,
						std::string glass3,std::string glass4,
					    double Length, double L,double T,double CM,double SA,
						int SolutionNo)
{
	// beta = �{��
	// t = ���ʒu(f=1)
	// e = �O�ʂƌ�ʂ̋���(f=1)
	// color1-3 = ��g��1�y�ѐF�����v�Z�g��2,3
	// glass1-4 = �O�ʂ̍\�������Y1,2, ��ʂ̍\�������Y3,4 �̍ގ�
	// Length = �O�ʂ��瑜���œ_�܂ł̒���(���̖������̂Ƃ��̑S��)(f=1)
	// L,T,CM,SA = �����W���ڕW
	// SolutionNo = �Ԃ����̔ԍ�(�����̉���������)
	//     SolutionNo>=1 : �����𖞂�����
	//     SolutionNo==0 : �����𖞂�����������ꍇ�́C�ł�SA5����������
	//                     �Ȃ��ꍇ�́C�ł�SA���ڕW�ɋ߂����iSA�ȊO�̏����͖������j
	const int SOLUTIONS=10; // ���̕ۑ��p�z��̃T�C�Y�i�����炭2�ŗǂ��͂�)
	const double qSTART=-2, qEND=2;
	int i,j;
	cLens1 x(6,3),result[SOLUTIONS+1];
	double phiA,phiB;
	double delta;
	double alA,alB,hA,hB, alpA,alpB,hpA,hpB, alpA_t;
	double s, gamma;
	double n1,n2,n3,n4, nu1,nu2,nu3,nu4;
	matrix<double> A(4,4),B(4,1),X(4,1);
	double phi1,phi2,phi3,phi4;
	double al1,al2,al3,al4;
	double q3, c1,c2,c3,c4, c21,c41, c1m,c2m,c3m,c4m,c21m,c41m;
	double A1,A2,A3,B3,A4,B4;
	double SA_A,CM_A,SA_B,CM_B, sa,sa_pre=0;
	int Solutions=0;
	double min=1e60;

	x.Set_d(3,e);
	x.Set_gname(1,glass1); x.Set_gname(2,glass2); 
	x.Set_gname(4,glass3); x.Set_gname(5,glass4);
	x.Set_color(1,color1); x.Set_color(2,color2); x.Set_color(3,color3);
	x.Set_t(t);

	// f=1�Ƃ��ăp���[�z�u�̐ݒ�
	if(e!=0){
		Telephoto(phiA,phiB,Length,e,1);
	}
	else{
		return cLens1(6,3);
	}

	// Length==e+1 �܂��� Length=1 �̋߂��Ƃ��C�O�ʂ܂��͌�ʂ̃p���[���Ȃ��Ȃ�̂Ŕr��
	if( fabs(phiA)<0.001 || fabs(phiB)<0.001 ) return cLens1(6,3);

	// ���̋����̐ݒ�
	if(beta==0) s=LN*10; else s=phiB*e-1+1/beta;
	x.Set_s(s);

	// �ߎ����������l
	if(s>=LN) {
		alA=0;   hA=1;
		alpA=-1; hpA=0;  // ���ʒu0�̂Ƃ�
		alpA_t=-1;       // ���ʒut�̂Ƃ�
	}
	else {
		delta=phiB*e;    // ������_
		alA=1/(s-delta); hA=s/(s-delta);
		alpA=-(s-delta)/s; hpA=0;                                // ���ʒu0�̂Ƃ�   
		if(fabs(t)>=LN ) alpA_t=0; else alpA_t=-(s-delta)/(s-t); // ���ʒut�̂Ƃ�
	}
	
	// �ڕW�l��t=0�ł̒l�Ɋ��Z����
	if(s>=LN) gamma=-t; else gamma=-(alpA-alpA_t)/alA;
	T=T-gamma*L;       // "�����Y�݌v�@" (4.20)
	CM=CM-gamma*SA;    // "�����Y�݌v�@" (4.14)

	// �F�����ڕWL,T���C�O�ʁC��ʂ̍\�������Y�̃p���[ phi1,2,3,4 �����肷��D
	alB=alA+hA*phiA;    hB=hA-alB*e;
	alpB=alpA+hpA*phiA; hpB=hpA-alpB*e;
	n1=Re(Index(glass1,x.wl[1]));
	n2=Re(Index(glass2,x.wl[1]));
	n3=Re(Index(glass3,x.wl[1]));
	n4=Re(Index(glass4,x.wl[1]));
	nu1=(n1-1)/( Re(Index(glass1,x.wl[x.short_cn]))-Re(Index(glass1,x.wl[x.long_cn])) );
	nu2=(n2-1)/( Re(Index(glass2,x.wl[x.short_cn]))-Re(Index(glass2,x.wl[x.long_cn])) );
	nu3=(n3-1)/( Re(Index(glass3,x.wl[x.short_cn]))-Re(Index(glass3,x.wl[x.long_cn])) );
	nu4=(n4-1)/( Re(Index(glass4,x.wl[x.short_cn]))-Re(Index(glass4,x.wl[x.long_cn])) );
	A[1][1]=1;          A[1][2]=1;          A[1][3]=0;          A[1][4]=0;
	A[2][1]=0;          A[2][2]=0;          A[2][3]=1;          A[2][4]=1;
	A[3][1]=hA*hA/nu1;  A[3][2]=hA*hA/nu2;  A[3][3]=hB*hB/nu3;  A[3][4]=hB*hB/nu4;
	A[4][1]=hA*hpA/nu1; A[4][2]=hA*hpA/nu2; A[4][3]=hB*hpB/nu3; A[4][4]=hB*hpB/nu4;
	B[1][1]=phiA;
	B[2][1]=phiB;
	B[3][1]=L;
	B[4][1]=T;
	if(rank(A)<4) return cLens1(6,3); else X=inv(A)*B;
	phi1=X[1][1];
	phi2=X[2][1];
	phi3=X[3][1];
	phi4=X[4][1];

	// �����Y3�̌`����qq3�ƑS�nCM�̖ڕW��背���Y�`��͌��肷��D
	// q3�𓮂����Ď��s����ɂ��SA�̖ڕW�𖞂�������T���D
	al1=alA; al2=al1+hA*phi1; al3=alB; al4=al3+hB*phi3;
	for(q3=qSTART; q3<=qEND; q3+=0.001) {
		c3=phi3*(q3+1)/2/(n3-1); 
		c4=c3-phi3/(n3-1);
		A3= ACoefficient(2,phi3,n3,al3,hB)*c3*c3 +ACoefficient(1,phi3,n3,al3,hB)*c3
		   +ACoefficient(0,phi3,n3,al3,hB);
		B3= BCoefficient(1,phi3,n3,al3,hB)*c3 +BCoefficient(0,phi3,n3,al3,hB);
		A4= ACoefficient(2,phi4,n4,al4,hB)*c4*c4 +ACoefficient(1,phi4,n4,al4,hB)*c4
		   +ACoefficient(0,phi4,n4,al4,hB);
		B4= BCoefficient(1,phi4,n4,al4,hB)*c4 +BCoefficient(0,phi4,n4,al4,hB);
		SA_B=hB*hB*hB*hB*(A3+A4);
		CM_B=hB*hB*hB*hpB*(A3+A4) +hB*hB*(B3+B4);
		CM_A=CM-CM_B;
		c1=( CM_A/hA/hA -BCoefficient(0,phi1,n1,al1,hA) -BCoefficient(0,phi2,n2,al2,hA) 
			+BCoefficient(1,phi2,n2,al2,hA)*phi1/(n1-1) )
		   / ( BCoefficient(1,phi1,n1,al1,hA) +BCoefficient(1,phi2,n2,al2,hA) );
		c2=c1-phi1/(n1-1);
		A1= ACoefficient(2,phi1,n1,al1,hA)*c1*c1 +ACoefficient(1,phi1,n1,al1,hA)*c1
		   +ACoefficient(0,phi1,n1,al1,hA);
		A2= ACoefficient(2,phi2,n2,al2,hA)*c2*c2 +ACoefficient(1,phi2,n2,al2,hA)*c2
		   +ACoefficient(0,phi2,n2,al2,hA);
		SA_A=hA*hA*hA*hA*(A1+A2);
		sa=SA_A+SA_B;

		c21=c2-phi2/(n2-1);
		c41=c4-phi4/(n4-1);
		
		if( q3!=qSTART && (SA-sa)*(SA-sa_pre)<0 ) {
			x.Set_r(1, c1==0  ? 0 : 1/c1);
			x.Set_r(2, c2==0  ? 0 : 1/c2);
			x.Set_r(3, c21==0 ? 0 : 1/c21);
			x.Set_r(4, c3==0  ? 0 : 1/c3);
			x.Set_r(5, c4==0  ? 0 : 1/c4);
			x.Set_r(6, c41==0 ? 0 : 1/c41);
			if(Solutions<SOLUTIONS) { ++Solutions; result[Solutions]=x; }
		}
		sa_pre=sa;

		if( fabs(sa-SA)<min ){  // result[0]�Ƃ��čł�SA���ڕW�ɋ߂�����ۑ����邽�߂̏���
			min=fabs(sa-SA);
			// result[0]=x; �Ƃ���Ɣ��ɏ������Ԃ�������
			c1m=c1;
			c2m=c2;
			c3m=c3;
			c4m=c4;
			c21m=c21;
			c41m=c41;
		}
	}

	// SA5�����������ɕ��ׂ� (SolutionNo=0�ȊO)
	for(i=1; i<=Solutions-1; ++i) for(j=i+1; j<=Solutions; ++j) {
		if( result[i].Coefficient("SA5")<result[j].Coefficient("SA5") ) {
			Swap(result[i],result[j]);
		}
	}

	// SolutionNo=0�ɑ΂���߂�l�̐ݒ�
	if(Solutions==0){
		x.Set_r(1, c1m==0  ? 0 : 1/c1m);
		x.Set_r(2, c2m==0  ? 0 : 1/c2m);
		x.Set_r(3, c21m==0 ? 0 : 1/c21m);
		x.Set_r(4, c3m==0  ? 0 : 1/c3m);
		x.Set_r(5, c4m==0  ? 0 : 1/c4m);
		x.Set_r(6, c41m==0 ? 0 : 1/c41m);
		result[0]=x;
	}
	else{
		result[0]=result[1];
	}

	if( 0<=SolutionNo && SolutionNo<=Solutions ) return result[SolutionNo]; else return cLens1(6,3);
}

double cLens1::DoubleThinDoubletSA5(double beta,double t,double e,
	                               std::string color1,std::string color2,std::string color3,
	                               std::string glass1,std::string glass2,
						           std::string glass3,std::string glass4,
					               double Length, double L,double T,double CM,double SA,
								   int SolutionNo)
{
	cLens1 x=DoubleThinDoublet(beta,t,e, color1,color2,color3, glass1,glass2,glass3,glass4,
	                          Length, L,T,CM,SA, SolutionNo);
	return x.Coefficient("SA5");
}
double cLens1::DoubleThinDoubletAS(double beta,double t,double e,
	                              std::string color1,std::string color2,std::string color3,
	                              std::string glass1,std::string glass2,
						          std::string glass3,std::string glass4,
					              double Length, double L,double T,double CM,double SA,
								  int SolutionNo)
{
	cLens1 x=DoubleThinDoublet(beta,t,e, color1,color2,color3, glass1,glass2,glass3,glass4,
	                          Length, L,T,CM,SA, SolutionNo);
	return x.Coefficient("AS");
}
double cLens1::DoubleThinDoubletPT(double beta,double t,double e,
	                              std::string color1,std::string color2,std::string color3,
	                              std::string glass1,std::string glass2,
						          std::string glass3,std::string glass4,
					              double Length, double L,double T,double CM,double SA,
								  int SolutionNo)
{
	cLens1 x=DoubleThinDoublet(beta,t,e, color1,color2,color3, glass1,glass2,glass3,glass4,
	                          Length, L,T,CM,SA, SolutionNo);
	return x.Coefficient("PT");
}
double cLens1::DoubleThinDoubletDS(double beta,double t,double e,
	                              std::string color1,std::string color2,std::string color3,
	                              std::string glass1,std::string glass2,
						          std::string glass3,std::string glass4,
					              double Length, double L,double T,double CM,double SA,
								  int SolutionNo)
{
	cLens1 x=DoubleThinDoublet(beta,t,e, color1,color2,color3, glass1,glass2,glass3,glass4,
	                          Length, L,T,CM,SA, SolutionNo);
	return x.Coefficient("DS");
}


void cLens1::DoubleThinlens(double& A1,double& A2,double& B2, double& f1,double& f2, double& SigmaPhi,
							double& beta1,double& beta2, double& DS,
							double B1,
						    double Length,double eLRatio,double beta,double t,double SA,double CM,double AS)
{
	// Length,eLRatio(=e/L)�Ō��܂�C2���̔������Y����Ȃ�n(f=1),
	//      L<1       -> �O��A�� ���B��
	//      L>1 e>L-1 -> �O��A�� ���B��
	//      L>1 e<L-1 -> �O��A�� ���B��
	// �ɂ����āC�O�ʂ�B2��^�����Ƃ��CSA,CM,AS�̏����𖞂����C
	// �O�ʂ�A1, ��ʂ�A2,B2, �y�ъe�����Y�̏œ_����,�{���ƑS�n�̃��ӂ��v�Z����D
	//
	// A1,B1,A2,B2 = ��1����ё�2�����Y��A,B. 
	//               ������ "�e�����Y��" f=1�Ƃ����Ƃ��̒l�Ƃ��� (��̓I�`��݌v�̖ڕW�Ƃ��ĕ�����Ղ�����)�D
	//               �ŗL�W��A0,B0��Ώۂɂ��Ă��ǂ����C���̖������łȂ������Y�ł́C
	//               ���ۂ̎��������ʂ�������Â炢�D
	// f1,f2 = �e�����Y�̏œ_����(�S�nf=1)
	// SigmaPhi = ����(�S�nf=1)
	// beta1,beta2 = �e�����Y�̔{��
	// DS = �S�n�c�ȌW��(�S�nf=1, ���K��(1))
	// Length = �O�ʂ��瑜���œ_�܂ł̒���(���̖������̂Ƃ��̑S��)(�S�nf=1)
	// eLRatio = e/Length, e = �O�ʂƌ�ʂ̋���(�S�nf=1)
	// beta = �S�n�{��
	// t = ���ʒu(�S�nf=1)
	// SA,CM,AS = �S�n�����W���ڕW(�S�nf=1, ���K��(1))
	cLens1 x(2,1);
	double e, p1,p2;
	matrix<double> A(3,3),B(3,1),X(3,1);
	double h1,h2, hp1,hp2;
	double a1,b1,a2,b2;  // f=1�ɋK�i������Ă��Ȃ�A,B

	e=Length*eLRatio;

	x.Set_asph_type(1,IDEAL);
	x.Set_asph_type(2,IDEAL);
	x.Set_d(1,e);

	// �p���[�z�u�̐ݒ�	
	if(e!=0 && e!=Length){
		p1=1-(Length-1)/e;    // "�����_" (4.3.2)
		p2=(1-p1)/(1-e*p1);   // "�����_" (4.3.1)
		f1= p1==0 ? 0 : 1/p1;
		f2= p2==0 ? 0 : 1/p2;
		SigmaPhi=p1+p2;
	}
	else{
		A1=A2=B2=DS=f1=f2=beta1=beta2=0;
		return;
	}
	// Length==e+1 �܂��� Length=1 �̋߂��Ƃ��C�O�ʂ܂��͌�ʂ̃p���[���Ȃ��Ȃ�̂Ŕr��
	if( fabs(p1)<0.001 || fabs(p2)<0.001 ){
		A1=A2=B2=DS=beta1=beta2=0;
		return;
	}
	x.fideal(1)=1/p1;
	x.fideal(2)=1/p2;

	// ����,���̐ݒ�
	x.SetM(beta);
	x.Set_t(t);

	// �e�����Y�̔{���̎擾
	beta1=x.M(1,1);
	beta2=x.M(2,2);
	
	// �ߎ������̐ݒ�
	x.NormalizeUnit="f=1";
	x.NormalizeType=1;
	x.CalcCoefficients();
	h1=x.h[1];
	h2=x.h[2];
	hp1=x.hp[1];
	hp2=x.hp[2];
	
	// "�����Y�݌v�@�h(4.30)���C
	//   SA = A1*h1^4                        +A2*h2^4
	//   CM = A1*h1^3*hp1   +B1*h1^2         +A2*h2^3*hp2   +B2*h2^2
	//   AS = A1*h1^2*hp1^2 +B1*2*h1*hp1 +p1 +A2*h2^2*hp2^2 +B2*2*h2*hp2 +p2
	
	b1=B1*p1*p1;    // ���FB1�͑O�ʂ�f=1�ł̒l���^�����Ă���
	A[1][1]=h1*h1*h1*h1;   A[1][2]=h2*h2*h2*h2;   A[1][3]=0;          // X[1][1]=a1
	A[2][1]=h1*h1*h1*hp1;  A[2][2]=h2*h2*h2*hp2;  A[2][3]=h2*h2;      // X[2][1]=a2
	A[3][1]=h1*h1*hp1*hp1; A[3][2]=h2*h2*hp2*hp2; A[3][3]=2*h2*hp2;   // X[3][1]=b2
	B[1][1]=SA;
	B[2][1]=CM-b1*h1*h1;
	B[3][1]=AS-b1*2*h1*hp1-p1-p2;

	if(rank(A)==3){
		X=inv(A)*B;
		a1=X[1][1];
		a2=X[2][1];
		b2=X[3][1];
		DS= h1*hp1*hp1*hp1*a1+3*hp1*hp1*b1+(hp1/h1)*(3*p1+p1/1.5)
		   +h2*hp2*hp2*hp2*a2+3*hp2*hp2*b2+(hp2/h2)*(3*p2+p2/1.5);  // "�����Y�݌v�@" (4.30)
		// f=1�ł̎����W���ɕϊ�����
		A1=a1/(p1*p1*p1);
		A2=a2/(p2*p2*p2);
		B2=b2/(p2*p2);
	}
	else{
		A1=A2=B2=DS=0;
		return;
	}
}
double cLens1::DoubleThinlensA1(double B1,
	                            double Length,double eLRatio,double beta,double t,double SA,double CM,double AS)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,B1,Length,eLRatio,beta,t,SA,CM,AS);
	return A1;
}
double cLens1::DoubleThinlensA2(double B1,
	                            double Length,double eLRatio,double beta,double t,double SA,double CM,double AS)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,B1,Length,eLRatio,beta,t,SA,CM,AS);
	return A2;
}
double cLens1::DoubleThinlensB2(double B1,
	                            double Length,double eLRatio,double beta,double t,double SA,double CM,double AS)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,B1,Length,eLRatio,beta,t,SA,CM,AS);
	return B2;
}
double cLens1::DoubleThinlensDS(double B1,
	                            double Length,double eLRatio,double beta,double t,double SA,double CM,double AS)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,B1,Length,eLRatio,beta,t,SA,CM,AS);
	return DS;
}
double cLens1::DoubleThinlensF1(double Length,double eLRatio)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,0,Length,eLRatio,0,0,0,0,0);
	return f1;
}
double cLens1::DoubleThinlensF2(double Length,double eLRatio)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,0,Length,eLRatio,0,0,0,0,0);
	return f2;
}
double cLens1::DoubleThinlensSigmaPhi(double Length,double eLRatio)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,0,Length,eLRatio,0,0,0,0,0);
	return SigmaPhi;
}
double cLens1::DoubleThinlensBeta1(double Length,double eLRatio,double beta)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,0,Length,eLRatio,beta,0,0,0,0);
	return beta1;
}
double cLens1::DoubleThinlensBeta2(double Length,double eLRatio,double beta)
{
	double A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS;
	DoubleThinlens(A1,A2,B2,f1,f2,SigmaPhi,beta1,beta2,DS,0,Length,eLRatio,beta,0,0,0,0);
	return beta2;
}


cLens1 cLens1::Triplet(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
                     std::string color1,std::string color2,std::string color3,
                     std::string glass1,std::string glass2,std::string glass3, 
                     double CM,double AS,double DS, double f,double d1,double d3,double d5) 
{
	cLens1 x(6,3);
	double phi1,phi2,phi3, e1,e2, s, al1,al2,al3,h1,h2,h3,hp1,hp2,hp3, n1,n2,n3;
	{
		// determin power distribution
		double A,B,C,D, x,x1,dx, f,g;
		A=-Epsilon*Length/(1+Epsilon)*Epsilon*Length;                              // (4.4.1)
		B= Epsilon*Length*(1-Kappa) -Epsilon*Length/(1+Epsilon)*Epsilon*Kappa
		  +Epsilon*Epsilon*Length*Length/(1+Epsilon)/(1+Epsilon)*SigmaPhi;         // (4.4.1)
		C=-Epsilon*Length/(1+Epsilon)*(2-Kappa)*SigmaPhi +(2-Kappa)*Epsilon*Kappa; // (4.4.1)
		D= (1-Kappa)*SigmaPhi +(1+Epsilon)*Kappa*Kappa/Length -1;                  // (4.4.1)
		x=2; dx=1000;
		while( dx>0.0001){
			f=A*x*x*x +B*x*x +C*x +D;
			g=3*A*x*x +2*B*x +C;
			x1=x-f/g; dx=fabs(x1-x); x=x1;
		}
		phi3=x;
		e1=Length/(1+Epsilon); e2=Epsilon*e1;             // (4.4.2)
		phi1=(Kappa+e2*phi3)/e1; phi2=SigmaPhi-phi1-phi3; //(4.4.3)
	}
	{
		// calc s, al,h,hp
		matrix<double> R1(2,2),T1(2,2),R2(2,2),T2(2,2),R3(2,2), A(2,2), X(2,1);
		double f,delta, bf, alp1,alp2;
		R1.a[1][1]=1;    R1.a[1][2]=0;		T1.a[1][1]=1;    T1.a[1][2]=-e1;
		R1.a[2][1]=phi1; R1.a[2][2]=1;		T1.a[2][1]=0;    T1.a[2][2]=1;
		R2.a[1][1]=1;    R2.a[1][2]=0;		T2.a[1][1]=1;    T2.a[1][2]=-e2;
		R2.a[2][1]=phi2; R2.a[2][2]=1;		T2.a[2][1]=0;    T2.a[2][2]=1;
		R3.a[1][1]=1;    R3.a[1][2]=0;
		R3.a[2][1]=phi3; R3.a[2][2]=1;
		A=R3*T2*R2*T1*R1;
		f=A.a[2][1]; delta=(1-A.a[2][2])/A.a[2][1];
		bf=-f+delta;
		s= beta==0 ? LN*10 : bf+f/beta;
		if(fabs(s)>=LN){ al1=0;           h1=1;           alp1=-1;           hp1=0;}
		else           { al1=1/(s-delta); h1=s/(s-delta); alp1=-(s-delta)/s; hp1=0;}
		X.a[1][1]=h1;
		X.a[2][1]=al1;
		X=T1*R1*X;
		h2 =X.a[1][1];
		al2=X.a[2][1];
		X=T2*R2*X;
		h3 =X.a[1][1];
		al3=X.a[2][1];
		X.a[1][1]=hp1;
		X.a[2][1]=alp1;
		X=T1*R1*X;
		hp2 =X.a[1][1];
		alp2=X.a[2][1];
		X=T2*R2*X;
		hp3 =X.a[1][1];
	}
	{
		// calc n
		n1=Re(Index(glass1,Wavelength(color1)));
		n2=Re(Index(glass2,Wavelength(color1)));
		n3=Re(Index(glass3,Wavelength(color1)));
	}
	{   
		// calc shapes
		double c1,c2,c3,c11,c21,c31, a,b,c;
		double DS3,DS32,DS31,DS30,DS22,DS21,DS20;
		double AS32,AS31,AS30,AS22,AS21,AS20,as,as_pre;
		double CM11,CM10,CM22,CM21,CM20,CM32,CM31,CM30;
		DS32=DSCoefficient(2,phi3,n3,al3,h3,hp3); 
		DS31=DSCoefficient(1,phi3,n3,al3,h3,hp3); 
		DS30=DSCoefficient(0,phi3,n3,al3,h3,hp3);
		DS22=DSCoefficient(2,phi2,n2,al2,h2,hp2); 
		DS21=DSCoefficient(1,phi2,n2,al2,h2,hp2); 
		DS20=DSCoefficient(0,phi2,n2,al2,h2,hp2);
		AS32=ASCoefficient(2,phi3,n3,al3,h3,hp3);
		AS31=ASCoefficient(1,phi3,n3,al3,h3,hp3);
		AS30=ASCoefficient(0,phi3,n3,al3,h3,hp3);
		AS22=ASCoefficient(2,phi2,n2,al2,h2,hp2);
		AS21=ASCoefficient(1,phi2,n2,al2,h2,hp2);
		AS20=ASCoefficient(0,phi2,n2,al2,h2,hp2);
		CM11=CMCoefficient(1,phi1,n1,al1,h1,hp1);
		CM10=CMCoefficient(0,phi1,n1,al1,h1,hp1);
		CM22=CMCoefficient(2,phi2,n2,al2,h2,hp2);
		CM21=CMCoefficient(1,phi2,n2,al2,h2,hp2);
		CM20=CMCoefficient(0,phi2,n2,al2,h2,hp2);
		CM32=CMCoefficient(2,phi3,n3,al3,h3,hp3);
		CM31=CMCoefficient(1,phi3,n3,al3,h3,hp3);
		CM30=CMCoefficient(0,phi3,n3,al3,h3,hp3);
		for(c3=0; c3<10; c3+=0.001) {
			DS3=DS32*c3*c3+DS31*c3+DS30;
			a=DS22; b=DS21; c=DS20+DS3-DS;
			if(b*b-4*a*c>=0) c2=( -b-sqrt(b*b-4*a*c) )/2/a; else return cLens1(6,3);
			as=phi1 +AS22*c2*c2+AS21*c2+AS20 +AS32*c3*c3+AS31*c3+AS30;
			if(c3!=0 && (as-AS)*(as_pre-AS)<0) {
				c1=(CM-CM10-CM22*c2*c2-CM21*c2-CM20-CM32*c3*c3-CM31*c3-CM30)/CM11;
				x.s=s; x.t=0;
				x.Set_gname(1,glass1); x.Set_gname(3,glass2); x.Set_gname(5,glass3);
				x.Set_color(1,color1); x.Set_color(2,color2); x.Set_color(3,color3);
				x.r(1)= c1==0 ? 0 : 1/c1; c11=c1-phi1/(n1-1); x.r(2)= c11==0 ? 0 : 1/c11;
				x.r(3)= c2==0 ? 0 : 1/c2; c21=c2-phi2/(n2-1); x.r(4)= c21==0 ? 0 : 1/c21;
				x.r(5)= c3==0 ? 0 : 1/c3; c31=c3-phi3/(n3-1); x.r(6)= c31==0 ? 0 : 1/c31;
				x.d(2)=e1; x.d(4)=e2;
				{
					//   thicken   
					double f1,f2,f3;
					x.AdjustFocalLength(f,0);
					f1=x.f(1,2); f2=x.f(3,4); f3=x.f(5,6);
					x.d(1)=d1; x.d(3)=d3; x.d(5)=d5;
					x.AdjustFocalLength(f1,1,2,1,0,0); 
					x.AdjustFocalLength(f2,3,4,1,0,0); 
					x.AdjustFocalLength(f3,5,6,1,0,0);
					x.d(2)+=x.Trimed(1,2).delta1(); x.d(2)-=x.Trimed(3,4).delta();
					x.d(4)+=x.Trimed(3,4).delta1(); x.d(4)-=x.Trimed(5,6).delta();
					x.SetM(beta);
					x.t=x.Trimed(1,2).delta();
				}
				return x;
			}
			as_pre=as;
		}
	}
	return cLens1(6,3);
}
double cLens1::TripletSA(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                    std::string color1,std::string color2,std::string color3,
					    std::string glass1,std::string glass2,std::string glass3,
					    double CM,double AS,double DS, double f,double d1,double d3,double d5)
{
	cLens1 x=Triplet(beta, SigmaPhi,Length,Kappa,Epsilon, 
	                color1,color2,color3, glass1,glass2,glass3, CM,AS,DS, f,d1,d3,d5);
	return x.Coefficient("SA");
}
double cLens1::TripletSA5(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                     std::string color1,std::string color2,std::string color3,
					     std::string glass1,std::string glass2,std::string glass3,
					     double CM,double AS,double DS, double f,double d1,double d3,double d5)
{
	cLens1 x=Triplet(beta, SigmaPhi,Length,Kappa,Epsilon, 
	                color1,color2,color3, glass1,glass2,glass3, CM,AS,DS, f,d1,d3,d5);
	return x.Coefficient("SA5");
}
double cLens1::TripletCM(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                    std::string color1,std::string color2,std::string color3,
					    std::string glass1,std::string glass2,std::string glass3,
					    double CM,double AS,double DS, double f,double d1,double d3,double d5)
{
	cLens1 x=Triplet(beta, SigmaPhi,Length,Kappa,Epsilon, 
	                color1,color2,color3, glass1,glass2,glass3, CM,AS,DS, f,d1,d3,d5);
	return x.Coefficient("CM");
}
double cLens1::TripletAS(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                    std::string color1,std::string color2,std::string color3,
					    std::string glass1,std::string glass2,std::string glass3,
					    double CM,double AS,double DS, double f,double d1,double d3,double d5)
{
	cLens1 x=Triplet(beta, SigmaPhi,Length,Kappa,Epsilon, 
	                color1,color2,color3, glass1,glass2,glass3, CM,AS,DS, f,d1,d3,d5);
	return x.Coefficient("AS");
}
double cLens1::TripletDS(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                    std::string color1,std::string color2,std::string color3,
					    std::string glass1,std::string glass2,std::string glass3,
					    double CM,double AS,double DS, double f,double d1,double d3,double d5)
{
	cLens1 x=Triplet(beta, SigmaPhi,Length,Kappa,Epsilon, 
	                color1,color2,color3, glass1,glass2,glass3, CM,AS,DS, f,d1,d3,d5);
	return x.Coefficient("DS");
}
double cLens1::TripletLC(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                    std::string color1,std::string color2,std::string color3,
					    std::string glass1,std::string glass2,std::string glass3,
					    double CM,double AS,double DS, double f,double d1,double d3,double d5)
{
	cLens1 x=Triplet(beta, SigmaPhi,Length,Kappa,Epsilon, 
	                color1,color2,color3, glass1,glass2,glass3, CM,AS,DS, f,d1,d3,d5);
	return x.Coefficient("LC");
}
double cLens1::TripletTC(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                    std::string color1,std::string color2,std::string color3,
					    std::string glass1,std::string glass2,std::string glass3,
					    double CM,double AS,double DS, double f,double d1,double d3,double d5)
{
	cLens1 x=Triplet(beta, SigmaPhi,Length,Kappa,Epsilon, 
	                color1,color2,color3, glass1,glass2,glass3, CM,AS,DS, f,d1,d3,d5);
	return x.Coefficient("TC");
}

int cLens1::ZernikeL(int term_no,int IsFringeOrder,int jBase){
	return cZernike::lNumber(term_no,IsFringeOrder,jBase);
}
int cLens1::ZernikeN(int term_no,int IsFringeOrder,int jBase){
	return cZernike::nNumber(term_no,IsFringeOrder,jBase);
}

double cLens1::GaussBeamTruncatedPower
       (double BeamPhiX,double BeamPhiY,double AperturePhiX,double AperturePhiY,
        double ApertureDX,double ApertureDY,int n)
{
	// 1/e~2���a��BeamPhiX,BeamPhiY�ł���K�E�X�r�[���́C
	// ���a AperturePhiX,AperturePhiY, �ΐS ApertureDx,ApertureDy
	// �ł���ȉ~�J���ɂ������铧�ߗ����V�~�����[�V��������D
	// (��]�Ώ̂̃r�[���ƊJ���C�ΐS�Ȃ��ȊO��exact�Ȍv�Z�͕s�\�̂悤�Ȃ̂ŁC
	//  �{�֐����쐬�����D 07/06/07 )
	double mx,my;

	if(BeamPhiX>0 && BeamPhiY>0){
		if(BeamPhiX==BeamPhiY && AperturePhiX==AperturePhiY && ApertureDX==0 && ApertureDY==0){
			return 1-exp(-2*AperturePhiY*AperturePhiY/BeamPhiY/BeamPhiY);
		}
		else{
			cLens1 lens(1,1);
			lens.nSpot= n>1 ? n : 300;
			lens.s=10000000000;
			lens.Afocal=1;
			lens.GaussianPhiX=lens.GaussianPhiY=1;  // �r�[�������a1�̉~�ƂȂ�悤�ȍ��W�Ōv�Z����
			lens.EPD=2;  // exact�Ȍv�Z��1/e^2�a��2�{�̌a�Ɋ܂܂�鋭�x��99.97%
			mx=1/BeamPhiX;
			my=1/BeamPhiY;
			lens.EAx(1)=AperturePhiX*mx;
			lens.EAy(1)=AperturePhiY*my;
			lens.EPx=ApertureDX*mx;
			lens.EPy=ApertureDY*my;
			return lens.Ton();
		}
	}
	else{
		return 0;
	}
}

std::string cLens1::ZoomCamTable(double fv,double mo,
                                 double th1,double m_th1,double th2,double m_th2,
                                 double th_start,double th_end,double th_step,
                                 int variator_leading)
{
	// fv    = �o���G�[�^�œ_����
	// mo    = �o���G�[�^�ȊO�̔{���i�܂��͏œ_�����j
	// th1   = ����J���ʒu
	// m_th1 = th1�ł̑S�n�{���i�܂��͏œ_�����j
	// th2   = ����J���ʒu
	// m_th2 = th2�ł̑S�n�{���i�܂��͏œ_�����j
	// th_start/th_end/th_step = �\���쐬����ŏ��̃J���ʒu/�Ō�̃J���ʒu/�J���ʒu�Ԋu
	// variator_leading = �o���G�[�^�̕����R���y���Z�[�^����̂Ƃ��^
	//
	// �ȏ���o���G�[�^�ړ���x�C�R���y���Z�[�^�ړ��ʂ��̕\���쐬����Dth=0��x=y=0�Ƃ���D
	//     x�� z=f/m ���邢�� z'=-f*m ���D
	//     y�� z'-z = -f*(m+1/m) ���D

	std::string s;
	char buf[100];
	int i;
	double th, a,m,m0, x,y;

	// ���o�I���Ԋu�̂��߁C log(m/m_th1) = a*(th-th1)+b �Ƃ���D
	// (tips: log(m)=a*(th-th1)+b �ł� m<0 �������Ȃ��j
	//   th=th1�Ƃ���ƁCb=log(m_th1/m_th1)=0
	//   th=th2�Ƃ���ƁClon(m_th2/m_th1) = a*(th2-th1)
	// ���Ca���߂�D
	a=log(m_th2/m_th1)/(th2-th1);
	// th=0(x=y=0)�ł̔{���́C
	m0=m_th1*exp(a*(-th1));

	//       ####.### ####.### ####.#####  ####.#####
	sprintf(buf,"   th       m     x(variator) y(compensator)\n"); s+=buf;
	
	for(i=0; (th=th_start+th_step*i)<=th_end*1.000001; i++){
		m=m_th1*exp(a*(th-th1));
		if(variator_leading){
			y=-fv*(m/mo+mo/m) +fv*(m0/mo+mo/m0);
			x=-fv*mo/m +fv*mo/m0;
		}
		else{
			y=fv*(m/mo+mo/m) -fv*(m0/mo+mo/m0);
			x=fv*m/mo -fv*m0/mo;
		}
		sprintf(buf,"%8.3f %8.3f %10.3f %10.3f\n", th,m,x,y); s+=buf;
	}

	return s;
}

void cLens1::Telephoto(double &p1,double &p2,double L,double e,double f){
	//  L : �O�ʂ���œ_�܂�
    //  e : �O��ʂ̊Ԋu
	//  f : �S�n�œ_����
	//    L<1     -> �O�ʓ� ��ʉ�
	//    1<L<e+1 -> �O�ʓ� ��ʓ�
	//    e+1<L   -> �O�ʉ� ��ʓ�

	if(e!=0 && f!=0){
		L/=f;
		e/=f;
		p1=1-(L-1)/e;        // "�����_" (4.3.2)
		p2=(1-p1)/(1-e*p1);  // "�����_" (4.3.1)
		p1/=f;
		p2/=f;
	}
	else{
		p1=p2=0;
	}
}
double cLens1::TelephotoF1(double L,double e,double f){
	double p1,dummy;

	Telephoto(p1,dummy,L,e,f);
	return p1==0 ? 0 : 1/p1;
}
double cLens1::TelephotoF2(double L,double e,double f){
	double dummy,p2;

	Telephoto(dummy,p2,L,e,f);
	return p2==0 ? 0 : 1/p2;
}

void cLens1::OffAxialConicAB(double &a,double &b,double t_deg,double LongR,double ShortR){
	// �ȉ~���̒����Z�����a LongR,ShortR, ����уI�t�A�L�V�����񎟋Ȗʂ�t���C
	// �I�t�A�L�V�����񎟋Ȗʂ�a,b�����߂�D
	// �v�Z���@�͊􉽊w�I�֌W�C
	//     2*LongR  = a+b
	//     (2*F)^2 = a^2 +b^2 -2ab*cos(2t)  (�]���藝, 2*F�͏œ_�ԋ����j
	// �ɂ��D
	// �����Ȃ��Ƃ��� a=b=0 �Ƃ���D
	double x1,x2,F; 

	x1=cos(2*t_deg*PI/180)+1;
	F=sqrt(LongR*LongR-ShortR*ShortR);             // �œ_�Ԃ̋����̔���
	x2=LongR*LongR*x1*x1 -2*x1*(LongR*LongR-F*F);  // ������
	if(x2>=0){
		a=( LongR*x1 +sqrt(x2) ) / x1;
		b=2*LongR-a;
	}
	else{
		a=b=0;
	}
}

void cLens1::MartinEq(double &Sph,double &Cyl,double Diopter,double th_deg,double N,int OriginalEq/*=0*/){
	// �������Y�̌�����ɎΓ��˂��������ɉ�����Sph�����Cyl�l(diopter)���v�Z����D
	// ������g�}�[�`��(Martin)�̎��h�ł̋ߎ��i���ˊp���������j���g��Ȃ��Ōv�Z����D
	double sn,Dt,Ds;

	sn=sin(th_deg*PI/180);

	if(OriginalEq==0){
		// �ȉ��́C�g�􉽌��w�i�O��a�v�j�h��(8.29)(8.30)����
		// �������Y�̕\�ʁC���ʂɘA�p���邱�Ƃœ������D
		Dt=Diopter * (N * sqrt(1-sn*sn/N/N) - sqrt(1-sn*sn)) / (N-1) / (1-sn*sn);
		Ds=Diopter * (N * sqrt(1-sn*sn/N/N) - sqrt(1-sn*sn)) / (N-1);
	}
	else{
		// �I���W�i���̃}�[�`���̎��ł́C��̎���sn�ɂ��ăe�[���[�W�J���C
		// sn��2��܂ł��Ƃ�ߎ������Ă���D
		Dt = Diopter * (1+(N*2+1)/(N*2)*sn*sn);
		Ds = Diopter * (1+1/(N*2)*sn*sn);
	}

	Sph= Dt>=Ds ? Dt : Ds;
	Cyl= Dt>=Ds ? Ds-Dt : Dt-Ds;  // �}�C�i�X�\��
}

double cLens1::MartinSph(double Diopter,double th_deg,double N){
	double sph,cyl;

	MartinEq(sph,cyl,Diopter,th_deg,N);
	return sph;
}

double cLens1::MartinCyl(double Diopter,double th_deg,double N){
	double sph,cyl;

	MartinEq(sph,cyl,Diopter,th_deg,N);
	return cyl;
}

void cLens1::MartinEqGeneral(double &Sph,double &Cyl,double &Axis_deg,double Sph0,double Cyl0,double Axis0_deg,double N,
							 double rox_deg,double roy_deg,int order,int inverse/*=0*/)
{
	// ���ܓx�� Sph0,Cyl0,Axis0_deg�C�ޗ��̋��ܗ���N�̔������Y���C
	// rox_deg,roy_deg �����X�����Ƃ��̋��ܓx Sph,Cyl,Axis_deg ��Ԃ��D
	// inv==true �̂Ƃ��́C
	// �X������Ԃŋ��ܓx�� Sph0,Cyl0,Axis0_deg �̃����Y�̌X���Ȃ���Ԃł̋��ܓx
	// Sph,Cyl,Axis_deg ��Ԃ��D

	vector<double> ez(0,0,1),v;
	matrix<double> Rx(3,3),Ry(3,3),Ez(3,1),V(3,1);
	double Axis0,Axis,rox,roy,phi,th,thg, A,B,C;
	
	Axis0=Axis0_deg*PI/180;
	rox=rox_deg*PI/180;
	roy=roy_deg*PI/180;

	Ez[1][1]=ez.x;
	Ez[2][1]=ez.y;
	Ez[3][1]=ez.z;

	Rx[1][1]=1; Rx[1][2]=0;        Rx[1][3]=0;          // x������̉�]�s��
	Rx[2][1]=0; Rx[2][2]=cos(rox); Rx[2][3]=-sin(rox);
	Rx[3][1]=0; Rx[3][2]=sin(rox); Rx[3][3]=cos(rox);

	Ry[1][1]=cos(roy);  Ry[1][2]=0; Ry[1][3]=sin(roy);  // y������̉�]�s��
	Ry[2][1]=0;         Ry[2][2]=1; Ry[2][3]=0;
	Ry[3][1]=-sin(roy); Ry[3][2]=0; Ry[3][3]=cos(roy);

	if(order==0){
		V=inv(Rx)*inv(Ry)*Ez;  // �������Y�̌����̕���
	}
	else{
		V=inv(Ry)*inv(Rx)*Ez;
	}

	v.x=V[1][1]; v.y=V[2][1]; v.z=V[3][1];

	th=acos(sProduct(v,ez));              // �X�����������Y�̌�����z���̂Ȃ��p
	phi= v.x==0 ? PI/2 : atan(v.y/v.x);   // �X�����������Y�̌�����xy�ʂւ̓��e��x������Ƃ���p�x

	Axis0-=phi;

	A=(-Sph0-Cyl0*cos(Axis0)*cos(Axis0))/2;
	B=-0.5*Cyl0*sin(Axis0*2);
	C=(-Sph0-Cyl0*sin(Axis0)*sin(Axis0))/2;

	thg=asin(sin(th)/N);  // �������Y�����ł̓��ˊp

	if(inverse==0){
		A=A/cos(th)/cos(th)/cos(thg);
		B=B/cos(th)/cos(thg);
		C=C/cos(thg);
	}
	else{
		A=A*cos(th)*cos(th)*cos(thg);
		B=B*cos(th)*cos(thg);
		C=C*cos(thg);
	}

	Cyl=-2*sqrt((A-C)*(A-C)+B*B);
	Sph=-(A+C)-Cyl/2;
	Axis= Cyl==0 ? 0 : 0.5*atan2(-2*B/Cyl,-2*(A-C)/Cyl);
	if(Axis<0) Axis+=PI;
	Axis_deg=Axis*180/PI;
}

/////////////////////////////////////////////////////////////////////////////////////////

cLens1::cLens1(int i,int j){
	k=i; cn=j;
	create();
	initialize();
}

cLens1::cLens1(const cLens1& x){
	cOptics::operator =(x);

	k=x.k; cn=x.cn;
	create();
	initialize(); // ref_sphere=0;
	assignment(x);
}

cLens1::~cLens1(){ erase(); }

void cLens1::Reset() { *this=cLens1(); }

cLens1& cLens1::operator=(const cLens1& x) { 
	// �N���X������Z�q�̒�`�̋L�q���Ȃ���΁C�N���X�̑���ɂ����ăf�t�H���g�őS�Ẵ����o����������D
	// �������C�L�q�����ꍇ�͖������Ȃ���΂��̃����o�͑������Ȃ��̂ŁC
	// ��������������o�� assignment() �̒��Ȃǂɑ���𖾎����邱�ƁD

	cOptics::operator =(x);  // �璷�ȋC�����邪�C��{�N���X�̃����o��������̂ɕK�v�ȗl�q�i20130205)

	if(mems<x.k || cn<x.cn){
		// k�܂���cn���������菬�����ꍇ�̂݁C�������̍Ċm�ۂ��s�Ȃ��D
		// ���������ă��������s�����Ȃ���΃����o�̃A�h���X�͕ς�Ȃ��D
		// ����́C�Ⴆ�� optimize() �̃R�[�h���Ȍ��ɂ���̂ɕK�v�ł���D
		// (optimize()�ł́Cpop(),*this=buf_ini �Ȃ�*this�ւ̑��������)
		erase();
		k=x.k; cn=x.cn;
		create();
	}
	else{
		k=x.k; cn=x.cn;
	}
	assignment(x);
	return *this;
}

bool operator==(const cLens1& a,const cLens1& b){
	int i,j;
	bool x;
	x=(a.k==b.k)&&(a.cn==b.cn);
	if (x==false) return x;
	for(j=1; j<=a.cn; j++) {
		x=x&&(a.color[j]==b.color[j])&&(a.colorweight[j]==b.colorweight[j]);
	}
	x=x&&(a.s==b.s)&&(a.t==b.t);
	for(i=0; i<=a.k+1; i++){
		x=x&&(a.surf[i]==b.surf[i]);
		x=x&&(a.med[i]==b.med[i]);
	}
	x=x&&(a.s1fix==b.s1fix);
	x=x&&(a.yObjectMax==b.yObjectMax)&&(a.xObjectMax==b.xObjectMax);
	x=x&&(a.stop==b.stop);
	x=x&&(a.EPD==b.EPD); x=x&&(a.EPDx==b.EPDx); x=x&&(a.EPy==b.EPy); x=x&&(a.EPx==b.EPx);
	x=x&&(a.nSpot==b.nSpot);
	x=x&&(a.NormalizeUnit==b.NormalizeUnit);
	x=x&&(a.NormalizeType==b.NormalizeType);
	return x;
}

int cLens1::push() { Stack.push(*this); return 1; }
int cLens1::pop()  { return Stack.pop(*this);     }
int cLens1::undo() { return Stack.undo(*this);    }
int cLens1::redo() { return Stack.redo(*this);    }

std::string cLens1::Get_Note(){ return Note; }
void cLens1::Set_Note(std::string val){ Note=val; }

double cLens1::Get_FRW(){ return FRW; }
void cLens1::Set_FRW(double value){ FRW=value; }

int cLens1::Get_k(){ return k; }
void cLens1::Set_k(int val){
	cLens1 buf=*this;
	int temp=cn;

	erase();
	k=val; cn=temp;
	create();
	assignment(buf);
	// ����assignment()�̑����intialize()���g���ƁC
	// shapes[LENSVIEW].Clear()�����s����邱�ƂɂȂ�D
	// ����ƁC"Lens.xls"�ɂ����ă����Y�f�[�^�̕ύX�̍ۂɁC
	// �{�֐������s�ɂ��shapes[LENSVIEW]���N���A����C
	// MakeLensView()�ɂ����ăE�C���h�E�̌Œ肪�ł��Ȃ��Ȃ��Ă��܂��D
}

int cLens1::Get_cn(){ return cn; }
void cLens1::Set_cn(int val){
	cLens1 buf=*this;
	int k0=k;

	erase();
	k=k0; cn=val;
	create();
	assignment(buf);
	// ���������� intialize() ���g���ƁC
	// shapes[LENSVIEW].Clear()�����s����邱�ƂɂȂ�D
	// ����ƁC"Lens.xls"�ɂ����ă����Y�f�[�^�̕ύX�̍ۂɁC
	// �{�֐������s�ɂ��shapes[LENSVIEW]���N���A����C
	// MakeLensView()�ɂ����ăE�C���h�E�̌Œ肪�ł��Ȃ��Ȃ��Ă��܂��D
}

std::string cLens1::Get_color(int j){
	if( 1<=j && j<=cn ) return color[j]; else return "";
}
void cLens1::Set_color(int j, std::string s){
	int i;
	double a, max,min,min2;

	if( 1<=j && j<=cn ){
		color[j]=s;
		wl[j]=Wl(j);
		for(i=0; i<=k; ++i) med[i].g.set_wl(j,wl[j]);
	}

	// �F�����v�Z���ɎQ�Ƃ���long_cn,short_cn,mid_cn��ݒ肷��D
	// Set_color(0,"")�Ƃ���΂��̐ݒ肾���s�����Ƃ��ł���D
	max=0; min=min2=1e30; long_cn=short_cn=mid_cn=1;

	for(j=1; j<=cn; ++j){
		if(wl[j]>max){
			max=wl[j];
			long_cn=j;
		}
		if(wl[j]<min){
			min=wl[j];
			short_cn=j;
		}
	}

	// mid_cn = 1 (cn<=2)
    //          long_cn �ł�short_cn�ł��Ȃ��� (cn==3)
	//          �g���̋t���̓���long_cn��short_cn�̒����ɍł��߂�j (cn>3)
	for(j=1; j<=cn; ++j){
		if(j!=long_cn && j!=short_cn){
			a=fabs( (1/wl[long_cn]/wl[long_cn]+1/wl[short_cn]/wl[short_cn])/2 - 1/wl[j]/wl[j] );
			if(a<min2){
				min2=a;
				mid_cn=j;
				// �R�[�V�[�̕��U�� n-1=A(1+B/��^2)���C
				// dn�� 1/��(1)^2-1/��(2)^2 �ɔ�Ⴗ��ƍl������D
				// �F�����W����dn�ɔ�Ⴗ��D
			}
		}
	}
}

double cLens1::Get_colorweight(int j){ if( 1<=j && j<=cn ) return colorweight[j]; else return 0; }
void   cLens1::Set_colorweight(int j,double value){ if( 1<=j && j<=cn ) colorweight[j]=value; }

double cLens1::Get_s(){ return s; }
void   cLens1::Set_s(double value){ s=value; }

double cLens1::Get_t(){ return t; }
void   cLens1::Set_t(double value){ t=value; }

double cLens1::Get_r(int i){ if( 0<=i && i<=k+1 ) return r(i); else return 0; }
void   cLens1::Set_r(int i,double value){ if( 0<=i && i<=k+1 ) r(i)=value; }

double cLens1::Get_rObj(){ return rObj(); }
void   cLens1::Set_rObj(double value){ rObj()=value; }

double cLens1::Get_rImage(){ return rImage(); }
void   cLens1::Set_rImage(double value){ rImage()=value; }

double cLens1::Get_Newton(int i){ if( 0<=i && i<=k ) return Newton(i); else return 0; }
void   cLens1::Set_Newton(int i,double value){ if( 0<=i && i<=k ) Newton(i)=value; }

double cLens1::Get_As0(int i){ if( 0<=i && i<=k ) return As0(i); else return 0; }
void   cLens1::Set_As0(int i,double value){ if( 0<=i && i<=k ) As0(i)=value; }

double cLens1::Get_As45(int i){ if( 0<=i && i<=k ) return As45(i); else return 0; }
void   cLens1::Set_As45(int i,double value){ if( 0<=i && i<=k ) As45(i)=value; }

int  cLens1::Get_rVariable(int i){ if( 1<=i && i<=k) return rVariable(i); else return 0; }
void cLens1::Set_rVariable(int i,int value){ if( 1<=i && i<=k) rVariable(i)=value; }

int  cLens1::Get_EAtype(int i){ if( 0<=i && i<=k ) return EAtype(i); else return 0; }
void cLens1::Set_EAtype(int i,int value){ if( 0<=i && i<=k ) EAtype(i)=value; }

double cLens1::Get_EAy(int i){ if( 0<=i && i<=k ) return EAy(i); else return 0; }
void   cLens1::Set_EAy(int i,double value){ if( 0<=i && i<=k ) EAy(i)=value; }

double cLens1::Get_EAx(int i){ if( 0<=i && i<=k ) return EAx(i); else return 0; }
void   cLens1::Set_EAx(int i,double value){ if( 0<=i && i<=k ) EAx(i)=value; }

double cLens1::Get_CHMy(int i){ if( 0<=i && i<=k ) return CHMy(i); else return 0; }
void   cLens1::Set_CHMy(int i,double value){ if( 0<=i && i<=k ) CHMy(i)=value; }

double cLens1::Get_CHMx(int i){ if( 0<=i && i<=k ) return CHMx(i); else return 0; }
void   cLens1::Set_CHMx(int i,double value){ if( 0<=i && i<=k ) CHMx(i)=value; }

double cLens1::Get_EAdy(int i){ if( 0<=i && i<=k ) return EAdy(i); else return 0; }
void   cLens1::Set_EAdy(int i,double value){ if( 0<=i && i<=k ) EAdy(i)=value; }

double cLens1::Get_EAdx(int i){ if( 0<=i && i<=k ) return EAdx(i); else return 0; }
void   cLens1::Set_EAdx(int i,double value){ if( 0<=i && i<=k ) EAdx(i)=value; }

int  cLens1::Get_asph_type(int i){ if( 0<=i && i<=k+1 ) return asph_type(i); else return 0; }
void cLens1::Set_asph_type(int i,int value) { if( 0<=i && i<=k+1 ) asph_type(i)=value; }

double cLens1::Get_kp(int i){ if( 0<=i && i<=k+1 ) return kp(i); else return 0; }
void   cLens1::Set_kp(int i,double value){ if( 0<=i && i<=k+1 ) kp(i)=value; }

double cLens1::Get_NormH(int i){ if( 0<=i && i<=k+1 ) return NormH(i); else return 0; }
void   cLens1::Set_NormH(int i,double value){ if( 0<=i && i<=k+1 ) NormH(i)=value; }

double cLens1::Get_a1(int i){ if( 0<=i && i<=k+1 ) return a1(i); else return 0; }
void   cLens1::Set_a1(int i,double value){ if( 0<=i && i<=k+1 ) a1(i)=value; }

double cLens1::Get_a2(int i){ if( 0<=i && i<=k+1 ) return a2(i); else return 0; }
void   cLens1::Set_a2(int i,double value){ if( 0<=i && i<=k+1 ) a2(i)=value; }

double cLens1::Get_a3(int i){ if( 0<=i && i<=k+1 ) return a3(i); else return 0; }
void   cLens1::Set_a3(int i,double value){ if( 0<=i && i<=k+1 ) a3(i)=value; }

double cLens1::Get_a4(int i){ if( 0<=i && i<=k+1 ) return a4(i); else return 0; }
void   cLens1::Set_a4(int i,double value){ if( 0<=i && i<=k+1 ) a4(i)=value; }

double cLens1::Get_a5(int i){ if( 0<=i && i<=k+1 ) return a5(i); else return 0; }
void   cLens1::Set_a5(int i,double value){ if( 0<=i && i<=k+1 ) a5(i)=value; }

double cLens1::Get_a6(int i){ if( 0<=i && i<=k+1 ) return a6(i); else return 0; }
void   cLens1::Set_a6(int i,double value){ if( 0<=i && i<=k+1 ) a6(i)=value; }

double cLens1::Get_a7(int i){ if( 0<=i && i<=k+1 ) return a7(i); else return 0; }
void   cLens1::Set_a7(int i,double value){ if( 0<=i && i<=k+1 ) a7(i)=value; }

double cLens1::Get_a8(int i){ if( 0<=i && i<=k+1 ) return a8(i); else return 0; }
void   cLens1::Set_a8(int i,double value){ if( 0<=i && i<=k+1 ) a8(i)=value; }

double cLens1::Get_a9(int i){ if( 0<=i && i<=k+1 ) return a9(i); else return 0; }
void   cLens1::Set_a9(int i,double value){ if( 0<=i && i<=k+1 ) a9(i)=value; }

double cLens1::Get_a10(int i){ if( 0<=i && i<=k+1 ) return a10(i); else return 0; }
void   cLens1::Set_a10(int i,double value){ if( 0<=i && i<=k+1 ) a10(i)=value; }

double cLens1::Get_a11(int i){ if( 0<=i && i<=k+1 ) return a11(i); else return 0; }
void   cLens1::Set_a11(int i,double value){ if( 0<=i && i<=k+1 ) a11(i)=value; }

double cLens1::Get_a12(int i){ if( 0<=i && i<=k+1 ) return a12(i); else return 0; }
void   cLens1::Set_a12(int i,double value){ if( 0<=i && i<=k+1 ) a12(i)=value; }

double cLens1::Get_a13(int i){ if( 0<=i && i<=k+1 ) return a13(i); else return 0; }
void   cLens1::Set_a13(int i,double value){ if( 0<=i && i<=k+1 ) a13(i)=value; }

double cLens1::Get_a14(int i){ if( 0<=i && i<=k+1 ) return a14(i); else return 0; }
void   cLens1::Set_a14(int i,double value){ if( 0<=i && i<=k+1 ) a14(i)=value; }

double cLens1::Get_a15(int i){ if( 0<=i && i<=k+1 ) return a15(i); else return 0; }
void   cLens1::Set_a15(int i,double value){ if( 0<=i && i<=k+1 ) a15(i)=value; }

double cLens1::Get_a16(int i){ if( 0<=i && i<=k+1 ) return a16(i); else return 0; }
void   cLens1::Set_a16(int i,double value){ if( 0<=i && i<=k+1 ) a16(i)=value; }

double cLens1::Get_a18(int i){ if( 0<=i && i<=k+1 ) return a18(i); else return 0; }
void   cLens1::Set_a18(int i,double value){ if( 0<=i && i<=k+1 ) a18(i)=value; }

double cLens1::Get_a20(int i){ if( 0<=i && i<=k+1 ) return a20(i); else return 0; }
void   cLens1::Set_a20(int i,double value){ if( 0<=i && i<=k+1 ) a20(i)=value; }

double cLens1::Get_b20(int i){ if( 0<=i && i<=k+1 ) return b(i,2,0); else return 0; }
void   cLens1::Set_b20(int i,double value){ if( 0<=i && i<=k+1 ) b(i,2,0)=value; }

double cLens1::Get_b11(int i){ if( 0<=i && i<=k+1 ) return b(i,1,1); else return 0; }
void   cLens1::Set_b11(int i,double value){ if( 0<=i && i<=k+1 ) b(i,1,1)=value; }

double cLens1::Get_b02(int i){ if( 0<=i && i<=k+1 ) return b(i,0,2); else return 0; }
void   cLens1::Set_b02(int i,double value){ if( 0<=i && i<=k+1 ) b(i,0,2)=value; }

double cLens1::Get_b30(int i){ if( 0<=i && i<=k+1 ) return b(i,3,0); else return 0; }
void   cLens1::Set_b30(int i,double value){ if( 0<=i && i<=k+1 ) b(i,3,0)=value; }

double cLens1::Get_b21(int i){ if( 0<=i && i<=k+1 ) return b(i,2,1); else return 0; }
void   cLens1::Set_b21(int i,double value){ if( 0<=i && i<=k+1 ) b(i,2,1)=value; }

double cLens1::Get_b12(int i){ if( 0<=i && i<=k+1 ) return b(i,1,2); else return 0; }
void   cLens1::Set_b12(int i,double value){ if( 0<=i && i<=k+1 ) b(i,1,2)=value; }

double cLens1::Get_b03(int i){ if( 0<=i && i<=k+1 ) return b(i,0,3); else return 0; }
void   cLens1::Set_b03(int i,double value){ if( 0<=i && i<=k+1 ) b(i,0,3)=value; }

double cLens1::Get_b40(int i){ if( 0<=i && i<=k+1 ) return b(i,4,0); else return 0; }
void   cLens1::Set_b40(int i,double value){ if( 0<=i && i<=k+1 ) b(i,4,0)=value; }

double cLens1::Get_b31(int i){ if( 0<=i && i<=k+1 ) return b(i,3,1); else return 0; }
void   cLens1::Set_b31(int i,double value){ if( 0<=i && i<=k+1 ) b(i,3,1)=value; }

double cLens1::Get_b22(int i){ if( 0<=i && i<=k+1 ) return b(i,2,2); else return 0; }
void   cLens1::Set_b22(int i,double value){ if( 0<=i && i<=k+1 ) b(i,2,2)=value; }

double cLens1::Get_b13(int i){ if( 0<=i && i<=k+1 ) return b(i,1,3); else return 0; }
void   cLens1::Set_b13(int i,double value){ if( 0<=i && i<=k+1 ) b(i,1,3)=value; }

double cLens1::Get_b04(int i){ if( 0<=i && i<=k+1 ) return b(i,0,4); else return 0; }
void   cLens1::Set_b04(int i,double value){ if( 0<=i && i<=k+1 ) b(i,0,4)=value; }

double cLens1::Get_b(int i,int m,int n){ if( 0<=i && i<=k+1 ) return b(i,m,n); else return 0; }
void   cLens1::Set_b(int i,int m,int n,double value){ if( 0<=i && i<=k+1 ) b(i,m,n)=value; }

std::string cLens1::Get_bTerms(int i){
	if( 0<=i && i<=k+1 ) return surf[i].b.GetTerms(); else return "";
}
void cLens1::Set_bTerms(int i,std::string s){
	if( 0<=i && i<=k+1 ) surf[i].b.SetTerms(s);
}

cZernike cLens1::Get_ZSurface(int i){
	return this->zernike(i);
}
void cLens1::Set_ZSurface(int i, cZernike zernike){
	this->zernike(i)=zernike;
}

std::string cLens1::Get_ZCoefficients(int i){
	if( 0<=i && i<=k+1) return surf[i].zernike.GetCoefficients(); else return "";
}
void cLens1::Set_ZCoefficients(int i,std::string s){
	if( 0<=i && i<=k+1 ) surf[i].zernike.SetCoefficients(s);
}

double cLens1::Get_ZernikeR0(int i){
	return zernike(i).GetR0();
}
void cLens1::Set_ZernikeR0(int i,double value){
	zernike(i).SetR0(value);
}

int cLens1::Get_ZernikeMaxOrder(int i){
	return zernike(i).GetMaxOrder();
}
void cLens1::Set_ZernikeMaxOrder(int i,int value){
	zernike(i).SetMaxOrder(value);
}

std::string cLens1::Get_DconTerms(int i){
	if( 0<=i && i<=k+1 ) return surf[i].Dcon.GetTerms(); else return "";
}
void cLens1::Set_DconTerms(int i,std::string s){
	if( 0<=i && i<=k+1 ) surf[i].Dcon.SetTerms(s);
}

std::string cLens1::Get_LCoefficients(int i){
	if( 0<=i && i<=k+1) return surf[i].legendre.GetCoefficients(); else return "";
}
void cLens1::Set_LCoefficients(int i,std::string s){
	if( 0<=i && i<=k+1 ) surf[i].legendre.SetCoefficients(s);
}

double cLens1::Get_LegendreR0(int i){
	return legendre(i).R0;
}
void cLens1::Set_LegendreR0(int i,double value){
	legendre(i).R0=value;
}

std::string cLens1::Get_SplineData(int i){
	if(surf[i].spline.NotNull()){
		return surf[i].spline.GetData();
	}
	else{
		return "";
	}
}
void cLens1::Set_SplineData(int i,std::string s){
	surf[i].spline.SetData(s);
}

double cLens1::Get_Apv(int i){ if( 0<=i && i<=k+1 ) return Apv(i); else return 0; }
void   cLens1::Set_Apv(int i,double value){ if( 0<=i && i<=k+1 ) Apv(i)=value; }

double cLens1::Get_sfx(int i){ if( 0<=i && i<=k+1 ) return sfx(i); else return 0; }
void   cLens1::Set_sfx(int i,double value){ if( 0<=i && i<=k+1 ) sfx(i)=value; }

double cLens1::Get_sfy(int i){ if( 0<=i && i<=k+1 ) return sfy(i); else return 0; }
void   cLens1::Set_sfy(int i,double value){ if( 0<=i && i<=k+1 ) sfy(i)=value; }

int  cLens1::Get_UserDefSurf(int i){ if( 0<=i && i<=k+1 ) return UserDefSurf(i); else return 0; }
void cLens1::Set_UserDefSurf(int i,int value){ if( 0<=i && i<=k+1 ) UserDefSurf(i)=value; }

int  cLens1::Get_cylinder(int i){ if( 0<=i && i<=k+1 ) return cylinder(i); else return 0; }
void cLens1::Set_cylinder(int i,int value){ if( 0<=i && i<=k+1 ) cylinder(i)=value; }

double cLens1::Get_ry(int i){ if( 0<=i && i<=k+1 ) return ry(i); else return 0; }
void   cLens1::Set_ry(int i,double value){ if( 0<=i && i<=k+1 ) ry(i)=value; }

double cLens1::Get_rx(int i){ if( 0<=i && i<=k+1 ) return rx(i); else return 0; }
void   cLens1::Set_rx(int i,double value){ if( 0<=i && i<=k+1 ) rx(i)=value; }
	
double cLens1::Get_kpy(int i){ if( 0<=i && i<=k+1 ) return kpy(i); else return 0; }
void   cLens1::Set_kpy(int i,double value){ if( 0<=i && i<=k+1 ) kpy(i)=value; }

double cLens1::Get_kpx(int i){ if( 0<=i && i<=k+1 ) return kpx(i); else return 0; }
void   cLens1::Set_kpx(int i,double value){ if( 0<=i && i<=k+1 ) kpx(i)=value; }

int    cLens1::Get_IsXToroid(int i){ if( 0<=i && i<=k+1 ) return IsXToroid(i); else return 0; }
void   cLens1::Set_IsXToroid(int i,int value){ if( 0<=i && i<=k+1 ) IsXToroid(i)=value; }

double cLens1::Get_fideal(int i){ if( 0<=i && i<=k+1 ) return fideal(i); else return 0; }
void   cLens1::Set_fideal(int i,double value){ if( 0<=i && i<=k+1 ) fideal(i)=value; }

double cLens1::Get_aCOA(int i){ if( 0<=i && i<=k+1 ) return aCOA(i); else return 0; }
void   cLens1::Set_aCOA(int i,double value){ if( 0<=i && i<=k+1 ) aCOA(i)=value; }

double cLens1::Get_bCOA(int i){ if( 0<=i && i<=k+1 ) return bCOA(i); else return 0; }
void   cLens1::Set_bCOA(int i,double value){ if( 0<=i && i<=k+1 ) bCOA(i)=value; }

double cLens1::Get_tCOA(int i){ if( 0<=i && i<=k+1 ) return tCOA(i); else return 0; }
void   cLens1::Set_tCOA(int i,double value){ if( 0<=i && i<=k+1 ) tCOA(i)=value; }

double cLens1::Get_SA0(int i){ if( 0<=i && i<=k+1 ) return SA0(i); else return 0; }
void   cLens1::Set_SA0(int i,double value){ if( 0<=i && i<=k+1 ) SA0(i)=value; }

double cLens1::Get_CM0(int i){ if( 0<=i && i<=k+1 ) return CM0(i); else return 0; }
void   cLens1::Set_CM0(int i,double value){ if( 0<=i && i<=k+1 ) CM0(i)=value; }

int  cLens1::Get_grating(int i){ if(0<=i && i<=k+1) return grating(i); else return 0; }
void cLens1::Set_grating(int i,int value){ if(0<=i && i<=k+1) grating(i)=value; }

int  cLens1::Get_difforder(int i){ if(0<=i && i<=k+1) return difforder(i); else return 0; }
void cLens1::Set_difforder(int i,int value){ if(0<=i && i<=k+1) difforder(i)=value; }

double cLens1::Get_gpitch(int i){ if(0<=i && i<=k+1) return gpitch(i); else return 0; }
void   cLens1::Set_gpitch(int i,double value){ if(0<=i && i<=k+1) gpitch(i)=value; }

double cLens1::Get_grx(int i){ if(0<=i && i<=k+1) return grx(i); else return 0; }
void   cLens1::Set_grx(int i,double value){ if(0<=i && i<=k+1) grx(i)=value; }

double cLens1::Get_gry(int i){ if(0<=i && i<=k+1) return gry(i); else return 0; }
void   cLens1::Set_gry(int i,double value){ if(0<=i && i<=k+1) gry(i)=value; }

double cLens1::Get_grz(int i){ if(0<=i && i<=k+1) return grz(i); else return 0; }
void   cLens1::Set_grz(int i,double value){ if(0<=i && i<=k+1) grz(i)=value; }

double cLens1::Get_Diffusion(int i){ if(0<=i && i<=k+1) return Diffusion(i); else return 0; }
void   cLens1::Set_Diffusion(int i,double value){ if(0<=i && i<=k+1) Diffusion(i)=value; }

int  cLens1::Get_Fish(int i){ if(0<=i && i<=k+1) return Fish(i); else return 0; }
void cLens1::Set_Fish(int i,int value){ if(0<=i && i<=k+1) Fish(i)=value; }

int  cLens1::Get_decenter_type(int i){ if( 0<=i && i<=k+1 ) return decenter_type(i); else return 0; }
void cLens1::Set_decenter_type(int i,int value){ if( 0<=i && i<=k+1 ) decenter_type(i)=value; }

double cLens1::Get_dx(int i){ if( 0<=i && i<=k+1 ) return dx(i); else return 0; }
void   cLens1::Set_dx(int i,double value){ if( 0<=i && i<=k+1 ) dx(i)=value; }

double cLens1::Get_dy(int i){ if( 0<=i && i<=k+1 ) return dy(i); else return 0; }
void   cLens1::Set_dy(int i,double value){ if( 0<=i && i<=k+1 ) dy(i)=value; }

double cLens1::Get_dz(int i){ if( 0<=i && i<=k+1 ) return dz(i); else return 0; }
void   cLens1::Set_dz(int i,double value){ if( 0<=i && i<=k+1 ) dz(i)=value; }

double cLens1::Get_rox(int i){ if( 0<=i && i<=k+1 ) return rox(i); else return 0; }
void   cLens1::Set_rox(int i,double value){ if( 0<=i && i<=k+1 ) rox(i)=value; }

double cLens1::Get_roy(int i){ if( 0<=i && i<=k+1 ) return roy(i); else return 0; }
void   cLens1::Set_roy(int i,double value){ if( 0<=i && i<=k+1 ) roy(i)=value; }

double cLens1::Get_roz(int i){ if( 0<=i && i<=k+1 ) return roz(i); else return 0; }
void   cLens1::Set_roz(int i,double value){ if( 0<=i && i<=k+1 ) roz(i)=value; }

int  cLens1::Get_order(int i){ if( 0<=i && i<=k+1 ) return order(i); else return 0; }
void cLens1::Set_order(int i,int value){ if( 0<=i && i<=k+1 ) order(i)=value; }

int  cLens1::Get_ret(int i){ if( 0<=i && i<=k+1 ) return ret(i); else return 0; }
void cLens1::Set_ret(int i,int value){ if( 0<=i && i<=k+1 ) ret(i)=value; }

double cLens1::Get_dx1(int i){ if( 0<=i && i<=k+1 ) return dx1(i); else return 0; }
void   cLens1::Set_dx1(int i,double value){ if( 0<=i && i<=k+1 ) dx1(i)=value; }

double cLens1::Get_dy1(int i){ if( 0<=i && i<=k+1 ) return dy1(i); else return 0; }
void   cLens1::Set_dy1(int i,double value){ if( 0<=i && i<=k+1 ) dy1(i)=value; }

double cLens1::Get_dz1(int i){ if( 0<=i && i<=k+1 ) return dz1(i); else return 0; }
void   cLens1::Set_dz1(int i,double value){ if( 0<=i && i<=k+1 ) dz1(i)=value; }

double cLens1::Get_rox1(int i){ if( 0<=i && i<=k+1 ) return rox1(i); else return 0; }
void   cLens1::Set_rox1(int i,double value){ if( 0<=i && i<=k+1 ) rox1(i)=value; }

double cLens1::Get_roy1(int i){ if( 0<=i && i<=k+1 ) return roy1(i); else return 0; }
void   cLens1::Set_roy1(int i,double value){ if( 0<=i && i<=k+1 ) roy1(i)=value; }

double cLens1::Get_roz1(int i){ if( 0<=i && i<=k+1 ) return roz1(i); else return 0; }
void   cLens1::Set_roz1(int i,double value){ if( 0<=i && i<=k+1 ) roz1(i)=value; }

int  cLens1::Get_order1(int i){ if( 0<=i && i<=k+1 ) return order1(i); else return 0; }
void cLens1::Set_order1(int i,int value){ if( 0<=i && i<=k+1 ) order1(i)=value; }

std::string cLens1::Get_CoatName(int i){ if( 1<=i && i<=k ) return CoatName(i); else return ""; }
void        cLens1::Set_CoatName(int i,std::string filename){ if( 1<=i && i<=k ) CoatName(i)=filename; }

int  cLens1::Get_CoatReverse(int i){ if( 1<=i && i<=k ) return CoatReverse(i); else return 0; }	
void cLens1::Set_CoatReverse(int i,int value){ if( 1<=i && i<=k ) CoatReverse(i)=value; }

std::string cLens1::Get_rem(int i){ if( 1<=i && i<=k ) return rem(i); else return ""; }
void        cLens1::Set_rem(int i,std::string value){ if( 1<=i && i<=k ) rem(i)=trim(value,0); }

double cLens1::Get_d(int i){ if( 0<=i && i<=k-1 ) return d(i); else return 0; }
void   cLens1::Set_d(int i,double value){ if( 0<=i && i<=k-1 ) d(i)=value; }

double cLens1::Get_delta_d(int i){ if( 0<=i && i<=k-1 ) return delta_d(i); else return 0; }
void   cLens1::Set_delta_d(int i,double value){ if( 0<=i && i<=k-1 ) delta_d(i)=value; }

int  cLens1::Get_dVariable(int i){ if( 1<=i && i<=k-1) return dVariable(i); else return 0; }
void cLens1::Set_dVariable(int i,int value){ if( 1<=i && i<=k-1) dVariable(i)=value; }

std::string cLens1::Get_gname(int i){
	if( 0<=i && i<=k+1) return gname(i); else return "";
}
void cLens1::Set_gname(int i,std::string s){
	if( 0<=i && i<=k+1){
		gname(i)= s=="" ? "1" : s;   // s����̂Ƃ��͋�C�Ƃ݂Ȃ�
	}
}

int  cLens1::Get_gVariable(int i){ if( 1<=i && i<=k-1) return gVariable(i); else return 0; }
void cLens1::Set_gVariable(int i,int value){ if( 1<=i && i<=k-1) gVariable(i)=value; }

double cLens1::Get_N(int i,int j){ if( 0<=i && i<=k+1 && 1<=j && j<=cn ) return N(i,j); else return 0; }

double cLens1::Get_s1fix(){ return s1fix; }
void   cLens1::Set_s1fix(double value){ s1fix=value; }

std::string cLens1::Get_AfocalMode(){
	switch(AfocalMode){
		case MIN:
			return "min";
		case RAD:
			return "rad";
		case FUNDUS:
			return "fundus";
		case VA:
			return "va";
	}
	return "";
}
void cLens1::Set_AfocalMode(std::string value){
	if(value=="min"    || value=="MIN"    ) AfocalMode=MIN;
	if(value=="rad"    || value=="RAD"    ) AfocalMode=RAD;
	if(value=="fundus" || value=="FUNDUS" ) AfocalMode=FUNDUS;
	if(value=="va"     || value=="VA"     ) AfocalMode=VA;
}

double cLens1::Get_yObjectMax(){ return yObjectMax; }
void   cLens1::Set_yObjectMax(double value){ yObjectMax=value; }

double cLens1::Get_xObjectMax(){ return xObjectMax; }
void   cLens1::Set_xObjectMax(double value){ xObjectMax=value; }

double cLens1::Get_yObjectMaxAng(){
	if(abs(s)<LN){
		return -atan(yObjectMax/(s-(abs(t)<LN ? t : 0)))*180/PI;		
	}
	else{
		return atan(yObjectMax)*180/PI;		
	}
}
void cLens1::Set_yObjectMaxAng(double value){
	if(abs(s)<LN){
		yObjectMax=-(s-(abs(t)<LN ? t : 0))*tan(value*PI/180);
	}
	else{
		yObjectMax=tan(value*PI/180);		
	}
}

double cLens1::Get_xObjectMaxAng(){
	if(abs(s)<LN){
		return -atan(xObjectMax/(s-(abs(t)<LN ? t : 0)))*180/PI;		
	}
	else{
		return atan(xObjectMax)*180/PI;		
	}
}
void cLens1::Set_xObjectMaxAng(double value){
	if(abs(s)<LN){
		xObjectMax=-(s-(abs(t)<LN ? t : 0))*tan(value*PI/180);		
	}
	else{
		xObjectMax=tan(value*PI/180);		
	}
}

int  cLens1::Get_stop(){ return stop; }
void cLens1::Set_stop(int i){ if(0<=i && i<=k) stop=i; }  // if(1<=i ...) �Ƃ����stop=0(�i��ʎw��Ȃ�)���ݒ�ł��Ȃ�

double cLens1::Get_EPD(){ return EPD; }
void   cLens1::Set_EPD(double value){ EPD=value; }
double cLens1::Get_EPDx(){ return EPDx; }
void   cLens1::Set_EPDx(double value){ EPDx=value; }
double cLens1::Get_EPy(){ return EPy; }
void   cLens1::Set_EPy(double value){ EPy=value; }
double cLens1::Get_EPx(){ return EPx; }
void   cLens1::Set_EPx(double value){ EPx=value; }

int  cLens1::Get_nSpot(){ return nSpot; }
void cLens1::Set_nSpot(int value){ if(value>=1) nSpot=value; }

double cLens1::Get_SourcePhiY(){ return SourcePhiY; }       
void   cLens1::Set_SourcePhiY(double value){ SourcePhiY=value; }
double cLens1::Get_SourcePhiX(){ return SourcePhiX; }       
void   cLens1::Set_SourcePhiX(double value){ SourcePhiX=value; }
double cLens1::Get_SourceAxis(){ return SourceAxis; }       
void   cLens1::Set_SourceAxis(double value){ SourceAxis=value; }
int    cLens1::Get_SourceAreaType(){ return SourceAreaType; }       
void   cLens1::Set_SourceAreaType(int value){ SourceAreaType=value; }

double cLens1::Get_GrinDs(){ return GrinDs; }
void   cLens1::Set_GrinDs(double value){ if(value>0) GrinDs=value; }
int    cLens1::Get_GrinToListStep(){ return GrinToListStep; }
void   cLens1::Set_GrinToListStep(int value){ if(value>0) GrinToListStep=value; }

int cLens1::Get_ExcludeVirtualRay(){ return ExcludeVirtualRay; }
void cLens1::Set_ExcludeVirtualRay(int value){ ExcludeVirtualRay=value; }
int cLens1::Get_ExcludeVirtualObject(){ return ExcludeVirtualObject; }
void cLens1::Set_ExcludeVirtualObject(int value){ ExcludeVirtualObject=value; }

int cLens1::SurfNo(std::string rem){
	// rem�����߂ƈ�v����Ƃ��C���̖ʔԍ���Ԃ��D
	// ����ȊO�̂Ƃ��́Catoi(rem.c_str()) ��Ԃ��D
	// ���F rem����������("0","1","2", .. "9","+","-",".")�Ŏn�܂�Ƃ��͒��߂ƈ�v���Ă��Ă�
	//      ���������łȂ������̒��O�܂ł̕����񂪕\�����l���Ԃ����D
	//      (���������āC���߂͐�����������n�߂Ȃ������ǂ��j
	this->is_numeric(rem);
	return atoi(rem.c_str());
}

double cLens1::AfocalImageUnit(){
	if(Afocal){
		switch(AfocalMode){
			case MIN:
				return (PI/180)/60; // 1min=(PI/180)/60rad
			case RAD:
				return 1;
			case FUNDUS:
				return 1.0/17;
			case VA:
				// �t��(���g��)�����͂ƂȂ�
				// ����1�ɑΉ����镪��\��2min =(PI/180)/60*2rad
				return (PI/180)/60*2;
		}
		return 1;
	}
	else{
		return 1;
	}
}

std::string cLens1::AfocalImageUnitStr(){
	if(Afocal){
		switch(AfocalMode){
			case MIN:
				return "'";
			case RAD:
				return "rad";
			case FUNDUS:
				return " rad/17";
			case VA:
				return "[2min]";
		}
		return "";
	}
	else{
		return "";
	}
}

std::string cLens1::AfocalFreqUnitStr(){
	if(Afocal){
		switch(AfocalMode){
			case MIN:
				return "/min";
			case RAD:
				return "/rad";
			case FUNDUS:
				return "[17/rad]";
			case VA:
				return "VA";
		}
		return "";
	}
	else{
		return "";
	}
}


double cLens1::tCalc() {
	// �i��ʔԍ�stop���t��ݒ肷��. EPphi�͕ς��Ȃ��D
	if(1<=stop && stop<=k){
		AMatrixCalc(1,stop,1);
		t=N(0,1)*(-A2[2]/A2[1]);
		return t;
	}
	else return 0;
}

double cLens1::EPDCalc() {
	// �i��ʔԍ�stop���EPD��ݒ肷��Dt�͕ς��Ȃ��D
	if(1<=stop && stop<=k){
		AMatrixCalc(1,stop,1);
		EPD=fabs( (-A2[2]*A2[3]/A2[1]+A2[4])*EAy(stop) );
		if(fabs(t)>=LN){
			EPD/=fabs(t);  // RayTrace() �ł� |t|>=LN �ɑ΂��鏈���ɑΉ�
		}
		EPDx=0;
		return EPD;
	}
	else return 0;
}

void cLens1::EPCalc(){
	// �i��ʔԍ�stop���t��EPD��ݒ肷��
	tCalc();
	EPDCalc();
}

double cLens1::Mstop(){
	// ���˓�/�i�� ��Ԃ�(��Βl)
	if(1<=stop && stop<=k){
		AMatrixCalc(1,stop,1);
		return fabs( -A2[2]*A2[3]/A2[1]+A2[4] );
	}
	else return 0;
}
double cLens1::Mstop1(){
	// �ˏo��/�i�� ��Ԃ�(��Βl)
	if(1<=stop && stop<=k){
		AMatrixCalc(stop,k,1);
		return fabs( 1/A2[4] );
	}
	else return 0;
}

double cLens1::ApodizationAmplitude(){
	// RayTrace()���s����̂ݗL���D
	// OPD_calc()�Ȃǂ�ref_sphere���RayTrace()�����s����̂ŁC
	// ref_sphere->ApodizationAmplitude()���ĂԂ��ƁD
	double a, x,y, rx,ry;

	if(ApodizationSurf==0){
		// ���˓����W�ɂ��
		x=this->x_pupil;
		y=this->y_pupil;
	}
	else if(1<=ApodizationSurf && ApodizationSurf<=k){
		// �ʍ��W�ɂ��
		x=this->x[ApodizationSurf];
		y=this->y[ApodizationSurf];
	}
	else{
		return 1;
	}

	a=1;
	if(GaussianPhiY!=0){
		// �U������=GaussianPhiY(X)��1/e�ɂȂ镪�z�D
		// ���Ȃ킿�C���x����=GaussianPhiY(X)��1/e^2�ɂȂ镪�z�D
		ry=GaussianPhiY;
		rx=GaussianPhiX;
		if(rx==0) rx=ry;
		a*=exp( -(y*y)/(ry*ry/4) -(x*x)/(rx*rx/4) );
	}
	// �ȉ��CLambert�ʂ̌v�Z�������֑g�ݍ��ޗ\��(06.09.14) /////////////////////

	return a;
}
double cLens1::ApodizationIntensity(){
	double a;
	a=ApodizationAmplitude();
	return a*a;
}

double cLens1::Wl(int j){
	// ��j�g����nm�P�ʂŕԂ�
	if(1<=j && j<<cn) return Wavelength(color[j]);
	else              return 0;
}

double cLens1::power(int i1, int i2,int j/*=1*/){
	if(1<=i1 && i1<=i2 && i2<=k){
		AMatrixCalc(i1,i2,j);
		return A2[3];
	}
	else return 0;
}
double cLens1::power(int j/*=1*/){
	return power(1,k,j);
}

double cLens1::dpower(int i1, int i2,int j/*=1*/){
	// power=1�ɋK�i������power�̐F������Ԃ��D
	// LC�̂悤�ɔ{���ɂ���ĕς�Ȃ����߁C�����Y���̂̐F������c�����₷���D
	// �����ނ� dp/p<<1 �ł���̂ŁC�����݌v�̃E�G�C�g��ݒ肵�₷���悤�Ƀp�[�Z���g�ŕԂ��D
	double p,dp;

	p=power(i1,i2,1);

	switch(cn){
	case 1:
		dp=0;
		break;
	case 2:
		dp= p==0 ? 0 : (power(i1,i2,2)-power(i1,i2,1))/p*100;
		break;
	default:
		// j=1����g��, j=2��{���g���[|�Z�g���[}, j=cn��{�Z�g���[|���g���[} �̕��т����҂���D
		dp= p==0 ? 0 : (power(i1,i2,cn)-power(i1,i2,2))/p*100;
		break;
	}
	return dp;
}
double cLens1::dpower(int j/*=1*/){
	return dpower(1,k,j);
}

double cLens1::f(int i1, int i2,int j/*=1*/){
	if(1<=i1 && i1<=i2 && i2<=k){
		AMatrixCalc(i1,i2,j);
		return 1/A2[3];
	}
	else return 0;
}
double cLens1::f(int j/*=1*/){
	return f(1,k,j);
}

double cLens1::bf(int i1,int i2,int j/*=1*/){
	if(1<=i1 && i1<=i2 && i2<=k){
		AMatrixCalc(i1,i2,j);
		return A2[3]!=0 ? N(i2,j)*A2[1]/A2[3] : 0;
	}
	else return 0;
}
double cLens1::bf(int j/*=1*/){
	return bf(1,k,j);
}

double cLens1::ff(int i1,int i2,int j/*=1*/){
	if(1<=i1 && i1<=i2 && i2<=k){
		AMatrixCalc(i1,i2,j);
		return A2[3]!=0 ? -N(i1-1,j)*A2[4]/A2[3] : 0;
	}
	else return 0;
}
double cLens1::ff(int j/*=1*/){
	return ff(1,k,j);
}

double cLens1::bfRatio(int i1,int i2){
	return (bf(i1,i2)/N(i2,1))/f(i1,i2);
}
double cLens1::bfRatio(){
	return bfRatio(1,k);
}

double cLens1::ffRatio(int i1,int i2){
	return (ff(i1,i2)/N(i1-1,1))/f(i1,i2);
}
double cLens1::ffRatio(){
	return ffRatio(1,k);
}

double cLens1::BFOverF2(int i1,int i2){
	// bf/f^2 ��Ԃ��DSLO,FC���̑Ε������Y�ŁC�Ε����ˋP�_����������Ƃ��C
	// ���ԑ��ʂ��Ε������Y�ʂ���ǂ̂��炢����Ă��邩��]������D
	double bf,f;

	bf=this->bf(i1,i2);
	f=this->f(i1,i2);
	return (bf/N(i2,1))/(f*f);
}
double cLens1::BFOverF2(){
	return BFOverF2(1,k);
}

double cLens1::FFOverF2(int i1,int i2){
	// ff/f^2 ��Ԃ��DSLO,FC���̑Ε������Y�ŁC�Ε����ˋP�_����������Ƃ��C
	// ���ԑ��ʂ��Ε������Y�ʂ���ǂ̂��炢����Ă��邩��]������D
	double ff,f;

	ff=this->ff(i1,i2);
	f=this->f(i1,i2);
	return (ff/N(i1-1,1))/(f*f);
}
double cLens1::FFOverF2(){
	return FFOverF2(1,k);
}

double cLens1::delta(int i1,int i2,int j/*=1*/){
	if(1<=i1 && i1<=i2 && i2<=k) return ff(i1,i2,j)+N(i1-1,j)*f(i1,i2,j);
	else return 0;
}
double cLens1::delta(int j/*=1*/){
	return delta(1,k,j);
}

double cLens1::delta1(int i1,int i2,int j/*=1*/){
	if(1<=i1 && i1<=i2 && i2<=k) return bf(i1,i2,j)-N(i2,j)*f(i1,i2,j);
	else return 0;
}
double cLens1::delta1(int j/*=1*/){
	return delta1(1,k,j);
}

double cLens1::g_hat(int j/*=1*/){
	return s-delta(j);
}

double cLens1::g1_hat(int j/*=1*/){
	return s1()-delta1(j);
}

double cLens1::g1(int j/*=1*/){
	// �ˏo������ߎ����ʂ܂ł̋���
	return stop==0 ? 0 : s1()-s1(vertex(stop),j);
}

double cLens1::nodal(int i1,int i2,int j/*=1*/){
	// �g�����Y�݌v�@�h (2.33b)
	if(1<=i1 && i1<=i2 && i2<=k && power(i1,i2,j)!=0){
		return ff(i1,i2,j)+fabs(N(i2,j))/sgn(N(i1-1,j))*f(i1,i2,j);
	}
	else return 0;
}
double cLens1::nodal(int j/*=1*/){
	return nodal(1,k,j);
}

double cLens1::nodal1(int i1,int i2,int j/*=1*/){
	// �g�����Y�݌v�@�h (2.33a)
	if(1<=i1 && i1<=i2 && i2<=k && power(i1,i2,j)!=0){
		return bf(i1,i2,j)-fabs(N(i1-1,j))/sgn(N(i2,j))*f(i1,i2,j);
	}
	else return 0;
}
double cLens1::nodal1(int j/*=1*/){
	return nodal1(1,k,j);
}

double cLens1::cc(int i){
	// i�ʋȗ����S(i�ʔ��� -1�{����_)�̕��̋�Ԃł̑��̑�1�ʂ���Ƃ���z�����ʒu�D
	// ���̓_��ʂ����i�ʂŔ��˂���Ɠ������H�ł��̓_�ɖ߂�D
	double h1,al1,h,al;
	
	if(1<=i && i<=k){
		h1=1;
		al1= r(i)==0 ? 0 : h1/r(i)*N(i,1);
		AMatrixCalc(1,i);
		InvertAMatrix();
		h =A2[1]*h1+A2[2]*al1;
		al=A2[3]*h1+A2[4]*al1;
		return al==0 ? LN : h/(al/N(0,1));
	}
	else{
		return 0;
	}
}

double cLens1::cc1(int i){
	// i�ʋȗ����S�̑���Ԃł̍ŏI�ʂ���Ƃ���ʒu�D
	double h1,al1,h,al;
	
	if(1<=i && i<=k){
		h=1;
		al= r(i)==0 ? 0 : h/r(i)*N(i-1,1);
		AMatrixCalc(i,k);
		h1 =A2[1]*h+A2[2]*al;
		al1=A2[3]*h+A2[4]*al;
		return al1==0 ? LN : h1/(al1/N(k,1));
	}
	else{
		return 0;
	}
}

double cLens1::vertex(int i,double si/*=0*/){
	// i�ʒ��_(i�ʔ��ˁE���܂̎啽��)�̕��̋�Ԃł̑��̑�1�ʂ���Ƃ���z�����ʒu�D
	// ���̖ʏ�̓_��ʂ����i�ʂŔ��˂���Ƃ��̓_�ɖ߂�D
	// si��i�ʂ̕��̋����i�������j
	double h1,al1,h,al;
	
	if(1<=i && i<=k){
		if(si==0){
			h1=0;
			al1=1;
		}
		else{
			h1=1;
			al1=h1/si*N(i-1,1);
		}
		AMatrixCalc(1,i);
		InvertAMatrix();
		h =A2[1]*h1+A2[2]*al1;
		al=A2[3]*h1+A2[4]*al1;
		return al==0 ? LN : h/(al/N(0,1));
	}
	else{
		return 0;
	}
}

double cLens1::vertex1(int i,double si/*=0*/){
	// i�ʒ��_�̑���Ԃł̍ŏI�ʂ���Ƃ���ʒu�D
	double h1,al1,h,al;
	
	if(1<=i && i<=k){
		if(si==0){
			h=0;
			al=1;
		}
		else{
			h=1;
			al=h/si*N(i-1,1);
		}
		AMatrixCalc(i,k);
		h1 =A2[1]*h+A2[2]*al;
		al1=A2[3]*h+A2[4]*al;
		return al1==0 ? LN : h1/(al1/N(k,1));
	}
	else{
		return 0;
	}
}

double cLens1::deadspace(int i1,int i2,int j/*=1*/){
	return TotalThickness(i1,i2)-delta(i1,i2)+delta1(i1,i2);
}
double cLens1::deadspace(int j){
	return deadspace(1,k,j);
}

double cLens1::M(double s,int i1,int i2,int j){
	double al,h, al_in,al_out;

	if( 1<=i1 && i1<=i2 && i2<=k ){
		h= s==0 ? 0 : 1;
		al= fabs(s)>=LN ? 0 : (s==0 ? 1 : N(0,j)*h/s);
		AMatrixCalc(1,i1-1,j);
		al_in=A2[3]*h+A2[4]*al;
		AMatrixCalc(1,i2,j);
		al_out=A2[3]*h+A2[4]*al;
		return al_in/al_out;
	}
	else return 0;
}
double cLens1::M(int i1,int i2,int j/*=1*/){
	return M(this->s,i1,i2,j);
}
double cLens1::M(int j/*=1*/){
	return M(this->s,1,k,j);
}

void cLens1::SetM(double value){
	if(fabs(value)==0) s=LN*10;
	else {
		s=ff()+f()*N(0,1)/value;
	}
}

double cLens1::Mpupil(int i1,int i2,int j/*=1*/){
	double al,h, al_in,al_out;
	if( 1<=i1 && i1<=i2 && i2<=k ){
		h= t==0 ? 0 : 1;
		al= fabs(t)>=LN ? 0 : (t==0 ? 1 : N(0,j)*h/t);
		AMatrixCalc(1,i1-1,j);
		al_in=A2[3]*h+A2[4]*al;
		AMatrixCalc(1,i2,j);
		al_out=A2[3]*h+A2[4]*al;
		return al_in/al_out;
	}
	else return 0;
}
double cLens1::Mpupil(int j/*=1*/){
	return Mpupil(1,k,j);
}

double cLens1::gamma(int i1, int i2,int j/*=1*/){
	// �g�����Y�݌v�@�h(2.40) afocal�n�̊p�{��
	if(1<=i1 && i1<=i2 && i2<=k){
		AMatrixCalc(i1,i2,j);
		return A2[1];
	}
	else return 0;
}
double cLens1::gamma(int j/*=1*/){
	return gamma(1,k,j);
}

double cLens1::si(double s,int i,int j){
	// ��i�ʑO�̕��̈ʒu�����߂�
	if( i<1 || this->k<i ) return 0;
	if(i==1){
		return s;
	}
	else{
		return s1i(s,i-1,j)-d(i-1);
	}
}

double cLens1::si(int i,int j/*=1*/){
	return si(this->s,i,j);
}

double cLens1::s1i(double s,int i,int j){
	// ��i�ʌ�̑��ʒu�����߂�
	double al,h, al1,h1;

	if( i<1 || this->k<i ) return 0;
	h= s==0 ? 0 : 1;
	al= fabs(s)>=LN ? 0 : (s==0 ? 1 : N(0,j)*h/s);
	AMatrixCalc(1,i,j);
	h1=A2[1]*h+A2[2]*al;
	al1=A2[3]*h+A2[4]*al;
	return al1==0 ? LN : N(i,j)*h1/al1;
}

double cLens1::s1i(int i,int j/*=1*/){
	return s1i(this->s,i,j);
}

double cLens1::s1(double s,int j){
	return s1i(s,this->k,j);
}

double cLens1::s1(int j/*=1*/){
	return s1(this->s,j);
}

double cLens1::vdiopter1i(double s,int i,int j){
	// ��i�ʌ�̒��_���ܗ�[diopter]�����߂�
	// SCA()�Ɠ������C�������̂Ƃ����Ƃ���D
	double s1;

	s1=s1i(s,i,j);
	return s1==0 ? 0 : 1000/s1*sgn(N(i,1));
}

double cLens1::vdiopter1i(int i,int j/*=1*/){
	return vdiopter1i(this->s,i,j);
}

double cLens1::LCPar(int from,int to){
	return Afocal ? 1000/s1(to)-1000/s1(from) : s1(to)-s1(from);
}

double cLens1::t1(int j/*=1*/){
	return s1(this->t,j);
}

double cLens1::ti(int i,int j/*=1*/){
	// ��i�ʑO�̓��ʒu�����߂�
	if( i<1 || this->k<i ) return 0;
	if(i==1){
		return t;
	}
	else{
		return t1i(i-1,j)-d(i-1);
	}
}

double cLens1::t1i(int i,int j/*=1*/){
	// ��i�ʌ�̓��ʒu�����߂�
	double al,h, al1,h1;

	if( i<1 || this->k<i ) return 0;
	h= t==0 ? 0 : 1;
	al= fabs(t)>=LN ? 0 : (t==0 ? 1 : N(0,j)*h/t);
	AMatrixCalc(1,i,j);
	h1=A2[1]*h+A2[2]*al;
	al1=A2[3]*h+A2[4]*al;
	return al1==0 ? LN : N(i,j)*h1/al1;
}

double cLens1::PupilToPupil(int i1/*=0*/,int i2/*=0*/){
	// ���˓�����ˏo���܂ł̋���
	if(i1==0 && i2==0){ i1=1; i2=k; };
	if( i1<1 || i2<i1 || k<i2 ) return 0;
	return -ti(i1) +TotalThickness(i1,i2) +t1i(i2);
}

double cLens1::ExitPupilZ(int findpupil){
	// �ŏI�ʂ���̓��˓��̈ʒu
	// t1�̋@�\�ɉ����Cfindpupil=true �̂Ƃ��͌����ǐՂ���J���i��ʂ����肷��D
	if(findpupil==0){
		return t1();
	}
	else{
		cLens1 x=*this;  // �f�[�^�ɕύX�������Ȃ��悤�R�s�[���쐬
		point dummy;

		x.stop=NotThruSurf(x.FindRay(dummy,0,0,"ymax",1,1));  // ��������𐧌�����ʔԍ�
		if(x.stop>0){
			x.tCalc();
			return x.t1();
		}
		else return 0;		
	}
}

double cLens1::ExitPupilDia(int findpupil){
	// ���˓��̒��a
	if(findpupil==0){
		return EPD*Mpupil();
	}
	else{
		cLens1 x=*this;  // �f�[�^�ɕύX�������Ȃ��悤�R�s�[���쐬
		point dummy;

		x.stop=NotThruSurf(x.FindRay(dummy,0,0,"ymax",1,1));  // ��������𐧌�����ʔԍ�
		if(x.stop>0){
			x.EPCalc();
			return x.EPD*x.Mpupil();
		}
		else return 0;	
	}
}

double cLens1::ThinT(int i1,int i2){
	// i1����i2�ʂ̕����n�𔖃����Y�ɕϊ������Ƃ���
	// ���˓��ʒu�i���Ȃ킿�C�������Y�ϊ��O�̎啽�ʂ�����˓��܂ł̋����j��Ԃ��D
	// ��C���̋����Ƃ���D
	return ( ti(i1,1)-delta(i1,i2,1) )/N(i1-1,1);
}

double cLens1::ThinT(){
	return ThinT(1,k);
}

double cLens1::ImageInfinity(){
	// ���̋���s1�̋t����Ԃ��D���Ȃ킿���_���œ_�ɂ��邩�ǂ�����]������D
	double s1;

	s1=this->s1(1);
	return s1==0 ? 0 : 1/s1;
}

double cLens1::TelecentricityObj(){
	// ���˓��ʒut�̋t����Ԃ��D
	cLens1 X;
	double t;

	X=*this;
	t=X.tCalc();
	return t==0 ? 0 : 1/t;
}

double cLens1::TotalThickness(int i1,int i2){
	// i1�ʂ���i2�ʂ܂ł̌���(���̖ʂ�0,���ʂ�k+1)
	double result;
	int i;
	if( 0<=i1 && i1<i2 && i2<=k+1 ){
		result=0;
		for(i=i1; i<=i2-1; i++){
			if(i==0)           result+= fabs(s)<LN ? -s : 0;
			if(1<=i && i<=k-1) result+= d(i);
			if(i==k)           result+= dk(0);
		}
		return result;
	}
	else return 0;	
}

double cLens1::TotalThickness(){
	return TotalThickness(1,k);
}

double cLens1::TotalInAirThickness(int i1,int i2){
	// i1�ʂ���i2�ʂ܂ł̋�C���Z����(���̖ʂ�0,���ʂ�k+1)
	double result;
	int i;
	if( 0<=i1 && i1<i2 && i2<=k+1 ){
		result=0;
		for(i=i1; i<=i2-1; i++){
			if(i==0)           result+= fabs(s)<LN ? -s/N(0,1) : 0;
			if(1<=i && i<=k-1) result+= d(i)/N(i,1);
			if(i==k)           result+= dk(0)/N(k,1);
		}
		return result;
	}
	else return 0;	
}

double cLens1::TotalInAirThickness(){
	return TotalInAirThickness(1,k);
}

double cLens1::ConjugateLength(){
	// ���̂��瑜�܂ł̋����Ds1fix=0�Ƃ���D
	cLens1 x;
	x=*this;
	x.s1fix=0;
	return x.TotalThickness(0,x.k+1);
}

double cLens1::TotalOpticalThickness(int i1,int i2){
	// i1�ʂ���i2�ʂ܂ł̌��w����(���̖ʂ�0,���ʂ�k+1)
	double result;
	int i;

	if( 0<=i1 && i1<i2 && i2<=k+1 ){
		result=0;
		for(i=i1; i<=i2-1; i++){
			if(i==0)           result+= fabs(s)<LN ? -N(0,1)*s : 0;
			if(1<=i && i<=k-1) result+= N(i,1)*d(i);
			if(i==k)           result+= N(k,1)*dk(0);
		}
		return result;
	}
	else return 0;	
}

double cLens1::TotalOpticalThickness(){
	return TotalOpticalThickness(1,k);
}

double cLens1::TotalOpticalThicknessDerivative2(int i1,int i2,int j){
	// ��j�g����i1�ʂ���i2�ʂ܂ł̌��w�����̔g��2�K����(���̖ʂ�0,���ʂ�k+1)
	// �P�ʂ� 1/nm (��: d(i)��mm�P��)
	double result;
	int i;
	if( 0<=i1 && i1<i2 && i2<=k+1 ){
		result=0;
		for(i=i1; i<=i2-1; i++){
			if(i==0)           result+= fabs(s)<LN ? -(s*1e6)*Re(IndexDerivative2(gname(0),Wl(j))) : 0;
			if(1<=i && i<=k-1) result+= (d(i)*1e6)*Re(IndexDerivative2(gname(i),Wl(j)));
			if(i==k)           result+= (dk(0)*1e6)*Re(IndexDerivative2(gname(k),Wl(j)));
		}
		return result;
	}
	else return 0;	
}

vector<double> cLens1::VertexGlobal(int i){
	// ��i�ʒ��_�̃O���[�o�����W��Ԃ�
	if(0<=i && i<=k+1){
		make_coordinate(0);
		return o[i];
	}
	else{
		return vector<double>(0,0,0);
	}
}

vector<double> cLens1::exGlobal(int i){
	// ��i�ʃ��[�J�����Wx���̕����x�N�g�����O���[�o�����W�ŕԂ�
	if(0<=i && i<=k+1){
		make_coordinate(0);
		return ex[i];
	}
	else{
		return vector<double>(0,0,0);
	}
}

vector<double> cLens1::eyGlobal(int i){
	// ��i�ʃ��[�J�����Wy���̕����x�N�g�����O���[�o�����W�ŕԂ�
	if(0<=i && i<=k+1){
		make_coordinate(0);
		return ey[i];
	}
	else{
		return vector<double>(0,0,0);
	}
}

vector<double> cLens1::ezGlobal(int i){
	// ��i�ʃ��[�J�����Wz���̕����x�N�g�����O���[�o�����W�ŕԂ�
	if(0<=i && i<=k+1){
		make_coordinate(0);
		return ez[i];
	}
	else{
		return vector<double>(0,0,0);
	}
}
vector<double> cLens1::IncidentPointGlobal(int i,double yObj,double xObj,std::string SetRay,int FindPupil){
	// ��i�ʂƌ����̌�_�̃O���[�o�����W��Ԃ�
	point p;
	vector<double> P;

	if(0<=i && i<=k+1){
		make_coordinate(0);
		FindRay(p,yObj,xObj,SetRay,FindPupil,1);
		RayTrace(yObj,xObj,p.y,p.x,0,1,0,0,0,0);
		P=vector<double>(this->x[i],this->y[i],this->z[i]);  // P�̃��[�J�����W
		P=o[i] +P.x*ex[i] +P.y*ey[i] +P.z*ez[i];             // P�̃O���[�o�����W
		return P;
	}
	else{
		return vector<double>(0,0,0);
	}
}

double cLens1::Clearance(int iPoint,double yObjPoint,double xObjPoint,std::string SetRayPoint,
						 int iLine,double yObjLine,double xObjLine,std::string SetRayLine,int FindPupil){
	// �O���[�o�����W�ōl����D
	// (1) ��������o�H(yObjPoint,xObjPoint,SetRayPoint�Ŏw��)��iPoint�ʂɓ��˂���_P
	// (2) ������̌����o�H(yObjLine,xObjLine,SetRayLine�Ŏw��)�� iLine�ʂւ̓��˓_P1, iLine+1�ʂւ̓��˓_P2
	//     ��ʂ钼��
	// �ɂ��āC(1)(2)�̍ŒZ������Ԃ��D���̕����́C
	//  P,P1,P2���O���[�o�����W��yz���ʂɓ��e�����Ƃ��C �x�N�g�� P1->P �� P1->P2 �̊O�ς�
	//  x�������������Ƃ��ɐ��Ƃ���D
	//  �Ⴆ�΁Coff-axial���w�n�̎����݌v�ŁC�~���[�Ƃ��̃~���[�ł͔��˂��Ȃ������������Ȃ����߂̕]���֐��Ƃ��Ďg���D
	vector<double> P,P1,P2,X;
	double distance;

	P =IncidentPointGlobal(iPoint ,yObjPoint,xObjPoint,SetRayPoint,FindPupil);
	P1=IncidentPointGlobal(iLine  ,yObjLine,xObjLine,SetRayLine,FindPupil);
	P2=IncidentPointGlobal(iLine+1,yObjLine,xObjLine,SetRayLine,FindPupil);
	X=NearestPointOnLine(P,P1,P2);  // P1��P2��ʂ钼����P����~�낵�������̑�
	distance=abs(X-P);              // �����̒���
	P.x=0; P1.x=0; P2.x=0;          // P,P1,P2���O���[�o�����W��yz���ʂɓ��e����
	return sgn(sProduct(vProduct(P-P1,P2-P1),vector<double>(1,0,0)))*distance;
}

double cLens1::GDD(int i1,int i2,int j){
	// i1�ʂ���i2�ʂ܂ł̌��w�����ɂ���j�g���ł�
	// �Q�x�����U (d/d��)(d��/d��) �� fs^2 �P�ʂŕԂ��D
	const double c=299.7925;   // �^������x[nm/fs]
	double wl;

	// �� = Nd/��
	// �� = 2��c/��
	// d��/d�� = (d��/d��)(d��/d��) = -{ (dN/d��)d*��-Nd }/c
	// (d/d��)(d��/d��) = {(d/d��)(d��/d��)}(d��/d��) = ��^3/(2��c^2){(d/d��)(dN/d��)}d
	wl=Wl(j);
	return wl*wl*wl/(2*PI*c*c)*TotalOpticalThicknessDerivative2(i1,i2,j);
}

double cLens1::ChirpedPulseWidth(int i1,int i2,int j,double dt0_fs){
	//  ��j�g���̃K�E�X�^�p���X��i1�ʂ���i2�ʂ܂ł̌��w�����ɂ��
	//  �Q�x�����U���󂯂���̃p���X����Ԃ��D
	//  dt0_fs �͍ŏ��̃p���X����t(fs)
	return cOptics::ChirpedPulseWidth(dt0_fs,GDD(i1,i2,j));
}

std::string cLens1::SurfaceSagTable(int i,double hStep){
	// x,y����̖ʃT�O�ʂ̕\���쐬����D�񋅖ʂ̐��}�ȂǂɎg���D
	std::string s;
	int k;
	double h,z;
	char buf[100];

	// "    y         z"
	// "yyyy.yyyyyzzzz.zzzzzzz"

	s+="    y         z\n";

	if( asph_type(i)==CONIC_OA){ // y�̐����Ŕ�Ώ̂Ȃ̂ŁC�������o��
		h=-phi_y(i)/2;
		z=surface_sag(i,h,0,0);
		sprintf(buf,"%10.5f%12.7f\n",h,z); s+=buf;

		k=ToInt(floor(phi_y(i)/2/hStep));

		for(k=k; k>0; k--){
			h=-hStep*k;
			z=surface_sag(i,h,0,0);
			sprintf(buf,"%10.5f%12.7f\n",h,z); s+=buf;
		}
	}

	h=0;
	while( h<phi_y(i)/2 ){
		h+=hStep;
		if( h>phi_y(i)/2 ) h=phi_y(i)/2;
		z=surface_sag(i,h,0,0);
		sprintf(buf,"%10.5f%12.7f\n",h,z); s+=buf;
	}

	if( asph_type(i)!=SPH && (asph_type(i)!=CONIC || cylinder(i)!=0) ){
		s+="\n";
		s+="    x         z\n";
		h=0;
		while( h<phi_x(i)/2 ){
			h+=hStep;
			if( h>phi_x(i)/2 ) h=phi_x(i)/2;
			z=surface_sag(i,0,h,0);
			sprintf(buf,"%10.5f%12.7f\n",h,z); s+=buf;
		}
	}

	return s;
}

double cLens1::SurfaceSagMax(int i,double hStep){
	// �O�a���ł̍ő�T�O��Ԃ��i���Ȃ炸�������̃T�O�ł͂Ȃ��j�D
	// y�����h_step�X�e�b�v�ŒT������D
	double y,sag,max;
	
	if(hStep==0) hStep=1;
	max=-1e30;
	for(y=hStep; fabs(y)<=fabs(phi_y(i)/2); y+=hStep){
		sag=surface_sag(i,y,0,0);
		if(sag>max) max=sag;
	}
	return max;
}

double cLens1::SurfaceSagMin(int i,double hStep){
	// �O�a���ł̍ŏ��T�O��Ԃ��D
	double y,sag,min;
	
	if(hStep==0) hStep=1;
	min=1e30;
	for(y=hStep; fabs(y)<=fabs(phi_y(i)/2); y+=hStep){
		sag=surface_sag(i,y,0,0);
		if(sag<min) min=sag;
	}
	return min;
}

double cLens1::SurfaceSlope(int i,double y,double x){
	// i�ʂ�(x,y)�ɂ�����@���ƌ������Ȃ��p��deg�ŕԂ��D
	// �����Y�̒[�ł̃R�[�g��������������̂ɍ쐬 (2015.08.10)
	vector<double> n;

	n=surface_normal(i,y,x,0);
	n=n/abs(n);
	return acos( sProduct(n,vector<double>(0,0,1)) )*180/PI;
}

double cLens1::SurfaceSlopeMax(int i,double hStep){
	// �O�a���ł̍ő�ʌ��z��deg�ŕԂ��D
	// y�����h_step�X�e�b�v�ŒT������D
	double y,slope,max;
	
	if(hStep==0) hStep=1;
	max=0;
	for(y=hStep; fabs(y)<=fabs(phi_y(i)/2); y+=hStep){
		slope=SurfaceSlope(i,y,0);
		if(slope>max) max=slope;
	}
	slope=SurfaceSlope(i,sgn(hStep)*fabs(phi_y(i)/2),0);
	if(slope>max) max=slope;
	return max;
}

std::string cLens1::SurfaceSlopeTable(int i,double hStep){
	// x,y����̖ʃX���[�v�̕\���쐬����D�񋅖ʃ����Y�̉��H���m�F�ȂǂɎg���D
	std::string s;
	double h,th;
	char buf[100];

	// "    y        ��(deg)"
	// "yyyy.yyyyyzzzz.zzzzz"

	s+="    y        ��(deg)\n";
	h=0;
	while( h<phi_y(i)/2 ){
		h+=hStep;
		if( h>phi_y(i)/2 ) h=phi_y(i)/2;
		th=SurfaceSlope(i,h,0);
		sprintf(buf,"%10.5f%10.5f\n",h,th);
		s+=buf;
	}

	if( asph_type(i)==SPH || (asph_type(i)==CONIC && cylinder(i)==0) ) return s;

	s+="\n";
	s+="    x        ��(deg)\n";
	h=0;
	while( h<phi_x(i)/2 ){
		h+=hStep;
		if( h>phi_x(i)/2 ) h=phi_x(i)/2;
		th=SurfaceSlope(i,0,h);
		sprintf(buf,"%10.5f%10.5f\n",h,th);
		s+=buf;
	}

	return s;
}

double cLens1::koba(int i){
	double z1,z2, ty,tx;
	if( 1<=i && i<=k-1 ){
		z1=surface_sag(i  ,phi_y(i)/2  ,0,0);
		z2=surface_sag(i+1,phi_y(i+1)/2,0,0);
		ty=d(i)>=0 ? -z1+z2+d(i) : z1-z2-d(i);

		z1=surface_sag(i  ,0,phi_x(i)/2  ,0);
		z2=surface_sag(i+1,0,phi_x(i+1)/2,0);
		tx=N(i,1)>=0 ? -z1+z2+d(i) : z1-z2-d(i);
		
		return Min(ty,tx);
	}
	else return 0;
}

double cLens1::ct(int i){
	// z���̐������ł͂Ȃ������̐i�s�����ɑ��钸�_�Ԋu
	// �Ⴆ�΁C�����݌v�̋��E������ d(i)>0 �łȂ� tc(i)>0 �Ƃ��ׂ��D
	return N(i,1)>0 ? d(i) : -d(i);
}

double cLens1::Steepness(int i){
	// |�O�`���a|/|�ȗ����a| ��Ԃ��D1�ȏ�͕s�\�Ȍ`��ł���i1�̂Ƃ������j�D
	// ���ʁC�񋅖ʂł�-1��Ԃ��D
	// 0.8�ȉ��ł��邱�Ƃ�����\�̖ڈ��ƂȂ�D
	double x;

	if( asph_type(i)==SPH && r(i)!=0 ){
		x=lens_phi(Max(EAx(i),EAy(i)))/2/fabs(r(i));
	}
	else{
		x=-1;
	}

	return x;
}

std::string cLens1::DefaultBC(double min_koba/*=1*/,double min_ct/*=0.5*/,double max_steepness/*=0.8*/,
							  double weight/*=1*/,int i1/*=0*/,int i2/*=0*/){
	// i1�ʂ���i2�ʂ̊Ԃ̃R�o���ƒ��S���� min_koba, min_ct �ȏ�ƂȂ�C
	// �ʋȗ����}�s�����Ȃ��isteepness = �����Y�O�`���a / �ʋȗ����a) �������C
	// optimize�֐���target�����Ƃ��Đݒ肷�镶����𐶐�����D
	//
	// weight�̃f�t�H���g��-1����1���D�܂����D-1�̂Ƃ��͎����݌v�őS�Ă�r,d��ϐ��ɂ��Ȃ��ƁC
	// Lagrange����搔�@�̕������������ɍœK�����~�܂��Ă��܂��D
	int i;
	char buf[100];
	std::string s;

	if(i1==0) i1=1;
	if(i2==0) i2=k;

	for(i=i1; i<=i2-1; i++){
		sprintf(buf,"(koba %d) > %g %g;\n",i,min_koba,weight); s+=buf;
	}
	for(i=i1; i<=i2-1; i++){
		sprintf(buf,"(ct %d) > %g %g;\n",i,min_ct,weight); s+=buf;
	}
	for(i=i1; i<=i2; i++){
		sprintf(buf,"(Steepness %d) < %g %g;\n",i,max_steepness,weight); s+=buf;
	}

	return s;
}

double cLens1::DAbsC(int i1,int i2){
	double c1,c2;
	c1= r(i1)==0 ? 0 : 1/r(i1);
	c2= r(i2)==0 ? 0 : 1/r(i2);
	return fabs(c1)-fabs(c2);
}

double cLens1::ShapeFactor(int i){
	// i�ʂ�i+1�ʂ̌`����qq��Ԃ��D
	double c1,c2;

	if( 1<=i && i<=k-1 ){
		c1=c(i);
		c2=c(i+1);
		return c1-c2==0 ? 0 : (c1+c2)/(c1-c2);
	}
	else{
		return 0;
	}
}

double cLens1::ZValue(int i,int i1){
	// ��i�ʂ�i+1�ʂ�Z�l���v�Z����D
	double ea,ea1;
	if( 1<=i && i<i1 && i1<=k ){
		switch(EAtype(i)){
		case 0:
		case 1:
			if(EAx(i)==0) ea=EAy(i);
			else          ea= EAy(i)<EAx(i) ? EAy(i) : EAx(i);
			break;
		default:
			ea=0;
			break;
		}
		switch(EAtype(i1)){
		case 0:
		case 1:
			if(EAx(i1)==0) ea1=EAy(i1);
			else           ea1= EAy(i1)<EAx(i1) ? EAy(i1) : EAx(i1);
			break;
		default:
			ea1=0;
			break;
		}
		return ZValue(r(i),r(i+1),ea,ea1);
	}
	else return 0;
}

double cLens1::PreformKoba(int i){
	// ���[���h�����Y�ł́C
	// �����̏ꍇ�v���t�H�[�������Y��R�͋ߎ�R���p������D
	// ���ʂł͒��S���瓖�Ă邽�߂�₻���菬����R���p������D
	// �������Y�̏ꍇ�͕s���D
	cLens1 X;

	X=*this;
	X.asph_type(i)=0;
	X.asph_type(i+1)=0;
	return X.koba(i);
}
/*
double cLens1::Inflection(int i){
	// y�����ɂ����ĕϋȂ̒��x������D
	// 0�܂��͕��̒l��Ԃ��C0�ł���Εϋȓ_�͂Ȃ����^���₷���D
	const int DIV=10;
	int j,maxj;
	double z[DIV+1];
	double ymax,y,max,min,sign;
	
	ymax=phi_y(i)/2;
	for(j=1; j<=DIV; j++){
		y=ymax/DIV*j;
		z[j]=surface_sag(i,y,0,0);
	}
	sign=sgn(z[1]);
	max=0;
	for(j=1; j<=DIV; j++){
		if( z[j]*sign>max ){
			max=z[j]*sign;  // max  : �ߎ��T�O�����𐳂Ƃ����Ƃ��̍ō�����
			maxj=j;         // maxj : �ō��n�_�imax�ɂȂ�ʒu�j
		}
	}
	min=0;
	for(j=maxj; j<=DIV; j++){
		if( z[j]*sign-max<min ){
			min=z[j]*sign-max;  // min: �ō��n�_�Ɖ��̊Ԃł̍Œፂ��
		}
	}
	return min;
}
*/
double cLens1::Inflection(int i){
	// 	�T�O��2�������̍ŏ���Ԃ��i�ߎ��ł�2�������𐳂Ƃ���j�D
	//  ���ꂪ0�ȏ�ł���Εϋȓ_���Ȃ��C�������₷���D
	const int DIV=10;
	int j;
	double z[DIV+1];
	double dy, ymax,y,min,sign, z2;
	
	ymax=phi_y(i)/2;
	dy=ymax/DIV;

	for(j=1; j<=DIV; j++){
		y=dy*j;
		z[j]=surface_sag(i,y,0,0);  // �e�ʒu�ł̃T�O��
	}

	sign=sgn(z[1]);   // �ߎ��ł̃T�O�̕���
	min=1e30;
	
	for(j=2; j<=DIV-1; j++){
		if( ( z2=( (z[j+1]-z[j])-(z[j]-z[j-1]) )/dy*sign ) < min ){  // z2 = 2�������W��
			min=z2;
		}
	}

	return min;
}

double cLens1::InflectionPoint(int i){
	// i�ʂ�y�����̕ϋȓ_�����߂�D�Ȃ��Ƃ���0��Ԃ��D
	const int DIV=10;
	int j;
	double ymax,y,y1,y2,sign;

	ymax=phi_y(i)/2;
	y1=ymax/DIV;
	sign=sgn(surface_normal(i,y1,0,0).y);

	y2=0;
	for(j=2; j<=DIV; j++){
		y=ymax/DIV*j;
		if( sign*surface_normal(i,y,0,0).y<0 ){
			y2=y;
			break;
		}
	}
	if(y2==0) return 0;

	do{
		if( sign*surface_normal(i,(y1+y2)/2,0,0).y<0 ){
			y2=(y1+y2)/2;
		}
		else{
			y1=(y1+y2)/2;
		}
	} while( fabs(y2-y1)>0.001 );

	return (y1+y2)/2;
}

double cLens1::AxialRadiusY(int i,double y){
	// Axial�ȗ����a��YZ�ʓ��ŋ��߂�
	//    Axial�ȗ����a = ��i�ʂ̓_(0,y)�Ƃ��̓_�̖@����Z���Ƃ̌�_�̋���
	vector<double> v;
	double z,za;

	z=surface_sag(i,y,0,0);      // i�ʏ�(0,y)�ɂ�����T�O
	v=surface_normal(i,y,0,0);   // (0,y)�̖@���̕����]��
	if(v.y==0){
		return 0;
	}
	else{
		za=z-(v.z/v.y)*y;  // �@����Z���Ƃ̌�_
		return sgn(za)*sqrt(y*y+(za-z)*(za-z));
	}
}

double cLens1::TangentialCurvature(int i,double y){
	// Tangential�ȗ����a��YZ�ʓ��ŋ��߂�
	vector<double> v;
	double cx,cxy,cy;
	matrix<double> T(3,3);

	v=surface_normal(i,y,0,0);   // (0,y)�̖@���̕����]��
	T=Tmatrix(v);
	surface_curvature(i,cx,cy,cxy,0,y,T,0);
	return cy;
}

double cLens1::TangentialCurvatureMax(int i,double hStep){
	// �O�a���ł�y�����ȗ����a�̐�Βl�̍ő��Ԃ�
	// y�����h_step�X�e�b�v�ŒT������D
	double y,cy,max;
	
	if(hStep==0) hStep=1;
	max=0;
	for(y=hStep; fabs(y)<=fabs(phi_y(i)/2); y+=hStep){
		cy=fabs(TangentialCurvature(i,y));
		if(cy>max) max=cy;
	}
	return max;
}

void cLens1::DzToD(int i){
	// ��i�ʂ�dz(i)�̌��ʂ�ʊԊud�ŕ\�����Cdz(i)��0�ɂ���D
	// dz(i)�͋ߎ��v�Z�Ŗ�������邪�C�{�֐��ɂ��dz(i)�̌��ʂ��ߎ��l�ɔ��f�ł���D
	int ii;

	if(decenter_type(i)==0) return;

	switch(decenter_type(i)){
	case 1:
		d(i-1)+=dz(i);
		d(i)+=dz1(i);
		for(ii=i+1; ii<=k; ++ii){
			if(ret(ii)==i && decenter_type(ii)==1){
				d(ii)-=dz(i);
			}
		}
		dz(i)=dz1(i)=0;
		break;
	case 2:
		d(i-1)+=dz(i);
		d(i)-=dz(i);
		dz(i)=0;
		break;
	case 3:
		d(i-1)+=dz(i);
		d(i)+=dz(i);
		for(ii=i+1; ii<=k; ++ii){
			if(ret(ii)==i && decenter_type(ii)==1){
				d(ii)-=dz(i)*2;
			}
		}
		dz(i)=0;
		break;
	}
}

void cLens1::RotateSurface(int i,double Sx,double Sy,double Sz,double Rx,double Ry,double Rz,double th){
	int j;
	vector<double> S,R, E,E1, N,N1,Np,Ns,Ns1, L,L1, K, Z,X;
	double thy,thx,thz;
	S=vector<double>(Sx,Sy,Sz);
	R=vector<double>(Rx,Ry,Rz);
	th=th*PI/180;
	R=R/abs(R);
	E=-S+sProduct(S,R)*R;					// (A.11.3)	
	for(j=1; j<=2; ++j){
		if(j==1) N=vector<double>(0,0,1);
		if(j==2) N=vector<double>(1,0,0);
		Np=sProduct(N,R)*R;                 // (A.11.6)
		Ns=N-Np;							// (A.11.6)
		L=E+Ns;								// (A.11.4)
		L1=L*cos(th)-vProduct(L,R)*sin(th); // (A.11.14)
		E1=E*cos(th)-vProduct(E,R)*sin(th); // (A.11.16)
		Ns1=L1-E1;							// (A.11.9)
		N1=Ns1+Np;							// (A.11.7)
		K=E1-E; 
		if(j==1) Z=N1;
		if(j==2) X=N1;
	}
	Set_dx(i,K.x); Set_dy(i,K.y); Set_dz(i,K.z);
	rotate_angle2(thx,thy,thz,Z.x,Z.y,Z.z,X.x,X.y);
	Set_rox(i,thx); Set_roy(i,thy), Set_roz(i,thz); Set_order(i,0);
}

void cLens1::RotateBlock(int i1,int i2,double Sx,double Sy,double Sz,double Rx,double Ry,double Rz,double th){
	// i1�ʂ���i2�ʂ܂ł��Ci1�ʋ��ܑO�ΐS�O���W���_����S�̓_��ʂ�C������R�̎��̎�����Ɖ�]����D
	// ����ɁCi2�ʋ��܌�ΐS��̌�����i1�ʋ��ܑO�ΐS�O�̌����̗��z�����ƈ�v������D
	int i;
	double x0,y0,z0,l,m,n;

	decenter_type(i1)=1;
	RotateSurface(i1,Sx,Sy,Sz,Rx,Ry,Rz,th);
	make_coordinate(0);
	x0=y0=z0=0; l=m=0; n=1;                   // i1�ʕΐS�O���̑�����
	transform_line(x0,y0,z0,l,m,n,i1,1,i1,0); // i1�ʒ��_���W
	for(i=i1; i<=i2-1; i++){
		image_line(x0,y0,z0,l,m,n,i);             // i�ʌ���
		transform_line(x0,y0,z0,l,m,n,i,0,i+1,0); // i+1�ʒ��_���W
	}
	image_line(x0,y0,z0,l,m,n,i2);  // i2�ʌ���
	decenter_type(i2)=1;
	order1(i2)=0;
	dx1(i2)=-z0*l/n+x0;
	dy1(i2)=-z0*m/n+y0;
	dz1(i2)=0;
	rotate_angle(rox1(i2),roy1(i2),l,m,n);
	roz1(i2)=0;
}

void cLens1::RotateBlockX(int i1,int i2,double th){
	// i1�ʂ���i2�ʂ܂ł��C���̕����̑O���ߓ_��ʂ�x�����̎��̎�����Ɖ�]����D
	// ����ɁCi2�ʋ��܌�ΐS��̌�����i1�ʋ��ܑO�ΐS�O�̌����̗��z�����ƈ�v������D
	RotateBlock(i1,i2, 0,0,nodal(i1,i2,1), 1,0,0, th);
}

void cLens1::RotateBlockY(int i1,int i2,double th){
	// i1�ʂ���i2�ʂ܂ł��C���̕����̑O���ߓ_��ʂ�y�����̎��̎�����Ɖ�]����D
	// ����ɁCi2�ʋ��܌�ΐS��̌�����i1�ʋ��ܑO�ΐS�O�̌����̗��z�����ƈ�v������D
	RotateBlock(i1,i2, 0,0,nodal(i1,i2,1), 0,1,0, th);
}

void cLens1::RotateBlockAroundPupil(int i1,int i2,int stop,double rox){
	// i1�ʂ���i2�ʂ܂ł��C���̕����̓��˓��𒆐S��rox��]����D
	// ���̂Ƃ��C
	//  �E���˓��͌�����Ɏc��
	//  �Ei1�ʒ��_��z�����ʒu���ێ�����
	// (�����ɂ͓��˓���z�����Ɉړ����邪�C�������������v�Z�\�t�g�ɓ��͂��Ղ��j
	double t;

	t=Trimed(i1,i2).tCalc();
	decenter_type(i1)=1;
	dy(i1)=t*sin(rox/180*PI);
	this->rox(i1)=rox;
	decenter_type(i2)=1;
	ret(i2)=i1;
}

void cLens1::RotateSurfaceXYZOrder
     (int i,double rox1,double roy1,double roz1,double rox2,double roy2,double roz2,int IsMirror){
	// rox1->roy1->roz1->rox2->roy2->roz2�̕ΐS�Ɠ�����rox,roy,roz��^����D
	// rox1,roy1,roz1�݂̂ł͕\�����ɂ����ꍇ���o���̂�rox2,roy2,roz2�������ɂ����D(07.09.21)
	//
	//	rox,roy�ŉ񂵂�z���ƁCrox/2,roy/2�ŉ񂵂�z���́C���ˌ����C���˖ʖ@��
	//	�̊֌W�Ƃ͌���Ȃ��D���������Ă��̊֐����Ȃ��Ɣ��˖ʖ@���̕ΐS�ʂ��킩��Ȃ��D
	//	�܂��CIsMirror=1�ŋ��߂�rox(i),roy(i)�ł�type3�ΐS�́C
	//  �ΐS��̌������ΐS�O�����Ɠ��˔��˂̊֌W�ɂȂ�Ƃ͌���Ȃ��D
	//  (�����ΐS��2��J��Ԃ��Ƃ��C�ez���͓��ˌ��C���˖ʖ@���C���ˌ��̊֌W�ɂȂ�Ƃ͌���Ȃ�)
	//  ���˖ʂ�type2�Ƃ��C�ˋ�ʂ��g���Ĕ��ˌ��������`����Ȃǂ��K�v�D
	//  �����ł́Crox1(i),roy1(i),roz1(i)�Ŕ��ˌ�����ݒ肵�Ă���D
	matrix<double> A(3,3);
	vector<double> Z(0,0,1), Z1;
	vector<double> X(1,0,0), X1;
	A=::Tmatrix(    0,    0,-roz2);
	A=::Tmatrix(    0,-roy2,    0)*A;
	A=::Tmatrix(-rox2,    0,    0)*A;
	A=::Tmatrix(    0,    0,-roz1)*A;
	A=::Tmatrix(    0,-roy1,    0)*A;
	A=::Tmatrix(-rox1,    0,    0)*A;
	Z1=A*Z;
	X1=A*X;

	if(IsMirror){
		this->decenter_type(i)=1;
		// ���ˌ�����ݒ�
		rotate_angle2(this->rox1(i),this->roy1(i),this->roz1(i),Z1.x,Z1.y,Z1.z,X1.x,X1.y);
		this->order1(i)=0;
		this->ret(i)=i;
		// �~���[�ʂ�ݒ�
		Z1=Z1+Z;
		X1=X1+X;
		rotate_angle2(this->rox(i),this->roy(i),this->roz(i),Z1.x,Z1.y,Z1.z,X1.x,X1.y);
		this->order(i)=0;
	}
	else{
		rotate_angle2(this->rox(i),this->roy(i),this->roz(i),Z1.x,Z1.y,Z1.z,X1.x,X1.y);
		this->order(i)=0;
	}
}

void cLens1::RotateSurfaceZYXOrder(double& roz,double& roy,double& rox,int IsMirror){
	matrix<double> A(3,3);
	vector<double> Z(0,0,1), Z1;
	vector<double> X(1,0,0), X1;
	A=::Tmatrix(-rox,-roy,-roz);
	Z1=A*Z;
	X1=A*X;

	if(IsMirror){
		Z1=Z1+Z;
		// IsMirror=1�ŋ��߂�rox,roy,roz���g����type3�ΐS���s���Ă�
		// �ΐS�OZ���ƕΐS��Z�������ˌ��C���ˌ��̊֌W�ɂȂ�Ƃ͌���Ȃ��D
	}
	rotate_angle2(rox,roy,roz,Z1.x,Z1.y,Z1.z,X1.x,X1.y);
}
void cLens1::RotateSurfaceZYXOrder(int i,double roz,double roy,double rox,int IsMirror){
	if(IsMirror){
		this->decenter_type(i)=1;
		// �~���[�ʂ�ݒ�
		RotateSurfaceZYXOrder(roz,roy,rox,1);
		this->rox(i)=rox;
		this->roy(i)=roy;
		this->roz(i)=roz;
		this->order(i)=0;
		this->ret(i)=i;
		// ���ˌ�����ݒ�
		RotateSurfaceZYXOrder(roz,roy,rox,0);
		this->rox1(i)=rox;
		this->roy1(i)=roy;
		this->roz1(i)=roz;
		this->order1(i)=0;		
	}
	else{
		RotateSurfaceZYXOrder(roz,roy,rox,0);
		this->rox(i)=rox;
		this->roy(i)=roy;
		this->roz(i)=roz;
		this->order(i)=0;
	}
}

double cLens1::Scan(std::string ScanXY,std::string RotateAxisXY,double th_deg,int GetData){
	int i, i1,i2;
	std::string s1,s2;

	if      (ScanXY=="x" || ScanXY=="X") { s1="x-scan"; s2="xscan"; }
	else if (ScanXY=="y" || ScanXY=="Y") { s1="y-scan"; s2="yscan"; }
	else    return 0;

	i1=i2=-1;
	for(i=1; i<=k; ++i){
		if( left(rem(i),(int)s1.length())==s1 || left(rem(i),(int)s2.length())==s2 ){i1=i; break;}
	}
	for(i=i1+1; i<=k; ++i){
		if( left(rem(i),(int)s1.length())==s1 || left(rem(i),(int)s2.length())==s2 ){i2=i; break;}
	}

	if(GetData){
		if(RotateAxisXY=="x" || RotateAxisXY=="X"){
			if(i1!=-1) th_deg=rox(i1); else th_deg=0;
		}
		if(RotateAxisXY=="y" || RotateAxisXY=="Y"){
			if(i1!=-1) th_deg=roy(i1); else th_deg=0;
		}
	}
	else{
		if(RotateAxisXY=="x" || RotateAxisXY=="X"){
			if(i1!=-1) rox(i1)=th_deg;
			if(i2!=-1) rox(i2)=th_deg;
		}
		if(RotateAxisXY=="y" || RotateAxisXY=="Y"){
			if(i1!=-1) roy(i1)=th_deg;
			if(i2!=-1) roy(i2)=th_deg;
		}
	}

	return th_deg;
}

void cLens1::xScan(std::string RotateAxisXY,double th_deg){
	// rem(i)�� "x-scan" �܂��� "xscan" �Ŏn�܂�ʂ�th_deg��]������D
	// RotateAxisXY="x" �̂Ƃ���]����x���C"y"�̂Ƃ���y���D
	Scan("x",RotateAxisXY,th_deg,0);
}
void cLens1::xScan(double th_deg){
	Scan("x","y",th_deg,0);   // x�����ł͑����̏ꍇ�C�~���[�̉�]����y��
}

double cLens1::xScan(std::string RotateAxisXY){
	// rem(i)�� "x-scan" �܂��� "xscan" �Ŏn�܂�ʂ���]�p(deg)�𓾂�D
	// RotateAxisXY="x" �̂Ƃ���]����x���C"y"�̂Ƃ���y���D
	return Scan("x",RotateAxisXY,0,1);
}
double cLens1::xScan(){
	return Scan("x","y",0,1);
}

void cLens1::yScan(std::string RotateAxisXY,double th_deg){
	// rem(i)�� "y-scan" �܂��� "yscan" �Ŏn�܂�ʂ�th_deg��]������D
	// RotateAxisXY="x" �̂Ƃ���]����x���C"y"�̂Ƃ���y���D
	Scan("y",RotateAxisXY,th_deg,0);
}
void cLens1::yScan(double th_deg){
	Scan("y","x",th_deg,0);    // y�����ł͑����̏ꍇ�C�~���[�̉�]����x��
}

double cLens1::yScan(std::string RotateAxisXY){
	// rem(i)�� "y-scan" �܂��� "yscan" �Ŏn�܂�ʂ���]�p(deg)�𓾂�D
	// RotateAxisXY="x" �̂Ƃ���]����x���C"y"�̂Ƃ���y���D
	return Scan("y",RotateAxisXY,0,1);
}
double cLens1::yScan(){
	return Scan("y","x",0,1);
}

void cLens1::zScan(double dz){
	// rem(i)�� "zscan" �܂��� "z-scan" �Ŏn�܂�ʂ� this->dz(i)=dz �Ƃ���D
	int i;
	std::string s1,s2;
	
	s1="zscan";
	s2="z-scan";
	for(i=1; i<=k; ++i){
		if( left(rem(i),(int)s1.length())==s1 || left(rem(i),(int)s2.length())==s2 ){
			this->dz(i)=dz;
		}
	}
}
double cLens1::zScan(){
	// rem(i)�� "zscan" �܂��� "z-scan" �Ŏn�܂�ʂ�dz(i)��Ԃ��D
	int i;
	std::string s1,s2;
	double dz=0;
	
	s1="zscan";
	s2="z-scan";
	for(i=1; i<=k; ++i){
		if( left(rem(i),(int)s1.length())==s1 || left(rem(i),(int)s2.length())==s2 ){
			dz=this->dz(i);
			break;
		}
	}

	return dz;
}

void cLens1::InverseDecenter(double& dx,double& dy,double& dz,double& rox,double& roy,double& roz,int& order){
	// �ΐS(dx,dy,dz),rox��],roy��],roz��]�����ɖ߂��ΐS�����߂�D
	dx=-dx;
	dy=-dy;
	dz=-dz;
	rox=-rox;
	roy=-roy;
	roz=-roz;
	order= order==0 ? 1 : 0;
}

void cLens1::TransformRetDecenter(int i){
	// type1�ΐS�ɂ����āCret(i),dx1,dy1,dz1,rox1,roy1,roz1 ���C
	// ret(i)���g��Ȃ��ŁCdx1,dy1,dz1,rox1,roy1,roz1�ŕ\������D
	vector<double> x,y,z,o;
	double dx,dy,dz, rox,roy,roz;
	
	if( i<1 || k<i ) return;
	if(decenter_type(i)==1  && 0<ret(i) && ret(i)<=i ){
		// ��]�p�̌v�Z
		// �ʍ��W�n���猩��ret�ΐS�Cdx1,dy1,dz1,rox1,roy1,roz1�ΐS���x,y���W�����v�Z����
		x=vector<double>(1,0,0);
		z=vector<double>(0,0,1);
		decenter_out_rev(i,x.x,x.y,x.z,0);
		decenter_out_rev(i,z.x,z.y,z.z,0);
		rotate_angle2(rox,roy,roz,z.x,z.y,z.z,x.x,x.y);
		
		// ���s�ړ��ʂ̌v�Z
		// �ʍ��W�n���猩��ret�ΐS�Cdx1,dy1,dz1,rox1,roy1,roz1�ΐS��̍��W�����_���v�Z����
		o=vector<double>(0,0,0);
		decenter_out_rev(i,o.x,o.y,o.z,1);
		dx=o.x; dy=o.y; dz=o.z;
		
		// �ȏ���Cret(i)�ΐS���g��Ȃ��f�[�^�ɏ���������
		this->dx1(i)=dx; this->dy1(i)=dy; this->dz1(i)=dz;
		this->rox1(i)=rox; this->roy1(i)=roy; this->roz1(i)=roz;
		this->order1(i)=0;
		this->ret(i)=0;
	}
}
void cLens1::TransformRetDecenter(){
	// ���w�n�S�̂ɂ���TransformRetDecenter���s���D
	int i;
	for(i=1; i<=k; i++) TransformRetDecenter(i);
}

void cLens1::DeleteDecenter(){
	// �ΐS���폜����
	int i;
	for(i=0; i<=k+1; i++) decenter_type(i)=0;
}

bool cLens1::IsRotationallySymmetric(){
	int i; bool x=true;
	for(i=1; i<=k; ++i) {
		x=( x && ( (decenter_type(i)==0) || (dy(i)==0 && dx(i)==0 && roy(i)==0 && rox(i)==0) ) );
	}
	return x;
}
bool cLens1::IsYAxisSymmetric(){
	int i; bool x=true;
	for(i=1; i<=k; ++i){
		x=( x && ( (decenter_type(i)==0) || (dx(i)==0 && roy(i)==0) ) ); 
	}
	return x;
}
bool cLens1::IsXAxisSymmetric(){
	int i; bool x=true;
	for(i=1; i<=k; ++i){
		x=( x && ( (decenter_type(i)==0) || (dy(i)==0 && rox(i)==0) ) ); 
	}
	return x;
}

void cLens1::Scale(double m,int i1,int i2,int WithD,int WithEA){
	int i;
	int count=0;
	if( i1<0 || k+1<i2 || i2<i1 ) return;
	if( m==0 ) return;
	if( i1==0 && i2==0) { i1=0; i2=k+1; };  // i1=i2=0�̂Ƃ��͕��̖ʂ��瑜�ʂ܂�
	
	for(i=i1; i<=i2; i++) surf[i].scale(m,WithEA,1);
	if(WithD){
		for(i=i1; i<=i2-1; i++) med[i].scale(m);
	}
	if( i1==0 && i2==k+1 ){
		// ���̖ʂ��瑜�ʂ܂ł̂Ƃ��͕��̋����C���ʒu�����ς���
		scale_condition(m);
	}
}

void cLens1::Scale(double m){
	Scale(m,0,0,1,1);
}

void cLens1::AdjustFocalLength(double fl,int i1,int i2,
							   int KeepD/*=1*/,int ByAllSystem/*=0*/,int WithEA/*=0*/){
	// �f�t�H���g�����̂Ƃ���r�݂̂𒲐�����D
	double cc,cc1,f0;
	int count=0;
	cLens1 buf;

	if( i1==0 && i2==0) { i1=1; i2=k; };
	if( i1<1 || k<i2 || i2<i1 ) return;
	if( fl==0 || power(i1,i2)==0 ) return;
	
	if( KeepD ){
		buf=*this;
		do{
			f0=f(i1,i2);
			
			cc=fl/f0;
			if(ByAllSystem) Scale(cc,0,0,0,WithEA); else Scale(cc,i1,i2,0,WithEA);
			
			cc1=fabs(1/cc)*( (fabs(cc)-1)*(fl-f0)/(f(i1,i2)-f0) +1 );                 // (1)
			if(ByAllSystem) Scale(cc1,0,0,0,WithEA); else Scale(cc1,i1,i2,0,WithEA);  // (2)

			count++;
			if(count>100) { *this=buf; return; } 
			
			// �s (1)(2) �͈ꌩ�Ȃ��Ă��悳���������C�Ȃ��ƐU�����Ė������[�v�ƂȂ邱�Ƃ�����D
			// ����ł��������Ȃ����Ƃ͂���̂Łi�����̌��������Y��fl�ƍŏ���f�̕������قȂ�Ƃ��j�C
			// �J�E���^�[�ɂ��ł��؂�͕K�v�D
		} while( fabs( f(i1,i2)/fl-1 )>0.00000001 );
	}
	else{
		cc=fl/f(i1,i2);
		if(ByAllSystem){
			Scale(cc,0,0,1,WithEA);
		}
		else{
			Scale(cc,i1,i2,1,WithEA);
		}
	}
}

void cLens1::AdjustFocalLength(double fl,int KeepD/*=1*/){
	AdjustFocalLength(fl,1,k,KeepD,1,0);
}

int cLens1::cBend(int i1,int i2,double dc){
	int i;
	if(1<=i1 && i1<=i2 && i2<=k) {
		for(i=i1; i<=i2; ++i){
			if(i==i2 && N(i2,1)-N(i2-1,1)!=0) dc=dc*(N(i1-1,1)-N(i2-1,1))/(N(i2,1)-N(i2-1,1));
			if(r(i)==0) r(i)=1/dc;
			else        r(i)=1/( 1/r(i)+dc );
		}
		return 1;
	}
	else return 0;
}

double cLens1::qValue(int i){
	// ��i,i+1�ʂɂ��āC�`����q q={r(i+1)+r(i)}/{r(i+1)-r(i)} ���v�Z����D
	double c1,c2;
	if(1<=i && i<=k-1) {
		c1=c(i);
		c2=c(i+1);
		return c1==c2 ? 0 : (c1+c2)/(c1-c2);
	}
	else return 0;
}

int cLens1::qBend(int i,double q){
	// ��i,i+1�ʂɂ��āC
	// p=c(i)-c(i+1) ��ۂ����܂� {r(i+1)+r(i)}/{r(i+1)-r(i)}=q �ƂȂ�悤�x���f�B���O����D
	double c1,c2,p;
	if(1<=i && i<=k-1) {
		c1=c(i);
		c2=c(i+1);
		p=c1-c2;
		c1=p*(q+1)/2;
		c2=p*(q-1)/2;
		r(i)  = c1==0 ? 0 : 1/c1;
		r(i+1)= c2==0 ? 0 : 1/c2;
		return 1;
	}
	else return 0;
}

void cLens1::Add2ndOrderTerm(int i,double a2){
	//  �񋅖ʂ̕\���C
	//     z=ch^2/[1+sqrt(1-(Q+1)(ch)^2)] +a4*h^4 +a6*h^6 +a8*h^8 +a10*h^10
	// ��a2*h^2���������Ƃ��C10���܂ł̓W�J,
	//     z=A2*h^2 +A4*h^4 .... +A10*h^10
	// ���������Ȃ�悤��, c,a4,a6,a8,a10��ύX����D
	// �͌^��Ŋp���ʂ�2���̍����g���Ă�����̂����邪�C�C
	// �񋅖ʂ̕\����2���̍���ǉ������r(i)���ߎ�r�Ɠ������Ȃ��Ȃ�v���O�����ύX�ӏ����������߁C
	// �{�֐���ǉ������D(06/06/23)
	double c,c1;
	c =this->c(i);
	c1=c+2*a2;

	r(i)= c1==0 ? 0 : 1/c1;
	a4(i) -=(kp(i)+1)*(c1*c1*c1-c*c*c)/8 *pow(NormH(i),4);
	a6(i) -=(kp(i)+1)*(kp(i)+1)*(c1*c1*c1*c1*c1-c*c*c*c*c)/16 *pow(NormH(i),6);
	a8(i) -=(kp(i)+1)*(kp(i)+1)*(kp(i)+1)*(c1*c1*c1*c1*c1*c1*c1-c*c*c*c*c*c*c)/128 *pow(NormH(i),8);
	a10(i)-=(kp(i)+1)*(kp(i)+1)*(kp(i)+1)*(kp(i)+1)*(c1*c1*c1*c1*c1*c1*c1*c1*c1-c*c*c*c*c*c*c*c*c)/256
		    *pow(NormH(i),10);
}

void cLens1::ToThinLens(int i1,int i2){
	if( i1<1 || i1>=i2 || i2>k ) return;
	int ii;
	cLens1 x=Trimed(i1,i2);
	double fl=x.f();
	double del=x.delta(), del1=x.delta1();
	for(ii=1; ii<=x.k-1; ii++){ x.d(ii)=0; }
	x.AdjustFocalLength(fl,0);
	for(ii=i1; ii<=i2; ii++){ surf[ii]=x.surf[ii-i1+1]; }
	for(ii=i1; ii<=i2-1; ii++){ d(ii)=0; }
	if(i1==1){
		if(fabs(s)<LN) s=s-del;
		if(fabs(t)<LN) t=t-del;
	}
	else{
		d(i1-1)+=del;
	}
	if(i2<k) d(i2)-=del1;
	if(i2==k && s1fix!=0) s1fix-=del1;
}

void cLens1::ToIdealLens(int i1,int i2){
	double f;
	f=this->f(i1,i2,1);
	ToThinLens(i1,i2);
	r(i1)=0;
	asph_type(i1)=IDEAL; fideal(i1)=f;
	Delete(i1+1,i2);
}

double cLens1::ToAplanaticSurf(int i){
	// i�ʂ��A�v���i�e�B�b�N�ȖʂƂ���(SA=CM=AS=0)�D�ȍ~�̋ߎ��֌W�͕ω�����D
	double r;
	
	if(1<=i && i<=k){
		r=si(i)*N(i-1,1)/(N(i,1)+N(i-1,1));
		Set_r(i,r);
		return r;
	}
	else{
		return 0;
	}
}

double cLens1::ToConcentricSurf(int i){
	// i�ʂ𕨓_�ɑ΂��ē��S�ȖʂƂ���(SA=CM=0)�D�ȍ~�̋ߎ��֌W�͕ω�����D
	double r;

	r=si(i);
	Set_r(i,r);
	return r;
}

double cLens1::Get_e(int pre_i1,int pre_i2,int post_i1,int post_i2){
	// pre_i1�ʂ���pre_i2�ʂ܂ł̑O�Q��post_i1�ʂ���post_i2�ʂ܂ł̌�Q��
	// ��_�ԋ��������ܗ��Ŋ��������Z�Ԋu�𓾂�D
	int i;
	double del1,del,e;

	del1= pre_i1==pre_i2   ? 0 : delta1(pre_i1,pre_i2)/N(pre_i2,1);
	del = post_i1==post_i2 ? 0 : delta(post_i1,post_i2)/N(post_i1-1,1);

	e=0;
	for(i=pre_i2; i<=post_i1-1; ++i){
		if(i==0){
			e+= -s/N(0,1);
		}
		else if(i==k){
			e+= dk(0)/N(k,1);
		}
		else{
			e+= d(i)/N(i,1);
		}
	}
	
	return -del1+e+del;
}

void cLens1::Set_e(int pre_i1,int pre_i2,int post_i1,int post_i2,double value){
	// pre_i1�ʂ���pre_i2�ʂ܂ł̑O�Q��post_i1�ʂ���post_i2�ʂ܂ł̌�Q��
	// ��_�ԋ��������ܗ��Ŋ��������Z�Ԋu��value�ɐݒ肷��D
	// d(pre_i2)�𒲐����邱�Ƃɂ��D
	double a;

	a= value-Get_e(pre_i1,pre_i2,post_i1,post_i2);

	if(pre_i2==0){
		s-=a*N(0,1);
	}
	else{
		d(pre_i2)+=a*N(pre_i2,1);
	}
}

double cLens1::e(int i1,int i2){
	// i1-1�ʂ�i1�ʂ���i2�ʂ܂ł̕����̎�_�܂ł̋��������ܗ��Ŋ������l�𓾂�D
	return Get_e(i1-1,i1-1,i1,i2);
}
double cLens1::e1(int i1,int i2){
	// i1�ʂ���i2�ʂ܂ł̕����̎�_����i2+1�ʂ܂ł̋��������ܗ��Ŋ������l�𓾂�D
	return Get_e(i1,i2,i2+1,i2+1);
}

void cLens1::TransformACoefficients(int i,double newNormH){
	a4(i)*=pow(newNormH/NormH(i),4);
	a6(i)*=pow(newNormH/NormH(i),6);
	a8(i)*=pow(newNormH/NormH(i),8);
	a10(i)*=pow(newNormH/NormH(i),10);
	a12(i)*=pow(newNormH/NormH(i),12);
	a14(i)*=pow(newNormH/NormH(i),14);
	a16(i)*=pow(newNormH/NormH(i),16);
	a18(i)*=pow(newNormH/NormH(i),18);
	a20(i)*=pow(newNormH/NormH(i),20);
	NormH(i)=newNormH;
}

void cLens1::ReduceAsphTerms(int i,int max_order,double phi/*=0*/,double h_step/*=1*/){
	// ��i�ʂ̔񋅖ʍ��̍ő原����max_order�Ƃ���D
	// �T�O�ʂ̃T���v�����O�𒼌a phi�̓����ŃX�e�b�v h_step�ōs���C�ŏ����@�ŌW�����Đݒ肷��D
	int terms,j;
	double h,z,z0;
	cFitting X;

	terms=max_order/2-1;
	X.SetNumberOfTerms(terms);
	if(phi==0) phi=ea_max(i);

	for(j=1; j<=terms; ++j) X.dimensionSet(j,j*2+2,0);

	for(h=h_step; h<=phi/2; h+=h_step){
		z0=c(i)*h*h/( 1+sqrt(1-(kp(i)+1)*c(i)*c(i)*h*h) );
		z=surface_sag(i,h,0,0)-z0;
		X.AddData(h,0,z);
	}

	X.CalcCoefficients();
	a4(i)=X.coefficientGet(4,0,0);
	a6(i)=X.coefficientGet(6,0,0);
	a8(i)=X.coefficientGet(8,0,0);
	a10(i)=X.coefficientGet(10,0,0);
	a12(i)=X.coefficientGet(12,0,0);
	a14(i)=X.coefficientGet(14,0,0);
	a16(i)=X.coefficientGet(16,0,0);
	a18(i)=X.coefficientGet(18,0,0);
	a20(i)=X.coefficientGet(20,0,0);
}

double cLens1::SetSpheroid(int i,char majoraxis,int convex,double rshort,double ftof){
	// ��i�ʂ���]�ȉ~�̖ʂƂ���D
	//   majoraxis : ����(��]�Ώ̎�)�̕��� x,y,z
	//   convex    : �ʖ�(z�̕��̕����Ɍ������ē�)�͔�0,���ʂ�0
	//   rshort    : �Z����(���a)
	//   ftof      : �œ_�Ԃ̋���
	//   �Z�������ɂ��L���a��ݒ肷��D
	double a,bb,b;

	a=fabs(rshort);       // �Z����(���a)
	bb=a*a+ftof*ftof/4;   // ������(���a)�̓��
	b=sqrt(bb);

	if(majoraxis=='x' || majoraxis=='X'){
		asph_type(i)=ANAMO;
		r(i)=0;
		ry(i)= convex==0 ? -a : a;
		rx(i)= convex==0 ? -bb/a : bb/a;
		kpy(i)=0;
		kpx(i)=bb/a/a-1;
		EAtype(i)=0;
		EAy(i)=lens_ea(a*2);
		EAx(i)=lens_ea(b*2);
	}
	if(majoraxis=='y' || majoraxis=='Y'){
		asph_type(i)=ANAMO;
		r(i)=0;
		rx(i)= convex==0 ? -a : a;
		ry(i)= convex==0 ? -bb/a : bb/a;
		kpx(i)=0;
		kpy(i)=bb/a/a-1;
		EAtype(i)=0;
		EAy(i)=lens_ea(b*2);
		EAx(i)=lens_ea(a*2);
	}
	if(majoraxis=='z' || majoraxis=='Z'){
		asph_type(i)=CONIC;
		r(i)= convex==0 ? -a*a/sqrt(bb) : a*a/sqrt(bb);
		kp(i)=a*a/bb-1;
		a4(i)=a6(i)=a8(i)=a10(i)=a12(i)=a14(i)=a16(i)=0;
		EAtype(i)=0;
		EAy(i)=lens_ea(a*2);
		EAx(i)=lens_ea(a*2);
	}
	
	return sqrt(bb)-ftof/2; // �œ_���璷�����_�܂ł̋���
}

void cLens1::SetSpheroid2(int i,int i1,int i2){
	// ��i�ʂ̌`����i1�ʂƑ�i2�ʂ̊e���_���œ_�Ƃ����]�ȉ~�̂ɂ���D
	// ��]�ȉ~�͖̂{�֐����s�O��i�ʃ��[�J�����W���_��ʂ���̂Ƃ���D
	// ��]�ȉ~�̂͒Z����[�̋O�Տ�̓_�𒸓_�Ƃ��ĕ\������D
	// (��) �I�v�g�X���L�pSLO  ��1�� �팟�ᓵ�C��2�� �ȉ~���C��3�� �X�L���i�~���[

	double a,b,ff;
	vector<double> A,B,O, v;

	make_coordinate(0);

	b=(abs(o[i2]-o[i])+abs(o[i]-o[i1]))/2;   // �������a = �e�œ_����ȉ~�܂ł̋����̘a
	ff=abs(o[i2]-o[i1]);                     // �œ_�Ԃ̋���
	a=sqrt(b*b-ff*ff/4);                     // �Z�����a

	O=(o[i1]+o[i2])/2.0;     // �ȉ~�̒��S
	B=o[i2]-O;
	B=B/abs(B);              // �ȉ~�̒��S����o[i2]�ɋ߂��������_�Ɍ������P�ʃx�N�g��

	A=o[i]-O;
	A=A-sProduct(A,B)*B;
	A=A/abs(A);              // �ȉ~�̒��S���璸�_�i�Z����[�̋O�Ղ�o[i],o[i1],o[i2]
	                         // �Ō��܂镽�ʂ�2�̌�_�̂����Co[i]�ɋ߂����́j�Ɍ������P�ʃx�N�g��
	                         // ���Ȃ킿i�ʒ��_�͒Z���̐�[�Ƃ���D
	
	o[i]=O+a*A;                             // i��(�ȉ~��)�̒��_��ݒ�
	ez[i]= sProduct(ez0[i],A)>0 ? A : -A;   // i�ʃ��[�J��z������
	ey[i]= sProduct(ey0[i],B)>0 ? B : -B;   // i�ʃ��[�J��y������
	ex[i]= vProduct(ey[i],ez[i]);           // i�ʃ��[�J��x������

	// o[i],ez[i],ey[i],ex[i]����i�ʂ̕ΐS�ʂ����߂�
	decenter_type(i)=1;

	v=o[i]-o0[i];
	dx(i)=sProduct(v,ex0[i]); dy(i)=sProduct(v,ey0[i]); dz(i)=sProduct(v,ez0[i]);
	rotate_angle3(rox(i),roy(i),roz(i),
	              sProduct(ez[i],ex0[i]),sProduct(ez[i],ey0[i]),sProduct(ez[i],ez0[i]),
				  sProduct(ey[i],ex0[i]),sProduct(ey[i],ey0[i]));

	v=o0[i+1]-d(i)*ez0[i+1]-o[i];
	dx1(i)=sProduct(v,ex[i]); dy1(i)=sProduct(v,ey[i]); dz1(i)=sProduct(v,ez[i]);
	rotate_angle3(rox1(i),roy1(i),roz1(i),
	              sProduct(ez0[i+1],ex[i]),sProduct(ez0[i+1],ey[i]),sProduct(ez0[i+1],ez[i]),
				  sProduct(ey0[i+1],ex[i]),sProduct(ey0[i+1],ey[i]));

	// �ʌ`����߂�
	SetSpheroid(i,'y',sProduct(A,ez[i])<0 ? 1:0,a,ff);
}

void cLens1::SetSpheroid3(int i,double z1,double z2){
	// ��i�ʂ̌`���(0,0,z1),(0,0,z2)���œ_�Ƃ����]�ȉ~�̂ɂ���D
	//    �����ŁC(0,0,z1)��i�ʂ̕ΐS�O���[�J�����W�C
	//    (0,0,z2)�͔��ˋ��܌�̕ΐS������������̃��[�J�����W(*1)�ɂ����̂Ƃ���D
	//        (*1)�F���Ȃ킿���_�́Ci+1�ʂ̕ΐS�O���[�J�����W�� (0,0,-d(i)) �ƂȂ�D
	// ��]�ȉ~�͖̂{�֐����s�O��i�ʃ��[�J�����W���_��ʂ���̂Ƃ���D
	// ��]�ȉ~�̂͒Z����[�̋O�Տ�̓_�𒸓_�Ƃ��ĕ\������D

	double a,b,ff;
	vector<double> A,B,O, v,v1,v2;

	make_coordinate(0);

	v1=o0[i]+ez0[i]*z1;
	v2=o0[i+1]+ez0[i+1]*(z2-d(i));

	b=(abs(v2-o[i])+abs(o[i]-v1))/2;   // �������a = �e�œ_����ȉ~�܂ł̋����̘a
	ff=abs(v2-v1);                     // �œ_�Ԃ̋���
	a=sqrt(b*b-ff*ff/4);               // �Z�����a

	O=(v1+v2)/2.0;     // �ȉ~�̒��S
	B=v2-O;
	B=B/abs(B);        // �ȉ~�̒��S����o[i2]�ɋ߂��������_�Ɍ������P�ʃx�N�g��

	A=o[i]-O;
	A=A-sProduct(A,B)*B;
	A=A/abs(A);              // �ȉ~�̒��S���璸�_�i�Z����[�̋O�Ղ�o[i],o[i1],o[i2]
	                         // �Ō��܂镽�ʂ�2�̌�_�̂����Co[i]�ɋ߂����́j�Ɍ������P�ʃx�N�g��
	                         // ���Ȃ킿i�ʒ��_�͒Z���̐�[�Ƃ���D
	
	o[i]=O+a*A;                             // i��(�ȉ~��)�̒��_��ݒ�
	ez[i]= sProduct(ez0[i],A)>0 ? A : -A;   // i�ʃ��[�J��z������
	ey[i]= sProduct(ey0[i],B)>0 ? B : -B;   // i�ʃ��[�J��y������
	ex[i]= vProduct(ey[i],ez[i]);           // i�ʃ��[�J��x������

	// o[i],ez[i],ey[i],ex[i]����i�ʂ̕ΐS�ʂ����߂�
	decenter_type(i)=1;

	v=o[i]-o0[i];
	dx(i)=sProduct(v,ex0[i]); dy(i)=sProduct(v,ey0[i]); dz(i)=sProduct(v,ez0[i]);
	rotate_angle3(rox(i),roy(i),roz(i),
	              sProduct(ez[i],ex0[i]),sProduct(ez[i],ey0[i]),sProduct(ez[i],ez0[i]),
				  sProduct(ey[i],ex0[i]),sProduct(ey[i],ey0[i]));

	v=o0[i+1]-d(i)*ez0[i+1]-o[i];
	dx1(i)=sProduct(v,ex[i]); dy1(i)=sProduct(v,ey[i]); dz1(i)=sProduct(v,ez[i]);
	rotate_angle3(rox1(i),roy1(i),roz1(i),
	              sProduct(ez0[i+1],ex[i]),sProduct(ez0[i+1],ey[i]),sProduct(ez0[i+1],ez[i]),
				  sProduct(ey0[i+1],ex[i]),sProduct(ey0[i+1],ey[i]));

	// �ʌ`����߂�
	SetSpheroid(i,'y',sProduct(A,ez[i])<0 ? 1:0,a,ff);
}


void cLens1::SetHyperboloid(int i,double f1,double f2){
	// i     : i�ʂ�o�Ȗʂɂ���D
	// f1,f2 : i�ʃ��[�J�����W���_����e�œ_�܂ł̋����D�œ_��z����ɂ���Ƃ���D
	
	// �o�Ȗʂ̕�������
	//    (z/a)^2 - (y/b)^2 =1
	//        a : z=0 (i�ʃ��[�J�����W�ł͂Ȃ�) �̓_����o�Ȗʒ��_�܂ł̋���
	//        b : b^2=c^2-a^2,  c�͏œ_�Ԃ̋����̔����D
	// �����z�ɂ��ĉ����ƁC
	//     z = �}a �}(a/2/b^2)y^2 ...
	// ����� ���ʃT�O�̓W�J�C
	//     z = (1/2/r)y^2 + ...
	// �Ɣ�r����΁C�ߎ��ȗ����a r = �} b^2/a �ƂȂ�D
	// r�̕����́C f1+f2=0 �̂Ƃ� r=�� (���ʋ��ɂ�錋��)
	//             f1+f2>0 �̂Ƃ�(���̑��̏œ_�̕�������) r<0
	//             f1+f2<0 �̂Ƃ�(���̑��̏œ_�̕�������) r>0
	// �R�[�j�b�N�W��kp�́C
	//     kp = -(b/a)^2 -1  (cLens1::surface_sag �̒��߂��Q�Ɓj
	
	double a,b,c;

	if(f1*f2>=0) return;  // �o�Ȗʂ�2�œ_�͎��Ƌ��ł���̂ŁCf1*f2�͕��ƂȂ�D

	c=(fabs(f1)+fabs(f2))/2;  // c���œ_�ԋ����̔����Ƃ���D
	a= fabs(f1) < fabs(f2) ? c-fabs(f1) : c-fabs(f2);  // ���m
	b= sqrt(c*c-a*a);                                  // ���m

	if(f1+f2>0)      r(i)=-b*b/a;
	else if(f1+f2<0) r(i)= b*b/a;
	else             r(i)= 0;      // f1+f2==0

	asph_type(i)=1;
	kp(i)= -b*b/a/a -1;
	a4(i)=a6(i)=a8(i)=a10(i)=a12(i)=a14(i)=a16(i)=0;
	EAy(i)=fabs(f1)+fabs(f2);   //  �b�� 2016.02.18
}

void cLens1::SetHyperboloid2(int i,int i1,int i2){
	// ��i�ʂ̌`����i1�ʂƑ�i2�ʂ̊e���_���œ_�Ƃ���o�Ȗʂɂ���D
	// �o�Ȗʂ͖{�֐����s�O�̑�i�ʃ��[�J�����W�̌��_��ʂ���̂Ƃ���D

	double a,b,c, a0,b0, f1,f2;
	vector<double> P,v;

	make_coordinate(0);

	a=abs(o[i1]-o[i]);   // i1�ʌ��_��i�ʌ��_�̊Ԋu
	b=abs(o[i2]-o[i]);   // i2�ʌ��_��i�ʌ��_�̊Ԋu
	c=abs(o[i1]-o[i2]);  // i1�ʌ��_��i2�ʌ��_�̊Ԋu = �œ_�Ԃ̋���
	P=(o[i2]-o[i1])/abs(o[i2]-o[i1]);   // i1�ʌ��_����i2�ʌ��_�֌������P�ʃx�N�g��

	// o[i1]��o[i2]�����Ԑ��Ƒo�Ȗʂ̌����_��P�Ƃ��C
	// o[i1]��P�̊Ԋu��a0, o[i2]��P�̊Ԋu��b0�Ƃ���ƁC
	// �A�������� a0-b0=a-b, a0+b0=c ���C
	a0 = ( a-b+c )/2;
	b0 = a0-a+b;

	o[i]=(b0*o[i1]+a0*o[i2])/(a0+b0); // o[i]�����(o[i1],o[i2]) �� a0:b0 �����_(�o�Ȗʂ̒��_)�Ɉڂ�
	if(sProduct(ez0[i],P)>0){         // P (o[i1]��o[i2]) ��ez0[i]�ɉ����Ă��邩�D
		ez[i]=P;
		f1=b0;     // f1,f2�͏œ_��z���W�i�����t���j
		f2=-a0;
	}
	else{
		ez[i]=-P;
		f1=a0;
		f2=-b0;
	}
	ez[i]= sProduct(ez0[i],P)>0 ? P : -P;      // �Ȃ�ׂ��ΐS�ʂ����Ȃ�����l��ez[i]��ݒ�
	ey[i]=ey0[i]-ez[i]*sProduct(ey0[i],ez[i]); // ey0[i]��ez[i]���܂ޖʓ��ŁCey0[i]�ɉ���������ez[i]�Ɛ����ɂ���D
	ey[i]=ey[i]/abs(ey[i]);
	ex[i]=vProduct(ey[i],ez[i]);               // ez[i]��ey[i]�����܂�΁Cex[i]�����܂�D

	// o[i],ez[i],ey[i],ex[i]����i�ʂ̕ΐS�ʂ����߂�
	decenter_type(i)=1;

	v=o[i]-o0[i];
	dx(i)=sProduct(v,ex0[i]); dy(i)=sProduct(v,ey0[i]); dz(i)=sProduct(v,ez0[i]);
	rotate_angle3(rox(i),roy(i),roz(i),
	              sProduct(ez[i],ex0[i]),sProduct(ez[i],ey0[i]),sProduct(ez[i],ez0[i]),
				  sProduct(ey[i],ex0[i]),sProduct(ey[i],ey0[i]));

	v=o0[i+1]-d(i)*ez0[i+1]-o[i];
	dx1(i)=sProduct(v,ex[i]); dy1(i)=sProduct(v,ey[i]); dz1(i)=sProduct(v,ez[i]);
	rotate_angle3(rox1(i),roy1(i),roz1(i),
	              sProduct(ez0[i+1],ex[i]),sProduct(ez0[i+1],ey[i]),sProduct(ez0[i+1],ez[i]),
				  sProduct(ey0[i+1],ex[i]),sProduct(ey0[i+1],ey[i]));

	// �ʌ`����߂�
	SetHyperboloid(i,f1,f2);
}

void cLens1::SetParaboloid(int i,double f){
	// i : i�ʂ�����ʂɂ���D
	// f : i�ʃ��[�J�����W���_����œ_�܂ł̋���(��������)�D�œ_��z����ɂ���Ƃ���D
	
	r(i)=f*2;
	asph_type(i)=1;
	kp(i)=-1;
	a4(i)=a6(i)=a8(i)=a10(i)=a12(i)=a14(i)=a16(i)=0;
	EAy(i)=fabs(f*2);   //  �b�� 2016.02.19
}

void cLens1::SetParaboloid2(int i,int i_f,int i_inf){
	// ��i�ʂ̌`���i_f�ʒ��_���œ_�Ƃ��C
	// ����(��]�Ώ̎�)��i�ʒ��_��i_inf�ʒ��_�����Ԑ��ƕ��s�ȕ����ʂƂ���D
	// �����ʂ͖{�֐����s�O��i�ʒ��_��ʂ���̂Ƃ���D
	// ���͓���邪�C�ȗ����a�̑傫�������̗p����D

	// �����ʂ̐����F
	//    �œ_F����ʏ�̓_P�܂ł̋�����a�Ƃ���D
	//    P��������ƕ��s��a�����i�񂾈ʒu(����2����)��
	//    �����ɉ�����F����̋�����|2f|�ł���Df�͏œ_�����D

	double a, f;
	vector<double> P,Q, v;

	make_coordinate(0);

	P=(o[i_inf]-o[i])/abs(o[i_inf]-o[i]);  // �����ʂ̌��������P�ʃx�N�g��
	Q=(o[i]-o[i_f])/abs(o[i]-o[i_f]);      // �œ_��������ʏ�̓_�֌������P�ʃx�N�g��
	a=abs(o[i]-o[i_f]);                    // �œ_�ƕ����ʏ�̓_�̊Ԋu

	f=( a*fabs(sProduct(Q,P)) + a )/2;  // �œ_�����̐�Βl�D��L�g�����ʂ̐����h���D
	                                    // (-a*fabs(sProduct(Q,P))+a)/2 �����������������̗p����D

	o[i]=o[i_f]+f*sgn(sProduct(Q,P))*P;   // �����ʂ̒��_�̈ʒu

	ez[i]= sProduct(ez0[i],P)>0 ? P : -P;      // �Ȃ�ׂ��ΐS�ʂ����Ȃ�����l��ez[i]��ݒ�
	ey[i]=ey0[i]-ez[i]*sProduct(ey0[i],ez[i]); // ey0[i]��ez[i]���܂ޖʓ��ŁCey0[i]�ɉ���������ez[i]�Ɛ����ɂ���D
	ey[i]=ey[i]/abs(ey[i]);
	ex[i]=vProduct(ey[i],ez[i]);               // ez[i]��ey[i]�����܂�΁Cex[i]�����܂�D

	if(sProduct(ez[i],o[i_f]-o[i])>0) f=f; else f=-f;

	// o[i],ez[i],ey[i],ex[i]����i�ʂ̕ΐS�ʂ����߂�
	decenter_type(i)=1;

	v=o[i]-o0[i];
	dx(i)=sProduct(v,ex0[i]); dy(i)=sProduct(v,ey0[i]); dz(i)=sProduct(v,ez0[i]);
	rotate_angle3(rox(i),roy(i),roz(i),
	              sProduct(ez[i],ex0[i]),sProduct(ez[i],ey0[i]),sProduct(ez[i],ez0[i]),
				  sProduct(ey[i],ex0[i]),sProduct(ey[i],ey0[i]));

	v=o0[i+1]-d(i)*ez0[i+1]-o[i];
	dx1(i)=sProduct(v,ex[i]); dy1(i)=sProduct(v,ey[i]); dz1(i)=sProduct(v,ez[i]);
	rotate_angle3(rox1(i),roy1(i),roz1(i),
	              sProduct(ez0[i+1],ex[i]),sProduct(ez0[i+1],ey[i]),sProduct(ez0[i+1],ez[i]),
				  sProduct(ey0[i+1],ex[i]),sProduct(ey0[i+1],ey[i]));

	// �ʌ`����߂�
	SetParaboloid(i,f);
}

void cLens1::Chamfer(int i,double Wx,double Wy,double Cx,double Cy,double Offset){
	// Wx x Wy �̋�`�ɑ΂���x�����̕�Cx�Cy�����̕�Cy �̖ʎ���
	// ��i�ʂ̗L���͈͂Ɏ{���DWx x Wy ���L���͈͂Ƃ͌���Ȃ��D
	// �ʎ�蕔���̌������Ւf����悤�ȗL���͈͂�����2�ʂ��i�ʂ̒���ɑ}�����邱�Ƃɂ��D
	// Offset��0�ȊO�̕����w�肷��΁C�L���͈͖͂ʎ�肩�炻�̕������I�t�Z�b�g�����D
	// ����ɂ��Ⴆ�΃R�[�g�̔�L�����𔽉f�ł���D

	double th, a,b,l,m, p;
	
	if(i<1 || k<i) return;
	if(Cx==0 && Cy==0) return;
	if(Cx==0 && Cy!=0) Cx=Cy;
	if(Cx!=0 && Cy==0) Cy=Cx;

	if(rem(i+1)!="chamfer" || rem(i+2)!="chamfer"){
		Add(i,2,0);
	}
	rem(i+1)=rem(i+2)="chamfer";
	decenter_type(i+1)=2;
	decenter_type(i+2)=2;
	switch(decenter_type(i)){
		case 1:
		case 2:
		case 3:
			decenter_type(i)=1;
			decenter_type(i+2)=1;
			ret(i+2)=i;
			dx1(i+2)=dx1(i); dx1(i)=0;
			dy1(i+2)=dy1(i); dy1(i)=0;
			dz1(i+2)=dz1(i); dz1(i)=0;
			rox1(i+2)=rox1(i); rox1(i)=0;
			roy1(i+2)=roy1(i); roy1(i)=0;
			roz1(i+2)=roz1(i); roz1(i)=0;
			break;
	}

	th=atan(Cy/Cx)*180/PI;

	// �E��̖ʎ��ƈ�v���钼���́C
	//   �_(a,b)=(Wx/2-Cx,Wy/2) ��ʂ�C�X���� -Cy/Cx �ł���D
	//   ���Ȃ킿�C�����̕����]��(l,m)=(Cx/��(Cx^2+Cy^2), -Cy/��(Cx^2+Cy^2)) �ƂȂ�D
	a=Wx/2-Cx;
	b=Wy/2;
	l= Cx/sqrt(Cx*Cx+Cy*Cy);
	m=-Cy/sqrt(Cx*Cx+Cy*Cy);
	// "���`�㐔(���c���s)" �� "7.3 �����ƕ���" �� "��3" ���C���_�ƒ����̋���p�����܂�D
	p=sqrt(a*a +b*b -(-l*a-m*b)*(-l*a-m*b));
	p=p-Offset;
	p=p*2;

	roz(i+1)=th;
	roz(i+2)=-th;
	EAy(i+1)=EAy(i+2)=p;                 // EAy�ɂ��ʎ�蕔���̌������Ղ�D
	EAx(i+1)=EAx(i+2)=sqrt(Wx*Wx+Wy*Wy); // �Ίp���ȏ゠���EAx���������Ղ邱�Ƃ͂Ȃ��D
	EAtype(i+1)=EAtype(i+2)=1;           // ��`�ɐݒ�
}

double cLens1::DconSetRn(int i,double val){
	Dcon(i).SetRn(val);
	return Dcon(i).Rn;
}

double cLens1::DconSetRn(int i){
	return DconSetRn(i,EAy(i)/2);
}

void cLens1::DconToPSeries(int i){
	if(asph_type(i)==CONIC){
		a4(i)+=Dcon(i).a4() *pow(NormH(i),4);
		a6(i)+=Dcon(i).a6() *pow(NormH(i),6);
		a8(i)+=Dcon(i).a8() *pow(NormH(i),8);
		a10(i)+=Dcon(i).a10() *pow(NormH(i),10);
		a12(i)+=Dcon(i).a12() *pow(NormH(i),12);
		a14(i)+=Dcon(i).a14() *pow(NormH(i),14);
		a16(i)+=Dcon(i).a16() *pow(NormH(i),16);
		a18(i)+=Dcon(i).a18() *pow(NormH(i),18);
		a20(i)+=Dcon(i).a20() *pow(NormH(i),20);
	}
	else{
		a4(i)=Dcon(i).a4() *pow(NormH(i),4);
		a6(i)=Dcon(i).a6() *pow(NormH(i),6);
		a8(i)=Dcon(i).a8() *pow(NormH(i),8);
		a10(i)=Dcon(i).a10() *pow(NormH(i),10);
		a12(i)=Dcon(i).a12() *pow(NormH(i),12);
		a14(i)=Dcon(i).a14() *pow(NormH(i),14);
		a16(i)=Dcon(i).a16() *pow(NormH(i),16);
		a18(i)=Dcon(i).a18() *pow(NormH(i),18);
		a20(i)=Dcon(i).a20() *pow(NormH(i),20);
	}
	Dcon(i)=cDcon();
}

void cLens1::PSeriesToDcon(int i){
	const int M=8;
	double a[M+1];

	if(DconRn(i)==0) DconSetRn(i);
	a[0]=a4(i) /pow(NormH(i),4);
	a[1]=a6(i) /pow(NormH(i),6);
	a[2]=a8(i) /pow(NormH(i),8);
	a[3]=a10(i) /pow(NormH(i),10);
	a[4]=a12(i) /pow(NormH(i),12);
	a[5]=a14(i) /pow(NormH(i),14);
	a[6]=a16(i) /pow(NormH(i),16);
	a[7]=a18(i) /pow(NormH(i),18);
	a[8]=a20(i) /pow(NormH(i),20);
	Dcon(i).PSeriesToDcon(20,a);
	a4(i)=a6(i)=a8(i)=a10(i)=a12(i)=a14(i)=a16(i)=a18(i)=a20(i)=0;
}

void cLens1::PSeriesToFreeForm(int i){
	b(i, 4,0)+= a4(i)/pow(NormH(i), 4); b(i,0, 4)+= a4(i)/pow(NormH(i), 4);
	b(i, 6,0)+= a6(i)/pow(NormH(i), 6); b(i,0, 6)+= a6(i)/pow(NormH(i), 6);
	b(i, 8,0)+= a8(i)/pow(NormH(i), 8); b(i,0, 8)+= a8(i)/pow(NormH(i), 8);
	b(i,10,0)+=a10(i)/pow(NormH(i),10); b(i,0,10)+=a10(i)/pow(NormH(i),10);
	b(i,12,0)+=a12(i)/pow(NormH(i),12); b(i,0,12)+=a12(i)/pow(NormH(i),12);
	b(i,14,0)+=a14(i)/pow(NormH(i),14); b(i,0,14)+=a14(i)/pow(NormH(i),14);
	b(i,16,0)+=a16(i)/pow(NormH(i),16); b(i,0,16)+=a16(i)/pow(NormH(i),16);
	b(i,18,0)+=a18(i)/pow(NormH(i),18); b(i,0,18)+=a18(i)/pow(NormH(i),18);
	b(i,20,0)+=a20(i)/pow(NormH(i),20); b(i,0,20)+=a20(i)/pow(NormH(i),20);
	a4(i)=a6(i)=a8(i)=a10(i)=a12(i)=a14(i)=a16(i)=a18(i)=a20(i)=0;
}

void cLens1::DconToFreeForm(int i){
	DconToPSeries(i);
	PSeriesToFreeForm(i);
}

void cLens1::LegendreToFreeForm(int i){
	int p,q;
	int LegendreMaxOrder=10;

	for(p=0; p<=LegendreMaxOrder; ++p) for(q=0; q<=LegendreMaxOrder; ++q){
		b(i,p,q)+=legendre(i).Z(1,1,p,q);
	}

	legendre(i)=cModLegendre();
}

void cLens1::SetSpline(int i,int n,double step){
	// ��i�ʂɁC�_��n�C�_�Ԋustep �ŃT�O��0�ɏ����������X�v���C��������ǉ�����
	int j;

	surf[i].spline.SetN(n);
	for(j=1; j<=n; ++j){
		surf[i].spline.SetX(j,step*(j-1));
		surf[i].spline.SetY(j,0);
	}
	if(asph_type(i)==0) asph_type(i)=CONIC;
}

int cLens1::SplineN(int i){
	// ��i�ʂ̃X�v���C���f�[�^����Ԃ�
	if(0<=i && i<=k+1){
		return is_spline_surface(i) ? surf[i].spline.GetN() : 0;
	}
	else{
		return 0;
	}
}

void cLens1::SplineDoubleN(int i){
	if(0<=i && i<=k+1){
		if(is_spline_surface(i)) surf[i].spline.DoubleN();
	}
}

void cLens1::SplineSetHStep(int i,double dh){
	if(0<=i && i<=k+1){
		if(is_spline_surface(i)) surf[i].spline.SetXStep(dh);
	}
}

double cLens1::dc(int i,double dz){
	// ��i�ʂŁC�L���͈͂̒[�ŃT�O��dz�����ω������郢c��Ԃ��D
	double h,u,v,u1,v1,z1,p;

	h=ea_max(i)/2;                       //  �L���͈͂̒[(���a)	
	p=1-(kp(i)+1)*c(i)*c(i)*h*h;

	if(p>0){
		u=c(i)*h*h;          
		v=1+sqrt(p);                    // z=u/v
		u1=h*h;                         // u1=du/dc
		v1=-(kp(i)+1)*c(i)*h*h/sqrt(p); // v1=dv/dc
		z1=(u1*v-u*v1)/v/v;             // z1=dz/dc
	}
	else{
		z1=h*h;
	}

	return dz/z1;       // dz=(dz/dc)dc
}


double cLens1::dkp(int i,double dz){
	// ��i�ʂŁC�L���͈͂̒[�ŃT�O��dz�����ω������郢kp��Ԃ��D
	double h,u,v,v1,z1,p;

	h=ea_max(i)/2;                       //  �L���͈͂̒[(���a)	
	p=1-(kp(i)+1)*c(i)*c(i)*h*h;

	
	if(p>0){
		u=c(i)*h*h;                   
		v=1+sqrt(p);                     // z=u/v
		v1=-c(i)*c(i)*h*h/2/sqrt(p);     // v1=dv/dkp
		z1=-u*v1/v/v;                    // z1=dz/dkp
	}
	else{
		// p=1-(kp+1)cchh=0 �� z=chh ���
		// z = h/��(kp+1)
		// h��萔�Ƃ���kp�Ŕ�������΁C
		z1=-h/2/pow(kp(i)+1,1.5);        // z1=dz/dkp
	}

	return dz/z1;                        // dz=(dz/dkp)dkp
}

double cLens1::da(int i,int n,double dz){
	// ��i�ʂŁCa*h^n = dz �̂Ƃ��� a�����߂� (h�͗L���͈͔��a)
	// �����݌v�Ŕ񋅖ʌW�� a4(i),a6(i),a8(i) ... ��ω�������X�e�b�v�Ƃ��ėp����D
	return dz/pow(ea_max(i)/2/NormH(i),n);
}

cLens1 cLens1::SwapObjPupil(){
	// ���̂Ɠ�����������D
	// ���ӁF2����s���Ă����ɖ߂�Ƃ͌���Ȃ��D
	//       ��j���s��Cs1fix�͕K��0�ɂȂ�D
	double t1fix,EPD_p,yObjectMax_p;
	int Afocal_p;

	t1fix=0;
	Afocal_p=abs(t1())>1000 ? 1 : 0;
	EPD_p=yObjectMax*2;
	yObjectMax_p=EPD/2;
	
	Swap(s,t);
	Swap(s1fix,t1fix);
	Swap(Afocal,Afocal_p);
	Swap(EPD,EPD_p); EPDx=0;
	Swap(yObjectMax,yObjectMax_p);

	return *this;
}

cLens1 cLens1::Deleted(int i1,int i2){
	int i, di;

	if(this->k<i2) i2=this->k;  // ��i1�ʈȍ~�S�č폜����Ƃ�, i2��this->k�𒲂ׂȂ��Ă��傫�����ł悢�D
	if( i1<1 || i2<i1 ) return *this;
	di=i2-i1+1;
	cLens1 x(this->k-di,cn);
	x.assignment(*this);
	for(i=i1; i<=x.k+1; i++){
		x.surf[i]=surf[i+di];
		if(x.ret(i)!=0) x.ret(i)-=di;
		x.med[i-1]=med[i+di-1];        // �폜���̉E���̔}�����g���D
	}

	if(i2==this->k){
		// �ŏI�ʂ܂ō폜����Ƃ��̂݁C����Ԕ}���͍폜���̍����̔}�����g���D�������Ȃ��ƁC
		//   �E�͌^�Ⴊ�Ō�ɕt�����n�ł͍폜��ɑ���Ԃ��Ɏq�̂̋��ܗ��ɂȂ��Ă��܂��D
		//   �E�~���[�n�ł͍폜��ɑ���Ԃ̋��ܗ������O�̋��ܗ��ƈٕ����ƂȂ��Ă��܂����Ƃ�����D
		x.med[x.k]=x.med[x.k+1]=med[i1-1];
	}

	if(1<i1 && i1<this->k){
		for(i=i1-1; i<=i2-1 && i<this->k; i++){
			// ���w�n�S����s�ςɂ���
			x.d(i1-1)+=d(i);
		}
	}

	if(i1==1){
		// ��1�ʂ��܂�ō폜����ꍇ�Cs,t���C������D
		double s,t;
		s=this->si(this->s,i2+1,1);
		t=this->si(this->t,i2+1,1);
		x.s=s;
		x.t=t;
	}

	if(i2==k){
		// �ŏI�ʂ��܂�ō폜����ꍇ
		s1fix=0;
		rImage()=0;
	}

	if(x.stop>i2) x.stop-=i2-i1+1;
	return x;
}

void cLens1::Delete(int i1,int i2){
	*this=Deleted(i1,i2);
}

void cLens1::DeleteBetween(int i){
	// ��i+1�ʂ���C���� rem(i)==rem(i1) �ƂȂ�i1�ʂ܂ł��폜����D
	//    �� �������Crem(i1)��ۑ����C�폜�� rem(i)�ɑ������D
	//       rem(i1)��������ƁCi1�ʂ��܂� DoubleTurn �ŁCsurf_no(rem(i1)) �������ɂł��Ȃ����߁D
	//       
	// ���̂Ƃ��C�����̔��f��͑O����v�Ƃ���i Turn �� DoubleTurn ���s�ɂ�� <2>,<3> .. �̕t���ɑΉ����邽�߁j
	// i�ʂ�i1�ʂ̊Ԃɔ��˖ʂ�����ꍇ�Ci�ʂ𔽎˖ʂƂ���n�����������D
	// �Ⴆ�΁CSLO�őΕ������Y�e�ʂ̔��˂���������̂Ɏg���D
	int ii,i1;
	std::string s;

	if(i<=0 || i>=k+1) return;
	i1=0;
	for(ii=i+1; ii<=k; ++ii){
		if(forward_match(rem(i),rem(ii))){
			i1=ii;
			break;
		}
	}
	if(i1==0) return;
	s=rem(i1);
	Delete(i+1,i1);
	rem(i)=s;  
	if(decenter_type(i)!=0) decenter_type(i)=2;


}

cLens1 cLens1::Trimed(int i1,int i2){
	cLens1 x=*this;

	if(1<=i1 && i1<=i2 && i2<=x.k){
		x.Delete(i2+1,x.k);
		x.med[i2]=this->med[i2];  // �c�������̒���̔}�����g���D
		x.Delete(1,i1-1);         // ���Fs,t, ret ��Delete()�ŏC�������D
		x.surf[0]=x.surf[x.k+1]=surface();  // ���̖ʁC���ʂ̕ΐS���͏�������
	}
	return x;
}

void cLens1::Trim(int i1,int i2){
	*this=Trimed(i1,i2);
}

cLens1 cLens1::Added(int pre_surf,int n,double z){
	// ��pre_surf�ʂ̌��z�̈ʒu��n�ʑ}������
	int i;

	if(pre_surf<0 || k<pre_surf) return *this;
	if(n<=0) return *this;
	
	cLens1 x(k+n,cn);
	x.assignment(*this);
	for(i=x.k+1; i>=pre_surf+n+1; --i){
		x.surf[i]=x.surf[i-n];
		if(x.surf[i].ret!=0 && x.surf[i].ret>pre_surf) x.surf[i].ret+=n;
	}
	for(i=x.k; i>=pre_surf+n+1; --i){
		x.med[i]=x.med[i-n];
	}
	for(i=pre_surf+1; i<=pre_surf+n; ++i){
		x.surf[i]=surface();
		x.EAy(i)=Max( fabs(x.EAy(pre_surf)),fabs(x.EAy(pre_surf+n+1)) ) *2; // *2�ɂ�����������}����
		x.med[i]=medium();
		x.gname(i)=gname(pre_surf);
	}
	if(pre_surf==0){
		if(fabs(s)<LN) x.s=-z;
	}
	else{
		x.d(pre_surf)=z;
	}
	for(i=pre_surf+1; i<=pre_surf+n-1; ++i) x.d(i)=0;   // �}�������̖ʊԊu�̓f�t�H���g��0
	if(pre_surf==k){
		if( s1fix!=0 && !Afocal ) x.s1fix=s1fix-z;
	}
	else{
		if(pre_surf==0){
			if(fabs(x.s)<LN){
				x.d(n)=-s-z;
			}
		}
		else{
			x.d(pre_surf+n)=d(pre_surf)-z;
		}
	}

	if(x.stop>pre_surf) x.stop+=n;
	return x;
}

void cLens1::Add(int pre_surf,int n,double z){
	*this=Added(pre_surf,n,z);
}

void cLens1::ToBlock(int i1,int i2){
	// i1�ʂ���i2�ʂ܂ł̕������u���b�N������D
	// i1�ʂ̑O��i2�ʂ̌�ɉˋ�ʂ��}�������D
	if(1<=i1 && i2<=k){
		Add(i2,1,0);
		if(i1==1) Add(0,1,-s); else Add(i1-1,1,d(i1-1));

		decenter_type(i2+2)=1;
		ret(i2+2)=i1;
		rem(i2+2)="block_out";
		decenter_type(i1)=1;
		rem(i1)="block_in";
	}
}

void cLens1::BreakCement(int i,std::string AdhesiveName/*="1"*/){
	// ��i�ʂ�����0�̋�C�w������2�ʂɕ�������i�ڍ��𔍂����j�D
	Add(i,1,0);
	surf[i+1]=surf[i];
	if(AdhesiveName=="") AdhesiveName="1";
	gname(i)= Re(Index(gname(i-1),wl[1]))>0 ? AdhesiveName : "-"+AdhesiveName;
}

cLens1 cLens1::Reversed(int z_reverse){
	// �ŏI�ʂ���t�ǐՂ��邽�߂̌��w�f�[�^�����.
	// ���Ȃ킿�C�ʔԍ����ŏI�ʂ���1,2, ... �ƂȂ�D
	int i; double temp;
	double dx,dy,dz,rox,roy,roz, dx1,dy1,dz1,rox1,roy1,roz1;
	int order,order1;
	cLens1 x0=*this, x=*this;
	x0.TransformRetDecenter();
	if(z_reverse){
		// z�����𔽓]����D
		for(i=0; i<=x.k+1; ++i) x.surf[i]=x0.surf[k+1-i].reversed();
		for(i=0; i<=x.k; ++i) x.med[i]=x0.med[k-i];

		for(i=0; i<=x.k+1; ++i){
			switch(x.decenter_type(i)){
			case 1:
				dx=-x.dx1(i);
				dy= x.dy1(i);
				dz=-x.dz1(i);
				rox=-x.rox1(i);
				roy= x.roy1(i);
				roz=-x.roz1(i);
				order=x.order1(i);
				dx1=-x.dx(i);
				dy1= x.dy(i);
				dz1=-x.dz(i);
				rox1=-x.rox(i);
				roy1= x.roy(i);
				roz1=-x.roz(i);
				order1=x.order(i);
				InverseDecenter(dx,dy,dz,rox,roy,roz,order);
				x.dx(i)=dx;
				x.dy(i)=dy;
				x.dz(i)=dz;
				x.rox(i)=rox;
				x.roy(i)=roy;
				x.roz(i)=roz;
				x.order(i)=order;
				InverseDecenter(dx1,dy1,dz1,rox1,roy1,roz1,order1);
				x.dx1(i)=dx1;
				x.dy1(i)=dy1;
				x.dz1(i)=dz1;
				x.rox1(i)=rox1;
				x.roy1(i)=roy1;
				x.roz1(i)=roz1;
				x.order1(i)=order1;
				break;
			case 2:
				x.dx(i)=-x.dx(i);
				x.dz(i)=-x.dz(i);
				x.rox(i)=-x.rox(i);
				x.roz(i)=-x.roz(i);
				break;
			case 3:
				InverseDecenter(x.dx(i),x.dy(i),x.dz(i),x.rox(i),x.roy(i),x.roz(i),x.order(i));
				x.dx(i)=-x.dx(i);
				x.dz(i)=-x.dz(i);
				x.rox(i)=-x.rox(i);
				x.roz(i)=-x.roz(i);
				break;
			}
		}
		x.s= s1fix ? -s1fix : -s1();
		x.s= fabs(x.s)>=LN ? LN*10 : x.s;
		temp=s; s=t; x.t=-s1(); s=temp; 
		x.t= fabs(x.t)>=LN ? LN*10 : x.t;
		x.EPD=fabs(EPD*Mpupil());
		x.EPDx=fabs(EPDx*Mpupil());
		x.EPx=EPx*Mpupil();
		x.EPy=EPy*Mpupil();
		if(x.stop!=0) x.stop=k-stop+1;
		x.s1fix=0;
		x.Afocal= fabs(s)>=LN ? 1 : 0;
	}
	else{
		// z�����𔽓]���Ȃ��D���ܗ��𔽓]����D
		for(i=0; i<=x.k+1; ++i) x.surf[i]=x0.surf[k+1-i];
		for(i=0; i<=x.k; ++i) x.med[i]=x0.med[k-i].reversed();
		for(i=0; i<=x.k+1; ++i){
			switch(x.decenter_type(i)){
			case 1:
				dx=x.dx1(i); dy=x.dy1(i); dz=x.dz1(i); rox=x.rox1(i); roy=x.roy1(i); roz=x.roz1(i);
				order=x.order1(i);
				dx1=x.dx(i); dy1=x.dy(i); dz1=x.dz(i); rox1=x.rox(i); roy1=x.roy(i); roz1=x.roz(i);
				order1=x.order(i);
				InverseDecenter(dx,dy,dz,rox,roy,roz,order);
				x.dx(i)=dx;
				x.dy(i)=dy;
				x.dz(i)=dz;
				x.rox(i)=rox;
				x.roy(i)=roy;
				x.roz(i)=roz;
				x.order(i)=order;
				InverseDecenter(dx1,dy1,dz1,rox1,roy1,roz1,order1);
				x.dx1(i)=dx1;
				x.dy1(i)=dy1;
				x.dz1(i)=dz1;
				x.rox1(i)=rox1;
				x.roy1(i)=roy1;
				x.roz1(i)=roz1;
				x.order1(i)=order1;
				break;
			case 2:
				break;
			case 3:
				InverseDecenter(x.dx(i),x.dy(i),x.dz(i),x.rox(i),x.roy(i),x.roz(i),x.order(i));
				break;
			
			}
		}
		x.s= s1fix ? s1fix : s1();
		x.s= fabs(x.s)>=LN ? LN*10 : x.s;
		temp=s; s=t; x.t=s1(); s=temp; 
		x.t= fabs(x.t)>=LN ? LN*10 : x.t;
		x.EPD=fabs(EPD*Mpupil());
		x.EPDx=fabs(EPDx*Mpupil());
		x.EPx=EPx*Mpupil();
		x.EPy=EPy*Mpupil();
		if(x.stop!=0) x.stop=k-stop+1;
		x.s1fix=0;
		x.Afocal= fabs(s)>=LN ? 1 : 0;
	}
	return x;
}

void cLens1::Reverse(){
	*this=Reversed(1);
}

cLens1 cLens1::zReversed(int ToMirrorImg){
	// z�������𔽓]����D
	// Reverse()�ł͖ʔԍ��𔽓]���邪zReverse()�ł͔��]���Ȃ��D
	// ToMirroImg���^�̂Ƃ��C�n�������Ƃ���D�n���~���[�Ő܂�Ȃ���Ƃ��Ɏg����D
	int i;
	cLens1 x=*this;

	for(i=0; i<=x.k+1; ++i) x.surf[i]=surf[i].reversed();
	for(i=0; i<=x.k; ++i) x.med[i]=med[i].reversed();
	for(i=0; i<=x.k+1; ++i){
		switch(x.decenter_type(i)){
		case 1:
		case 2:
		case 3:
			if(ToMirrorImg){
				x.dx(i)= dx(i);
				x.dy(i)= dy(i);
				x.dz(i)=-dz(i);
				x.rox(i)=-rox(i);
				x.roy(i)=-roy(i);
				x.roz(i)= roz(i);
				x.dx1(i)= dx1(i);
				x.dy1(i)= dy1(i);
				x.dz1(i)=-dz1(i);
				x.rox1(i)=-rox1(i);
				x.roy1(i)=-roy1(i);
				x.roz1(i)= roz1(i);
			}
			else{
				x.dx(i)=-dx(i);
				x.dy(i)= dy(i);
				x.dz(i)=-dz(i);
				x.rox(i)=-rox(i);
				x.roy(i)= roy(i);
				x.roz(i)=-roz(i);
				x.dx1(i)=-dx1(i);
				x.dy1(i)= dy1(i);
				x.dz1(i)=-dz1(i);
				x.rox1(i)=-rox1(i);
				x.roy1(i)= roy1(i);
				x.roz1(i)=-roz1(i);
			}
			break;
		}
	}
	x.s=-s;
	x.t=-t;
	x.s1fix=-s1fix;
	
	return x;
}

void cLens1::zReverse(int ToMirrorImg){
	*this=zReversed(0);
}

void cLens1::Replace(int i1,int i2,cLens1 X,int i1o/*=0*/,int i2o/*=0*/,int AdjustD/*=1*/){
	// i1�ʂ���i2�ʂ̊Ԃ̕�����X��i1o�ʂ���i2o�ʂ܂ł̕����Œu������D
	// i1�ʂ̑O�C�����i2�ʂ̌�̔}���͕ς��Ȃ��D
	// AdjustD���^�̂Ƃ��́C�u���O��Œu�������̎�_�ʒu�������ɂȂ�悤�ɁC�u�������O��̖ʊԊu�𒲐�����D
	int i;
	cLens1 x;
	medium med1,med2;

	if(i1<1 || i2>k || i1>i2) return;
	
	if(i1o!=0 || i2o!=0){
		X.Trim(i1o,i2o);
	}

	med1=med[i1-1];
	med2=med[i2];
	x=*this;
	
	x.Delete(i1,i2);
	x.Add(i1-1,X.k,0);
	
	for(i=1; i<=X.k; i++){
		x.surf[i1-1+i]=X.surf[i];
		if( X.surf[i].decenter_type==1 && X.surf[i].ret!=0 ){
			x.surf[i1-1+i].ret+=i1-1;
		}
	}
	for(i=1; i<=X.k-1; i++){
		x.med[i1-1+i]=X.med[i];
	}

	x.med[i1-1]=med1;
	x.med[i1+X.k-1]=med2;
	
	if(AdjustD){
		if(i1==1){
			x.s+=x.delta(i1,i1+X.k-1)-this->delta(i1,i2);
		}
		else{
			x.d(i1-1)-=x.delta(i1,i1+X.k-1)-this->delta(i1,i2);
		}

		if(i2==k){
			if(s1fix!=0) x.s1fix+=x.delta1(i1,i1+X.k-1)-this->delta1(i1,i2);
		}
		else{
			x.d(i1+X.k-1)+=x.delta1(i1,i1+X.k-1)-this->delta1(i1,i2);
		}
	}

	*this=x;
}

int cLens1::Replace(int i1,int i2,std::string filename,int i1o/*=0*/,int i2o/*=0*/,int AdjustD/*=1*/){
	cLens1 x;
	std::ifstream from(filename.c_str());
	if(from) {
		from >> x;
		Replace(i1,i2,x,i1o,i2o,AdjustD);
		return 1;
	}
	else return 0;
}

void cLens1::Merge(cLens1 X,double d,int Auto/*=1*/){
	// *this��X����������D*this�̑�k��(�ŏI��)��X�̑�1�ʂ̋�����d�Ƃ���D
	// Auto���^�̂Ƃ��C����Ԃ̋��ܗ������̏ꍇ�ɂ͓K�؂Ǝv���鎩������������D
	// (���Ȃ݂� "auto" �͗\���j
	cLens1 x=*this;

	if(Auto!=0 && x.N(x.k,1)<0){  // *this�̑���Ԃ̋����������̂Ƃ�
		// �E����Ԃ�X�̕���Ԃ̋��ܗ�����v������D
		// �EX�͕���Ԃ̋��ܗ��𐳂Ƃ��ċL�q���Ă���Ɛ��肷��D
		X.zReverse(0);
		d=-d;
	}

	x.Add(x.k,1,d);
	x.Replace(x.k,x.k,X);
	x.d(this->k)=d;   // Replace()�͒u�������������̗�����d����_�ʒu�ŕ␳����̂ŁC���ɖ߂��D
	x.s1fix=X.s1fix;
	x.surf[x.k+1]=X.surf[X.k+1];
	x.med[x.k]=X.med[X.k];
	x.Afocal=X.Afocal;
	*this=x;
}

int cLens1::Merge(std::string filename,double d){
	cLens1 x;
	std::ifstream from(filename.c_str());
	if(from) {
		from >> x;
		Merge(x,d);
		return 1;
	}
	else return 0;
}

cLens1 cLens1::Bended(int i){
	// ��i�ʂŌ��w�n��܂�Ȃ���D
	cLens1 x,x1;
	if( 1<=i && i<=k ){
		x =Deleted(i+1,k);
		x1=Deleted(1,i-1).zReversed(1);
		x.Merge(x1,0,0);
		switch(x.decenter_type(i)){
		case 1:
		case 3:
			x.decenter_type(i)=2;
			x.dx1(i)=x.dy1(i)=x.dz1(i)=x.rox1(i)=x.roy1(i)=x.roz1(i)=0;
			x.ret(i)=0;
			break;
		}
		x.Delete(i+1,i+1);
		x.stop= stop>i ? 0 : stop;
		return x;
	}
	else return *this;		
}

void cLens1::Bend(int i) {
	*this=Bended(i);
}

cLens1 cLens1::Turned(int i){
	// ��i�ʂŔ��˂��t�s����D
	int ii;
	cLens1 x(i,cn),x_rev;

	if( 1<=i && i<=k ){
		x.assignment(*this);
		x_rev=x;
		x_rev=x_rev.Reversed(0);
		for(ii=1; ii<=x_rev.k; ++ii){ if(x_rev.rem(ii)!="") x_rev.rem(ii)+="<2>"; }
		x.Merge(x_rev,0,0);
		switch(x.decenter_type(i)){
		case 1:
		case 3:
			x.decenter_type(i)=2;
			x.dx1(i)=x.dy1(i)=x.dz1(i)=x.rox1(i)=x.roy1(i)=x.roz1(i)=0;
			x.ret(i)=0;
			break;
		}
		x.Delete(i+1,i+1);
		x.stop= stop>i ? 0 : stop;
		return x;
	}
	else return *this;		
}

void cLens1::Turn(int i) {
	*this=Turned(i);
}

cLens1 cLens1::DoubleTurned(int i1,int i2) {
	// i1�ʂŔ��˂�����i2��(i2<i1)�Ŕ��˂��Ăё����֌�����.
	int i,ii;
	cLens1 x(i1,cn),temp;
	cLens1 x0=*this;
//	x0.TransformRetDecenter();   
//  �� ���̍s����(2017.07.10) �X�L���i����ret��������Ǝ�C�����ʓ|�Ȃ̂ŏ���. �������Cret�������ʂ��Q�Ƃ��Ă���Ɛ������Ȃ����ʂƂȂ�\������D
	
	if(1<=i1 && i1<=k && 1<=i2 && i2<i1){
		x.assignment(*this);                  // x= 1�`i2�`i1
		
		temp=x0.Trimed(i2,i1).Reversed(0);
		for(ii=1; ii<=temp.k; ++ii){ if(temp.rem(ii)!="") temp.rem(ii)+="<2>"; }
		x.Merge(temp,0,0);                    // x= 1�`i2�`i1 & i1�`i2
		switch(x.decenter_type(i1)){
		case 1:
		case 3:
			x.decenter_type(i1)=2;
			x.dx1(i1)=x.dy1(i1)=x.dz1(i1)=x.rox1(i1)=x.roy1(i1)=x.roz1(i1)=0;
			x.ret(i1)=0;
			break;
		}
		x.Delete(i1+1,i1+1);                  // x= 1�`i2�`i1 & i1-1�`i2
		
		temp=x0.Trimed(i2,x0.k);
		for(ii=1; ii<=1+i1-i2; ++ii){ if(temp.rem(ii)!="") temp.rem(ii)+="<3>"; }
		x.Merge(temp,0,0);                    // x= 1�`i2�`i1 & i1-1�`i2 & i2�`k
		switch(x.decenter_type(i=i1+(i1-i2))){
		case 1:
		case 3:
			x.decenter_type(i)=2;
			x.dx1(i)=x.dy1(i)=x.dz1(i)=x.rox1(i)=x.roy1(i)=x.roz1(i)=0;
			x.ret(i)=0;
			break;
		}
		x.Delete(i+1,i+1);                    // x=1�`i2�`i1 & i1-1�`i2 & i2+1�`k

		x.stop= i1<stop ? stop+2*(i1-i2) : stop;
		return x;
	}
	else return *this;
}

void cLens1::DoubleTurn(int i1,int i2) {
	*this=DoubleTurned(i1,i2);
}

cLens1 cLens1::TripleTurned(int i1,int i2,int i3) {
	// i1��,i2��(i2<i1),i3��(i2<i3)�̏��ɔ��˂����̑��֌������D
	int i,ii;
	cLens1 x(i1,cn),temp;
	cLens1 x0=*this;
//	x0.TransformRetDecenter();
//  �� ���̍s����(2017.07.10) �X�L���i����ret��������Ǝ�C�����ʓ|�Ȃ̂ŏ���. �������Cret�������ʂ��Q�Ƃ��Ă���Ɛ������Ȃ����ʂƂȂ�\������D

	if(1<=i1 && i1<=k && 1<=i2 && i2<i1 && i2<i3 && i3<=k){
		x.assignment(*this);          // x= 1�`i2�`i1
		
		temp=x0.Trimed(i2,i1).Reversed(0);
		for(ii=1; ii<=temp.k; ++ii){ if(temp.rem(ii)!="") temp.rem(ii)+="<2>"; }
		x.Merge(temp,0,0);            // x= 1�`i2�`i1 & i1�`i2
		switch(x.decenter_type(i1)){
		case 1:
		case 3:
			x.decenter_type(i1)=2;
			x.dx1(i1)=x.dy1(i1)=x.dz1(i1)=x.rox1(i1)=x.roy1(i1)=x.roz1(i1)=0;
			x.ret(i1)=0;
			break;
		}
		x.Delete(i1+1,i1+1);          // x= 1�`i2�`i1 & i1-1�`i2
		
		temp=x0.Trimed(i2,i3);
		for(ii=1; ii<=1+i3-i2; ++ii){ if(temp.rem(ii)!="") temp.rem(ii)+="<3>"; }
		x.Merge(temp,0,0);             // x= 1�`i2�`i1 & i1-1�`i2 & i2�`i3
		switch(x.decenter_type(i=i1+(i1-i2))){
		case 1:
		case 3:
			x.decenter_type(i)=2;
			x.dx1(i)=x.dy1(i)=x.dz1(i)=x.rox1(i)=x.roy1(i)=x.roz1(i)=0;
			x.ret(i)=0;
			break;
		}
		x.Delete(i+1,i+1);             // x= 1�`i2�`i1 & i1-1�`i2 & i2+1�`i3

		temp=x0.Trimed(1,i3).Reversed(0);
		for(ii=1; ii<=temp.k; ++ii){ if(temp.rem(ii)!="") temp.rem(ii)+="<4>"; }
		x.Merge(temp,0,0);             // x= 1�`i2�`i1 & i1-1�`i2 & i2+1�`i3 & i3�`1
		switch(x.decenter_type(i=i1+(i1-i2)+(i3-i2))){
		case 1:
		case 3:
			x.decenter_type(i)=2;
			x.dx1(i)=x.dy1(i)=x.dz1(i)=x.rox1(i)=x.roy1(i)=x.roz1(i)=0;
			x.ret(i)=0;
			break;
		}
		x.Delete(i+1,i+1);              // x= 1�`i2�`i1 & i1-1�`i2 & i2+1�`i3 & i3-1�`1
		
		if(i1<stop && stop<=i3){
			x.stop=stop+(i1-i2)+(i3-i2);
		}
		else if(i1<stop && i3<stop){
			x.stop=1;
		}

		return x;
	}
	else return *this;
}

void cLens1::TripleTurn(int i1,int i2,int i3) {
	*this=TripleTurned(i1,i2,i3);
}

int cLens1::Zoom(int i1,int i2,int i3,double MorF) {
	if( !( 1<=i1 && i1<i2 && i2<i3 && i3<k ) ) return 0;
	matrix<double> H(2,1),H1(2,1);
	matrix<double> S1(2,2),S2(2,2),S3(2,2),S4(2,2);
	matrix<double> P1(2,2),P2(2,2),P3(2,2),P4(2,2),P5(2,2),P6(2,2);
	matrix<double> T(2,2),I(2,2);
	double al1;
	double x,x1,x2,y,y1,y2,t, a1,a2,a3,a4,b1,b2,b3, q1,q2,q3;

	if( fabs(s)>=LN ){    // (2.5.2)
		H.a[1][1]=1;
		H.a[2][1]=0;
		al1=1/MorF;
		H1.a[1][1]=s1(1)*al1/N(k,1);
		H1.a[2][1]=al1;
	}
	else if(s!=0){        // (2.5.3)
		H.a[1][1]=1;
		H.a[2][1]=N(0,1)/s;
		al1=H.a[2][1]/MorF;
		H1.a[1][1]=s1(1)*al1/N(k,1);
		H1.a[2][1]=al1;
	}
	else{                // s==0�̂Ƃ�
		H.a[1][1]=0;
		H.a[2][1]=1;
		al1=H.a[2][1]/MorF;
		H1.a[1][1]=s1(1)*al1/N(k,1);
		H1.a[2][1]=al1;
	}
	
	T.a[1][1]=0; T.a[1][2]=1;  // (2.5.5)
	T.a[2][1]=0; T.a[2][2]=0;
	I.a[1][1]=1; I.a[1][2]=0;  // (2.5.5)
	I.a[2][1]=0; I.a[2][2]=1;
	AMatrixCalc(1,i1);  
	S1.a[1][1]=A2[1]; S1.a[1][2]=A2[2];
	S1.a[2][1]=A2[3]; S1.a[2][2]=A2[4];
	AMatrixCalc(i1+1,i2);
	S2.a[1][1]=A2[1]; S2.a[1][2]=A2[2];
	S2.a[2][1]=A2[3]; S2.a[2][2]=A2[4];
	AMatrixCalc(i2+1,i3);
	S3.a[1][1]=A2[1]; S3.a[1][2]=A2[2];
	S3.a[2][1]=A2[3]; S3.a[2][2]=A2[4];
	AMatrixCalc(i3+1,k);
	S4.a[1][1]=A2[1]; S4.a[1][2]=A2[2];
	S4.a[2][1]=A2[3]; S4.a[2][2]=A2[4];
	P1=inv(S3)*inv(S4)*H1;              // (2.5.10)
	P2=inv(S3)*T*inv(S4)*H1/N(i3,1);    // (2.5.12)
	P3=S2*S1*H;                         // (2.5.13)
	P4=S2*T*S1*H/N(i1,1);               // (2.5.15)
	P5=T*S2*S1*H/N(i2,1);               // (2.5.16)
	P6=T*S2*T*S1*H/N(i1,1)/N(i2,1);     // (2.5.17)
	x=d(i1);
	y=d(i2);
	t=d(i1)+d(i2)+d(i3);
	a1=P1.a[1][1]+t*P2.a[1][1]-P3.a[1][1];  // (2.5.19)
	a2=P4.a[1][1]-P2.a[1][1];               // (2.5.19)
	a3=P5.a[1][1]-P2.a[1][1];               // (2.5.19)
	a4=-P6.a[1][1];                         // (2.5.19)
	b1=P1.a[2][1]+t*P2.a[2][1]-P3.a[2][1];  // (2.5.19)
	b2=P4.a[2][1]-P2.a[2][1];               // (2.5.19)
	b3=-P2.a[2][1];                         // (2.5.19)
	q1=-a4*b2; q2=a2*b3-a3*b2-a4*b1; q3=a1*b3-a3*b1;   // (2.5.24)
	if( q2*q2-4*q1*q3<0 ) return 0;
	x1=( -q2+sqrt( q2*q2-4*q1*q3 ) )/2/q1;  // (2.5.23)
	x2=( -q2-sqrt( q2*q2-4*q1*q3 ) )/2/q1;  // (2.5.23)
	y1=-(b2/b3)*x1-(b1/b3);      // (2.5.22)
	y2=-(b2/b3)*x2-(b1/b3);      // (2.5.22)
	if( x1>0 && y1>0 && t-x1-y1>0 ) { d(i1)=x1; d(i2)=y1; d(i3)=t-x1-y1; return 1; }
	if( x2>0 && y2>0 && t-x2-y2>0 ) { d(i1)=x2; d(i2)=y2; d(i3)=t-x2-y2; return 1; }
	return 0;
}

void cLens1::DispersionToZero(int i1,int i2){
	// ��i1�ʂ�i2�ʂ̊Ԃ̔}���̕��U��0�ɂ���D
	// �e�����̐F�����ւ̊�^�𒲂ׂ�Ƃ��ȂǂɎg���D
	int i;

	if( 0<=i1 && i1<i2 && i2<=k+1 ){
		for(i=i1; i<i2; ++i){
			Set_gname(i,str(N(i,1)));
		}
	}
}

void cLens1::rRound(int n){
	if(n>0){
		for(int i=1; i<=k; i++) r(i)=Round(r(i),-n);
	}
}
void cLens1::dRound(int n){
	if(n>0){
		for(int i=1; i<=k-1; i++) d(i)=Round(d(i),-n);
	}
}
void cLens1::RoundA(int n){
	// �񋅖ʌW�� a4,a6,... ��L������n���Ɋۂ߂�
	if(n>0){
		for(int i=1; i<=k; i++){
			a4(i)=Round(a4(i),n);
			a6(i)=Round(a6(i),n);
			a8(i)=Round(a8(i),n);
			a10(i)=Round(a10(i),n);
			a12(i)=Round(a12(i),n);
			a14(i)=Round(a14(i),n);
			a16(i)=Round(a16(i),n);
			a18(i)=Round(a18(i),n);
			a20(i)=Round(a20(i),n);
		}
	}
}

std::string cLens1::FileName(){
	return remove_path(filename);
}

std::string cLens1::LensData(int colors/*=1*/){
	int i,j; unsigned int len;
	std::string s; double sgn;
	char buf[1000];
	std::string a,format;
	if(colors<0 || colors>cn) colors=cn;
	// i : 3d
	// r : 11.3f
	// * : 1
	// d : 11.3f
	// N : 9.5f
	// one space
	// glassname : -[len+2]s(-:left filled)
	// EA,�� : 8.2fx6.2f��6.2fx6.2f�� (32��)
	//         8.2f��6.2f��           (18��, �~�␳���`�݂̂̂Ƃ��j

	//                     r           d     N(..)  [len+3    ]   EA              ��      
	//           iiirrrrrrr.rrr*ddddddd.dddnnn.nnnnn [(len+2)g]eeeee.eexeee.ee��pppp.ppxppp.pp��
	//           iiirrrrrrr.rrr*ddddddd.dddnnn.nnnnn [(len+2)g]eeeee.ee��       pppp.pp��
	//                     r           d     N(..)                EA       ��
	//           iiirrrrrrr.rrr*ddddddd.dddnnn.nnnnn [(len+2)g]eeeee.ee��pppp.pp��

	len=0;
	for(i=0; i<=k; i++){ // gname(i)�̍ő咷
		if(gname(i).length()>len) len=(unsigned int)gname(i).length();
	}

	sprintf(buf,"          r           d   "); a+=buf;
	for(j=1; j<=colors; j++){
		sprintf(buf,"%-9s", ( "  N("+Get_color(j)+")" ).c_str() ); a+=buf;
	}
	if(ea_is_defined()){
		a+=spc(len+3);
		if(ea_all_regular()){
			sprintf(buf,"   EA       ��"); a+=buf;
		}
		else{
			sprintf(buf,"   EA              ��"); a+=buf;
		}
	}
	sprintf(buf,"\n"); a+=buf;

	for(i=0; i<=k; i++){
		if(i==0) { sprintf(buf,"   "); a+=buf; }
		else{ sprintf(buf,"%3d",i); a+=buf; }

		if(i==0){ 
			sprintf(buf,"            "); a+=buf; 
		}
		else if(r(i)==0){ 
			sprintf(buf,"       ��  "); a+=buf;
			if( is_aspheric(i) ) a+="*"; else a+=" ";
		}
		else{ 
			sprintf(buf,"%11.3f",r(i)); a+=buf; 
			if( is_aspheric(i) ) a+="*"; else a+=" ";
		}

		if(i==0 || i==k){ sprintf(buf,"           "); a+=buf; }
		else{ sprintf(buf,"%11.3f",d(i)); a+=buf; }
		
		for(j=1; j<=colors; j++){ sprintf(buf,"%9.5f",N(i,j)); a+=buf; }

		sprintf(buf," "); a+=buf;
		
		s=Get_gname(i);
		if( s[0]=='-' ){ sgn=-1; s=s.erase(0,1); }

		format="%-"+trim(str(len+2),1)+"s";
		if( s!="1" ){ sprintf(buf,format.c_str(), s.c_str()); a+=buf; }
		else{ a+=spc(len+2); }

		if( i!=0 && ea_is_defined() ){
			// �L���a
			sprintf(buf,"%8.2f",EAy(i)); a+=buf;
			if(EAx(i)!=0){
				sprintf(buf,"x%6.2f",EAx(i)); a+=buf;
				if(EAtype(i)==1){ sprintf(buf,"��"); a+=buf; }
				else            { a+="  "; }
			}
			else{
				if(ea_all_regular()){
					if(EAtype(i)==1){ a+="��"; }
					else            { a+="  "; }
				}
				else{
					if(EAtype(i)==1){ a+="��       "; }
					else            { a+="         "; }
				}
			}
			// ��������܂񂾃����Y�a
			sprintf(buf,"%7.2f",phi_y(i)); a+=buf;
			if(EAx(i)!=0){
				sprintf(buf,"x%6.2f",phi_x(i)); a+=buf;
				if(EAtype(i)==1){ sprintf(buf,"��"); a+=buf; }
				else            { a+="  "; }
			}
			else{
				if(ea_all_regular()){
					if(EAtype(i)==1){ a+="��"; }
					else            { a+="  "; }
				}
				else{
					if(EAtype(i)==1){ a+="��       "; }
					else            { a+="         "; }
				}
			}

		}
		
		a+=" "+rem(i);

		sprintf(buf,"\n"); a+=buf;
	}

	if(r(0)  !=0){ sprintf(buf," rObj=%g\n"  , r(0)  ); a+=buf; }
	if(r(k+1)!=0){ sprintf(buf," rImage=%g\n", r(k+1)); a+=buf; }

	for(i=1; i<=k; i++){
		if(CHMy(i)>0 || CHMx(i)>0){
			sprintf(buf,"no.%2d ChamferY=%g ChamferX=%g\n", i,CHMy(i),CHMx(i));
			a+=buf;
		}
	}

	for(i=1; i<=k; i++){
		if(EAdy(i)!=0 || EAdx(i)!=0){
			sprintf(buf,"no.%2d EAdy=%g EAdx=%g\n", i,EAdy(i),EAdx(i));
			a+=buf;
		}
	}

	for(i=0; i<=k+1; i++){
		if(i==0){
			s="obj  ";
		}
		else if(i==k+1){
			s="img  ";
		}
		else{
			sprintf(buf,"no.%2d",i); s=buf;
		}
		switch(asph_type(i)){
			case CONIC:
				a+=s;
				sprintf(buf," aspherical kp=%g NormH=%g\n",kp(i),NormH(i));
				a+=buf;
				sprintf(buf,"        a4=%.15g a6=%.15g a8=%.15g a10=%.15g a12=%.15g\n"
				        ,a4(i),a6(i),a8(i),a10(i),a12(i));
				a+=buf;
				sprintf(buf,"        a14=%.15g a16=%.15g a18=%.15g a20=%.15g\n"
				        ,a14(i),a16(i),a18(i),a20(i));
				a+=buf;
				if(a1(i)!=0 || a2(i)!=0 || a3(i)!=0 || a5(i)!=0 || a7(i)!=0 || a9(i)!=0 || a11(i)!=0 || a13(i)!=0 || a15(i)!=0){
					sprintf(buf,"        a1=%.15g a2=%.15g a3=%.15g a5=%.15g a7=%.15g\n"
						    ,a1(i),a2(i),a3(i),a5(i),a7(i));
					a+=buf;
					sprintf(buf,"        a9=%.15g a11=%.15g a13=%.15g a15=%.15g\n"
						    ,a9(i),a11(i),a13(i),a15(i));
					a+=buf;
				}
				break;
			case TOROID:
				a+=s;
				if(IsXToroid(i)) sprintf(buf," xtoroidal ry=%g rx=%g kpx=%g\n",ry(i),rx(i),kpx(i));
				else             sprintf(buf," ytoroidal ry=%g rx=%g kpy=%g\n",ry(i),rx(i),kpy(i));
				a+=buf;
				break;
			case ANAMO:
				a+=s;
				sprintf(buf," biconic ry=%g rx=%g kpy=%g kpx=%g\n",ry(i),rx(i),kpy(i),kpx(i));
				a+=buf;
				break;
			case IDEAL:
				a+=s;
				sprintf(buf," ideal lens f=%g SA0=%g CM0=%g\n",fideal(i),SA0(i),CM0(i));
				a+=buf;
				break;
			case CONIC_OA:
				a+=s;
				sprintf(buf," off-axial conic a=%g b=%g t=%g\n",aCOA(i),bCOA(i),tCOA(i));
				a+=buf;
				break;
		}
		switch(asph_type(i)){
			case CONIC:
			case IDEAL:
				if(cylinder(i)==GENY) a+="      cylinder(generator Y-axis)\n";
				if(cylinder(i)==GENX) a+="      cylinder(generator X-axis)\n";
				break;
		}
		if(is_freeform_surface(i)){
			int m,n,N;
			for(N=0; N<=b_size(i); N++){
				for(m=0; m<=N; m++){
					n=N-m;
					if(b(i,m,n)!=0){
						sprintf(buf,"      b(%d,%d)=%g\n",m,n,b(i,m,n)); a+=buf;
					}
				}
			}
		}
		if(is_zernike_surface(i)){
			sprintf(buf,"      zernike= %s\n", Get_ZCoefficients(i).c_str());
			a+=buf;
		}
		if(is_Dcon_surface(i)){
			sprintf(buf,"      Dcon Rn=%g A4=%g A6=%g A8=%g A10=%g A12=%g A14=%g A16=%g\n",
					DconRn(i),DconA4(i),DconA6(i),DconA8(i),DconA10(i),DconA12(i),DconA14(i),DconA16(i));
			a+=buf;
		}
		if(is_legendre_surface(i)){
			sprintf(buf,"      legendre= %s\n", Get_LCoefficients(i).c_str());
			a+=buf;
		}
		if(is_spline_surface(i)){
			sprintf(buf,"      spline= %s\n", Get_SplineData(i).c_str());
			a+=buf;
		}
		if(is_periodic_surface(i)){
			sprintf(buf,"      Apv=%g sfx=%g sfy=%g\n", Apv(i),sfx(i),sfy(i));
			a+=buf;
		}
		if(is_userdef_surface(i)){
			sprintf(buf,"      UserDefSurf\n");
			a+=buf;
		}
		if(grating(i)){
			sprintf(buf,"no.%2d grating",i);
			a+=buf;
			sprintf(buf," order=%d pitch=%g grx=%g gry=%g grz=%g\n",
				    difforder(i),gpitch(i),grx(i),gry(i),grz(i));
			a+=buf;
		}
		if(Diffusion(i)){
			sprintf(buf,"no.%2d diffusion\n",i);
			a+=buf;
		}
		if(Fish(i)){
			sprintf(buf,"no.%2d fish surface\n",i);
			a+=buf;
		}
	}

	for(i=0; i<=k+1; i++){
		if(i==0){
			s="obj  ";
		}
		else if(i==k+1){
			s="img  ";
		}
		else{
			sprintf(buf,"no.%2d",i); s=buf;
		}
		if(decenter_type(i)!=0){
			a+=s;
			switch(decenter_type(i)){
				case 1:
					sprintf(buf," decenter type%d",decenter_type(i));
					a+=buf;
					sprintf(buf," dx=%g dy=%g dz=%g rox=%g roy=%g roz=%g order=%d  ret=%d\n",
						dx(i),dy(i),dz(i),rox(i),roy(i),roz(i),order(i),ret(i));
					a+=buf;
					sprintf(buf,"                    ");
					a+=buf;
					sprintf(buf," dx1=%g dy1=%g dz1=%g rox1=%g roy1=%g roz1=%g order1=%d\n",
							dx1(i),dy1(i),dz1(i),rox1(i),roy1(i),roz1(i),order1(i));
					a+=buf;
					break;
				case 2:
				case 3:
					sprintf(buf," decenter type%d",decenter_type(i));
					a+=buf;
					sprintf(buf," dx=%g dy=%g dz=%g rox=%g roy=%g roz=%g order=%d\n",
							dx(i),dy(i),dz(i),rox(i),roy(i),roz(i),order(i));
					a+=buf;
					break;
			}
		}
	}

	for(i=0; i<=k+1; i++){
		if(IsGrin(i)!=0){
			sprintf(buf, "no.%2d Selfoc no=%g ��A=%g phi=%g\n", i,fabs(N(i,1)),rA(i,1),GrinPhi(i));
			a+=buf;
		}
	}

	for(i=1; i<=k; i++){
		if(CoatName(i)!=""){
			sprintf(buf, "no.%2d CoatName=%s reverse=%d\n", i,CoatName(i).c_str(),CoatReverse(i));
			a+=buf;
		}
	}

	return a;
}

std::string cLens1::LensDataZemax(){
	// ZEMAX�`���e�L�X�g�f�[�^���쐬����D
	// EPI�����́g���}�⏕�p�G�N�Z���h�ɑΉ����邽�߂ɍ쐬 (2015.02.25)
	std::string a;
	char buf[1000];
	int i;
	
	sprintf(buf,"\tType\tComment\tCurvature\tThickness\tGlass\tSemi-Diameter\n");
	a+=buf;
	for(i=1; i<=k; i++){
		sprintf(buf,"%2d\tSTANDARD\t%s\t%.15g\t%.8g\t%s\t%.8g\n",
			    i,rem(i).c_str(),c(i),d(i),MaterialName(gname(i)).c_str(),EAy(i)/2);
		        // ��) .c_str()��t���Ȃ��ƃG���[��������
		a+=buf;
	}
	return a;
}

std::string cLens1::ParaxialValues(){
	char buf[1000];
	std::string a;

	if(fabs(s)>=LN) { sprintf(buf,"s=�� "); a+=buf; }
	else            { sprintf(buf,"s=%g ", s); a+=buf; }

	if(fabs(t)>=LN) { sprintf(buf,"t=�� "); a+=buf; }
	else            { sprintf(buf,"t=%g ", t); a+=buf; }

	if(Afocal)        { sprintf(buf,"s1fix=�� "); a+=buf; }
	else if(s1fix!=0) { sprintf(buf,"s1fix=%g ", s1fix); a+=buf; }

	sprintf(buf,"��d=%g \n", TotalThickness());
	a+=buf;
	sprintf(buf,"f/n=%.5f bf=%.5f ff=%.5f  ��=%.5f ��'=%.5f  HH'=%.5f\n",
	        f(),bf(),ff(),delta(),delta1(),deadspace());
	a+=buf;
	sprintf(buf,"s'=%.5f ��=%g\n", s1(),M());
	a+=buf;

	return a;
}

std::string cLens1::Coefficients() {
	int i;
	char buf[1000];
	std::string a;

	CalcCoefficients();
	sprintf(buf,"normalize=%d", NormalizeType); a+=buf;
	if(NormalizeUnit.substr(0,2)=="f=") {
		double x=atof( NormalizeUnit.substr(2,NormalizeUnit.length()-2).c_str() );
		sprintf(buf,"  unit=f/%g\n", x==0 ? 1:x); a+=buf;
	}
	else {
		double x=atof(NormalizeUnit.c_str());
		sprintf(buf,"  unit=%g\n", x==0 ? 1:x); a+=buf;
	}
	sprintf(buf,"         SA      SA5     CM      AS      DS      PT      L       T\n");
	a+=buf;
	for(i=1; i<=RowsNumber(); i++){
		sprintf(buf,"%5s%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
			strRow(i).c_str(),
			Coefficient("SA" ,1,0,0,i1Row(i),i2Row(i)) / pow(10,SAe) ,
			Coefficient("SA5",1,0,0,i1Row(i),i2Row(i)) / pow(10,SA5e),
			Coefficient("CM" ,1,0,0,i1Row(i),i2Row(i)) / pow(10,CMe) ,
			Coefficient("AS" ,1,0,0,i1Row(i),i2Row(i)) / pow(10,ASe) ,
			Coefficient("DS" ,1,0,0,i1Row(i),i2Row(i)) / pow(10,DSe) ,
			Coefficient("PT" ,1,0,0,i1Row(i),i2Row(i)) / pow(10,PTe) ,
			Coefficient("LC" ,1,0,0,i1Row(i),i2Row(i)) / pow(10,LCe) ,
			Coefficient("TC" ,1,0,0,i1Row(i),i2Row(i)) / pow(10,TCe)  );
		a+=buf;
	}
	sprintf(buf,"   ��%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n",
		Coefficient("SA" ) / pow(10,SAe) ,
		Coefficient("SA5") / pow(10,SA5e),
		Coefficient("CM" ) / pow(10,CMe) ,
		Coefficient("AS" ) / pow(10,ASe) ,
		Coefficient("DS" ) / pow(10,DSe) ,
		Coefficient("PT" ) / pow(10,PTe) ,
		Coefficient("LC" ) / pow(10,LCe) ,
		Coefficient("TC" ) / pow(10,TCe)  );
	a+=buf;
	sprintf(buf,"         e%3g    e%3g    e%3g    e%3g    e%3g    e%3g    e%3g    e%3g\n",
		SAe,SA5e,CMe,ASe,DSe,PTe,LCe,TCe);
	a+=buf;
	return a;
}

std::string cLens1::alpha_h_table() {
	int i;
	char buf[1000];
	std::string s;
	CalcCoefficients();
	sprintf(buf,"normalize=%d", NormalizeType); s+=buf;
	if(NormalizeUnit.substr(0,2)=="f="){
		double x=atof( NormalizeUnit.substr(2,NormalizeUnit.length()-2).c_str() );
		sprintf(buf,"  unit=f/%g\n", x==0 ? 1:x); s+=buf;
	}
	else{
		double x=atof(NormalizeUnit.c_str());
		sprintf(buf,"  unit=%g\n", x==0 ? 1:x); s+=buf;
	}
    //           ### #####.###### #####.###### #####.###### #####.######  #####.#### #####.#### #####.####   @@@..
	sprintf(buf,"        ��'           h           ��p'          hp           h*h        h*hp      hp*hp\n");
	s+=buf;
	sprintf(buf,"    %12.6f              %12.6f \n", u[1]*N(0,1),up[1]*N(0,1)); s+=buf;
	for(i=1; i<=k; i++){
		sprintf(buf,"%3d %12.6f %12.6f %12.6f %12.6f  %10.4f %10.4f %10.4f   %s\n",
		        i,u[i+1]*N(i,1),h[i],up[i+1]*N(i,1),hp[i],h[i]*h[i],h[i]*hp[i],hp[i]*hp[i],rem(i).c_str());
		s+=buf;
	}
	return s;
}

std::string cLens1::Focallengths(int i1,int i2,int ForDrawing){
	char buf[10000];
	std::string s;

	if(ForDrawing){
		cLens1 x=Trimed(i1,i2);
		x.Set_color(1,"e");
		// �O��̔}������C�Ƃ���D(�ڍ������Y�̕��i�̐��}�ɕK�v�j
		x.Set_gname(0,"1"); x.Set_gname(x.k,"1");
		s+=x.LensData(1);
		s+="    e-line\n";
		sprintf(buf,"f  =%.3f\n",  x.f() ); s+=buf;
		sprintf(buf,"bf'=%.3f\n",  x.bf()); s+=buf;
		sprintf(buf,"bf =%.3f\n", -x.ff()); s+=buf;

	}
	else{
		if(i1==1 && i2==k){
			sprintf(buf,"f (%s)=%.8g",          color[1].c_str(),f(i1,i2)        ); s+=buf;
			sprintf(buf," (%gdiopter)\n",       1000/f(i1,i2)                    ); s+=buf; 
			sprintf(buf,"bf(%s)=%.8g\n",        color[1].c_str(),bf(i1,i2)       ); s+=buf;
			sprintf(buf,"ff(%s)=%.8g\n",        color[1].c_str(),ff(i1,i2)       ); s+=buf;
			sprintf(buf,"��'(%s)=%.8g\n",       color[1].c_str(),delta1(i1,i2)   ); s+=buf;
			sprintf(buf,"�� (%s)=%.8g\n",       color[1].c_str(),delta(i1,i2)    ); s+=buf;
			sprintf(buf,"deadspace(%s)=%.8g\n", color[1].c_str(),deadspace(i1,i2)); s+=buf;
		}
		else{
			sprintf(buf,"f (%s,%d-%d)=%.8g",          color[1].c_str(),i1,i2,f(i1,i2)        ); s+=buf;
			sprintf(buf," (%gdiopter)\n",             1000/f(i1,i2)                          ); s+=buf;
			sprintf(buf,"bf(%s,%d-%d)=%.8g\n",        color[1].c_str(),i1,i2,bf(i1,i2)       ); s+=buf;
			sprintf(buf,"ff(%s,%d-%d)=%.8g\n",        color[1].c_str(),i1,i2,ff(i1,i2)       ); s+=buf;
			sprintf(buf,"��'(%s,%d-%d)=%.8g\n",       color[1].c_str(),i1,i2,delta1(i1,i2)   ); s+=buf;
			sprintf(buf,"�� (%s,%d-%d)=%.8g\n",       color[1].c_str(),i1,i2,delta(i1,i2)    ); s+=buf;
			sprintf(buf,"deadspace(%s,%d-%d)=%.8g\n", color[1].c_str(),i1,i2,deadspace(i1,i2)); s+=buf;
		}
	}
	return s;
}

std::ostream& operator<<(std::ostream& to,cLens1& x){
	int file_ver=138;
	int i,j;
	to<<file_ver<<std::endl;
	to<<x.k<<' '<<x.cn<<std::endl;
	to<<x.var<<' '<<x.var1<<' '<<x.var2<<' '<<x.var3<<' '<<x.var4<<' '<<x.var5<<std::endl;
	to<<blank_filled(x.Note)<<std::endl;
	to<<x.FRW<<std::endl;
	to<<x.s<<' '<<x.t<<' '<<x.s1fix<<' '<<x.Afocal<<' '<<x.AfocalMode<<' '<<x.AfocalRotateEye<<std::endl;
	to<<x.StopDominate<<std::endl;
	to<<x.yObjectMax<<' '<<x.xObjectMax<<' ';
	to<<x.stop<<' '<<x.EPD<<' '<<x.EPDx<<' '<<x.EPy<<' '<<x.EPx<<std::endl;
	to<<x.SourcePhiY<<' '<<x.SourcePhiX<<' '<<x.SourceAxis<<' '<<x.SourceAreaType<<std::endl;
	to<<x.ExcludeVirtualRay<<' '<<x.ExcludeVirtualObject<<std::endl;
	for(j=1; j<=x.cn; j++) {
		to<<x.color[j]<<' '<<x.colorweight[j]<<std::endl;
	}
	for(i=0; i<=x.k+1; i++){
		to<<x.surf[i]<<std::endl;
		to<<x.med[i]<<std::endl;
	}
	return to;
}

std::istream& operator>>(std::istream& from,cLens1& x){
	int i,j, file_ver;
	
	from >> file_ver; 
	surface::file_ver=file_ver;
	medium::file_ver=file_ver;

	// file_ver��ύX����Ƃ��́Csurface,medium�̃t�����h�֐� operator>> �ւ�case����ǉ����邱��
	//                           �{�֐������� "�ǂݍ��݌�̏���" �ւ�case����ǉ����邱��
	//                           operator<< �� file_ver ���ύX���邱�ƁD

	switch(file_ver){
		case 138:
			// medium��gVariable��ǉ� 200622
		case 137:
			// surface��EAdy,EAdx��ǉ� 200508
		case 136:
			// surface��a9,a11,a13,a15��ǉ� 200423
		case 135:
			// surface��spline��ǉ� 200123
		case 134:
			// var,var1,..var5��ǉ� 190729
			// (�������Ȃ��ƁC�����݌v�������l(=0)����J�n���邱�ƂɂȂ�C�]���֐����傫�����Ď������Ȃ����Ƃ�����D�j
			from>>i>>j; x=cLens1(i,j);
			from>>x.var>>x.var1>>x.var2>>x.var3>>x.var4>>x.var5;
			from>>x.Note; x.Note=inv_blank_filled(x.Note);
			from>>x.FRW;
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode>>x.AfocalRotateEye;
			from>>x.StopDominate;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPDx>>x.EPy>>x.EPx;
			from>>x.SourcePhiY>>x.SourcePhiX>>x.SourceAxis>>x.SourceAreaType;
			from>>x.ExcludeVirtualRay>>x.ExcludeVirtualObject;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 133:
			// EPphi��EPD�ɉ��́CEPDx��ǉ� 180919
			from>>i>>j; x=cLens1(i,j);
			from>>x.Note; x.Note=inv_blank_filled(x.Note);
			from>>x.FRW;
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode>>x.AfocalRotateEye;
			from>>x.StopDominate;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPDx>>x.EPy>>x.EPx;
			from>>x.SourcePhiY>>x.SourcePhiX>>x.SourceAxis>>x.SourceAreaType;
			from>>x.ExcludeVirtualRay>>x.ExcludeVirtualObject;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 132:
			// aCOA,bCOA�����ւ� 180329
		case 131:
			// surface��aCOA,bCOA,tCOA��ǉ� 180329
		case 130:
			// surface��a1,a2,a3,a5,a7��ǉ� 180123
		case 129:
			// surface��CHMy,CHMx��ǉ� 171006
		case 128:
			// surface��legendre��ǉ� 170526
		case 127:
			// surface��NormH��ǉ� 160610
		case 126:
			// surface��a18,a20��ǉ� 160317
		case 125:
			// surface��Fish��ǉ� 160223
		case 124:
			// AfocalRotateEye��ǉ� 150925
			from>>i>>j; x=cLens1(i,j);
			from>>x.Note; x.Note=inv_blank_filled(x.Note);
			from>>x.FRW;
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode>>x.AfocalRotateEye;
			from>>x.StopDominate;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			from>>x.SourcePhiY>>x.SourcePhiX>>x.SourceAxis>>x.SourceAreaType;
			from>>x.ExcludeVirtualRay>>x.ExcludeVirtualObject;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 123:
			// surface��a16��ǉ� 150421
		case 122:
			// surface��a12,a14��ǉ� 140619
		case 121:
			// StopDominate��ǉ� 130820
			from>>i>>j; x=cLens1(i,j);
			from>>x.Note; x.Note=inv_blank_filled(x.Note);
			from>>x.FRW;
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode>>x.StopDominate;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			from>>x.SourcePhiY>>x.SourcePhiX>>x.SourceAxis>>x.SourceAreaType;
			from>>x.ExcludeVirtualRay>>x.ExcludeVirtualObject;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 120:
			// surface��Dcon��ǉ� 130317
		case 119:
			// sufface��CoatName,CoatReverse��ǉ� 130213
		case 118:
			// sufface��SA0,CM0��ǉ� 120912
		case 117:
			// surface�Ŏ��R�Ȗʂ̕\����ύX 120906
		case 116:
			// surface��IsUserDefSurf��ǉ� 120605
		case 115:
			// surface��Apv,sfx,sfy��ǉ� 120410
		case 114:
			// Note��ǉ� 120406
			from>>i>>j; x=cLens1(i,j);
			from>>x.Note; x.Note=inv_blank_filled(x.Note);
			from>>x.FRW;
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			from>>x.SourcePhiY>>x.SourcePhiX>>x.SourceAxis>>x.SourceAreaType;
			from>>x.ExcludeVirtualRay>>x.ExcludeVirtualObject;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 113:
			// surface��IsXToroid��ǉ� 120316
		case 112:
			// surface��kpy,kpx��ǉ� 111220	
		case 111:
			// surface��zernike��ǉ� 110601
		case 110:
			// SourcePhiY����ǉ� 100805
			from>>x.FRW;  // �����ł̓R���X�g���N�^�̑O��FRW��ݒ肵�Ă��邪�ԈႢ�D
			              // �ݒ肵�Ă��R���X�g���N�^�ŏ���������Ă��܂��Dcase 114�ł͏C�����Ă���D
			from>>i>>j; x=cLens1(i,j);
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			from>>x.SourcePhiY>>x.SourcePhiX>>x.SourceAxis>>x.SourceAreaType;
			from>>x.ExcludeVirtualRay>>x.ExcludeVirtualObject;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 109:
			// medium�̃����ogname�ŋ󔒂����e���� 100412
		case 108:
			// surface�̃����o��order,order1��ǉ� 100202
		case 107:
			// ExcludeVirtualRay,ExcludeVirtualObject�ǉ� 091019
			from>>x.FRW;
			from>>i>>j; x=cLens1(i,j);
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			from>>x.ExcludeVirtualRay>>x.ExcludeVirtualObject;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 106:
			// FRW�ǉ� 090108
			from>>x.FRW;
			from>>i>>j; x=cLens1(i,j);
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 105:
			// Afocal�ǉ� 080912
			from>>i>>j; x=cLens1(i,j);
			from>>x.s>>x.t>>x.s1fix>>x.Afocal>>x.AfocalMode;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		case 104:
			// surface�̃����o��b20,...b04��ǉ� 080328
		case 103:
			// surface�̃����o��rem��ǉ� 080312
		case 102:
			// surface�̃����o��Diffusion��ǉ� 070723
		case 101:
			// surface�̃����o��dx1,dy1,dz1,rox1,roy1,roz1��ǉ� 070416
		case 100: 
			// AfocalMode�ǉ� 060807
			from>>i>>j; x=cLens1(i,j);
			from>>x.s>>x.t>>x.s1fix>>x.AfocalMode;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
		default:
			// version�̂Ȃ����`���ɑΉ�����D
			// �����ɂ͖ʐ�=100�`���݂̃o�[�W������ �̂Ƃ��듮��ƂȂ邪�������̂悤�ȃf�[�^�͂Ȃ��D
			i=file_ver;
			from>>j;
			x=cLens1(i,j);
			from>>x.s>>x.t>>x.s1fix;
			from>>x.yObjectMax>>x.xObjectMax;
			from>>x.stop>>x.EPD>>x.EPy>>x.EPx;
			for(j=1; j<=x.cn; j++) {
				from>>x.color[j]>>x.colorweight[j];
				x.Set_color(j,x.color[j]);
			}
			for(i=0; i<=x.k+1; i++){
				from>>x.surf[i];
				from>>x.med[i];
			}
			break;
	}

	// �ǂݍ��݌�̏���
	if(file_ver<105){
		// ver.105�ɔ������� 080912
		if( fabs(x.s1fix)>=LN ){ x.Afocal=1; x.s1fix=0; }
		x.AfocalMode+=1; // enum��0����n�܂��Ă����̂�1����n�߂�D���C���I���R�D
	}

	return from;
}

int cLens1::open(std::string filename){
	std::ifstream from(filename.c_str());
	if(from) {
		from >> *this;
		this->filename=filename;
		return 1;
	}
	else return 0;
}

int cLens1::save(std::string filename){
	std::ofstream to(filename.c_str());
	if(to) {
		to.precision(10);  // �f�t�H���g��6��
		to << *this;
		this->filename=filename;
		return 1;
	}
	else return 0;
}

int cLens1::SaveAsRTPFile(std::string filename){
	int i;
	const char C=',', Q='"';
	std::ofstream to(filename.c_str());
	if(to){
		to <<k <<C; 
		if(cn>3) to <<3 <<C; else to <<cn <<C;
		to <<s <<C <<t <<C <<s1fix <<std::endl;
		switch(cn){
		case  1: to <<Q <<color[1] <<Q <<C <<Q <<"C'"     <<Q <<C <<Q <<"F'"     <<Q <<std::endl; break;
		case  2: to <<Q <<color[1] <<Q <<C <<Q <<color[2] <<Q <<C <<Q <<"F'"     <<Q <<std::endl; break;
		default: to <<Q <<color[1] <<Q <<C <<Q <<color[2] <<Q <<C <<Q <<color[3] <<Q <<std::endl;
		}
		for(i=0; i<=k+1; ++i){
				if(i==0){
				to <<r(i) <<C <<d(i) <<C <<Q <<"" <<Q <<C 
				   <<EAtype(i) <<C <<EAy(i) <<C <<EAx(i)<<std::endl;
			}
			else{
				to <<r(i) <<C <<d(i) <<C <<Q <<gname(i-1) <<Q <<C 
				   <<EAtype(i) <<C <<EAy(i) <<C <<EAx(i)<<std::endl;
			}
			to <<asph_type(i) <<C <<kp(i) <<C <<0 <<C 
			   <<a4(i) <<C <<a6(i) <<C <<a8(i) <<C <<a10(i) <<std::endl;
			to <<cylinder(i) <<C <<fresnel(i) <<C <<rbase(i) <<C 
			   <<ry(i) <<C <<rx(i) <<C <<coneangle(i) <<C <<fideal(i)<<std::endl;
			to <<decenter_type(i) <<C <<dy(i) <<C <<dx(i) <<C <<dz(i) <<C
				<<roy(i) <<C <<rox(i) <<C <<roz(i) <<std::endl;
		}
		return 1;
	}
	else return 0;

}
	
int cLens1::Get_NormalizeType(){ return NormalizeType; }
void cLens1::Set_NormalizeType(int value){ 
	if(value==1 || value==2 || value==3) NormalizeType=value; 
}

std::string cLens1::Get_NormalizeUnit(){ return NormalizeUnit; }
void cLens1::Set_NormalizeUnit(std::string value){ NormalizeUnit=value; }

void cLens1::CalcCoefficients(int j2/*=0*/){
	// j2!=0 �̂Ƃ��́C�F�����W����j2�g���Ƒ�1�g���̊ԂŌv�Z����D
	// j2==0 �̂Ƃ��́C�Œ��g���[�ƍŒZ�g���[�̊ԂŌv�Z����D

	// static cLens1 x, *px;  
	// if (x==*this && px==this) return; else {x=*this; px=this;} ////////////////////////////////////
	// static cLens1 x;  
	// if (x==*this) return; else x=*this;  (�s�\���j
	//
	// static �L�[���[�h��t���Đ錾�������\�b�h�ƕϐ��ɂ̓R�s�[�� 1 �����Ȃ��A
	// �N���X�̂��ׂẴI�u�W�F�N�g��������Q�Ƃ��邱�ƂɂȂ�܂��B
	// �N���X�̃C���X�^���X�������������Ă��A
	// static ���\�b�h�� static �ϐ��̃R�s�[�� 1 �����쐬����܂���B
	// (ex.)
	// cLens1 x,y;
	// x.open("file");
	// y=x;
	// x.CalcCoefficients();
	// y.CalcCoefficients(); .... �v�Z����Ȃ�

	int i,j, ii;
	double dN,dN1;
	double g,g_hat;
	double c;
	double unit;
	double al,al1,alp,alp1;
	double A4;

	// ������
	for(i=0; i<=k+1; i++){
		SA[i]=CM[i]=AS[i]=DS[i]=PT[i]=LC[i]=TC[i]=LC2[i]=TC2[i]=0; 
		SAP[i]=CMP[i]=ASP[i]=DSP[i]=LCP[i]=0;
		SA5[i]=0; 
		CM41Z[i]=CM41[i]=CM41ALL[i]=0; 
		SA32F[i]=SA32Z[i]=SA32[i]=SA32ALL[i]=0;
		CM23[i]=CM23P[i]=CM23Z[i]=CM23ALL[i]=0;
		AS5[i]=SG5[i]=0;
		DS5[i]=0;
		PRE[i]=DSE1[i]=DSE2[i]=ASE[i]=PTE[i]=CME[i]=0;
	}

	// �G���[����i�s��l������邽�߁C�������Ɏ����W�������������Ă����j
	for(i=0; i<=k; i++) if(N(i,1)==0) return;  // N�Ŋ��邱�Ƃ������̂ŁC�G���[�������D
	if( power()==0 ) return;                   // power()=0�͂��Ȃ킿A[3]=0�ł���Cf()���ŃG���[�����D
	if( fabs(s)>=LN && fabs(t)>=LN ) return;
	if( NormalizeType!=1 && NormalizeType!=2 && NormalizeType!=3 ) return;

	double *cv=new double [k+2];
	double *hD=new double [k+2]; double *hDp=new double [k+2];
	double *P=new double [k+2];
	double *O=new double [k+2];
	
	SAt=CMt=ASt=DSt=PTt=LCt=TCt=LC2t=TC2t=0; 
	SAPt=CMPt=ASPt=DSPt=LCPt=0;
	SA5t=0; 
	CM41Zt=CM41t=CM41ALLt=0; 
	SA32Ft=SA32Zt=SA32t=SA32ALLt=0;
	CM23t=CM23Pt=CM23Zt=CM23ALLt=0;
	AS5t=SG5t=0;
	DS5t=0;
	PREt=DSE1t=DSE2t=ASEt=PTEt=CMEt=0;

	SAe=CMe=ASe=DSe=PTe=LCe=TCe=LC2e=TC2e=0; 
	SAPe=CMPe=ASPe=DSPe=LCPe=0;
	SA5e=0; 
	CM41Ze=CM41e=CM41ALLe=0; 
	SA32Fe=SA32Ze=SA32e=SA32ALLe=0;
	CM23e=CM23Pe=CM23Ze=CM23ALLe=0;
	AS5e=SG5e=0;
	DS5e=0;
	PREe=DSE1e=DSE2e=ASEe=PTEe=CMEe=0;

	g=s-t;
	g_hat=this->g_hat(1);
	
	// �����̒P�ʂ̐ݒ�
	if(NormalizeUnit.substr(0,2)=="f=") {
		double x=atof( NormalizeUnit.substr(2,NormalizeUnit.length()-2).c_str() );
		unit= x==0 ? f() : f()/x;
	}
	else {
		unit=atof(NormalizeUnit.c_str());
	}
	if( unit==0 ) c=1; else c=1/unit;

	for(i=0; i<=k+1; i++){
		if( r(i)==0 ) cv[i]=0; else cv[i]=1/(r(i)*c);
	}

	// �e�ʂ̃p���[
	for(i=1; i<=k; i++){
		phi[i]=(N(i,1)-N(i-1,1))*cv[i];
		if(is_ideallens(i)) phi[i]= fideal(i)==0 ? 0 : 1/(fideal(i)*c);
	}
		
	//  marginal ray  //
	switch(NormalizeType){
	case 1:
		if( fabs(s)>=LN ){
			u[1]=0; h[1]=1;
		}
		else{
			u[1]=1/(g_hat*c); h[1]=s/g_hat; 
		}
		break;
	case 2:
		if( fabs(s)>=LN ){
			u[1]=0; h[1]=f()*c; 
		}
		else{
			u[1]=M()/N(0,1); h[1]=(s*c)*M()/N(0,1); 
		}
		break;
	case 3:
		if( fabs(t)>=LN ){
			u[1]=-1/(f()*c)/N(0,1); h[1]=-(s/f())/N(0,1);
		}
		else{
			u[1]=1/Mpupil()/(g*c); h[1]=(s/g)/Mpupil();
		}
		break;
	}
	for(i=1; i<=k; i++){
		u[i+1]=( phi[i]*h[i]+N(i-1,1)*u[i] )/N(i,1);
		if(i!=k) h[i+1]=h[i]-(d(i)*c)*u[i+1];
	}
	for(i=1; i<=k; i++){
		hQ[i]=h[i]*N(i-1,1)*cv[i]-N(i-1,1)*u[i];
		hD[i]=u[i+1]/N(i,1)-u[i]/N(i-1,1);
	}

	//  principal ray  //
	switch(NormalizeType){
	case 1:
		if( fabs(s)>=LN && fabs(t)<LN ){
			up[1]=-1/N(0,1); hp[1]=-(t*c)/N(0,1);
		}
		if( fabs(s)<LN && fabs(t)>=LN ){
			up[1]=0; hp[1]=(g_hat*c)/N(0,1); 
		}
		if( fabs(s)<LN && fabs(t)<LN ){
			up[1]=-g_hat/g/N(0,1); hp[1]=up[1]*(t*c);
		}
		break;
	case 2:
		if( fabs(s)>=LN && fabs(t)<LN ){
			up[1]=-1/N(0,1)/(f()*c); hp[1]=-t/N(0,1)/f(); 
		}
		if( fabs(s)<LN && fabs(t)>=LN ){
			up[1]=0; hp[1]=1/M(); 
		}
		if( fabs(s)<LN && fabs(t)<LN ){
			up[1]=-1/(g*c)/M(); hp[1]=-t/g/M();
		}
		break;
	case 3:
		if(fabs(t)>=LN){
			up[1]=0; hp[1]=-f()*c;
		}
		else{
			up[1]=-Mpupil()/N(0,1); hp[1]=-Mpupil()*(t*c/N(0,1));
		}
		break;
	}
	for(i=1; i<=k; i++){
		up[i+1]=( phi[i]*hp[i]+N(i-1,1)*up[i] )/N(i,1);
		if(i!=k) hp[i+1]=hp[i]-(d(i)*c)*up[i+1];
	}
	for(i=1; i<=k; i++){
		hQp[i]=hp[i]*N(i-1,1)*cv[i]-N(i-1,1)*up[i];
		hDp[i]=up[i+1]/N(i,1)-up[i]/N(i-1,1);
	}

	// "�����Y�݌v�@" (4.5)����psi,omega
	for(i=1; i<=k; i++){		
		if(asph_type(i)==CONIC){
			A4=a4(i)/NormH(i)/NormH(i)/NormH(i)/NormH(i);
			P[i]=(N(i,1)-N(i-1,1))*( kp(i)*cv[i]*cv[i]*cv[i] + 8*(A4/c/c/c) );
			O[i]=3*(N(i,1)-N(i-1,1))
				 *( (2*kp(i)+kp(i)*kp(i))*cv[i]*cv[i]*cv[i]*cv[i]*cv[i] +16*(a6(i)/c/c/c/c/c) );
		}
		else{
			P[i]=O[i]=0;
		}
	}

	//  3rd order and chromatic coefficients  //
	for(i=1; i<=k; i++){
		if(asph_type(i)==IDEAL){
			// ���z�����Y�̎����W���͓����W������v�Z����D
			matrix<double> X(3,1),A,Y;

			X[1][1]=SA0(i);
			X[2][1]=CM0(i);
			X[3][1]=1;
			A=CMatrix(i);
			Y=A*X;           // �g�����_�h(3.5.5a)��
			SA[i]=Y[1][1];
			CM[i]=Y[2][1];
			AS[i]=Y[3][1];
			PT[i]=Y[4][1];
			DS[i]=Y[5][1];
			SAP[i]=Y[6][1];
			CMP[i]=ASP[i]=DSP[i]=LC[i]=TC[i]=LC2[i]=TC2[i]=0;
		}
		else{
			SA[i]=hQ[i]*hQ[i]*hD[i]*h[i]                  +h[i]*h[i]*h[i]*h[i]*P[i];
			CM[i]=hQ[i]*hQp[i]*hD[i]*h[i]                 +h[i]*h[i]*h[i]*hp[i]*P[i];
			AS[i]=hQp[i]*hQp[i]*hD[i]*h[i]                +h[i]*h[i]*hp[i]*hp[i]*P[i];
			DS[i]=hQp[i]*hQp[i]*hD[i]*hp[i]+hQp[i]*hDp[i] +h[i]*hp[i]*hp[i]*hp[i]*P[i];
			PT[i]=-( 1/N(i,1)-1/N(i-1,1) )*cv[i];
			SAP[i]=hp[i]*hQp[i]*hQp[i]*hDp[i]             +hp[i]*hp[i]*hp[i]*hp[i]*P[i];
			CMP[i]=hp[i]*hQ[i]*hQp[i]*hDp[i]              +h[i]*hp[i]*hp[i]*hp[i]*P[i];
			ASP[i]=hp[i]*hQ[i]*hQ[i]*hDp[i]               +h[i]*h[i]*hp[i]*hp[i]*P[i];
			DSP[i]=hQ[i]*hQ[i]*hDp[i]*h[i]-hQ[i]*hD[i]    +h[i]*h[i]*h[i]*hp[i]*P[i];
			
			if(1<=j2 && j2<=cn){
				dN=N(i-1,j2)-N(i-1,1); dN1=N(i,j2)-N(i,1);
			}
			else{
				dN=N(i-1,short_cn)-N(i-1,long_cn); dN1=N(i,short_cn)-N(i,long_cn);
			}
			LC[i]=h[i]*hQ[i]*( dN1/N(i,1)-dN/N(i-1,1) );
			TC[i]=h[i]*hQp[i]*( dN1/N(i,1)-dN/N(i-1,1) );
			LCP[i]=hp[i]*hQp[i]*( dN1/N(i,1)-dN/N(i-1,1) );
			
			dN=N(i-1,mid_cn)-N(i-1,long_cn); dN1=N(i,mid_cn)-N(i,long_cn);
			LC2[i]=h[i]*hQ[i]*( dN1/N(i,1)-dN/N(i-1,1) );
			TC2[i]=h[i]*hQp[i]*( dN1/N(i,1)-dN/N(i-1,1) );
		}
	}

	if(rObj()!=0){                // ���̖ʌW��
		PT[0]=-cv[0]/N(0,1);
		DS[0]= u[1]==0 ? 0 : -(up[1]/u[1])*cv[0]/N(0,1);
	}

	if(rImage()!=0){              // ���ʌW��
		PT[k+1]=cv[k+1]/N(k,1);
		DS[k+1]= u[k+1]==0 ? 0 : (up[k+1]/u[k+1])*cv[k+1]/N(k,1);
	}

	
	//  5th order coefficients  //
	for(i=1; i<=k; i++){

		// "�����Y�݌v�@" (4.8a)
		SA5[i]=3*SA[i]*( h[i]*cv[i]*( h[i]*cv[i]-2*u[i] ) - hQ[i]*hD[i] );
		CM41Z[i]=h[i]*hQ[i]*cv[i]*hD[i]*( u[i+1]-3*u[i] );
		CM41[i]=3*CM[i]*( h[i]*cv[i]*(h[i]*cv[i]-2*u[i]) -hQ[i]*hD[i] )-CM41Z[i];
		SA32Z[i]=h[i]*hQp[i]*cv[i]*hD[i]*( u[i+1]-3*u[i] ) +h[i]*cv[i]*cv[i]*hD[i];
		SA32F[i]=3*AS[i]*( h[i]*cv[i]*(h[i]*cv[i]-2*u[i])-hQ[i]*hD[i] )
			-h[i]*hQp[i]*cv[i]*hD[i]*(u[i+1]-3*u[i]) - SA32Z[i];
		SA32[i]=SA32F[i] +hQ[i]*cv[i]*( 1/N(i,1)*hD[i]+h[i]/N(i-1,1)*PT[i] );
		DS5[i]=3*h[i]*(hQp[i]*hQp[i]*hDp[i]+hp[i]*hp[i]*hp[i]*P[i])
			*(hp[i]*cv[i]*(hp[i]*cv[i]-2*up[i])-hQp[i]*hDp[i])
			-2*hp[i]*hQp[i]*hQp[i]*cv[i]*(hDp[i]/N(i,1)+hp[i]*PT[i]/N(i-1,1))
			-hp[i]*hp[i]*hQp[i]*PT[i]*( 2*hQp[i]*PT[i]-hp[i]*cv[i]*cv[i] );
		SG5[i]=3*h[i]*CMP[i]*cv[i]*(hp[i]*cv[i]-2*up[i]) -3*h[i]*hQ[i]*hQp[i]*hQp[i]*hDp[i]*hDp[i]
			-2*hQ[i]*hQp[i]*hp[i]*cv[i]*( hDp[i]/N(i,1)+PT[i]*hp[i]/N(i-1,1) )
			-hQ[i]*hp[i]*hp[i]*PT[i]*(2*hQp[i]*PT[i]-hp[i]*cv[i]*cv[i]);
		AS5[i]=SG5[i] +h[i]*hQp[i]*cv[i]*hDp[i]*(up[i+1]-3*up[i])
			-hQp[i]*cv[i]*(hDp[i]/N(i,1)+hp[i]*PT[i]/N(i-1,1))
			-hp[i]*PT[i]*(2*hQp[i]*PT[i]-hp[i]*cv[i]*cv[i]);
		CM23[i]=3*h[i]*ASP[i]*cv[i]*(hp[i]*cv[i]-2*up[i])
			-3*h[i]*hQ[i]*hQ[i]*hQp[i]*hDp[i]*hDp[i]
			+h[i]*hQ[i]*cv[i]*hDp[i]*(up[i+1]-3*up[i])
			-hQ[i]*cv[i]*(hp[i]*hQ[i]+h[i]*hQp[i])*(hDp[i]/N(i,1)+hp[i]*PT[i]/N(i-1,1))
			-h[i]*hp[i]*hQ[i]*PT[i]*( 2*hQp[i]*PT[i]-hp[i]*cv[i]*cv[i] );
		CM23P[i]=CM23[i]
			+h[i]*hQ[i]*cv[i]*hDp[i]*(up[i+1]-3*up[i])
			-hQ[i]*cv[i]*( hDp[i]/N(i,1)+hp[i]*PT[i]/N(i-1,1) )
			-h[i]*cv[i]*PT[i]*(up[i+1]-3*up[i]);
		CM23Z[i]=h[i]*hQ[i]*cv[i]*hDp[i]*(up[i+1]-3*up[i])
			-h[i]*cv[i]*PT[i]*(up[i+1]-3*up[i]) +PT[i]*PT[i];
		
		// (4.8b)
		if(P[i]!=0){
			SA5[i]+=((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQ[i]*hQ[i] +6*h[i]*hQ[i]*PT[i] 
				-9*h[i]*h[i]*cv[i]*cv[i])*h[i]*h[i]*h[i]*h[i]*P[i];
			CM41Z[i]+=((2/N(i,1)/N(i,1))*hQ[i] +h[i]*PT[i]
				-(2/N(i,1))*h[i]*cv[i])*h[i]*h[i]*h[i]*P[i];
			CM41[i]+=(((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQ[i]*hQ[i]
				+6*h[i]*hQ[i]*PT[i] -9*h[i]*h[i]*cv[i]*cv[i])*h[i]*h[i]*h[i]*hp[i]
				+((1/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQ[i]+h[i]*PT[i])*h[i]*h[i]*h[i])*P[i];
			SA32Z[i]+=(((2/N(i,1)/N(i,1))*hQ[i] +h[i]*PT[i] -(2/N(i,1))*h[i]*cv[i])*h[i]*h[i]*hp[i]
				-h[i]*h[i]/N(i,1)/N(i,1))*P[i];
			SA32F[i]+=(((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQ[i]*hQ[i] +6*h[i]*hQ[i]*PT[i]
				-9*h[i]*h[i]*cv[i]*cv[i])*h[i]*h[i]*hp[i]*hp[i]
				+2*((1/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQ[i]+h[i]*PT[i])*h[i]*h[i]*hp[i])*P[i];
			SA32[i]+=(((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQ[i]*hQ[i]+6*h[i]*hQ[i]*PT[i]
				-9*h[i]*h[i]*cv[i]*cv[i])*h[i]*h[i]*hp[i]*hp[i]
				+2*((1/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQ[i]+h[i]*PT[i])*h[i]*h[i]*hp[i]
				+(1/N(i,1))*(1/N(i,1)+1/N(i-1,1))*h[i]*h[i])*P[i];
			DS5[i]+=(((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]*hQp[i]+6*hp[i]*hQp[i]*PT[i]
				-9*hp[i]*hp[i]*cv[i]*cv[i])*h[i]*hp[i]*hp[i]*hp[i]
				+(-(1/N(i,1)/N(i-1,1))*hQp[i]-(2/N(i,1))*hp[i]*cv[i])*hp[i]*hp[i]*hp[i])*P[i];
			SG5[i]+=(((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]*hQp[i]+6*hp[i]*hQp[i]*PT[i]
				-9*hp[i]*hp[i]*cv[i]*cv[i])*h[i]*h[i]*hp[i]*hp[i]
				-2*((1/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]+hp[i]*PT[i])*h[i]*hp[i]
				+(1/N(i,1)/N(i-1,1))*hp[i]*hp[i])*P[i];
			AS5[i]+=(((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]*hQp[i]+6*hp[i]*hQp[i]*PT[i]
				-9*hp[i]*hp[i]*cv[i]*cv[i])*h[i]*h[i]*hp[i]*hp[i]
				-((1/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]+hp[i]*PT[i])*h[i]*hp[i]*hp[i]
				+(-(1/N(i,1)/N(i-1,1))*hQp[i]-(2/N(i,1))*hp[i]*cv[i])*h[i]*hp[i]*hp[i])*P[i];
			CM23[i]+=(((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]*hQp[i]+6*hp[i]*hQp[i]*PT[i]
				-9*hp[i]*hp[i]*cv[i]*cv[i])*h[i]*h[i]*h[i]*hp[i]
				-3*((1/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]+hp[i]*PT[i])*h[i]*h[i]*hp[i]
				+(1/N(i,1))*(1/N(i,1)+1/N(i-1,1))*h[i]*hp[i])*P[i];
			CM23P[i]+=(((3/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]*hQp[i]+6*hp[i]*hQp[i]*PT[i]
				-9*hp[i]*hp[i]*cv[i]*cv[i])*h[i]*h[i]*h[i]*hp[i]
				-2*((1/N(i,1))*(2/N(i,1)+1/N(i-1,1))*hQp[i]+hp[i]*PT[i])*h[i]*h[i]*hp[i]
				+(-(1/N(i,1)/N(i-1,1))*hQp[i]-(2/N(i,1))*hp[i]*cv[i])*h[i]*h[i]*hp[i])*P[i];
			CM23Z[i]+=((2/N(i,1)/N(i,1))*hQp[i]+hp[i]*PT[i]-(2/N(i,1))*hp[i]*cv[i])*h[i]*h[i]*hp[i]*P[i];
		}

		// (4.8c)
		if(O[i]!=0){
			SA5[i]+=h[i]*h[i]*h[i]*h[i]*h[i]*h[i]*O[i];
			SA32[i]+=h[i]*h[i]*h[i]*h[i]*hp[i]*hp[i]*O[i];
			CM41Z[i]+=0;
			CM23Z[i]+=0;
			CM41[i]+=h[i]*h[i]*h[i]*h[i]*h[i]*hp[i]*O[i];
			CM23[i]+=h[i]*h[i]*h[i]*hp[i]*hp[i]*hp[i]*O[i];
			SA32Z[i]+=0;
			AS5[i]+=h[i]*h[i]*hp[i]*hp[i]*hp[i]*hp[i]*O[i];
			SA32F[i]+=h[i]*h[i]*h[i]*h[i]*hp[i]*hp[i]*O[i];
			SG5[i]+=h[i]*h[i]*hp[i]*hp[i]*hp[i]*hp[i]*O[i];
			CM23P[i]+=h[i]*h[i]*h[i]*hp[i]*hp[i]*hp[i]*O[i];
			DS5[i]+=h[i]*hp[i]*hp[i]*hp[i]*hp[i]*hp[i]*O[i];
		}
		
		// (4.9)
		for(j=1; j<=i-1; j++){
			SA5[i] += 6*SA[i]*DSP[j] -6*CM[i]*SA[j];
			CM41Z[i] += -SA[i]*(2*ASP[j]-PT[j]) +2*CM[i]*(DSP[j]+CM[j]) 
				-(2*AS[i]-PT[i])*SA[j];
			CM41[i] += SA[i]*(4*ASP[j]+PT[j]) +2*CM[i]*(DSP[j]-2*CM[j]) 
				-(2*AS[i]+PT[i])*SA[j];
			SA32Z[i] += -SA[i]*CMP[j] +CM[i]*(AS[j]+2*PT[j]) 
				+2*PT[i]*CM[j] +AS[i]*DSP[j] -DS[i]*SA[j];
			SA32F[i] += 2*SA[i]*CMP[j] +2*CM[i]*(2*ASP[j]-AS[j]) -2*(2*AS[i]+PT[i])*CM[j];
			CM23P[i] += 4*CM[i]*CMP[j] +2*AS[i]*(ASP[j]-2*AS[j]) -2*DS[i]*CM[j];
			SA32[i] += 3*SA[i]*CMP[j] +CM[i]*(2*ASP[j]-3*AS[j]-PT[j]) 
				+(AS[i]+PT[i])*(DSP[j]-2*CM[j]) -DS[i]*SA[j];
			CM23Z[i] += -2*CM[i]*CMP[j] +AS[i]*(2*ASP[j]+2*AS[j]+3*PT[j])
				+PT[i]*(3*AS[j]+2*PT[j]) -2*DS[i]*CM[j];
			CM23[i] += SA[i]*SAP[j] +CM[i]*(4*CMP[j]-DS[j]) 
				+AS[i]*(ASP[j]-4*AS[j]-2*PT[j])
				+PT[i]*(ASP[j]-3*AS[j]-PT[j]) -DS[i]*CM[j];
			AS5[i] += 2*CM[i]*SAP[j] +2*AS[i]*(2*CMP[j]-DS[j])
				+PT[i]*CMP[j] -DS[i]*(4*AS[j]+PT[j]);
			SG5[i] += 4*CM[i]*SAP[j] +2*(AS[i]+PT[i])*(CMP[j]-2*DS[j])
				-2*DS[i]*(AS[j]+PT[j]);
			DS5[i] += 2*(3*AS[i]+PT[i])*SAP[j] -6*DS[i]*DS[j];
		}
		CM41ALL[i]=5*CM41[i] +CM41Z[i];                  // (4.4)
		SA32ALL[i]=SA32[i] +4*SA32F[i] +2*SA32Z[i];      // (4.4)
		CM23ALL[i]=3*CM23[i] +2*CM23P[i] +CM23Z[i];      // (4.4)
	}

	// �ΐS�����W�� //
	for(i=1; i<=k; i++){	
		al=N(i-1,1)*u[i];   al1=N(i,1)*u[i+1];
		alp=N(i-1,1)*up[i]; alp1=N(i,1)*up[i+1];

		PRE[i]=-2*(al1-al);
		DSE1[i]=-al*DS[i]+alp*AS[i];
		DSE2[i]=-alp*PT[i];
		ASE[i]=-al*AS[i]+alp*CM[i];
		PTE[i]=-al*PT[i];
		CME[i]=-al*CM[i]+alp*SA[i];

		for(ii=i+1; ii<=k; ii++){
			DSE1[i]+=(al1-al)*DS[ii]-(alp1-alp)*AS[ii];
			DSE2[i]+=(alp1-alp)*PT[ii];
			ASE[i] +=(al1-al)*AS[ii]-(alp1-alp)*CM[ii];
			PTE[i] +=(al1-al)*PT[ii];
			CME[i] +=(al1-al)*CM[ii]-(alp1-alp)*SA[ii];
		}
	}

	//  total of coefficients  //
	for(i=0; i<=k+1; i++){
		SAt+=SA[i]; CMt+=CM[i]; ASt+=AS[i]; DSt+=DS[i]; PTt+=PT[i];
		LCt+=LC[i]; TCt+=TC[i]; LC2t+=LC2[i]; TC2t+=TC2[i];
		SAPt+=SAP[i]; CMPt+=CMP[i]; ASPt+=ASP[i]; DSPt+=DSP[i]; LCPt+=LCP[i];
		SA5t+=SA5[i]; 
		CM41Zt+=CM41Z[i]; CM41t+=CM41[i]; CM41ALLt+=CM41ALL[i];
		SA32Ft+=SA32F[i]; SA32Zt+=SA32Z[i]; SA32t+=SA32[i]; SA32ALLt+=SA32ALL[i];
		CM23t+=CM23[i]; CM23Pt+=CM23P[i]; CM23Zt+=CM23Z[i]; CM23ALLt+=CM23ALL[i];
		AS5t+=AS5[i]; SG5t+=SG5[i];
		DS5t+=DS5[i];
		PREt+=PRE[i]; DSE1t+=DSE1[i]; DSE2t+=DSE2[i]; ASEt+=ASE[i]; PTEt+=PTE[i]; CMEt+=CME[i];
	}

	//  exponential  //
	SAe=exponent(SA);
	CMe=exponent(CM);
	ASe=exponent(AS);
	DSe=exponent(DS);
	PTe=exponent(PT);
	LCe=exponent(LC);
	TCe=exponent(TC);
	LC2e=exponent(LC2);
	TC2e=exponent(TC2);
	SA5e=exponent(SA5);
	SAPe=exponent(SAP);
	CMPe=exponent(CMP);
	ASPe=exponent(ASP);
	DSPe=exponent(DSP);
	LCPe=exponent(LCP);
	CM41Ze=exponent(CM41Z);
	CM41e=exponent(CM41);
	CM41ALLe=exponent(CM41ALL);
	SA32Fe=exponent(SA32F);
	SA32Ze=exponent(SA32Z);
	SA32e=exponent(SA32);
	SA32ALLe=exponent(SA32ALL);
	CM23e=exponent(CM23);
	CM23Pe=exponent(CM23P);
	CM23Ze=exponent(CM23Z);
	CM23ALLe=exponent(CM23ALL);
	AS5e=exponent(AS5);
	SG5e=exponent(SG5);
	DS5e=exponent(DS5);
	PREe=exponent(PRE);
	DSE1e=exponent(DSE1);
	DSE2e=exponent(DSE2);
	ASEe=exponent(ASE);
	PTEe=exponent(PTE);
	CMEe=exponent(CME);

	delete [] cv;
	delete [] hQ;  delete [] hD;
	delete [] hQp; delete [] hDp;
	delete [] P;
	delete [] O;
}

double cLens1::H(int i){
	if(i<1 || k<i) return 0;
	CalcCoefficients();
	return h[i];
}

double cLens1::Hp(int i){
	if(i<1 || k<i) return 0;
	CalcCoefficients();
	return hp[i];
}

double cLens1::U(int i){
	if(i<1 || k+1<i) return 0;
	CalcCoefficients();
	return u[i];
}

double cLens1::Up(int i){
	if(i<1 || k+1<i) return 0;
	CalcCoefficients();
	return up[i];
}

double cLens1::Al(int i){
	if(i<1 || k+1<i) return 0;
	CalcCoefficients();
	return u[i]*N(i-1,1);
}

double cLens1::Alp(int i){
	if(i<1 || k+1<i) return 0;
	CalcCoefficients();
	return up[i]*N(i-1,1);
}

double cLens1::Coefficient(std::string name,int sum_rms,int exponent,int total,int i1,int i2){
	// name     : �����W����
	// sum_rms  : 1�̂Ƃ������a�C2�̂Ƃ�root mean square.
	// exponent : �^�̂Ƃ��C�e�ʌW���ƑS�n�a�𓯈�w���̎w���`���ŕ\������Ƃ�
	//            �����_�ȏオ�ꌅ�ƂȂ�w��
	// total    : sum,rms��S�n�Ōv�Z
	// i2,i2    : sum,rms��i1�ʂ���i2�ʂŌv�Z
	//
	// sum_rms�̐^�l��exponent�̐^�l�ɗD�悷��D
	// total�̐^�l��i1,i2�ɗD�悷��D
	CalcCoefficients();

	if(total){
		i1=0;
		i2=k+1;
	}

	if(sum_rms==1){
		if(name=="hQ"      || name=="hq"     ) return sum(hQ,i1,i2);
		if(name=="hQp"     || name=="hqp"    ) return sum(hQp,i1,i2);
		if(name=="SA"      || name=="sa"     ) return sum(SA,i1,i2);
		if(name=="CM"      || name=="cm"     ) return sum(CM,i1,i2);
		if(name=="AS"      || name=="as"     ) return sum(AS,i1,i2);
		if(name=="DS"      || name=="ds"     ) return sum(DS,i1,i2);
		if(name=="PT"      || name=="pt"     ) return sum(PT,i1,i2);
		if(name=="LC"      || name=="lc"     ) return sum(LC,i1,i2);
		if(name=="TC"      || name=="tc"     ) return sum(TC,i1,i2);
		if(name=="LC2"     || name=="lc2"    ) return sum(LC2,i1,i2);
		if(name=="TC2"     || name=="tc2"    ) return sum(TC2,i1,i2);
		if(name=="SAP"     || name=="sap"    ) return sum(SAP,i1,i2);
		if(name=="CMP"     || name=="cmp"    ) return sum(CMP,i1,i2);
		if(name=="ASP"     || name=="asp"    ) return sum(ASP,i1,i2);
		if(name=="DSP"     || name=="dsp"    ) return sum(DSP,i1,i2);
		if(name=="LCP"     || name=="lcp"    ) return sum(LCP,i1,i2);
		if(name=="SA5"     || name=="sa5"    ) return sum(SA5,i1,i2);
		if(name=="CM41"    || name=="cm41"   ) return sum(CM41,i1,i2);
		if(name=="CM41Z"   || name=="cm41z"  ) return sum(CM41Z,i1,i2);
		if(name=="CM41ALL" || name=="cm41all") return sum(CM41ALL,i1,i2);
		if(name=="SA32"    || name=="sa32"   ) return sum(SA32,i1,i2);
		if(name=="SA32F"   || name=="sa32f"  ) return sum(SA32F,i1,i2);
		if(name=="SA32Z"   || name=="sa32z"  ) return sum(SA32Z,i1,i2);
		if(name=="SA32ALL" || name=="sa32all") return sum(SA32ALL,i1,i2);
		if(name=="CM23"    || name=="cm23"   ) return sum(CM23,i1,i2);
		if(name=="CM23P"   || name=="cm23p"  ) return sum(CM23P,i1,i2);
		if(name=="CM23Z"   || name=="cm23z"  ) return sum(CM23Z,i1,i2);
		if(name=="CM23ALL" || name=="cm23all") return sum(CM23ALL,i1,i2);
		if(name=="AS5"     || name=="as5"    ) return sum(AS5,i1,i2);
		if(name=="SG5"     || name=="sg5"    ) return sum(SG5,i1,i2);
		if(name=="DS5"     || name=="ds5"    ) return sum(DS5,i1,i2);
		if(name=="PRE"     || name=="pre"    ) return sum(PRE,i1,i2);
		if(name=="DSE1"    || name=="dse1"   ) return sum(DSE1,i1,i2);
		if(name=="DSE2"    || name=="dse2"   ) return sum(DSE2,i1,i2);
		if(name=="ASE"     || name=="ase"    ) return sum(ASE,i1,i2);
		if(name=="PTE"     || name=="pte"    ) return sum(PTE,i1,i2);
		if(name=="CME"     || name=="cme"    ) return sum(CME,i1,i2);
		else return 0;
	}
	if(sum_rms==2){
		if(name=="hQ"      || name=="hq"     ) return rms(hQ,i1,i2);
		if(name=="hQp"     || name=="hqp"    ) return rms(hQp,i1,i2);
		if(name=="SA"      || name=="sa"     ) return rms(SA,i1,i2);
		if(name=="CM"      || name=="cm"     ) return rms(CM,i1,i2);
		if(name=="AS"      || name=="as"     ) return rms(AS,i1,i2);
		if(name=="DS"      || name=="ds"     ) return rms(DS,i1,i2);
		if(name=="PT"      || name=="pt"     ) return rms(PT,i1,i2);
		if(name=="LC"      || name=="lc"     ) return rms(LC,i1,i2);
		if(name=="TC"      || name=="tc"     ) return rms(TC,i1,i2);
		if(name=="LC2"     || name=="lc2"    ) return rms(LC2,i1,i2);
		if(name=="TC2"     || name=="tc2"    ) return rms(TC2,i1,i2);
		if(name=="SAP"     || name=="sap"    ) return rms(SAP,i1,i2);
		if(name=="CMP"     || name=="cmp"    ) return rms(CMP,i1,i2);
		if(name=="ASP"     || name=="asp"    ) return rms(ASP,i1,i2);
		if(name=="DSP"     || name=="dsp"    ) return rms(DSP,i1,i2);
		if(name=="LCP"     || name=="lcp"    ) return rms(LCP,i1,i2);
		if(name=="SA5"     || name=="sa5"    ) return rms(SA5,i1,i2);
		if(name=="CM41"    || name=="cm41"   ) return rms(CM41,i1,i2);
		if(name=="CM41Z"   || name=="cm41z"  ) return rms(CM41Z,i1,i2);
		if(name=="CM41ALL" || name=="cm41all") return rms(CM41ALL,i1,i2);
		if(name=="SA32"    || name=="sa32"   ) return rms(SA32,i1,i2);
		if(name=="SA32F"   || name=="sa32f"  ) return rms(SA32F,i1,i2);
		if(name=="SA32Z"   || name=="sa32z"  ) return rms(SA32Z,i1,i2);
		if(name=="SA32ALL" || name=="sa32all") return rms(SA32ALL,i1,i2);
		if(name=="CM23"    || name=="cm23"   ) return rms(CM23,i1,i2);
		if(name=="CM23P"   || name=="cm23p"  ) return rms(CM23P,i1,i2);
		if(name=="CM23Z"   || name=="cm23z"  ) return rms(CM23Z,i1,i2);
		if(name=="CM23ALL" || name=="cm23all") return rms(CM23ALL,i1,i2);
		if(name=="AS5"     || name=="as5"    ) return rms(AS5,i1,i2);
		if(name=="SG5"     || name=="sg5"    ) return rms(SG5,i1,i2);
		if(name=="DS5"     || name=="ds5"    ) return rms(DS5,i1,i2);
		if(name=="PRE"     || name=="pre"    ) return rms(PRE,i1,i2);
		if(name=="DSE1"    || name=="dse1"   ) return rms(DSE1,i1,i2);
		if(name=="DSE2"    || name=="dse2"   ) return rms(DSE2,i1,i2);
		if(name=="ASE"     || name=="ase"    ) return rms(ASE,i1,i2);
		if(name=="PTE"     || name=="pte"    ) return rms(PTE,i1,i2);
		if(name=="CME"     || name=="cme"    ) return rms(CME,i1,i2);
		else return 0;
	}
	if(exponent){
		if(name=="SA"      || name=="sa"     ) return SAe;
		if(name=="CM"      || name=="cm"     ) return CMe;
		if(name=="AS"      || name=="as"     ) return ASe;
		if(name=="DS"      || name=="ds"     ) return DSe;
		if(name=="PT"      || name=="pt"     ) return PTe;
		if(name=="LC"      || name=="lc"     ) return LCe;
		if(name=="TC"      || name=="tc"     ) return TCe;
		if(name=="LC2"     || name=="lc2"    ) return LC2e;
		if(name=="TC2"     || name=="tc2"    ) return TC2e;
		if(name=="SAP"     || name=="sap"    ) return SAPe;
		if(name=="CMP"     || name=="cmp"    ) return CMPe;
		if(name=="ASP"     || name=="asp"    ) return ASPe;
		if(name=="DSP"     || name=="dsp"    ) return DSPe;
		if(name=="LCP"     || name=="lcp"    ) return LCPe;
		if(name=="SA5"     || name=="sa5"    ) return SA5e;
		if(name=="CM41"    || name=="cm41"   ) return CM41e;
		if(name=="CM41Z"   || name=="cm41z"  ) return CM41Ze;
		if(name=="CM41ALL" || name=="cm41all") return CM41ALLe;
		if(name=="SA32"    || name=="sa32"   ) return SA32e;
		if(name=="SA32F"   || name=="sa32f"  ) return SA32Fe;
		if(name=="SA32Z"   || name=="sa32z"  ) return SA32Ze;
		if(name=="SA32ALL" || name=="sa32all") return SA32ALLe;
		if(name=="CM23"    || name=="cm23"   ) return CM23e;
		if(name=="CM23P"   || name=="cm23p"  ) return CM23Pe;
		if(name=="CM23Z"   || name=="cm23z"  ) return CM23Ze;
		if(name=="CM23ALL" || name=="cm23all") return CM23ALLe;
		if(name=="AS5"     || name=="as5"    ) return AS5e;
		if(name=="SG5"     || name=="sg5"    ) return SG5e;
		if(name=="DS5"     || name=="ds5"    ) return DS5e;
		if(name=="PRE"     || name=="pre"    ) return PREe;
		if(name=="DSE1"    || name=="dse1"   ) return DSE1e;
		if(name=="DSE2"    || name=="dse2"   ) return DSE2e;
		if(name=="ASE"     || name=="ase"    ) return ASEe;
		if(name=="PTE"     || name=="pte"    ) return PTEe;
		if(name=="CME"     || name=="cme"    ) return CMEe;
		else return 0;
	}
	else{
		return 0;
	}
}

double cLens1::Coefficient(std::string name){
	return Coefficient(name,1,0,1,0,0);
}


matrix<double> cLens1::CMatrix(int i){
	// i�ʂ̓����s�� ("�����_(����)" (3.5.5b) �ɂ��)
	// �y���Ӂz
	//    ���O�ł� this->phi,h,hp,u,up �̐ݒ�iCalcCoefficien()�̎��s�Ȃǂɂ��j���K�v
	//    CalcCoefficients()�ɂ��ꍇ�͌��ʂ� NormalizeUnit,NormalizeType �ɍ��E�����
	// �y�Q�l�z
	//    �{�֐��̏������e��cMatrix(i1,i2)�Ɋ܂܂��āCcMatrix(i,i)�Ƃ��Ă��悢���C
	//    cMatrix(i1,i2)�͓�����CalcCoefficients()���Ă�ł��邽�߁CCalcCoefficients()�̓�������
	//    cMatrix(i1,i2)���ĂԂƎ��ȎQ�ƂɂȂ�G���[����������D

	matrix<double> A(6,3);
	double p,N,N1,h,hp,a,ap,P;
	
	if( i<1 || i>k ) return A;

	p=this->phi[i];
	N =this->N(i-1,1);
	N1=this->N(i,1);
	h =this->h[i];
	hp=this->hp[i];
	a =N*this->u[i];
	ap=N*this->up[i];

	if(is_ideallens(i) && N*N1>0){  // ���z�����Y�Ŕ��˖ʂłȂ��Ƃ�
		P=p/N/1.52;                 // �����n�����ܗ�N�̔}�����ɔz�u���ꂽ���ܗ�1.52�̔������Y�Q�Ɖ���
	}
	else{
		P=p/N/N1;                   // ���z�����Y�Ŕ��˖ʂł���Ƃ�(�P�ꔽ�˖ʂƉ���)���܂�
	}
	
	A[1][1]=h*h*h*h*p*p*p;
	A[1][2]=-4*a*h*h*h*p*p;
	A[1][3]=-a*h*h*h*p*p/N1/N1 +3*a*a*h*h*p/N1/N1 +2*a*a*h*h*P +(1/N1/N1-1/N/N)*a*a*a*h;

	A[2][1]=h*h*h*hp*p*p*p;
	A[2][2]=-h*h*p*p*(4*a*hp-1);
	A[2][3]=-a*h*h*hp*p*p/N1/N1 +a*h*p*(3*a*hp-2)/N1/N1
	        +a*h*P*(2*a*hp-1) +(1/N1/N1-1/N/N)*a*a*ap*h;

	A[3][1]=h*h*hp*hp*p*p*p;
	A[3][2]=-2*h*hp*p*p*(2*ap*h+1);
	A[3][3]=-a*h*hp*hp*p*p/N1/N1 +ap*h*p*(3*a*hp-1)/N1/N1
	        +2*a*ap*h*hp*P +(1/N1/N1-1/N/N)*a*ap*ap*h;

	A[4][1]=0;
	A[4][2]=0;
	A[4][3]=P;

	A[5][1]=h*hp*hp*hp*p*p*p;
	A[5][2]=-hp*hp*p*p*(4*ap*h+1);
	A[5][3]=-a*hp*hp*hp*p*p/N1/N1 +3*ap*ap*h*hp*p/N1/N1 +ap*hp*P*(2*a*hp-1)
	        +(1/N1/N1-1/N/N)*ap*ap*ap*h;

	A[6][1]=hp*hp*hp*hp*p*p*p;
	A[6][2]=-4*ap*hp*hp*hp*p*p;
	A[6][3]=-ap*hp*hp*hp*p*p/N1/N1 +3*ap*ap*hp*hp*p/N1/N1 +2*ap*ap*hp*hp*P
	        +(1/N1/N1-1/N/N)*ap*ap*ap*hp;

	return A;
}

matrix<double> cLens1::CMatrix(int i1,int i2){
	// i1�ʂ���i2�ʂ܂ł̓����s�� 
	
	cLens1 x;

	x=*this;
	if( i1<1 || i1>i2 || i2>x.k ) return matrix<double>(6,3);
	x.ToIdealLens(i1,i2);
	x.CalcCoefficients();
	return x.CMatrix(i1);
}

double cLens1::CMatrix(int i1,int i2,int i,int j){
	matrix<double> A;
	A=CMatrix(i1,i2);
	return A[i][j];
}

std::string cLens1::CMatrixStr(int i1,int i2){
	return str(CMatrix(i1,i2));
}

double cLens1::SA0(int i1,int i2){
	// �����n�ŗL�W�� A0
	cLens1 x;
	
	x=*this;
	x.Trim(i1,i2);
	x.ToThinLens(1,x.k);
	x.Set_NormalizeUnit("f=1");
	x.Set_NormalizeType(1);
	x.s=LN;
	x.t=0;
	return x.Coefficient("SA");
}
double cLens1::SA0(){
	return SA0(1,k);
}

double cLens1::CM0(int i1,int i2){
	// �����n�ŗL�W�� B0
	cLens1 x;
	
	x=*this;
	x.Trim(i1,i2);
	x.ToThinLens(1,x.k);
	x.Set_NormalizeUnit("f=1");
	x.Set_NormalizeType(1);
	x.s=LN;
	x.t=0;
	return x.Coefficient("CM");
}
double cLens1::CM0(){
	return CM0(1,k);
}

double cLens1::cDerivative(std::string command,int i,double dc){
	// �����W��df/dc���v�Z����
	if(1<=i && i<=k){
		return Derivative(command,&c(i),dc);
	}
	else return 0;
}

double cLens1::dDerivative(std::string command,int i,double dd){
	// �����W��df/dd���v�Z����
	if(1<=i && i<=k){
		return Derivative(command,&d(i),dd);
	}
	else return 0;
}

void cLens1::ToScheimpflugImagePlane(){
	// �ߎ����_�iScheimpflug���_�j�ɂ����
	// ���� ax+by+cz+d=0 (k+1�ʕΐS�O���W) �����߁C
	// ���w�n�̑��ʁik+1�ʁj��ΐS�ɂ���Ă���Ɉ�v������D
	int i,i0;
	double a,b,c,d;

	decenter_type(k+1)=0;  // ���ʂ̕ΐS���Ȃ����Ă����D��������ƕΐS�ʂ̌v�Z���y�D
	s1fix=0;               // s1fix����������D�������Ȃ��ƁC
	                       // �{�֐���OptimizeS1Fix()�̌��݌J��Ԃ���
	                       // s1fix�Ƒ��ʂ�dz�̘a�͕ۂ������̂̔z�����ς��Ă����C
	                       // s1fix�����ł͑��ʂ̈ʒu���킩��ɂ����Ȃ��Ă��܂��D
	make_coordinate(0);

	if(fabs(s)<LN){
		a=b=d=0; c=1;  // ����z=0
		for(i=1; i<=k; i++){
			transform_plane(a,b,c,d,i-1,0,i,0);
			image_plane(a,b,c,d,i);	
		}
		transform_plane(a,b,c,d,k,0,k+1,0);
	}
	else{
		i0=0;
		for(i=1; i<=k; i++){
			if( power(i,i,1)!=0 ){
				i0=i;   // �ŏ��̃p���[��0�łȂ���
				a=b=0;
				c=1;
				d=-bf(i,i,1); // �ŏ��Ƀp���[�����ʂ̑����œ_��
				break;
			}
		}
		if(i0==0) return;  // �S�Ă̖ʂŃp���[��0
		for(i=i0+1; i<=k; i++){
			transform_plane(a,b,c,d,i-1,0,i,0);
			image_plane(a,b,c,d,i);	
		}
		transform_plane(a,b,c,d,k,0,k+1,0);
	}

	if(c<0){
		// �Ȃ�ׂ��ȒP�ȕΐS�ƂȂ�悤��(a,b,c)�𑜖ʍ��Wz�������Ƃ̓��ς����ɂȂ�����Ɍ�����
		a*=-1;
		b*=-1;
		c*=-1;
		d*=-1;
	}

	if(c!=0){
		dx(k+1)=0;
		dy(k+1)=0;
		dz(k+1)=-d/c;
	}
	else if(b!=0){
		dx(k+1)=0;
		dy(k+1)=-d/b;
		dz(k+1)=0;	
	}
	else if(a!=0){
		dx(k+1)=-d/a;
		dy(k+1)=0;
		dz(k+1)=0;
	}
	else{
		return;
		// a=b=c=0 �͖��������ʂƍl������D
	}

	decenter_type(k+1)=2;
	order(k+1)=0;
	rotate_angle(rox(k+1),roy(k+1),a,b,c);
	roz(k+1)=0;
}

double cLens1::ScheimpflugImagePlaneTiltX(){
	cLens1 X;

	X=*this;
	X.ToScheimpflugImagePlane();
	return X.rox(k+1);
}

double cLens1::ScheimpflugImagePlaneTiltY(){
	cLens1 X;

	X=*this;
	X.ToScheimpflugImagePlane();
	return X.roy(k+1);
}

double cLens1::ImagePlaneTiltX(int i,double yObj,double xObj,
							  std::string SetRay,int FindPupil,double yPupil,double xPupil,double dyObj/*=0*/){
	// ������ߖT�ł̑��ʂ�yz�ʓ��ł̌X�����v�Z����D��i�ʂ̑����Ōv�Z���C�{�����悶�čŏI���ʂł̒l��Ԃ��D
	// ���w�n��yz�ʂɊւ��đΏ̂�z�肵�Ă���D
	vector<double> Ro,R,Q1o,Q1,S1o,S1, Ez,Ey, dQ1,dQ1k;
	double dy,dz, m;
	matrix<double> T;

	if(i<=0 || k<i) i=k;  // i�������Ȓl�̂Ƃ��͍ŏI�ʂƂ���
	
	SingleRayTrace(yObj,xObj,SetRay,FindPupil,yPupil,xPupil,0,0,1);
	Ro=vector<double>(x[i],y[i],z[i]);
	Q1o=vector<double>(X1[i],Y1[i],Z1[i]);
	S1o=Ro+Q1o*(-ha1[i].y/dQa1[i].y);  // ��������_�ʒu�ii�ʒ��_���W�j

	T=Tmatrix(vector<double>(X1[i],Y1[i],Z1[i]));  // �o�ˌ������W�ɕϊ�����s��
	dQ1=T*dQa1[i];
	T=Tmatrix(vector<double>(X1[k],Y[k],Z[k]));    // �ŏI�ʂ̏o�ˌ������W�ɕϊ�����s��
	dQ1k=T*dQa1[k];
	m=N(i,1)*dQ1.y/(N(k,1)*dQ1k.y);                // i+1�ʁ`�ŏI�� �ɂ�������ɉ�����y�����̉��{��

	if(dyObj==0) dyObj=0.01;
	SingleRayTrace(yObj+dyObj,xObj,SetRay,FindPupil,yPupil,xPupil,0,0,1);
	R=vector<double>(x[i],y[i],z[i]);
	Q1=vector<double>(X1[i],Y1[i],Z1[i]);
	S1=R+Q1*(-ha1[i].y/dQa1[i].y);    // �ߖT�������_�ʒu�ii�ʒ��_���W�j
	
	Ez=Q1o*sgn(sProduct(Q1o,vector<double>(0,0,1))); // ������ɉ�����z�����𐳂Ƃ���P�ʃx�N�g��
	Ey=vProduct(Ez,vector<double>(1,0,0));           // �q�ߖʓ��Ŏ�����ɐ�����y�����𐳂Ƃ���P�ʃx�N�g��
	
	dz=sProduct(S1-S1o,Ez);
	dy=sProduct(S1-S1o,Ey);

	dz*=m*m;   // �ŏI���ʂɕϊ�����
	dy*=m;

	return atan(dz/dy)*180/PI;
}

double cLens1::ImagePlaneTiltX(int i){
	return ImagePlaneTiltX(i,0,0,"",0,0,0);
}

double cLens1::OffAxialMirrorComa(int i){
	// ��i�ʂł���~���[�ʂɂ�鎋�쒆���̃R�}�������v�Z����Di=0�̂Ƃ��͑S�n�̘a��Ԃ��D
	// �E���w�n��yz�ʂɊւ��đΏ̂Ƃ���D
	// �E�~���[�ʂ̕ΐS��decenter_type(i)=3�ł�����̂Ƃ���D
	// �E�X�̖ʂɂ��Ă̓����ʂɖ�������3���R�}�����W������v�Z����D
	int i1,i2;
	matrix<double> T;
	vector<double> H,dQ,dQ1,dQ1k;
	double dummy, s,c, al,al1,alp,hQ,hQp,hD,cm,dy, sum;

	RayTrace(0,0,0,0,0,1,0,0,0,1);  // ���̍�0�̎�����i������j�ɉ����Ĕ�_�������ǐՂ��s��

	if(i==0){
		i1=1; i2=k;
	}
	else if(1<=i && i<=k){
		i1=i2=i;
	}
	else{
		return 0;
	}

	sum=0;
	for(i=i1; i<=i2; ++i){
		T=Tmatrix(vector<double>(X[i],Y[i],Z[i]));     // ���ˌ������W�ɕϊ�����s��
		H=T*ha[i]; dQ=T*dQa[i];

		s= dQ.y==0 ? LN : H.y*cos(rox(i)*PI/180)/dQ.y; // ���ˌ����ɉ��������̋���
		s*=cos(rox(i)*PI/180);                         // �ʖ@���ɉ��������̋���
		ParaxialExpansion(i,dummy,dummy,c);            // y�����ߎ��ȗ����a
		
		al=N(i-1,1)/s;    alp=-1;                      // "�����Y�݌v�@" (4.1)
		hQ=N(i-1,1)*c-al; hQp=1;                       // "�����Y�݌v�@" (4.5)  (h=1,hp=0)
		al1=al+(N(i,1)-N(i-1,1))*c;                    // "�����Y�݌v�@" (2.12)
		hD=al1/N(i,1)/N(i,1)-al/N(i-1,1)/N(i-1,1);     // "�����Y�݌v�@" (4.5)
		cm=hQ*hQp*hD;                                  // "�����Y�݌v�@" (4.6a)

		dy=-1.0/(al1*2)*cm*N(i-1,1)*tan(rox(i)*PI/180)*(H.y*H.y)*3;  // �R�}����, "�����Y�݌v�@" (4.4)
		dy*=cos(rox(i)*PI/180);                                      // ������ɐ����Ȑ���
		
		T=Tmatrix(vector<double>(X1[i],Y1[i],Z1[i]));  // �o�ˌ������W�ɕϊ�����s��
		dQ1=T*dQa1[i];
		T=Tmatrix(vector<double>(X1[k],Y[k],Z[k]));    // �ŏI�ʂ̏o�ˌ������W�ɕϊ�����s��
		dQ1k=T*dQa1[k];
		dy*=N(i,1)*dQ1.y/(N(k,1)*dQ1k.y);              // ���{�����|���čŏI����Ԃł̒l�ɂ���

		sum+=dy;
	}

	return sum;
}


int cLens1::RayTrace(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j,
                    int EA_enable,int mask_enable,int IsLambert,int as_trace, 
					int E_trace,const vector<complex>& E0, int MakeCoordinate){
	// MakeCoordinate=0 �Ƃ���� make_coordinate() �����s���Ȃ��D�������̌��ʂ�����D
	// �������C�{�֐��̒��O�Ŏ��s���Ă����K�v������D
	int i, errorcode;
	double o;
	double x0,y0,z0, x,y,z, X,Y,Z, X1,Y1,Z1;
	vector<double> ha,dQa, hb,dQb;
	vector<complex> E;

	if(s==0 && t==0) return INVALID;

	if(MakeCoordinate) make_coordinate(defocus);

	gname(k+1)=gname(k); // ���ʂł͋��܂��Ȃ��Ƃ���

	GrinRay_i.RemoveAll();
	GrinRay.RemoveAll();

	this->y_pupil=yPupil; this->x_pupil=xPupil;  // Apodization�v�Z�p�ɕۑ�
	yPupil+=EPy; xPupil+=EPx;

	if(fabs(s)<LN){
		z0=surface_sag(0,yObj,xObj,0); y0=yObj; x0=xObj;

		// x0,y0,z0�𕨑̖ʒ��_���W�����1�ʕΐS�O�̍��W�֕ϊ�����D
		// ���˓��͂��̍��W����ɂ���D
		transform(x0,y0,z0,0,0,1,1,1);

		if(fabs(t)<LN){
			// Y=(yPupil-yObj)/(t-s-z0); X=(xPupil-xObj)/(t-s-z0);
			// o=sqrt(1+Y*Y+X*X);
			// Z=1/o; Y=Y/o; X=X/o;
			//
			// ����ł�t-s-z0<0�̂Ƃ�Q=(X,Y,Z)�ƌ����̐i�s��������v���Ȃ��̂ŕύX(06.06.14)�D
			//   �Ⴆ�΁C���̂���1�ʂ��E(z>0�̕���)�ɂ���Ƃ��g�ʎ��������]���Ă��܂��D
			// X=xPupil-xObj; Y=yPupil-yObj; Z=t-s-z0;
			// o=sqrt(X*X+Y*Y+Z*Z);
			// X/=o; Y/=o; Z/=o;
			//
			// 06.10.18 "*sgn(Z)*sgn(N(0,j))" ��ǉ��D
			//   Q��Z�����̕����͕���ԋ��ܗ��̕����ƈ�v������D
			//   �������Ȃ��ƁC�Ⴆ�Δ��˂̂Ȃ��n��s>0,t=0�̂Ƃ��C
			//   ���̂����1�ʂ̊Ԃ�Q���t�ɂȂ�C���H���̌v�Z�ɕs�s�����N����D
			X=xPupil-x0; Y=yPupil-y0; Z=t-z0;
			o=sqrt(X*X+Y*Y+Z*Z)*sgn(Z)*sgn(N(0,j));
			if( fabs(o)<1e-6 ) return INVALID; 
			// if( o==0 )�ł͕��̖ʂɕΐS������ꍇ�C���S��0�ɂȂ炸��ŃG���[���������邱�Ƃ�����D
			// (09.10.07)
			o=1/o;  X*=o; Y*=o; Z*=o;
		}
		else{
			// if t=inf, Y/Z=yPupil and X/Z=xPupil.
			o=1/sqrt(1+yPupil*yPupil+xPupil*xPupil);
			Z=o; Y=yPupil*o; X=xPupil*o;
		}

		// ���̖ʒ��_���W�ɖ߂�
		transform(x0,y0,z0,1,1,0,0,1);
		transform(X,Y,Z,1,1,0,0,0);

		if(as_trace){
			matrix<double> T;
			ha=hb=vector<double>(0,0,0);
			dQa=vector<double>(0,1,0); dQb=vector<double>(1,0,0);
			T=inv( Tmatrix(vector<double>(X,Y,Z)) );
			dQa=T*dQa; dQb=T*dQb;                     // ���̖ʒ��_���W�
		}

		if(E_trace){
			E=E0;    // ���̖ʒ��_���W�
		}

		optical_path[0]=0;
	}
	else {
		// ���������̂̕ΐS�͂Ƃ肠�����l���Ȃ� (09.02.12)
		d(0)=0;
		z0=0; y0=yPupil-t*yObj; x0=xPupil-t*xObj;
		o=1/sqrt(1+yObj*yObj+xObj*xObj);
		Z=o*sgn(N(0,j)); Y=yObj*o; X=xObj*o;

		if(as_trace){
			matrix<double> T;
			ha=vector<double>(0,1,0); hb=vector<double>(1,0,0);
			dQa=dQb=vector<double>(0,0,0);
			T=inv( Tmatrix(vector<double>(X,Y,Z)) );
			ha=T*ha; hb=T*hb;
			dQa=T*dQa; dQb=T*dQb;
		}

		if(E_trace){
			E=E0;
		}

		{
			double k,xr,yr,zr;
			k=X*x0+Y*y0+Z*z0;
			xr=x0-k*X; yr=y0-k*Y; zr=z0-k*Z;
			// (xr,yr,zr): ���_��ʂ�����i�s����(X,Y,Z)�ɐ����ȕ��ʂƌ����̌�_
			optical_path[0]=-N(0,j)*sgn( (xr-x0)*X+(yr-y0)*Y+(zr-z0)*Z )
			                *sqrt( (xr-x0)*(xr-x0)+(yr-y0)*(yr-y0)+(zr-z0)*(zr-z0) );
		}
	}

	if(IsLambert) {
		if( Random(0,1,0)>Z*Z*Z*Z ) return INVALID;
	}

	d(k)=dk(defocus);

	// ���̖ʂɊւ���ʂ�ۑ�����D��1�ʈȍ~�� raytrace() ���ŕۑ������D			
	this->x[0]=x0; this->y[0]=y0; this->z[0]=z0;
	this->X[0]=this->X1[0]=X; 
	this->Y[0]=this->Y1[0]=Y; 
	this->Z[0]=this->Z1[0]=Z;

	this->xi1[0]=cosine(vector<double>(this->X1[0],this->Y1[0],this->Z1[0]),
			            surface_normal(0,this->y[0],this->x[0],1));

	if(E_trace){
		this->E[0]=this->E1[0]=E;
	}

	if(as_trace){
		this->ha[0]=this->ha1[0]=ha;
		this->hb[0]=this->hb1[0]=hb;
		this->dQa[0]=this->dQa1[0]=dQa;
		this->dQb[0]=this->dQb1[0]=dQb;
	}

	for(i=1; i<=k+1; ++i){
		errorcode=raytrace(i,j, x0,y0,z0,X,Y,Z, x,y,z,X1,Y1,Z1, EA_enable,mask_enable,
			               as_trace,ha,dQa,hb,dQb, E_trace,E);			
		
		if(errorcode!=0){
			if(Diffusion(i-1)==2 && i>=2){
				// Diffusion(i)==2�̂Ƃ��C
				// �������̂��ߊg�U�ʂ́g���̖ʁh��ʉ߂���܂ŌJ��Ԃ��D
				// (�����ʉ߂����Ȃ��ʂł���Ζ������[�v�ƂȂ�̂Œ���)
				// (�g�U�ʏ�̊e�_���猩���Ƃ��Ɏ��̖ʂ̌��������قȂ�ꍇ��
				//  �s���R�ȏd�݂Â��ƂȂ�̂Œ��Ӂj
				x0=this->x[i-2]; y0=this->y[i-2]; z0=this->z[i-2];
				X=this->X1[i-2]; Y=this->Y1[i-2]; Z=this->Z1[i-2];
				i-=2;				
			}
			else{
				return errorcode;
			}
		}
		else{
			x0=x; y0=y; z0=z; X=X1; Y=Y1; Z=Z1;
		}
	}

	return 0;
}

int cLens1::NotThruSurf(int errorcode){
	return errorcode%NOT_THRU;
}

std::string cLens1::ErrorMessage(int errorcode){
	std::string s;
	errorcode=errorcode-NotThruSurf(errorcode);
	switch(errorcode){
		case NOT_THRU: 
			s="not thru";
			break;
		case NOT_INTERSECT: 
			s="not intersect";
			break;
		case CLIPPED: 
			s="clipped";
			break;
		case TOTAL_REF: 
			s="total reflection";
			break;
		case INVALID:
			s="invalid";
			break;
		case GRIN_OUT:
			s="grin out";
			break;
	}
	return s;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int cLens1::FindAThruRay(double& yPupil,double& xPupil,double yObj,double xObj,
						 int j,int priority,int makecoordinate)
{
	const double TAN=10;
	int i,ii, ix,iy;
	point ymax,ymin,xmax,xmin, p;
	vector<double> v, po, y1max,y1min,x1max,x1min;
	bool is_mask, b;
	static list<point> cache;
	const int CACHE_SIZE=300;
	clock_t start,now;
	vector<complex> dummy;

	cache.SetMaxSize(CACHE_SIZE);

	if(makecoordinate!=0) make_coordinate(0);
	
	if(EPD==0){
		if( RayTrace(yObj,xObj,0,0,0,j,1,0,0,0,0,dummy,0)==0 ){
			yPupil=0; xPupil=0;
			return 1;
		}
	}
	else{
		for(ix=0; ix<=30; ++ix) for(iy=0; iy<=30; ++iy){
			// ���˓��̋ߖT�ł܂����݂�
			// <note>
			//   EPCalc()�Ȃǂɂ��Ct��EPD�ɐݒ肵�Ă������Ƃɂ��C
			//   ��1�ʂł̌������L���a�ɑ΂��Ĕ��ɍׂ��ꍇ�ł�
			//   ���̒i�K�œ��ߌ����𔭌��ł���ꍇ������D
			p.x=EPD/2*ix;
			p.y=EPD/2*iy;

			if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ){
				yPupil=p.y; xPupil=p.x;
				return 1;
			}

			if(ix!=0){
				p.x=-EPD/2*ix;
				p.y= EPD/2*iy;
				if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ){
					yPupil=p.y; xPupil=p.x;
					return 1;
				}
			}

			if(iy!=0){
				p.x= EPD/2*ix;
				p.y=-EPD/2*iy;
				if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ){
					yPupil=p.y; xPupil=p.x;
					return 1;
				}
			}

			if(ix!=0 && iy!=0){
				p.x=-EPD/2*ix;
				p.y=-EPD/2*iy;
				if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ){
					yPupil=p.y; xPupil=p.x;
					return 1;
				}
			}
		}
	}
	
	for(i=1; i<=cache.GetSize(); i++){
		// ���ɃL���b�V�����̓_�Ŏ��݂�
		cache.GetData(p,i);
		if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ){
			yPupil=p.y; xPupil=p.x;
			return 1;
		}
	}
	
	for(ii=1; ii<=3; ii++){
		// �T���͈͂̉��ɑΉ�����4�_�����߂�D
		// �����𕨓_�������Ƃ��ē��˓��ʂɓ��e���T���͈͂̉��Ƃ���D
		// <note>
		//      �����̓_��P���ɑ�1�ʐڕ��ʏ��EAy(1),EAx(1)�̑傫���Ƃ���ƁC
		//    ��1�ʂ̃T�O���傫���Ƃ��C�������ʂɓ�����ʒu�Ɛڕ��ʂɓ�����ʒu�̍����̍���
		//    �傫�����߁C�m���ɒT�����Ƃ��ł��Ȃ������D
		//      �܂��C��1�ʏ��EAy(1),EAz(1)�ɑΉ�����_(�܂�T�O���l������)�Ƃ��Ă�
		//    �m���ɒT�����Ƃ��ł��Ȃ����D���Ƃ��΁C���_���v���`�h�����O�C��1�ʂ�
		//    �p���ʂŕ��_�������Ƃ��p�����_�͂��̒T���͈͊O�ƂȂ��Ă��܂����Ƃ�����D
		//      �����Ŏ��̂悤�ɁCiii=1,2,3���J��Ԃ����@���l�����D (10.05.21)
		if(EAy(1)>0 && EAx(1)==0){
			y1max.x=0;                 y1max.y= EAy(1)/2+EAdy(1);
			y1min.x=0;                 y1min.y=-EAy(1)/2+EAdy(1);
			x1max.x= EAy(1)/2+EAdx(1); x1max.y=0; 
			x1min.x=-EAy(1)/2+EAdx(1); x1min.y=0;
			switch(ii){
			case 1:
				// max,min���ʏ�
				y1max.z=surface_sag(1,y1max.y,y1max.x,0);
				y1min.z=surface_sag(1,y1min.y,y1min.x,0);
				x1max.z=surface_sag(1,x1max.y,x1max.x,0); 
				x1min.z=surface_sag(1,x1min.y,x1min.x,0);
				break;
			case 2:
				// max�͖ʏ�Cmin�͐ڕ��ʏ�
				y1max.z=surface_sag(1,y1max.y,y1max.x,0);
				y1min.z=0;
				x1max.z=surface_sag(1,x1max.y,x1max.x,0); 
				x1min.z=0;
				break;
			case 3:
				// max�͐ڕ��ʏ�Cmin�͖ʏ�
				y1max.z=0;
				y1min.z=surface_sag(1,y1min.y,y1min.x,0);
				x1max.z=0;
				x1min.z=surface_sag(1,x1min.y,x1min.x,0);
				break;
			}
			is_mask=false;
		}
		else if(EAy(1)<0 && EAx(1)==0){
			is_mask=true;
			ii=3;
		}
		else if(EAy(1)>0 && EAx(1)>0){
			y1max.x=0;                 y1max.y= EAy(1)/2+EAdy(1);
			y1min.x=0;                 y1min.y=-EAy(1)/2+EAdy(1);
			x1max.x= EAx(1)/2+EAdx(1); x1max.y=0;
			x1min.x=-EAx(1)/2+EAdx(1); x1min.y=0;
			switch(ii){
			case 1:
				y1max.z=surface_sag(1,y1max.y,y1max.x,0); 
				y1min.z=surface_sag(1,y1min.y,y1min.x,0);
				x1max.z=surface_sag(1,x1max.y,x1max.x,0);
				x1min.z=surface_sag(1,x1min.y,x1min.x,0);
				break;
			case 2:
				y1max.z=surface_sag(1,y1max.y,y1max.x,0); 
				y1min.z=0;
				x1max.z=surface_sag(1,x1max.y,x1max.x,0);
				x1min.z=0;
				break;
			case 3:
				y1max.z=0;
				y1min.z=surface_sag(1,y1min.y,y1min.x,0);
				x1max.z=0;
				x1min.z=surface_sag(1,x1min.y,x1min.x,0);
				break;
			}
			is_mask=false;
		}
		else if(EAy(1)<0 && EAx(1)<0){
			is_mask=true;
			ii=3;
		}
		else return 0;

		// ���˓��͑�1�ʕΐS�O���W�Ɠ����ł��邽�߁C�ΐS�O���W��ɕϊ�����D
		transform(y1max.x,y1max.y,y1max.z,1,0,1,1,1);
		transform(y1min.x,y1min.y,y1min.z,1,0,1,1,1);
		transform(x1max.x,x1max.y,x1max.z,1,0,1,1,1);
		transform(x1min.x,x1min.y,x1min.z,1,0,1,1,1);

		// ���̍��Wpo���1�ʕΐS�O���W��ŕ\���D
		po.y=yObj;
		po.x=xObj;
		po.z=surface_sag(0,po.y,po.x,0);
		transform(po.x,po.y,po.z,0,0,1,1,1);

		// y1max,... ��z���W�����̂ƈ�v���邩�ǂ����D
		b= (y1max.z==po.z) || (y1min.z==po.z) || (x1max.z==po.z) || (x1min.z==po.z);
		
		if(s==0 && t==0){
			return 0;
		}
		else if(fabs(s)<LN && fabs(t)<LN){
			if(b || is_mask){
				ymax.x=-EPx;           ymax.y=po.y+t*TAN-EPy; 
				ymin.x=-EPx;           ymin.y=po.y-t*TAN-EPy;
				xmax.x=po.x+t*TAN-EPx; xmax.y=-EPy;
				xmin.x=po.x-t*TAN-EPx; xmin.y=-EPy;
			}
			else{
				v=y1max+(t-y1max.z)/(y1max-po).z*(y1max-po); // y1max��po�ɂ�钼�������ƌ����ʒu
				ymax.x=v.x-EPx; ymax.y=v.y-EPy;
				v=y1min+(t-y1min.z)/(y1min-po).z*(y1min-po);
				ymin.x=v.x-EPx; ymin.y=v.y-EPy;
				v=x1max+(t-x1max.z)/(x1max-po).z*(x1max-po);
				xmax.x=v.x-EPx; xmax.y=v.y-EPy;
				v=x1min+(t-x1min.z)/(x1min-po).z*(x1min-po);
				xmin.x=v.x-EPx; xmin.y=v.y-EPy;
			}
		}
		else if(fabs(s)<LN && fabs(t)>=LN){
			if(b || is_mask){
				ymax.x=-EPx;     ymax.y= TAN-EPy;
				ymin.x=-EPx;     ymin.y=-TAN-EPy;
				xmax.x= TAN-EPx; xmax.y=-EPy;
				ymin.x=-TAN-EPx; xmin.y=-EPy;
			}
			else{
				v=(y1max-po)/(y1max-po).z;
				ymax.x=v.x-EPx; ymax.y=v.y-EPy;
				v=(y1min-po)/(y1min-po).z;
				ymin.x=v.x-EPx; ymin.y=v.y-EPy;
				v=(x1max-po)/(x1max-po).z;
				xmax.x=v.x-EPx; xmax.y=v.y-EPy;
				v=(x1min-po)/(x1min-po).z;
				xmin.x=v.x-EPx; xmin.y=v.y-EPy;
			}
		}
		else if(fabs(s)>=LN && fabs(t)<LN){
			if(is_mask){
				ymax.x=-EPx;                            ymax.y=y1max.y*10+(t-y1max.z)*po.y-EPy; 
				ymin.x=-EPx;                            ymin.y=y1min.y*10+(t-y1min.z)*po.y-EPy;
				xmax.x=x1max.x*10+(t-x1max.z)*po.x-EPx; xmax.y=-EPy;
				xmin.x=x1min.x*10+(t-x1min.z)*po.x-EPx; xmin.y=-EPy;			
			}
			else{
				ymax.x=-EPx;                         ymax.y=y1max.y+(t-y1max.z)*po.y-EPy; 
				ymin.x=-EPx;                         ymin.y=y1min.y+(t-y1min.z)*po.y-EPy;
				xmax.x=x1max.x+(t-x1max.z)*po.x-EPx; xmax.y=-EPy;
				xmin.x=x1min.x+(t-x1min.z)*po.x-EPx; xmin.y=-EPy;
			}
		}
		else return 0;

		switch(priority) {
		case 1:
			p=(ymax+ymin)/2;
			for(i=1; i<=10000; ++i) {			
				if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ){
					yPupil=p.y; xPupil=p.x;
					cache.AddTail(p);
					return 1;
				}
				p=ymin+(ymax-ymin)*Random(0,1,0);
			}
			break;
		case 2:
			p=(xmax+xmin)/2;
			for(i=1; i<=10000; ++i) {			
				if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ){
					yPupil=p.y; xPupil=p.x;
					cache.AddTail(p);
					return 1;
				}
				p=xmin+(xmax-xmin)*Random(0,1,0);
			}
			break;
		default:
			p=(ymax+ymin+xmax+xmin)/4;
			start=clock();
			for(i=1; i<=10000000; ++i){
				// 06.12.04 ���s�񐔂�10,000,000�ɐݒ�D
				// ����ȉ��ł͒ʉߌ���������̂ɔ����ł��Ȃ��ꍇ�����������߁D
				// ���������̏ꍇ�͎��Ԃ����Ȃ肩�������D
				// 
				// 080404 �ʐ����������w�n�ł͎��Ԃ������肷����̂Ő�������TIMEOUT=10sec��ݒ�D
				// 090723 �ꗥ10�b�ł́C�ʐ������Ȃ��Ƃ��ʉߌ������Ȃ����f������̂Ɏ��Ԃ�
				//        �����肷���邽�߁CTIMEOUT��p�~���Cthis->FindAThruRayTimeOut
				//        �ɂ�蒲���ł���悤�ɂ����D
				if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ) {
					yPupil=p.y; xPupil=p.x;
					cache.AddTail(p);
					return 1;
				}

				p.x=(xmin+(xmax-xmin)*Random(0,1,0)).x;
				p.y=(ymin+(ymax-ymin)*Random(0,1,0)).y;
				// ��: p=(xmin+(xmax-xmin)*Random(0,1,0))+(ymin+(ymax-ymin)*Random(0,1,0)) �ł͂Ȃ�

				now=clock();
				if( (double)(now-start)/CLOCKS_PER_SEC > this->FindAThruRayTimeOut ) return 0;
			}
			break;
		}
	}

	return 0;
}

int cLens1::FindMarginalRay(double& yPupil,double& xPupil,double azimuth,
                           double yObj,double xObj,double yPupilIni,double xPupilIni,int j,int makecoordinate){
	const double K=1.5, TAN=10;
	double ystep,xstep, ymax,ymin,xmax,xmin, y1max,y1min,x1max,x1min, sag, yp,xp, sn,cs,ra;
	point p1,p2,p;
	vector<complex> dummy;

	if(makecoordinate) make_coordinate(0);
	
	if(EPD>0){
		// ��1�ʂ̌a�������ɑ΂��đ傫���C���C
		// �Ӑ}������˓��̊O���ɂ����������߂��铇�i�����ɂ��j������Ƃ��C
		// else�ȍ~�̂悤�ɑ�1�ʂ̌a����ɂ����̂ł́C���̓��ɂ��C
		// �Ӑ}������˓������߂��Ȃ����Ƃ�����D
		// EPD���K�؂Ȓl�ł���ۏ؂͂Ȃ����C
		// EPD>0�̂Ƃ��͑�1�ʂ̌a�ł͂Ȃ��CEPD���ymax,ymin,...�����߂�D 2013.08.23
		ymax= EPD/2; xmax= EPDx==0 ?  EPD/2: EPDx/2;
		ymin=-EPD/2; xmin= EPDx==0 ? -EPD/2:-EPDx/2;
	}
	else{
		// ��1�ʂ������Ռ��̏ꍇ���Ƃ肠�����Ռ����̌a���y1max, ... �����߂�D
		// y1max, ... �̒��ł����T���Ȃ��̂ł͂Ȃ��̂ł���ł悢�D
		// �����Օ��͖�������

		if(EAy(1)!=0 && EAx(1)==0){
			y1max=dy(1)+EAdy(1)+EAy(1)*K/2; y1min=dy(1)+EAdy(1)-EAy(1)*K/2;
			x1max=dx(1)+EAdx(1)+EAy(1)*K/2; x1min=dx(1)+EAdx(1)-EAy(1)*K/2;
		}
		else if(EAy(1)!=0 && EAx(1)!=0){
			y1max=dy(1)+EAdy(1)+EAy(1)*K/2; y1min=dy(1)+EAdy(1)-EAy(1)*K/2;
			x1max=dx(1)+EAdx(1)+EAx(1)*K/2; x1min=dx(1)+EAdx(1)-EAx(1)*K/2;
		}
		else return 0;
		
		if(s==0 && t==0) return 0;
		else if(fabs(s)<LN && fabs(t)<LN){
			sag=surface_sag(0,yObj,xObj,0);
			if(s+sag!=0){
				ymax=yObj+(y1max-yObj)*(t-s-sag)/(-s-sag)-EPy; 
				ymin=yObj+(y1min-yObj)*(t-s-sag)/(-s-sag)-EPy;
				xmax=xObj+(x1max-xObj)*(t-s-sag)/(-s-sag)-EPx; 
				xmin=xObj+(x1min-xObj)*(t-s-sag)/(-s-sag)-EPx;
			}
			else{
				ymax=yObj+t*TAN-EPy; 
				ymin=yObj-t*TAN-EPy;
				xmax=xObj+t*TAN-EPx; 
				xmin=xObj-t*TAN-EPx;
			}
		}
		else if(fabs(s)<LN && fabs(t)>=LN){
			sag=surface_sag(0,yObj,xObj,0);
			if(s+sag!=0){
				ymax=(y1max-yObj)/(-s-sag)-EPy; 
				ymin=(y1min-yObj)/(-s-sag)-EPy;
				xmax=(x1max-xObj)/(-s-sag)-EPx; 
				xmin=(x1min-xObj)/(-s-sag)-EPx;
			}
			else{
				ymax=TAN;
				ymin=-TAN;
				xmax=TAN;
				ymin=-TAN;
			}
		}
		else if(fabs(s)>=LN && fabs(t)<LN){
			ymax=y1max+t*yObj-EPy; 
			ymin=y1min+t*yObj-EPy;
			xmax=x1max+t*xObj-EPx; 
			xmin=x1min+t*xObj-EPx;
		}
		else return 0;
	}

	// 081210�폜 -> if( RayTrace(yObj,xObj,yPupilIni,xPupilIni,0,j,1,0,0,0)!=0 ) return 0;
	//   ���˓��̓����Ɏ����ɂ��Õ��̓������邱�Ƃ����邽�߁C
	//   (���˓����W�ɑ΂��邠�镔���̌������̕ω����C�L���a�O�ŋɒl�����Ȃ�)
	//   (yPupilIni,xPupilIni)�̌��������ۂɌn��ʉ߂ł��邱�Ƃ͗v�����Ȃ��D
	//     (ex) FindPupilEdge()�ł́Cp1=(pymax+pymin)/2 �����ۂɌn��ʉ߂ł��Ȃ������ƂȂ邱�Ƃ�����C
	//          ���̂Ƃ��������[�v���N����D

	// �������ʂ�_p1
	p1=point(xPupilIni,yPupilIni);
	
	// �������ʂ�Ȃ��_p2���m��
	yp=yPupilIni; xp=xPupilIni;
	ystep=fabs(ymax-ymin)/4;
	xstep=fabs(xmax-xmin)/4;

	// x,y���a��xstep,ysetp�̑ȉ~��azimuth�����̔��a��ra�ł���D
	//   (ra cos��)^2/xstep^2 + (ra sin��)^2/ystep^2 =1 ���D
	if(azimuth==0 || azimuth==360){
		sn=0; cs=1;
	}
	else if(azimuth==90){
		sn=1; cs=0;
	}
	else if(azimuth==180){
		sn=0; cs=-1;
	}
	else if(azimuth==270){
		sn=-1; cs=0;
	}
	else{
		sn=sin(azimuth*PI/180); cs=cos(azimuth*PI/180);
	}

	ra=sqrt(1/(cs*cs/xstep/xstep+sn*sn/ystep/ystep));
	ystep=ra*sn;
	xstep=ra*cs;

	do{
		yp+=ystep;
		xp+=xstep;
	} while( RayTrace(yObj,xObj,yp,xp,0,j,1,0,0,0,0,dummy,0)==0 );
	p2=point(xp,yp);

	// �����I�Ɏ������������߂�
	do{
		const double M=0.95;
		// p1 : �����ʉ߉\
		// p2 : �����ʉߕs�\
		//   �����ɂ����˓������ɈÕ��̓������邱�Ƃ����邽�߁C
		//   p��P��p1,p2�̒��_�ł͂Ȃ��Cp2���ɂ��āC
		//   p2�����˓��O���璼�ړ��˓��Õ��Ɉڂ邱�Ƃ��ɗ͖h���D(081210����)
		//
		// M=0.9��0.95�ɕύX. ��L�G���[�͌��������C
		// �Ⴆ�΃X���[�t�H�[�J�XMTF�̑��x�͏����x���Ȃ�D(100115)
		p=(1.0-M)*p1+M*p2;   // |p1-p| : |p2-p| = M : 1-M
		if( RayTrace(yObj,xObj,p.y,p.x,0,j,1,0,0,0,0,dummy,0)==0 ) p1=p; else p2=p;
	} while( distance(p1,p2)>0.000001 );

	yPupil=p1.y; xPupil=p1.x;
	// yPupil=p.y; xPupil=p.x; �Ƃ���Ƌ��ʂƌ����Ȃ��Ȃǂ�(yPupil,xPupil)�̌������ʂ�Ȃ����Ƃ�����D
	return RayTrace(yObj,xObj,p2.y,p2.x,0,j,1,0,0,0,0,dummy,0); // �ǂ��̖ʂ���������Ԃ�
}

int cLens1::FindMarginalRay2(double& yPupil,double& xPupil,int direction,
                           double yObj,double xObj,double yPupilIni,double xPupilIni,int j,int makecoordinate){
	double azimuth;

	switch(direction){
		case 1: azimuth= 90; break;
		case 2: azimuth=270; break;
		case 3: azimuth=  0; break;
		case 4: azimuth=180; break;
		default: return 0;
	}
	return FindMarginalRay(yPupil,xPupil,azimuth,yObj,xObj,yPupilIni,xPupilIni,j,makecoordinate);
}

int cLens1::FindPupilEdge(point& ymax,point& ymin, point& xmax,point& xmin, 
                     int& ymax_i,int& ymin_i, int& xmax_i,int& xmin_i,
					 double yObj,double xObj,int j,int findpupil,int makecoordinate) {
	point p1,p2, pymax,pymin, pxmax,pxmin;
	int found;

	// FindMarginalRay()�Œ����Օ��͖�������邽�߁C�{�֐��ł����������D

	if(findpupil && ea_is_defined()){
		if(makecoordinate) make_coordinate(0);
		found=0;

		if(yObj==0 && xObj==0){
			found=FindAThruRay(p2.y,p2.x,yObj,xObj,j,1,0);
			if(found==0){
				found=FindAThruRay(p2.y,p2.x,yObj,xObj,j,2,0);
			}
		}
		else if(yObj!=0 && xObj==0){
			found=FindAThruRay(p2.y,p2.x,yObj,xObj,j,1,0);
		}
		else if(yObj==0 && xObj!=0){
			found=FindAThruRay(p2.y,p2.x,yObj,xObj,j,2,0);
		}
		
		if(found==0){
			if( FindAThruRay(p2.y,p2.x,yObj,xObj,j,3,0)==0 ) return 0;
		}

		do {
			p1=p2;
			ymax_i=FindMarginalRay2(pymax.y,pymax.x,1,yObj,xObj,p1.y,p1.x,j,0);
			ymin_i=FindMarginalRay2(pymin.y,pymin.x,2,yObj,xObj,p1.y,p1.x,j,0);
			p1=(pymax+pymin)/2;
			xmax_i=FindMarginalRay2(pxmax.y,pxmax.x,3,yObj,xObj,p1.y,p1.x,j,0);
			xmin_i=FindMarginalRay2(pxmin.y,pxmin.x,4,yObj,xObj,p1.y,p1.x,j,0);
			p2=(pxmax+pxmin)/2;
		} while( distance(p1,p2)>0.0001 );
		ymax=pymax; ymin=pymin; xmax=pxmax; xmin=pxmin;
		return 1;
	}
	else{
		ymax=point(0,EPD/2); ymin=point(0,-EPD/2);
		xmax=point((EPDx==0 ? EPD:EPDx)/2,0); xmin=point((EPDx==0 ? -EPD:-EPDx)/2,0);
		ymax_i=ymin_i=xmax_i=xmin_i=0;
		return 1;
	}
}
int cLens1::FindPupilEdge(point& ymax,point& ymin, point& xmax,point& xmin, 
					 double yObj,double xObj,int j,int findpupil,int makecoordinate) {
	int dummy;
	return FindPupilEdge(ymax,ymin,xmax,xmin,dummy,dummy,dummy,dummy,yObj,xObj,j,findpupil,makecoordinate);
}
int cLens1::FindPupilEdge(point& ymax,point& ymin, point& xmax,point& xmin,
						  double yObj,double xObj,int j,int makecoordinate){
	return FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,j,1,makecoordinate);
}

int cLens1::FindRay(point &p,double yObj,double xObj,std::string kind,int findpupil,int j){
	point ymax,ymin,xmax,xmin;
	int ymax_i,ymin_i,xmax_i,xmin_i;
	// findpupil=true�ŁCkind�� ymax,ymin,xmax,xmin�̂Ƃ��C�����𐧌�����ʂ��������̐�����Ԃ��D
	// kind�� principal,ymax,ymin,xmax,xmin�̂ǂꂩ�ŁC�����̐ݒ�Ɏ��s�����Ƃ�0��Ԃ��D
	// ���̑���-1��Ԃ��D����kind==""�̂Ƃ���yPupil,xPupil��ς��Ȃ���-1��Ԃ��D

	if(kind==""){
		return -1;
	}
	if(kind=="principal"){ 
		if(FindPupilEdge(ymax,ymin,xmax,xmin,ymax_i,ymin_i,xmax_i,xmin_i,yObj,xObj,j,findpupil,1)){
			p=(ymax+ymin)/2;
			return -1;
		}
		else return 0;
	}
	if(kind=="ymax"){ 
		if(FindPupilEdge(ymax,ymin,xmax,xmin,ymax_i,ymin_i,xmax_i,xmin_i,yObj,xObj,j,findpupil,1)){
			p=ymax;
			return findpupil!=0 ? ymax_i : -1;
		}
		else return 0;
	}
	if(kind=="ymin"){ 
		if(FindPupilEdge(ymax,ymin,xmax,xmin,ymax_i,ymin_i,xmax_i,xmin_i,yObj,xObj,j,findpupil,1)){
			p=ymin;
			return findpupil!=0 ? ymin_i : -1;
		}
		else return 0;
	}
	if(kind=="xmax"){ 
		if(FindPupilEdge(ymax,ymin,xmax,xmin,ymax_i,ymin_i,xmax_i,xmin_i,yObj,xObj,j,findpupil,1)){
			p=xmax;
			return findpupil!=0 ? xmax_i : -1;
		}
		else return 0;
	}
	if(kind=="xmin"){ 
		if(FindPupilEdge(ymax,ymin,xmax,xmin,ymax_i,ymin_i,xmax_i,xmin_i,yObj,xObj,j,findpupil,1)){
			p=xmin;
			return findpupil!=0 ? xmin_i : -1;
		}
		else return 0;
	}
	return -1;
}

double cLens1::OptimizedDefocus(double yObj,double xObj,int FindPupil,int j){
	// �g�ʎ���rms���ŏ��ɂ��邽�߂�Zernike�W�J���ʍ�(n=2,l=0)�̒����ʂ��f�t�H�[�J�X�ʂ��v�Z����D
	// ��j�g���̂ݍl���D
	// ���F���l���������̂ɂ���ɂ͂ǂ�������悢�����܂̂Ƃ��땪����Ȃ�(07.10.24).
	matrix<double> X,Y,W,P,A;
	matrix<int> E;
	cZernike zernike;
	double dW,dW0, fno_y,fno_x,fno, defocus, costh;

	zernike.SetNormalize(0); // �P�ʉ~����Ł}1�ɂȂ鐳�K��
	MakeOPDMap(X,Y,W,P,A,E, zernike, yObj,xObj,0,j,j,FindPupil,0,3,0,0,0,0,"");
	dW =zernike.GetC( zernike.jNumber(0,2) )*2/1000; // ���ʊ֐�(n=2,l=0)�̍��፷��2�Ȃ̂�2�{����
	MakeOPDMap(X,Y,W,P,A,E, zernike, yObj,xObj,0,j,j,FindPupil,0,3,1,0,0,0,"");
	dW0=zernike.GetC( zernike.jNumber(0,2) )*2/1000;
	fno_for_unit_pupil(fno_y,fno_x,yObj,xObj,0,j,1);
	fno=fno_y/(zernike.GetR0()*2);

	if(Afocal){
		// "�����Y�݌v�̂��߂̔g�ʌ��w(����)" (3.14)���
		// dW<0(������)�̂Ƃ��K�v��defocus�𐳂Ƃ���D
		defocus=-1000*8*(dW-dW0)*fno*fno*fabs(N(k,j));
	}
	else{
		// �ȑO�� defocus=8*(dW-dW0)*fno*fno*N(k,j); �Ƃ��Ă����D
		// �������C�g�􉽌��w(�O��a�v)�h��(7.7)���Ɣ�ׂ��
		//   1-cos��=��^2/2  (�Ƃ͌����̊J���p)
		// �̋ߎ����g���Ă��邱�Ƃ��킩��D
		// �����ł��̋ߎ����g��Ȃ��悤�ɕύX�����D
		// �ŏ��̑��ʂ��œK���ʂ���傫������Ă���Ɛ��x���o�Ȃ��ꍇ�����邪�C
		// ���̕ύX�ŉ��P����邱�Ƃ��������D
		// �܂��C�ߎ����̂܂܂� |dW-dW0| ���������Ȃ�܂Ń��[�v���s�Ȃ�
		// ���Ƃ������������������Ԃ�������̂ł�߂��D
		// (091023)
		costh=sqrt(1-1.0/4.0/N(k,j)/N(k,j)/fno/fno);
		defocus=(dW-dW0)/N(k,j)/(1-costh);   // "�􉽌��w(�O��)" (7.7)
	}

	return defocus;
}

double cLens1::OptimizedDefocus(int FindPupil){
	// ���㕨�_�C��1�g���ɑ΂���œK�f�t�H�[�J�X��
	return OptimizedDefocus(0,0,FindPupil,1);
}

double cLens1::OptimizeS1fix(double yObj,double xObj,int FindPupil,int j){
	return this->s1fix= dk(0)+OptimizedDefocus(yObj,xObj,FindPupil,j);
}

double cLens1::OptimizeS1fix(int FindPupil){
	return OptimizeS1fix(0,0,FindPupil,1);
}

std::string cLens1::SingleRayTrace(double yObj,double xObj,
								   std::string SetRay,int findpupil,double yPupil,double xPupil,
								   double defocus,int lastsurf,int j,
								   double pol_phi/*=90*/,double pol_ratio/*=0*/){
	// pol_phi   : �Ό���A��X���ƂȂ��p�x
	// pol_ratio : �Ό���A�����̋��xIa�ƕΌ���B�����̋��xIb�̔� Ib/Ia,  Ib<Ia�ł���K�v�͂Ȃ��D
	// pol_phi, pol_ratio �͕��̖ʂ̍��W�ŗ^����D
	//
	// �Q�l : �~�Ό�(pro_ratio==1)��pol_phi�͖{���s��ł��邽�߁C�s����Ȍ��ʂ��o�͂����D
	// �Q�l : ���̎p��(��]�Ȃ�)�����邽�߂� ha[i].y,ha[i].x ���o�͂���悤�ɂ����D 2019.04.02
	int i;
	char buf[1000];
	const std::string NO="No ray go through\n";
	std::string s;
	point pupil,image;
	int limit_surf, errorcode;
	double sph,cyl,axis, mx,my;
	vector<complex> E0; double Ax,Ay,Phi;

	cPolarization::Analyze(Ax,Ay,Phi,1,sqrt(pol_ratio),pol_phi);
	E0.x=Ax;
	E0.y=complex(Ay*cos(Phi*PI/180),Ay*sin(Phi*PI/180));
	E0.z=0;

	pupil=point(xPupil,yPupil);
	limit_surf=FindRay(pupil,yObj,xObj,SetRay,findpupil,j);
	if( limit_surf==0 ) return NO;
	if( errorcode=RayTrace(yObj,xObj,pupil.y,pupil.x,defocus,j,1,1,0,1,1,E0) ) {
		sprintf(buf,"Can't go through : surf=%d(%s)\n", NotThruSurf(errorcode),ErrorMessage(errorcode).c_str());
		s+=buf;
		lastsurf=NotThruSurf(errorcode)-1;
	}

	// 1<=lastsurf<=k�̂Ƃ���lastsurf�ʏo�ˌ�ɂ���sph,cyl�����v�Z����C����ȊO�ł͑��ʂŌv�Z����D
	if(lastsurf<1 || k<lastsurf) lastsurf=0;
	SCA(sph,cyl,axis,yObj,xObj,pupil.y,pupil.x,defocus,lastsurf==0 ? k+1:lastsurf,j);
	mx=Mx(yObj,xObj,pupil.y,pupil.x,1,lastsurf==0 ? k:lastsurf,j);
	my=My(yObj,xObj,pupil.y,pupil.x,1,lastsurf==0 ? k:lastsurf,j);

    //           ###### #####.##### #####.##### ##.##### ##.##### ####.### ####.### ###.### ###.###   ###.# #.###  ##.#### ##.####
	sprintf(buf,"          y           x          Y1       X1      i[deg]   i'[deg]  dW/dz  dz(��m)/F  pol�� Ib/Ia  ha.y    ha.x\n");
	s+=buf;

	sprintf(buf,"object %11.5f %11.5f %8.5f %8.5f          %8.3f                   %5.1f %5.3f  %7.4f %7.4f\n",
	        yObj,xObj,Y1[0],X1[0],ArcCos(fabs(xi1[0]))*180/PI,EllipsePhi(0),EllipseIntensityRatio(0),ha[0].y,ha[0].x);
	s+=buf;
	sprintf(buf," pupil %11.5f %11.5f\n",pupil.y,pupil.x);
	s+=buf;

	for(i=1; i<=(lastsurf==0 ? k:lastsurf); i++){
		sprintf(buf,"%6d %11.5f %11.5f %8.5f %8.5f %8.3f %8.3f %7.3f %7.3f   %5.1f %5.3f  %7.4f %7.4f  %s\n", 
		        i,y[i],x[i],Y1[i],X1[i],ArcCos(fabs(xi[i]))*180/PI,ArcCos(fabs(xi1[i]))*180/PI,
		        DWDSag(i,yObj,xObj,SetRay,yPupil,xPupil),
				DSagDFringe(i,yObj,xObj,SetRay,yPupil,xPupil),
				EllipsePhi(i),EllipseIntensityRatio(i),  // ��i�ʂ̓��ˑ��̕Ό����
				ha[i].y,ha[i].x, rem(i).c_str());
		s+=buf;
	}

	if(lastsurf==0){
		ImageHeight(image,point(xObj,yObj),pupil,defocus,j,0,1,0);
		sprintf(buf," image %11.5f %11.5f                   %8.3f                            %5.1f %5.3f  %7.4f %7.4f\n",
		        image.y,image.x,ArcCos(fabs(xi[k+1]))*180/PI,EllipsePhi(k+1),EllipseIntensityRatio(k+1),ha[k+1].y,ha[k+1].x);
		s+=buf;
	}

	sprintf(buf,"sph=%.6f cyl=%.6f (s+c=%.6f ave=%.6f) axis=%.6f\n", sph,cyl,sph+cyl,sph+cyl/2,axis);
	s+=buf;
	sprintf(buf,"Mx=%.6f My=%.6f\n", mx,my);
	s+=buf;
	if(limit_surf>0){
		sprintf(buf,"%s limit surf=%d(%s)\n", SetRay.c_str(),NotThruSurf(limit_surf),ErrorMessage(limit_surf).c_str());
		s+=buf;
	}

	if(lastsurf==0){
		sprintf(buf,"image intensity=%g\n", sqabs(E[k+1].x)+sqabs(E[k+1].y)+sqabs(E[k+1].z) );
		s+=buf;
	}
	return s;
}

std::string cLens1::UsedRange(char y_or_x,double HObjMax,double HObjMin){
	// ���̍�HObjMin����HObjMax�ɑΉ�����e�ʂ̌����ʉߔ͈͂����߂�i���J�݌v�p�Ȃǁj�D
	int i;
	std::string s;
	char buf[1000];
	point p;
	std::string setraymax,setraymin;
	double yobjmax,yobjmin,xobjmax,xobjmin;
	ray_data ray1(k+1),ray2(k+1),ray3(k+1),ray4(k+1);
	list<double> hmax,hmin;

	if(y_or_x=='y'){
		setraymax="ymax";
		setraymin="ymin";
		yobjmax=HObjMax;
		yobjmin=HObjMin;
		xobjmax=xobjmin=0;
	}
	else if(y_or_x=='x'){
		setraymax="xmax";
		setraymin="xmin";
		yobjmax=yobjmin=0;
		xobjmax=HObjMax;
		xobjmin=HObjMin;
	}
	else{
		return "argument error\n";
	}

	if(FindRay(p,yobjmax,xobjmax,setraymax,1,1)==0) return "can't find ray\n";
	RayTrace(yobjmax,xobjmax,p.y,p.x,0,1,0,0,0,0);
	for(i=0; i<=k+1; ++i){
		ray1.x[i]=this->x[i]; ray1.y[i]=this->y[i];
	}

	if(FindRay(p,yobjmin,xobjmin,setraymax,1,1)==0) return "can't find ray\n";
	RayTrace(yobjmin,xobjmin,p.y,p.x,0,1,0,0,0,0);
	for(i=0; i<=k+1; ++i){
		ray2.x[i]=this->x[i]; ray2.y[i]=this->y[i];
	}

	if(FindRay(p,yobjmax,xobjmax,setraymin,1,1)==0) return "can't find ray\n";
	RayTrace(yobjmax,xobjmax,p.y,p.x,0,1,0,0,0,0);
	for(i=0; i<=k+1; ++i){
		ray3.x[i]=this->x[i]; ray3.y[i]=this->y[i];
	}

	if(FindRay(p,yobjmin,xobjmin,setraymin,1,1)==0) return "can't find ray\n";
	RayTrace(yobjmin,xobjmin,p.y,p.x,0,1,0,0,0,0);
	for(i=0; i<=k+1; ++i){
		ray4.x[i]=this->x[i]; ray4.y[i]=this->y[i];
	}

	for(i=0; i<=k+1; ++i){
		if(y_or_x=='y'){
			hmax.AddTail( Max(ray1.y[i],ray2.y[i],ray3.y[i],ray4.y[i]) );
			hmin.AddTail( Min(ray1.y[i],ray2.y[i],ray3.y[i],ray4.y[i]) );
		}
		else if(y_or_x=='x'){
			hmax.AddTail( Max(ray1.x[i],ray2.x[i],ray3.x[i],ray4.x[i]) );
			hmin.AddTail( Min(ray1.x[i],ray2.x[i],ray3.x[i],ray4.x[i]) );
		}
	}

	if(y_or_x=='y'){
		//          "    ## #####.### #####.###"
		sprintf(buf,"          ymax      ymin\n"); s+=buf;
	}
	else if(y_or_x=='x'){
		sprintf(buf,"          xmax      xmin\n"); s+=buf;
	}

	sprintf(buf,"object %9.3f %9.3f\n",hmax[1],hmin[1]); s+=buf;
	for(i=1; i<=k; i++){
		sprintf(buf,"    %2d %9.3f %9.3f  %s\n",i,hmax[i+1],hmin[i+1],rem(i).c_str()); s+=buf;
	}
	sprintf(buf,"image  %9.3f %9.3f\n",hmax[k+2],hmin[k+2]); s+=buf;

	return s;
}

std::string cLens1::yUsedRange(double yObjMax,double yObjMin){
	return UsedRange('y',yObjMax,yObjMin);
}

std::string cLens1::xUsedRange(double xObjMax,double xObjMin){
	return UsedRange('x',xObjMax,xObjMin);
}

vector<double> cLens1::RayPos(int i,double yObj,double xObj,
                              std::string SetRay,double yPupil,double xPupil,int j){
	point pupil;
	pupil=point(xPupil,yPupil);
	int errcode;

	// SetRay�� "principal","ymax", ... �ƈ�v���C�����̐ݒ�Ɏ��s�����Ƃ���0��Ԃ��D
	if(FindRay(pupil,yObj,xObj,SetRay,1,j)==0) return vector<double>(0,0,0);
	errcode=RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,0,0,0);  // EA�͖���
	if(errcode==0 || NotThruSurf(errcode)>i){
		if( 1<=i && i<=k+1 ){
			return vector<double>(x[i],y[i],z[i]);
		}
		else{
			return vector<double>(0,0,0);
		}
	}
	else{
		return vector<double>(0,0,0);
	}
}

int cLens1::RayPos(double* x,double* y,double* z,
				   int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j){
	vector<double> v;

	v=RayPos(i,yObj,xObj,SetRay,yPupil,xPupil,j);
	*x=v.x; *y=v.y; *z=v.z;
	return v==vector<double>(0,0,0) ? 0 : 1;
}

double cLens1::RayPosX(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j){
	return RayPos(i,yObj,xObj,SetRay,yPupil,xPupil,j).x;
}
double cLens1::RayPosY(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j){
	return RayPos(i,yObj,xObj,SetRay,yPupil,xPupil,j).y;
}
double cLens1::RayPosZ(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j){
	return RayPos(i,yObj,xObj,SetRay,yPupil,xPupil,j).z;
}

point cLens1::RayHeight(int i,double yObj,double xObj,
                       std::string SetRay,double yPupil,double xPupil,int j){
	vector<double> v;
	v=RayPos(i,yObj,xObj,SetRay,yPupil,xPupil,j);
	return point(v.x,v.y);
}

void cLens1::Footprint(double &ymax,double &ymin,double &xmax,double &xmin,
                     double &ydia,double &xdia,double &ycenter,double &xcenter,
                     double yObj,double xObj,int i,int j)
{
	point pymax,pymin,pxmax,pxmin;
	pymax=RayHeight(i,yObj,xObj,"ymax",0,0,j);
	pymin=RayHeight(i,yObj,xObj,"ymin",0,0,j);
	pxmax=RayHeight(i,yObj,xObj,"xmax",0,0,j);
	pxmin=RayHeight(i,yObj,xObj,"xmin",0,0,j);
	ymax=pymax.y;
	ymin=pymin.y;
	xmax=pxmax.x;
	xmin=pxmin.x;
	ydia=fabs(ymax-ymin);
	xdia=fabs(xmax-xmin);
	ycenter=(ymax+ymin)/2;
	xcenter=(xmax+xmin)/2;
}
double cLens1::FootprintYdia(double yObj,double xObj,int i,int j){
	double ymax,ymin,xmax,xmin,ydia,xdia,ycenter,xcenter;
	Footprint(ymax,ymin,xmax,xmin,ydia,xdia,ycenter,xcenter,yObj,xObj,i,j);
	return ydia;
}
double cLens1::FootprintXdia(double yObj,double xObj,int i,int j){
	double ymax,ymin,xmax,xmin,ydia,xdia,ycenter,xcenter;
	Footprint(ymax,ymin,xmax,xmin,ydia,xdia,ycenter,xcenter,yObj,xObj,i,j);
	return xdia;
}

double cLens1::IncidentAngle(int i,double yObj,double xObj,
                            std::string SetRay,double yPupil,double xPupil,int j,int InAir/*=0*/){
	// InAir���^�̂Ƃ��͋�C���Z�ł̓��ˊp��o��Ԃ��D sin��o = Nsin��  �� ��o = asin(Nsin��)
	int errcode;
	point pupil(xPupil,yPupil);
	double th;

	FindRay(pupil,yObj,xObj,SetRay,1,j);
	if( errcode=RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,0,0,0) ){  // �������ʂ�Ȃ��Ƃ�
		if( NotThruSurf(errcode)<=i ) return 0;                     // ���Ci�ʈȑO�ł�����Ƃ�
	}
	if(1<=i && i<=k+1){
		th=ArcCos(fabs(xi[i]));   // ArcCos�֐��͐��l�v�Z�덷�ɂ��1���킸���ɒ����������ɑΉ����Ă���
		if(InAir){
			return asin( fabs(N(i-1,1))*sin(th) )*180/PI;  // ��C���Z�p�x
		}
		else{
			return th*180/PI;  // ���p�x
		}
	}
	else return 0;
}
double cLens1::ExitAngle(int i,double yObj,double xObj,
                        std::string SetRay,double yPupil,double xPupil,int j,int InAir/*=0*/){
	// InAir���^�̂Ƃ��͋�C���Z�ł̏o�ˊp��o��Ԃ��D sin��o = Nsin�� �� ��o = asin(Nsin��)
	int errcode;
	double th;
	point pupil(xPupil,yPupil);

	FindRay(pupil,yObj,xObj,SetRay,1,j);
	if( errcode=RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,0,0,0) ){  // �������ʂ�Ȃ��Ƃ�
		if( NotThruSurf(errcode)<=i ) return 0;                     // ���Ci�ʈȑO�ł�����Ƃ�
	}
	if( 0<=i && i<=k ){
		th=ArcCos(fabs(xi1[i]));    // ArcCos�֐��͐��l�v�Z�덷�ɂ��1���킸���ɒ����������ɑΉ����Ă���D
		if(InAir){
			return asin( fabs(N(i,1))*sin(th) )*180/PI;  // ��C���Z�p�x
		}
		else{
			return th*180/PI; // ���p�x
		}
	}
	else return 0;
}

double cLens1::IncidentAngleMax(double yObj,double xObj,
                                std::string SetRay,double yPupil,double xPupil,int j){
	// �e�ʂւ̓��ˊp(��C���Z)�̍ő�����߂�D
	// SR-1�̍L�p�A�_�v�^�݌v�p�ɍ쐬  2010.10.05
	// (���������̍ő���ˊp���������قǔ��˂𓦂����X���p�x���������Ă��ނƎv����)
	int i;
	double th,max;

	max=0;
	for(i=1; i<=k; i++){
		th=IncidentAngle(i,yObj,xObj,SetRay,yPupil,xPupil,j,1);
		if(fabs(th)>max) max=fabs(th);
	}
	return max;
}

double cLens1::IncidentAngleYFanMin(int i,double yObj,double xObj,int FindPupil,int j){
	// ��i�ʂŁCy-fan�̗��[�̓��ˌ���������Q1,Q2,�ʖ@��������E1,E2�Ƃ���Ƃ��C
	//   (Q1xE1,Q2xE2)>0�̂Ƃ��C��Βl�̏��������̓��ˊp(��C���Z)�𐳂̒l�ŕԂ��D
	//   (Q1xE1,Q2xE2)<0�̂Ƃ��C��Βl�̏��������̓��ˊp(��C���Z)�𕉂̒l�ŕԂ��i0�ł̘A�����j�D
	// ��F
	// SR-1�̍L�p�A�_�v�^(�����Y���X��������)�̊e�����Y�ʂł��̊֐��̒l�����ł���Ȃ�ׂ��傫�����Ƃ�
	// �摜�ɓ������˃m�C�Y���o�Ȃ����Ƃ�\���D
	// (Q1xE1,Q2xE2)�����Ƃ������Ƃ́Cy-fan�̒��ɓ��ˊp��0�ƂȂ����������C���̌����̔��˂�
	// �Ɩ������̌o�H�Ŏ���n�֖߂�D
	int errcode;
	point p1,p2;
	vector<double> Q1,E1,Q2,E2;
	double th1,th2;

	if( i<1 || k<i ) return 0;

	FindRay(p1,yObj,xObj,"ymax",FindPupil,j);
	if( errcode=RayTrace(yObj,xObj,p1.y,p1.x,0,j,0,0,0,0) ){  // �������ʂ�Ȃ��Ƃ�
		if( NotThruSurf(errcode)<=i ) return 0;               // ���Ci�ʈȑO�ł�����Ƃ�
	}
	Q1=vector<double>(X[i],Y[i],Z[i]);
	E1=surface_normal(i,y[i],x[i],1);
	th1=IncidentAngle(i,yObj,xObj,"",p1.y,p1.x,j,1);

	FindRay(p2,yObj,xObj,"ymin",FindPupil,j);
	if( errcode=RayTrace(yObj,xObj,p2.y,p2.x,0,j,0,0,0,0) ){
		if( NotThruSurf(errcode)<=i ) return 0;
	}
	Q2=vector<double>(X[i],Y[i],Z[i]);
	E2=surface_normal(i,y[i],x[i],1);
	th2=IncidentAngle(i,yObj,xObj,"",p2.y,p2.x,j,1);
	
	return sgn(sProduct(vProduct(Q1,E1),vProduct(Q2,E2))) >0 ? Min(th1,th2) : -Min(th1,th2);
}

double cLens1::IncidentAnglePrinMin(int i,double yObj1,double yObj2,int FindPupil,int j){
	// ��i�ʂŁC����(yObj1,0),(yObj2,0)�ɑ΂��������̕�����Q1,Q2,�ʖ@��������E1,E2�Ƃ���Ƃ��C
	//   (Q1xE1,Q2xE2)>0�̂Ƃ��C��Βl�̏��������̓��ˊp(��C���Z)�𐳂̒l�ŕԂ��D
	//   (Q1xE1,Q2xE2)<0�̂Ƃ��C��Βl�̏��������̓��ˊp(��C���Z)�𕉂̒l�ŕԂ��i0�ł̘A�����j�D
	// ��F
	// SR-1�̍L�p�A�_�v�^(�����Y���X��������)�̊e�����Y�ʂł��̊֐��̒l�����ł���Ȃ�ׂ��傫�����Ƃ�
	// �摜�ɓ������˃m�C�Y���o�Ȃ����Ƃ�\���D
	// (Q1xE1,Q2xE2)�����Ƃ������Ƃ́C(yObj1,0)��(yObj2,0)�̊Ԃ̕��_��
	// ���ˊp��0�ƂȂ�ꍇ������C���̂Ƃ����˂͏Ɩ������̌o�H�Ŏ���n�֖߂邱�ƂɂȂ�D
	int errcode;
	point p;
	vector<double> Q1,E1,Q2,E2;
	double th1,th2;

	if( i<1 || k<i ) return 0;

	FindRay(p,yObj1,0,"principal",FindPupil,j);
	if( errcode=RayTrace(yObj1,0,p.y,p.x,0,j,0,0,0,0) ){ // �������ʂ�Ȃ��Ƃ�
		if( NotThruSurf(errcode)<=i ) return 0;          // ���Ci�ʈȑO�ł�����Ƃ�
	}
	Q1=vector<double>(X[i],Y[i],Z[i]);
	E1=surface_normal(i,y[i],x[i],1);
	th1=IncidentAngle(i,yObj1,0,"",p.y,p.x,j,1);

	FindRay(p,yObj2,0,"FindPupil",FindPupil,j);
	if( errcode=RayTrace(yObj2,0,p.y,p.x,0,j,0,0,0,0) ){
		if( NotThruSurf(errcode)<=i ) return 0;
	}
	Q2=vector<double>(X[i],Y[i],Z[i]);
	E2=surface_normal(i,y[i],x[i],1);
	th2=IncidentAngle(i,yObj2,0,"",p.y,p.x,j,1);
	
	return sgn(sProduct(vProduct(Q1,E1),vProduct(Q2,E2))) >0 ? Min(th1,th2) : -Min(th1,th2);
}

double cLens1::IncidentAngleY(int i,double yObj,double xObj,
							  std::string SetRay,double yPupil,double xPupil,int j,int InAir/*=0*/){
    // IncidentAngle�֐��́C0���œ��֐���0�ł���C�œK���̕]���֐��Ƃ��Ďg���ɂ����D
	// ����C�{�֐��̖߂�l�ɂ͐���������C�܂����˖ʂ�YZ���ʓ��̂Ƃ��͐��m�ȓ��ˊp��Ԃ��D
	int errcode;
	point pupil(xPupil,yPupil);
	double th;
	vector<double> Q,E,ex;

	FindRay(pupil,yObj,xObj,SetRay,1,j);
	if( errcode=RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,0,0,0) ){  // �������ʂ�Ȃ��Ƃ�
		if( NotThruSurf(errcode)<=i ) return 0;                     // ���Ci�ʈȑO�ł�����Ƃ�
	}
	if(1<=i && i<=k+1){
		E=surface_normal(i,y[i],x[i],1);   E=E/abs(E);  // �@�������P�ʃx�N�g��
		Q=vector<double>(X[i],Y[i],Z[i]);  Q=Q/abs(Q);  // ���ˌ��������P�ʃx�N�g��
		ex=vector<double>(1,0,0);                       // x�����P�ʃx�N�g��
		th=asin( sProduct(vProduct(Q,E),ex) );
		if(InAir){
			return asin( fabs(N(i-1,1))*sin(th) )*180/PI;  // ��C���Z�p�x
		}
		else{
		 	return th*180/PI;  // ���p�x
		}
	}
	else return 0;

	return 0;
}

double cLens1::IncidentAngleX(int i,double yObj,double xObj,
							  std::string SetRay,double yPupil,double xPupil,int j,int InAir/*=0*/){
	// IncidenAngleY�֐��Ɠ��l�D���˖ʂ�XZ���ʓ��̂Ƃ��͐��m�ȓ��ˊp��Ԃ��D
	int errcode;
	point pupil(xPupil,yPupil);
	double th;
	vector<double> Q,E,ey;

	FindRay(pupil,yObj,xObj,SetRay,1,j);
	if( errcode=RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,0,0,0) ){  // �������ʂ�Ȃ��Ƃ�
		if( NotThruSurf(errcode)<=i ) return 0;                     // ���Ci�ʈȑO�ł�����Ƃ�
	}
	if(1<=i && i<=k+1){
		E=surface_normal(i,y[i],x[i],1);   E=E/abs(E);  // �@�������P�ʃx�N�g��
		Q=vector<double>(X[i],Y[i],Z[i]);  Q=Q/abs(Q);  // ���ˌ��������P�ʃx�N�g��
		ey=vector<double>(0,1,0);                       // y�����P�ʃx�N�g��
		th=asin( sProduct(vProduct(Q,E),ey) );
		if(InAir){
			return asin( fabs(N(i-1,1))*sin(th) )*180/PI;  // ��C���Z�p�x
		}
		else{
		 	return th*180/PI;  // ���p�x
		}
	}
	else return 0;

	return 0;
}

vector<double> cLens1::IncidentDirectionCosine(int i,double yObj,double xObj,
	                   std::string SetRay,int FindPupil,double yPupil,double xPupil,int j){
	point pupil;
	pupil=point(xPupil,yPupil);
	FindRay(pupil,yObj,xObj,SetRay,FindPupil,j);
	if( RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,0,0,0) ) return vector<double>(0,0,0);
	if( 1<=i && i<=k+1 ){
		return vector<double>(X[i],Y[i],Z[i]);
	}
	else{
		return vector<double>(0,0,0);
	}
}

vector<double> cLens1::ExitDirectionCosine(int i,double yObj,double xObj,
					   std::string SetRay,int FindPupil,double yPupil,double xPupil,int j){
	point pupil;
	pupil=point(xPupil,yPupil);
	FindRay(pupil,yObj,xObj,SetRay,FindPupil,j);
	if( RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,0,0,0) ) return vector<double>(0,0,0);
	if( 1<=i && i<=k+1 ){
		return vector<double>(X1[i],Y1[i],Z1[i]);
	}
	else{
		return vector<double>(0,0,0);
	}
}

double cLens1::ExitDirectionCosineX(int i,double yObj,double xObj,
               std::string SetRay,int FindPupil,double yPupil,double xPupil,int j){
	return ExitDirectionCosine(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j).x;
}
double cLens1::ExitDirectionCosineY(int i,double yObj,double xObj,
               std::string SetRay,int FindPupil,double yPupil,double xPupil,int j){
	return ExitDirectionCosine(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j).y;
}
double cLens1::ExitDirectionCosineZ(int i,double yObj,double xObj,
               std::string SetRay,int FindPupil,double yPupil,double xPupil,int j){
	return ExitDirectionCosine(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j).z;
}

double cLens1::ExitTanX(int i,double yObj,double xObj,
                        std::string SetRay,int FindPupil,double yPupil,double xPupil,int j){
	return ExitDirectionCosineX(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j)
	     / ExitDirectionCosineZ(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j);
}
double cLens1::ExitTanY(int i,double yObj,double xObj,
                        std::string SetRay,int FindPupil,double yPupil,double xPupil,int j){
	return ExitDirectionCosineY(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j)
	     / ExitDirectionCosineZ(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j);
}

double cLens1::DWDSag(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil){
	// �ʂ�z�����ψʁi�@�������ψʂł͂Ȃ��j�ɑ΂���g�ʎ����ω��̊������v�Z����D

	point p;
	vector<double> E,Q,Q1;
	const vector<double> ez(0,0,1);
	cLens1 ln=*this;  // �������Ȃ���x,y,z,X,Y,Z��ύX������̂�SingleRayTrace()�̒��Ŏg���ɂ����D

	if( 1<=i && i<=ln.k ){
		// �Q�l���F�����Y�H�w�̗��_�Ǝ��� 2.3 �ۓ������̌v�Z (���r�Y AOET 1984)
		// (2.12)�̉����̎��C
		//     dW = (df,g)( |N|(i,g) - |N'|(i',g) )
		//          dW   : �g�ʎ���
		//          df   : �ʂ̕ψ�(dx,dy,dz)
		//          g    : �ʖ@���i�����ł͒P�ʃx�N�g���Ƃ���j
		//          i,i' : ���������P�ʃx�N�g��
		//     �� i,i'�̌����������i�s�����Ɏ��Ƃ��́CN,N'��|N|,|N'|�Ƃ��ׂ��D
		//        �������Ȃ��ƁC�Ⴆ�΋�C���̔��˖ʂ�dW=0�ɂȂ��Ă��܂��D
		// ���,
		//     dW = dz(ez,E)( |N|(Q,E) - |N'||(Q',E)| )

		p=point(xPupil,yPupil);
		ln.FindRay(p,yObj,xObj,SetRay,1,1);
		if( ln.RayTrace(yObj,xObj,p.y,p.x,0,1,0,0,0,0) ) return 0;
		E =DirectionCosine(ln.surface_normal(i,ln.y[i],ln.x[i],0));
		Q =DirectionCosine(vector<double>(ln.X[i],ln.Y[i],ln.Z[i]));
		Q1=DirectionCosine(vector<double>(ln.X1[i],ln.Y1[i],ln.Z1[i]));
		return sProduct(E,ez)*( fabs(ln.N(i-1,1))*sProduct(Q,E)-fabs(ln.N(i,1))*sProduct(Q1,E) );
	}
	else{
		return 0;
	}
}

double cLens1::DSagDFringe(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil){
	// �g�ʎ�����Ȃ������z�����̖ʕψʂ���m�P�ʂŕԂ��D
	// �Ⴆ�΁CFT�̊��x���v�Z�ł���D
	double dwdsag;

	dwdsag=DWDSag(i,yObj,xObj,SetRay,yPupil,xPupil);
	return dwdsag==0 ? 0 : (Wl(1)/1000)/dwdsag;
}

vector<complex> cLens1::IncidentE(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j,
								  double pol_phi/*=90*/,double pol_ratio/*=0*/){
	// ��i�ʂɓ��˂���d���Ԃ��Di�ʒ��_���W�ɂ��D
	int errcode;
	point pupil(xPupil,yPupil);
	vector<complex> E0; double Ax,Ay,Phi;
	double opl;

	cPolarization::Analyze(Ax,Ay,Phi,1,sqrt(pol_ratio),pol_phi);
	E0.x=Ax;
	E0.y=complex(Ay*cos(Phi*PI/180),Ay*sin(Phi*PI/180));
	E0.z=0;

	FindRay(pupil,yObj,xObj,SetRay,1,j);
	if( errcode=RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,0,0,0,1,E0) ){  // �������ʂ�Ȃ��Ƃ�
		if( NotThruSurf(errcode)<=i ) return vector<complex>(0,0,0);     // ���Ci�ʈȑO�ł�����Ƃ�
	}
	if(i==0) i=k+1;  // i==0�̂Ƃ��͑��ʂƂ���
	if(1<=i && i<=k+1){
		opl=OPL(yObj,xObj,"",yPupil,xPupil,0,j,0,i,1);
		return E[i]*exp(complex(0,1)*2*PI*opl/(Wl(j)*1e-9));
	}
	else return vector<complex>(0,0,0);
}

double cLens1::IncidentAmplitude(std::string xyz,
								 int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j,
								 double pol_phi/*=90*/,double pol_ratio/*=0*/){
	// ��i�ʂɓ��˂���d��̐U����Ԃ��Di�ʒ��_���W�ɂ��D
	vector<complex> E=IncidentE(i,yObj,xObj,SetRay,yPupil,xPupil,j,pol_phi,pol_ratio);
	if     (xyz=="x" || xyz=="X") return abs(E.x);
	else if(xyz=="y" || xyz=="Y") return abs(E.y);
	else if(xyz=="z" || xyz=="Z") return abs(E.z);
	else                          return 0;
}

double cLens1::IncidentPhase(std::string xyz,
							 int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j,
							 double pol_phi/*=90*/,double pol_ratio/*=0*/){
	// ��i�ʂɓ��˂���d��̈ʑ���deg�ŕԂ��Di�ʒ��_���W�ɂ��D
	vector<complex> E=IncidentE(i,yObj,xObj,SetRay,yPupil,xPupil,j,pol_phi,pol_ratio);
	if     (xyz=="x" || xyz=="X") return arg(E.x)*180/PI;
	else if(xyz=="y" || xyz=="Y") return arg(E.y)*180/PI;
	else if(xyz=="z" || xyz=="Z") return arg(E.z)*180/PI;
	else                          return 0;
}

double cLens1::IncidentIntensity(int i,double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,int j,
								 double pol_phi/*=90*/,double pol_ratio/*=0*/){
	// ��i�ʂɓ��˂��鋭�x��Ԃ��D
	vector<complex> E=IncidentE(i,yObj,xObj,SetRay,yPupil,xPupil,j,pol_phi,pol_ratio);
	return sqabs(E.x)+sqabs(E.y)+sqabs(E.z);
}

double cLens1::FNumber(double yObj,double xObj,int findpupil){
	// �����L��F�i���o�����߂�i��1�g��)
	point ymax,ymin,xmax,xmin;
	vector<double> v1,v2;

	if(Afocal) return 0; // �A�t�H�[�J���̂Ƃ��͉������Ȃ�
	FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,1,findpupil,1);	
	if(yObj==0 && xObj==0){
		v1=ExitDirectionCosine(this->k,yObj,xObj,"",0,ymax.y,ymax.x,1);
		v2=ExitDirectionCosine(this->k,yObj,xObj,"",0,ymin.y,ymin.x,1);
	}
	else if(yObj!=0 && xObj==0){
		v1=ExitDirectionCosine(this->k,yObj,xObj,"",0,ymax.y,ymax.x,1);
		v2=ExitDirectionCosine(this->k,yObj,xObj,"",0,ymin.y,ymin.x,1);
	}
	else if(yObj==0 && xObj!=0){
		v1=ExitDirectionCosine(this->k,yObj,xObj,"",0,xmax.y,xmax.x,1);
		v2=ExitDirectionCosine(this->k,yObj,xObj,"",0,xmin.y,xmin.x,1);
	}
	else{
		// ���_��Y����܂���X����łȂ��Ƃ��͌v�Z���Ȃ�
		return 0;
	}
	// �L��F�i���o = ��������p�x�̔�����a�Ƃ���Ƃ�, 1/(2*sina)  
	//  (JOEM�����p�Z�p1988 "���w�n�̓���(�ғ�����)" ���)
	return 1/( 2*sin( acos(sProduct(v1,v2))/2 ) );
}

double cLens1::FNumberXYAve(double yObj,double xObj,int findpupil){
	// �����L��F�i���o�����߂�i��1�g��)
	// x,y�����̑��敽�ςƂ���i���邳�͋t���̓��ɔ��j
	point ymax,ymin,xmax,xmin;
	vector<double> vy1,vy2,vx1,vx2;
	double Fx,Fy;

	if(Afocal) return 0; // �A�t�H�[�J���̂Ƃ��͉������Ȃ�
	FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,1,findpupil,1);	
	vy1=ExitDirectionCosine(this->k,yObj,xObj,"",0,ymax.y,ymax.x,1);
	vy2=ExitDirectionCosine(this->k,yObj,xObj,"",0,ymin.y,ymin.x,1);
	vx1=ExitDirectionCosine(this->k,yObj,xObj,"",0,xmax.y,xmax.x,1);
	vx2=ExitDirectionCosine(this->k,yObj,xObj,"",0,xmin.y,xmin.x,1);
	// �L��F�i���o = ��������p�x�̔�����a�Ƃ���Ƃ�, 1/(2*sina)  
	//  (JOEM�����p�Z�p1988 "���w�n�̓���(�ғ�����)" ���)
	Fy=1/( 2*sin( acos(sProduct(vy1,vy2))/2 ) );
	Fx=1/( 2*sin( acos(sProduct(vx1,vx2))/2 ) );
	return sqrt(Fy*Fx);	
}

void cLens1::FNumberObj(double &FNumY,double &FNumX,double yObj,double xObj,int findpupil){
	//�����L��F�i���o�����߂�i��1�g��, Y������X����)
	point ymax,ymin,xmax,xmin;
	vector<double> v1,v2;

	if(fabs(s)>LN){
		FNumY=FNumX=0;
		return; // ���������̂̂Ƃ��͉������Ȃ�
	}

	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,1,findpupil,1)){
		v1=IncidentDirectionCosine(1,yObj,xObj,"",0,ymax.y,ymax.x,1);
		v2=IncidentDirectionCosine(1,yObj,xObj,"",0,ymin.y,ymin.x,1);
		// �L��F�i���o = ��������p�x�̔�����a�Ƃ���Ƃ�, 1/(2*sina)  
		//  (JOEM�����p�Z�p1988 "���w�n�̓���(�ғ�����)" ���)
		FNumY=1/( 2*sin( acos(sProduct(v1,v2))/2 ) );

		v1=IncidentDirectionCosine(1,yObj,xObj,"",0,xmax.y,xmax.x,1);
		v2=IncidentDirectionCosine(1,yObj,xObj,"",0,xmin.y,xmin.x,1);
		FNumX=1/( 2*sin( acos(sProduct(v1,v2))/2 ) );
	}
	else{
		FNumY=FNumX=0;
	}
}

double cLens1::FNumberObj(double yObj,double xObj,int findpupil){
	//�����L��F�i���o�����߂�i��1�g��)
	double FNum,FNumY,FNumX;

	if(fabs(s)>LN) return 0; // ���������̂̂Ƃ��͉������Ȃ�
	FNumberObj(FNumY,FNumX,yObj,xObj,findpupil);
	if(yObj==0 && xObj==0){
		FNum=FNumY;
	}
	else if(yObj!=0 && xObj==0){
		FNum=FNumY;
	}
	else if(yObj==0 && xObj!=0){
		FNum=FNumX;
	}
	else{
		// ���_��Y����܂���X����łȂ��Ƃ��͌v�Z���Ȃ�
		return 0;
	}

	return FNum;
}

double cLens1::FNumberObjXYAve(double yObj,double xObj,int findpupil){
	// �����L��F�i���o�����߂�i��1�g��)
	// x,y�����̑��敽�ςƂ���i���邳�͋t���̓��ɔ��j
	point ymax,ymin,xmax,xmin;
	vector<double> vy1,vy2,vx1,vx2;
	double Fy,Fx;

	if(fabs(s)>LN) return 0; // ���������̂̂Ƃ��͉������Ȃ�
	FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,1,findpupil,1);	
	vy1=IncidentDirectionCosine(1,yObj,xObj,"",0,ymax.y,ymax.x,1);
	vy2=IncidentDirectionCosine(1,yObj,xObj,"",0,ymin.y,ymin.x,1);
	vx1=IncidentDirectionCosine(1,yObj,xObj,"",0,xmax.y,xmax.x,1);
	vx2=IncidentDirectionCosine(1,yObj,xObj,"",0,xmin.y,xmin.x,1);
	// �L��F�i���o = ��������p�x�̔�����a�Ƃ���Ƃ�, 1/(2*sina)  
	//  (JOEM�����p�Z�p1988 "���w�n�̓���(�ғ�����)" ���)
	Fy=1/( 2*sin( acos(sProduct(vy1,vy2))/2 ) );
	Fx=1/( 2*sin( acos(sProduct(vx1,vx2))/2 ) );
	return sqrt(Fy*Fx);
}

double cLens1::FNumberParaxial(int findpupil){
	// �����L��F�i���o���ߎ����_�ŋ��߂�i��1�g��)�D
	// ���������̒ǐՂɂ��FNumber()�͎���(�Ⴆ�΍������ʎ���)�ɉe������Ă��܂��D
	cLens1 X;
	point p;
	int ymax_i,i;
	double h,NA;
	
	X=*this;
	if(findpupil){
		X.FindPupilEdge(p,p,p,p,ymax_i,i,i,i,0,0,1,findpupil,1);  // �J���i��̖ʔԍ������߂�
		if(ymax_i>0){
			X.Set_stop(NotThruSurf(ymax_i));
			X.EPCalc();
		}
		else{
			return 0;
		}
	}
	h=X.EPD/2;
	NA= M()==0 ? h/f() : sin(atan(h/(s-t)))/M();
	return fabs(1/(NA*2));  // M()==0�̂Ƃ��Cf()/(2*h)
}

double cLens1::NA(double yObj,double xObj,int findpupil){
	// ����NA�����߂�i��1�g��)
	double x;
	x=FNumber(yObj,xObj,findpupil);
	return x==0 ? 0 : N(k,1)/2/x;
}

double cLens1::NAXYAve(double yObj,double xObj,int findpupil){
	// ����NA�����߂�i��1�g��)
	// x,y�����̑��敽�ςƂ���i���邳�͓��ɔ��j
	double x;
	x=FNumberXYAve(yObj,xObj,findpupil);
	return x==0 ? 0 : N(k,1)/2/x;
}

void cLens1::NAObj(double &NAObjY,double &NAObjX,double yObj,double xObj,int findpupil){
	double FNumY,FNumX;

	FNumberObj(FNumY,FNumX,yObj,xObj,findpupil);
	NAObjY= FNumY==0 ? 0 : N(0,1)/2/FNumY;
	NAObjX= FNumX==0 ? 0 : N(0,1)/2/FNumX;
}

double cLens1::NAObj(double yObj,double xObj,int findpupil){
	// ����NA�����߂�i��1�g��)
	double x;
	x=FNumberObj(yObj,xObj,findpupil);
	return x==0 ? 0 : N(0,1)/2/x;
}

double cLens1::NAObjXYAve(double yObj,double xObj,int findpupil){
	// ����NA�����߂�i��1�g��)
	// x,y�����̑��敽�ςƂ���i���邳�͓��ɔ��j
	double x;
	x=FNumberObjXYAve(yObj,xObj,findpupil);
	return x==0 ? 0 : N(0,1)/2/x;
}

double cLens1::NAParaxial(int findpupil){
	// ����NA���ߎ����_�ŋ��߂�i��1�g��)�D
	double x;
	x=FNumberParaxial(findpupil);
	return x==0 ? 0 : N(k,1)/2/x;
}

void cLens1::DiffractObj(double yObj,double xObj,int findpupil){
	// �ʕ��̂̑傫���ɑ�1�g���̉�܂�ݒ肷��
	// ���`��̊����厲���X���Ă����(����X,Y���ɉ����ȉ~�łȂ���)���m�łȂ��D
	double NAy,NAx;

	NAObj(NAy,NAx,yObj,xObj,findpupil);

	SourcePhiY= NAy==0 ? 0 : 1.22*(Wl(1)/1000000)/NAy;  // �G�A���[�f�B�X�N��1�Ê̒��a
	SourcePhiX= NAx==0 ? 0 : 1.22*(Wl(1)/1000000)/NAx;
	SourceAxis=0;
	SourceAreaType=0;
}

void cLens1::DiffractObj(double yObj,double xObj,int FindPupil,int nSpot,int DeleteFrom){
	// �ʕ��̂̑傫���ɑ�1�g���̉�܂�ݒ肷��
	// ���`��̊����厲���X���Ă���ꍇ�ɂ��Ή�����D
	// (SLO�Ŋ􉽌��w�I�ɂ͍��_�őΕ����ˋP�_��}���Ă���̂ɁC
	//  �����ł͗}������Ă��Ȃ����̌����p�ɖ{�֐����쐬�j
	cLens1 X;
	point g;
	double th,phix,phiy, x,y, Fx,Fy,Axis,NAx,NAy;
	vector<double> xmax,xmin,ymax,ymin, o;

	X=*this;

	if(nSpot>0) X.nSpot=nSpot;
	if(DeleteFrom>1) X.Delete(DeleteFrom,X.k);   // ��DeleteFrom�ʈȍ~�͍폜�D
	X.make_coordinate(0);

	if(X.MakeSpot(yObj,xObj,FindPupil,0,1,1,1,0,0,0,0)==0){   // ��1�ʏ�̃t�b�g�v�����g�����
		this->SourcePhiX=0;
		this->SourcePhiY=0;
		this->SourceAxis=0;
		this->SourceAreaType=0;
		return;      // �t�b�g�v�����g���쐬�ł��Ȃ�(�X�|�b�g����0)�ꍇ�͏I��
	}

	g.x=X.SpotXGravityCenter();                        // x�d�S
	g.y=X.SpotYGravityCenter();                        // y�d�S
	th=X.SpotInertiaPrincipalAxis()*PI/180;            // �����厲�̕���(rad)
	phix=X.SpotPrincipalWx();                          // �����厲�����̃X�|�b�g���a
	phiy=X.SpotPrincipalWy();                          // �����厲���������̃X�|�b�g���a

	x=phix/2*cos(th);
	y=phix/2*sin(th);
	xmax=vector<double>(g.x+x, g.y+y, X.surface_sag(1,g.x+x,g.y+y,0));  // xmax,xmin�͎厲�������a��
	xmin=vector<double>(g.x-x, g.y-y, X.surface_sag(1,g.x-x,g.y-y,0));  // ���[�̈ʒu�x�N�g��

	x=-phiy/2*sin(th);
	y= phiy/2*cos(th);
	ymax=vector<double>(g.x+x, g.y+y, X.surface_sag(1,g.x+x,g.y+y,0));  // ymax,ymin�͎厲�������a��
	ymin=vector<double>(g.x-x, g.y-y, X.surface_sag(1,g.x-x,g.y-y,0));  // ���[�̈ʒu�x�N�g��

	X.transform(xmax,1,0,0,0,1);    // ��0��(���̖�)�̍��W�ɕϊ�
	X.transform(xmin,1,0,0,0,1);
	X.transform(ymax,1,0,0,0,1);
	X.transform(ymin,1,0,0,0,1);

	o=vector<double>(xObj,yObj,X.surface_sag(0,yObj,xObj,0));   // ���_�̍��W

	Fx=1/( 2*sin( acos(cosine(xmax-o,xmin-o))/2 ) );    // F�i���o�[
	Fy=1/( 2*sin( acos(cosine(ymax-o,ymin-o))/2 ) );

	NAx= Fx==0 ? 0 : X.N(0,1)/2/Fx;
	NAy= Fy==0 ? 0 : X.N(0,1)/2/Fy;

	x=xmax.x-xmin.x; 
	y=xmax.y-xmin.y;
	Axis= x==0 ? 90 : atan(y/x)*180/PI;  // ���� xmax�`xmin ��XY�ʂւ̓��e�̌X��

	this->SourcePhiX= NAx==0 ? 0 : 1.22*(Wl(1)/1000000)/NAx;
	this->SourcePhiY= NAy==0 ? 0 : 1.22*(Wl(1)/1000000)/NAy;
	this->SourceAxis=Axis;
	this->SourceAreaType=0;
}

double cLens1::FresnelNumber(int findpupil,double defocus){
	// �t���l���� = (2N/��){��(Z^2+A^2)-sgn(Z)*Z}  (��������0�Ƃ�������. �����ɉ������ˏo�����ʂƔg�ʂƂ̋���)
	//                  N= ����ԋ��ܗ� 
	//                 ��= �^�󒆔g��
	//                  Z= ���˓�����ˏo���܂ł̋���
	//                  A= ���˓����a
	//  ��Ԃ��D���̒l���������قǉ�����C�傫���قǋߎ���ƂȂ�D
	double wl,n,z,z0,a, frenum0,frenum1;

	wl=Wl(1)/1000000;             // ��1�g��(mm)
	n=N(k,1);                     // ����ԋ��ܗ�
	a=ExitPupilDia(findpupil)/2;  // �ˏo�����a
	z0=ExitPupilZ(findpupil);     // �ˏo���ʒu

	z=s1()-z0;
	frenum0=(2*n/wl)*(sqrt(z*z+a*a)-sgn(z)*z);  // �ߎ����ʂł̃t���l����

	if(Afocal){
		z= defocus==0 ? z=LN : 1000/defocus*sgn(n)-z0;
	}
	else{
		z=dk(defocus)-z0;
	}
	frenum1=(2*n/wl)*(sqrt(z*z+a*a)-sgn(z)*z);  // �w�葜�ʂł̃t���l����

	return fabs(frenum0-frenum1);
}

double cLens1::AngleOfView(double h){
	// ����p��deg�P�ʂŕԂ��D
	if(fabs(s)>=LN){
		return atan(h)*180/PI;
	}
	else{
		return atan(h/g_hat(1))*180/PI;
	}
}

double cLens1::LSA(double yPupil,double xPupil,double yPupilPrincipal,double xPupilPrincipal,int j,int AbeTheory/*=0*/){
	// ���F�X���b�g�����v�Ȃǎ���������W��(0,0)�Ƃ͌���Ȃ��DyPupilPrincipal,xPupilPrincipal�͕K�v�D
	double dz,dz0;
	double yp,zp,Yp,Zp, y,z,Y,Z;
	
	if(AbeTheory==0){
		// z�������_�͑�1�g��������̏œ_�D�ΐS�n���l����y�����̏œ_�Ƃ���D
		dz0=DeltaM(0,0,yPupilPrincipal,xPupilPrincipal,0,1);  // �����i�s��������

		if(yPupil==yPupilPrincipal && xPupil==xPupilPrincipal){
			dz=DeltaM(0,0,yPupil,xPupil,0,j);  // �����i�s��������
		}
		else{
			// �����(��j�g��)��ǐ�
			// ���j�ΐS�n�ł͑�1�g���Ƒ�j�g���Ŏ�������قȂ邽�߁C
			//     ��1�g���̎�����ł́C�Ⴂ���ˍ��ŋߖT�lDeltaM()�Ɏ������Ȃ��D
			RayTrace(0,0,yPupilPrincipal,xPupilPrincipal,0,j,0,0,0,0);
			yp=this->y[k+1];
			zp=this->z[k+1];
			Yp=this->Y[k+1];
			Zp=this->Z[k+1];
			// ������ǐ�
			RayTrace(0,0,yPupil,xPupil,0,j,0,0,0,0);
			y=this->y[k+1];
			z=this->z[k+1];
			Y=this->Y[k+1];
			Z=this->Z[k+1];
			// ������Ƃ������Ō����ʒudz�����߂�
			// yp+(dz-xp)Yp/Zp = y+(dz-x)Y/Z ���C
			if(Yp/Zp-Y/Z!=0){
				dz=(-yp+zp*Yp/Zp+y-z*Y/Z)/(Yp/Zp-Y/Z);
			}
			else{
				return 0;
			}

			if(Afocal) dz=1000/dz;
			dz*=sgn(Z); // �����i�s��������(������Z�̕��̕����ɐi�s����Ƃ��C�␳�̃A���_�[�^�I�[�o�[�𒼊ϓI�ɂ��邽�߁j
		}

		return dz-dz0;
	}
	else{
		return LSAAbe(yPupil,xPupil,yPupilPrincipal,xPupilPrincipal,j,AbeTheory);
	}
}

double cLens1::LSAAbe(double yPupil,double xPupil,double yPupilPrincipal,double xPupilPrincipal,int j,int order/*=3*/){
	// �����ǐՂł͂Ȃ��C�����W���ɂ�苅�ʏc�������v�Z����
	std::string normalizeunit0;
	int normalizetype0;
	point pupil;
	double g,g_hat,g1_hat,N1, R,dz;

	normalizeunit0=NormalizeUnit; normalizetype0=NormalizeType; // �ۑ�
	NormalizeUnit="1";   // �X�P�[���̋K�i���͂��Ȃ�
	NormalizeType=1;     // �g���K��1�h�Ƃ��C�g�����Y�݌v�@�h�� (4.22) ��p����D
	CalcCoefficients(j);  // �F�����͑�j�g�� vs ��1�g��

	pupil=point(yPupil-yPupilPrincipal,xPupil-xPupilPrincipal);

	if(fabs(s)<LN){
		g=s-t;
		g_hat=s-delta();
		pupil.y*=g_hat/g; pupil.x*=g_hat/g;  // �����W�͎啽�ʏ�Ƃ���
	}
	
	R=abs(pupil);
	g1_hat=this->g1_hat();
	N1=N(k,1);
	dz=-(g1_hat*g1_hat/2/N1)*( SAt*R*R +2*LCt );
	if(order==5) dz-=(g1_hat*g1_hat/2/N1)*(SA5t/4-(g1_hat/2/N1)*SAt*SAt)*R*R*R*R;
	//   �g�����Y�݌v�@�h(4.22)���D�F������(4.18) ����сg�����_�h(3.3.1) ���
	if(Afocal) dz=1000/dz;
	dz*=sgn(N1);  // �����i�s��������

	NormalizeUnit=normalizeunit0; NormalizeType=normalizetype0; // ���ɖ߂�
	return dz;
}

double cLens1::LSA(double yPupilNormalized,int findpupil,int j/*=1*/){
	point ymax,ymin,xmax,xmin, p0,p;

	if(FindPupilEdge(ymax,ymin,xmax,xmin,0,0,1,findpupil,1)){
		p0=(ymax+ymin)/2;
		p=p0+(ymax-p0)*yPupilNormalized;
		return LSA(p.y,p.x,p0.y,p0.x,j);
	}
	else{
		return 0;
	}
}

double cLens1::LSA(){
	return LSA(1,1,1);
}
double cLens1::LSA70(){
	return LSA(0.7,1,1);
}
double cLens1::LSA50(){
	return LSA(0.5,1,1);
}

double cLens1::LSA2nd(){
	return LSA(1,1,2);
}
double cLens1::LSA3rd(){
	return LSA(1,1,3);
}
double cLens1::DLSA2to3(){
	return LSA3rd()-LSA2nd();
}

double cLens1::LSAp(double yObjNormalized,int findfield,int j/*=1*/){
	cLens1 buf=*this;
	return buf.SwapObjPupil().LSA(yObjNormalized,findfield,j);
}

double cLens1::LSAp(){
	cLens1 buf=*this;
	return buf.SwapObjPupil().LSA();
}
double cLens1::LSAp70(){
	cLens1 buf=*this;
	return buf.SwapObjPupil().LSA70();
}
double cLens1::LSAp50(){
	cLens1 buf=*this;
	return buf.SwapObjPupil().LSA50();
}

double cLens1::LSAp2nd(double yObjNormalized,int findfield){
	cLens1 buf=*this;
	return buf.SwapObjPupil().LSA(yObjNormalized,findfield,2);
}
double cLens1::LSAp2nd(){
	cLens1 buf=*this;
	return buf.SwapObjPupil().LSA2nd();
}

std::string cLens1::LSAs(int findpupil,int colors,double weight/*=1*/){
	// �����̓��������ł̏c���ʎ�����]���֐��Ƃ���Coptimize()�̈���targets�𐶐�����D
	int j;
	std::string s;
	char buf[1000];
	double w;

	if(colors>cn) colors=cn;
	for(j=1; j<=colors; j++){
		w=colorweight[j]*weight;
		sprintf(buf,"(LSA 1    %d %d) 0 %g;\n", findpupil,j,w     ); s+=buf;
		sprintf(buf,"(LSA 0.95 %d %d) 0 %g;\n", findpupil,j,w*0.95); s+=buf;
		sprintf(buf,"(LSA 0.9  %d %d) 0 %g;\n", findpupil,j,w*0.9 ); s+=buf;
		sprintf(buf,"(LSA 0.85 %d %d) 0 %g;\n", findpupil,j,w*0.85); s+=buf;
		sprintf(buf,"(LSA 0.8  %d %d) 0 %g;\n", findpupil,j,w*0.8 ); s+=buf;
		sprintf(buf,"(LSA 0.75 %d %d) 0 %g;\n", findpupil,j,w*0.75); s+=buf;
		sprintf(buf,"(LSA 0.7  %d %d) 0 %g;\n", findpupil,j,w*0.7 ); s+=buf;
		sprintf(buf,"(LSA 0.65 %d %d) 0 %g;\n", findpupil,j,w*0.65); s+=buf;
		sprintf(buf,"(LSA 0.6  %d %d) 0 %g;\n", findpupil,j,w*0.6 ); s+=buf;
		sprintf(buf,"(LSA 0.55 %d %d) 0 %g;\n", findpupil,j,w*0.55); s+=buf;
		sprintf(buf,"(LSA 0.5  %d %d) 0 %g;\n", findpupil,j,w*0.5 ); s+=buf;
		sprintf(buf,"(LSA 0.4  %d %d) 0 %g;\n", findpupil,j,w*0.4 ); s+=buf;
		sprintf(buf,"(LSA 0.3  %d %d) 0 %g;\n", findpupil,j,w*0.3 ); s+=buf;
		sprintf(buf,"(LSA 0.2  %d %d) 0 %g;\n", findpupil,j,w*0.2 ); s+=buf;
	}
	return s;
}

std::string cLens1::LSAps(int findpupil,int colors,double weight/*=1*/){
	// �����̓��������ł̏c���ʎ�����]���֐��Ƃ���Coptimize()�̈���targets�𐶐�����D
	int j;
	std::string s;
	char buf[1000];
	double w;

	if(colors>cn) colors=cn;
	for(j=1; j<=colors; j++){
		w=colorweight[j]*weight;
		sprintf(buf,"(LSAp 1    %d %d) 0 %g;\n", findpupil,j,w     ); s+=buf;
		sprintf(buf,"(LSAp 0.95 %d %d) 0 %g;\n", findpupil,j,w*0.95); s+=buf;
		sprintf(buf,"(LSAp 0.9  %d %d) 0 %g;\n", findpupil,j,w*0.9 ); s+=buf;
		sprintf(buf,"(LSAp 0.85 %d %d) 0 %g;\n", findpupil,j,w*0.85); s+=buf;
		sprintf(buf,"(LSAp 0.8  %d %d) 0 %g;\n", findpupil,j,w*0.8 ); s+=buf;
		sprintf(buf,"(LSAp 0.75 %d %d) 0 %g;\n", findpupil,j,w*0.75); s+=buf;
		sprintf(buf,"(LSAp 0.7  %d %d) 0 %g;\n", findpupil,j,w*0.7 ); s+=buf;
		sprintf(buf,"(LSAp 0.65 %d %d) 0 %g;\n", findpupil,j,w*0.65); s+=buf;
		sprintf(buf,"(LSAp 0.6  %d %d) 0 %g;\n", findpupil,j,w*0.6 ); s+=buf;
		sprintf(buf,"(LSAp 0.55 %d %d) 0 %g;\n", findpupil,j,w*0.55); s+=buf;
		sprintf(buf,"(LSAp 0.5  %d %d) 0 %g;\n", findpupil,j,w*0.5 ); s+=buf;
		sprintf(buf,"(LSAp 0.4  %d %d) 0 %g;\n", findpupil,j,w*0.4 ); s+=buf;
		sprintf(buf,"(LSAp 0.3  %d %d) 0 %g;\n", findpupil,j,w*0.3 ); s+=buf;
		sprintf(buf,"(LSAp 0.2  %d %d) 0 %g;\n", findpupil,j,w*0.2 ); s+=buf;
	}
	return s;
}


double cLens1::SC(double yPupil,double xPupil){
/////////////////////////////////////////////////////////////	
	if( Afocal ) return 0;    // afocal�̂Ƃ�������
/////////////////////////////////////////////////////////////
	double R,cosu,sinu1,ghat1; 
	RayTrace(0,0,yPupil,xPupil,0,1,0,0,0,0);
	if( fabs(s)<LN ) R=(delta()-s)*Y[1]/Z[1]; else R=yPupil;
	cosu=Z[1];
	sinu1=-Y[k+1];
	ghat1=s1()-delta1();
	if(R==0 || sinu1==0){
		return 0;
	}
	else{
		return (R*cosu/sinu1-ghat1) *sgn(Z[k+1]);  // "�����Y�݌v�@" (3.47)
		                                           // �����i�s���������i�O���t�Ɉꏏ�ɕ\������LSA�ɍ��킹��j
	}
}

double cLens1::SC(double yPupilNormalized,int findpupil){
	point ymax,ymin,xmax,xmin, p0,p;

	if(FindPupilEdge(ymax,ymin,xmax,xmin,0,0,1,findpupil,1)){
		p0=(ymax+ymin)/2;
		p=p0+(ymax-p0)*yPupilNormalized;
		return SC(p.y,p.x);
	}
	else{
		return 0;
	}
}

double cLens1::SC(){
	point ymax,ymin,xmax,xmin;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,0,0,1,1)){
		return SC(ymax.y,ymax.x);
	}
	else{
		return 0;
	}
}
double cLens1::SC70(){
	point ymax,ymin,xmax,xmin, p;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,0,0,1,1)){
		p=(ymax+ymin)/2 + (ymax-(ymax+ymin)/2)*0.7;
		return SC(p.y,p.x);
	}
	else{
		return 0;
	}
}
double cLens1::SC50(){
	point ymax,ymin,xmax,xmin, p;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,0,0,1,1)){
		p=(ymax+ymin)/2 + (ymax-(ymax+ymin)/2)*0.5;
		return SC(p.y,p.x);
	}
	else{
		return 0;
	}
}

double cLens1::OSC(double yPupilNormalized,int findpupil){
	return SC(yPupilNormalized,findpupil)-LSA(yPupilNormalized,findpupil,1);
}

double cLens1::OSC(){
	// "�����Y�݌v�@" (3.46)��g1_hat/g_hat���ȗ���������
	return SC()-LSA();
}
double cLens1::OSC70(){
	return SC70()-LSA70();
}
double cLens1::OSC50(){
	return SC50()-LSA50();
}

double cLens1::OSCp(double yObjNormalized,int findfield){
	cLens1 buf=*this;
	return buf.SwapObjPupil().OSC(yObjNormalized,findfield);
}

double cLens1::OSCp(){
	cLens1 buf=*this;
	return buf.SwapObjPupil().OSC();
}

/*
double cLens1::DeltaM(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j){
	RayTrace(yObj,xObj,yPupil,xPupil,defocus,j,0,0);
	return fabs(s1fix)<LN ? hyxi1*Z[k+1]/uy1 : ( -uy1/hyxi1/Z[k+1]+defocus/1000 )*1000;
}
*/
double cLens1::DeltaM(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j,int AbeTheory/*=0*/){
	if(AbeTheory==0){
		RayTrace(yObj,xObj,yPupil,xPupil,defocus,j,0,0,0,1);
		if(Afocal){
			if(AfocalRotateEye){
				// Z�Ŋ���Ȃ��ق�����̉���⎋�x�]����ƍ���
				return (-dQa[k+1].y/ha[k+1].y-defocus/1000)*1000;
			}
			else{
				// �����i�s���������DZ���ɉ������ʂɂ��邽�߂�|Z|�Ŋ���D
				return (-dQa[k+1].y/ha[k+1].y/fabs(Z[k+1])-defocus/1000)*1000;
			}
		}
		else{
			return (-ha[k+1].y/dQa[k+1].y)*fabs(Z[k+1]);
		}
	}
	else{
		return DeltaMAbe(yObj,xObj,AbeTheory);
	}
}

double cLens1::DeltaMAbe(double yObj,double xObj,int order/*=3*/){
	// �����ǐՂł͂Ȃ��C�����W���ɂ��q�ߑ��ʂ��v�Z����
	std::string normalizeunit0;
	int normalizetype0;
	double g1_hat,N1, Y,tw,dz;

	normalizeunit0=NormalizeUnit; normalizetype0=NormalizeType; // �ۑ�
	NormalizeUnit="1";   // �X�P�[���̋K�i���͂��Ȃ�
	NormalizeType=1;     // �g���K��1�h�Ƃ��C�g�����Y�݌v�@�h�� (4.22) ��p����D
	CalcCoefficients();
	
	Y=sqrt(yObj*yObj+xObj*xObj);
	if(fabs(s)<LN) tw=N(0,1)*Y/g_hat(); else tw=N(0,1)*Y;

	g1_hat=this->g1_hat();
	N1=N(k,1);
	dz=-(g1_hat*g1_hat/2/N1)*( (3*ASt+PTt)*tw*tw );  // �g�����Y�݌v�@�h (4.22)
	if(order==5) dz-=(g1_hat*g1_hat/2/N1)*( (4*AS5t+SG5t)/4 -(g1_hat/2/N1)*(3*ASt+PTt)*(3*ASt+PTt) )
		             *tw*tw*tw*tw;
	if(Afocal) dz=1000/dz;

	dz*=sgn(N1);  // �����i�s�����𐳂Ƃ���

	NormalizeUnit=normalizeunit0; NormalizeType=normalizetype0; // ���ɖ߂�
	return dz;
}

double cLens1::DeltaM(double yObjNormalized,int findpupil,int j/*=1*/){
	point ymax,ymin,xmax,xmin, pupil;

	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax*yObjNormalized,0,1,findpupil,1)) {
		pupil=(ymax+ymin)/2;
		return DeltaM(yObjectMax*yObjNormalized,0,pupil.y,pupil.x,0,j);
	}
	else return 0;
}
double cLens1::DeltaM(){
	return DeltaM(1,1);
}
double cLens1::DeltaM70(){
	return DeltaM(0.7,1);
}
double cLens1::DeltaM50(){
	return DeltaM(0.5,1);
}
/*
double cLens1::DeltaS(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j){
	RayTrace(yObj,xObj,yPupil,xPupil,defocus,j,0,0);
	return fabs(s1fix)<LN ? hx*Z[k+1]/ux1 : ( -ux1/hx/Z[k+1]+defocus/1000 )*1000;
}
*/
double cLens1::DeltaS(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j,int AbeTheory/*=0*/){
	if(AbeTheory==0){
		RayTrace(yObj,xObj,yPupil,xPupil,defocus,j,0,0,0,1);
		if(Afocal){
			if(AfocalRotateEye){
				// Z�Ŋ���Ȃ��ق�����̉���⎋�x�]����ƍ���
				return ( -dQb[k+1].x/hb[k+1].x-defocus/1000 )*1000;
			}
			else{
				// �����i�s���������DZ���ɉ������ʂɂ��邽�߂�|Z|�Ŋ���D
				return ( -dQb[k+1].x/hb[k+1].x/fabs(Z[k+1])-defocus/1000 )*1000;
			}
		}
		else{
			return (-hb[k+1].x/dQb[k+1].x)*fabs(Z[k+1]);
		}
	}
	else{
		return DeltaSAbe(yObj,xObj,AbeTheory);
	}
}

double cLens1::DeltaSAbe(double yObj,double xObj,int order/*=3*/){
	// �����ǐՂł͂Ȃ��C�����W���ɂ��q�ߑ��ʂ��v�Z����
	std::string normalizeunit0;
	int normalizetype0;
	double g1_hat,N1, Y,tw,dz;

	normalizeunit0=NormalizeUnit; normalizetype0=NormalizeType; // �ۑ�
	NormalizeUnit="1";   // �X�P�[���̋K�i���͂��Ȃ�
	NormalizeType=1;     // �g���K��1�h�Ƃ��C�g�����Y�݌v�@�h�� (4.22) ��p����D
	CalcCoefficients();
	
	Y=sqrt(yObj*yObj+xObj*xObj);
	if(fabs(s)<LN) tw=N(0,1)*Y/g_hat(); else tw=N(0,1)*Y;

	g1_hat=this->g1_hat();
	N1=N(k,1);
	dz=-(g1_hat*g1_hat/2/N1)*( (ASt+PTt)*tw*tw );  // �g�����Y�݌v�@�h (4.22)
	if(order==5) dz-=(g1_hat*g1_hat/2/N1)*( SG5t/4 -(g1_hat/2/N1)*(ASt+PTt)*(ASt+PTt) )
		             *tw*tw*tw*tw;
	if(Afocal) dz=1000/dz;
	dz*=sgn(N1);  // �����i�s�����𐳂Ƃ���

	NormalizeUnit=normalizeunit0; NormalizeType=normalizetype0; // ���ɖ߂�
	return dz;
}

double cLens1::DeltaS(double yObjNormalized,int findpupil,int j/*=1*/){
	point ymax,ymin,xmax,xmin, pupil;

	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax*yObjNormalized,0,1,findpupil,1)) {
		pupil=(ymax+ymin)/2;
		return DeltaS(yObjectMax*yObjNormalized,0,pupil.y,pupil.x,0,j);
	}
	else return 0;
}
double cLens1::DeltaS(){
	return DeltaS(1,1);
}
double cLens1::DeltaS70(){
	return DeltaS(0.7,1);
}
double cLens1::DeltaS50(){
	return DeltaS(0.5,1);
}

double cLens1::DSToDM(double yObj,double xObj,double yPupil,double xPupil,double defocus,int j){
	// Cyl()���ƁC-(+)�\���Ȃ̂Œl��0�Ő܂�Ԃ��ƂȂ�̂ŁC�����݌v�ɕs�����D
	return DeltaM(yObj,xObj,yPupil,xPupil,defocus,j)-DeltaS(yObj,xObj,yPupil,xPupil,defocus,j);
}
double cLens1::DSToDM(double yObjNormalized,int findpupil){
	point ymax,ymin,xmax,xmin, pupil;

	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax*yObjNormalized,0,1,findpupil,1)) {
		pupil=(ymax+ymin)/2;
		return DSToDM(yObjectMax*yObjNormalized,0,pupil.y,pupil.x,0,1);
	}
	else return 0;
}
double cLens1::DSToDM(){
	return DeltaM()-DeltaS();
}
double cLens1::DSToDM70(){
	return DeltaM70()-DeltaS70();
}
double cLens1::DSToDM50(){
	return DeltaM50()-DeltaS50();
}

std::string cLens1::DSDMs(int findpupil,double s_weight/*=1*/,double m_weight/*=1*/,int cols/*=1*/,int fine_pitch/*=0*/){
	// �����̕��̍��ł̔�_������]���֐��Ƃ���Coptimize()�̈���targets�𐶐�����D
	int j;
	std::string s;
	char buf[1000];

	for(j=1; j<=cols; ++j){
		                  sprintf(buf,"(DeltaS 1     %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.975 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.95  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.925 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.9   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.875 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.85  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.825 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.8   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.775 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.75  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.725 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.7   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.675 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.65  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.625 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.6   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.575 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.55  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.525 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.5   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.475 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.45  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.425 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.4   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.375 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.35  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.325 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.3   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.275 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.25  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.225 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.2   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.175 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.15  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.125 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0.1   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.075 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.05  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaS 0.025 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaS 0     %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;

		                  sprintf(buf,"(DeltaM 1     %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.975 %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.95  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.925 %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.9   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.875 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.85  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.825 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.8   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.775 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.75  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.725 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.7   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.675 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.65  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.625 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.6   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.575 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.55  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.525 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.5   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.475 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.45  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.425 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.4   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.375 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.35  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.325 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.3   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.275 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.25  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.225 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.2   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.175 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.15  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.125 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0.1   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.075 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.05  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		if(fine_pitch!=0) sprintf(buf,"(DeltaM 0.025 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		                  sprintf(buf,"(DeltaM 0     %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;

		if(!IsXAxisSymmetric()){
			                  sprintf(buf,"(DeltaS -1     %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.975 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.95  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.925 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.9   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.875 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.85  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.825 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.8   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.775 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.75  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.725 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.7   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.675 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.65  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.625 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.6   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.575 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.55  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.525 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.5   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.475 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.45  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.425 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.4   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.375 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.35  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.325 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.3   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.275 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.25  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(Deltas -0.225 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.2   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.175 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.15  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.125 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaS -0.1   %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.075 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.05  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaS -0.025 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;

			                  sprintf(buf,"(DeltaM -1     %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.975 %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.95  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.925 %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.9   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.875 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.85  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.825 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.8   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.775 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.75  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.725 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.7   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.675 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.65  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.625 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.6   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.575 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.55  %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.525 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.5   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.475 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.45  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.425 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.4   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.375 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.35  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.325 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.3   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.275 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.25  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.225 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.2   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.175 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.15  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.125 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			                  sprintf(buf,"(DeltaM -0.1   %d %d) 0 %g;\n", findpupil,j,m_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.075 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.05  %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
			if(fine_pitch!=0) sprintf(buf,"(DeltaM -0.025 %d %d) 0 %g;\n", findpupil,j,s_weight); s+=buf;
		}
	}

	return s;
}

int cLens1::SCA(double& Sph,double& Cyl,double& Axis, 
                double yObj,double xObj,double yPupil,double xPupil,double defocus,int i,int j) {
	// 1<=i<=k�̂Ƃ��͑�i�ʁC����ȊO�ł͑��ʂɂ��Čv�Z����
	// �߂�l�F�����ǐՂ̐���=1�^���s=0
	matrix<double> A(3,3), T(3,3);
	vector<double> X,B, ha,hb,dQa,dQb;
	double a,b,c, sph,cyl,axis;

	Sph=Cyl=Axis=0;
	if( RayTrace(yObj,xObj,yPupil,xPupil,defocus,j,0,0,0,1) ) return 0;
	if(i<1 || i>k) i=k+1;
	T=Tmatrix(vector<double>(this->X1[i],this->Y1[i],this->Z1[i]));
	ha=T*this->ha1[i]; hb=T*this->hb1[i]; dQa=T*this->dQa1[i]; dQb=T*this->dQb1[i];
	if( !Afocal || i<k+1 ) {
		// h = -SdQ - C(dQ,A)A, A=(cos��,sin��) (��:cyl������)  ..�y���zS,C�̓p���[�ł͂Ȃ�������(dQ= .. �łȂ� h= .. �ɂȂ�)
		//     A�ɕ��s : h�a=-(S+C)dQ�a ->  t=-h�a/dQ�a=S+C
		//     A�ɐ��� : h��=-SdQ��     ->  t=-h��/dQ��=S
		// ����𐬕��ɕ����ď����ƁC
		//   hx = -SdQx - C(dQxcos��+dQysin��)cos�� = -(S+C(cos��)^2)dQx - Csin��cos��dQy
		//   hy = -SdQy - C(dQxcos��+dQysin��)sin�� = -Csin��cos��dQx - (S+C(sin��)^2)dQy
		// �ƂȂ�D
		//   a=-S-C(cos��)^2
		//   b=-Csin��cos��
		//   c=-S-C(sin��)^2
		// �Ƃ���΁C
		//   hx = a*dQx + b*dQy
		//   hy = b*dQx + c*dQy
		// �ƂȂ�C������a,b,c�����߂�D
		// ����ɁCC<0�Ƃ��� S,C,�Ƃ����肳���D
		A.a[1][1]=dQa.x; A.a[1][2]=dQa.y; A.a[1][3]=0;
		A.a[2][1]=0;     A.a[2][2]=dQa.x; A.a[2][3]=dQa.y;
		A.a[3][1]=dQb.x; A.a[3][2]=dQb.y; A.a[3][3]=0;
		B=vector<double>(ha.x,ha.y,hb.x);
		if(rank(A)<3){  // dQa.y<<1 �̂Ƃ����A�͐����łȂ�
			A.a[1][1]=dQa.x; A.a[1][2]=dQa.y; A.a[1][3]=0;
			A.a[2][1]=dQb.x; A.a[2][2]=dQb.y; A.a[2][3]=0;
			A.a[3][1]=0;     A.a[3][2]=dQb.x; A.a[3][3]=dQb.y;
			B=vector<double>(ha.x,hb.x,hb.y);
			// �����ɂ��ĉ�]�Ώ̂Ȍn�̏ꍇ�C����A��
			//   | 0 1 0 |
			//   | 1 0 0 |
			//   | 0 1 0 |  �̂悤�ȍs��ƂȂ邽�ߐ������͎ア�i1,3�s�����s�j�D
		}
		X=inv(A)*B;
		a=X.x; b=X.y; c=X.z;
		cyl=-sqrt( (a-c)*(a-c)+4*b*b );  // cyl<0�Ƃ���
		sph=-(a+c+cyl)/2;
		axis= cyl==0 ? 0 : 0.5*atan2(-2*b/cyl,-(a-c)/cyl);

		if(i==k+1){
			// ���ʂ̂Ƃ���Z���ɉ�����������Ԃ��D
			// (T���悶�Ă���,�܂�,��_�����ǐՂ�Z��������Ƃ��Ă���̂�
			//  sph,cyl�̕�����z�����ł͂Ȃ��CZ����(��������)�ɂ�邱�Ƃɒ��ӁD)
			Sph=sph*this->Z1[i];
			Cyl=cyl*this->Z1[i];
		}
		else{
			// ���ʈȊO�̂Ƃ��͑�i�ʂ���̌����ɉ�����������Ԃ�
			// �����ɉ�����decenter_type=3(�܂�Ȃ��~���[�Ȃ�)�ł��܂�Ӗ��̂Ȃ����ʂƂȂ�
			// defocus�͖���
			Sph=sph;
			Cyl=cyl;
		}

		Axis=axis*180/PI;
		
		if(Cyl>0){
			Sph=Sph+Cyl;
			Cyl=-Cyl;
			Axis= Axis>0 ? Axis-90 : Axis+90;
		}
	}
	else { 
		// if( Afocal && i==k+1 )
		//
		// dQ = -Sh - C(h,a)a, a=(cos��,sin��) (��:cyl������)
		// 
		// Sph�̕����͎������̂Ƃ����Ƃ���D
		// (���x�␳�����̌��w�@��ł͊�Ɏ�������������(����(���ዾ�v)�ŏœ_������)���Ƃɂ��킹���D
	    //      ��: SLO�Ŏ��x�␳�����̎��͎�����
		//      ��: �ڊᎋ�x�␳����(�����Y���������񂾏��)�̂Ƃ��͐ڊႩ�甭�U�����o�Ă���
		//  �t�ɖڂ̑����猩��ƁC�ڂ�����𕨑̂Ƃ��������������ďo���Ԃ͋ߎ��ł���C
		//  �}�C�i�X�����Y�ŕ␳����邩��C���ܗ͕͂��ƂȂ�C����Ƃ͋t�D
		//  �܂��C�ڊᎋ�x�ŃX���b�g�����v�̑��ʘp�Ȃ𑪒肷��Ƃ��C���������o�Ă���ΐڊ����������
		//  (���x-����)�̂ŁC���̑��茋�ʂɂȂ�C������t�D)
		//
		// ���� i=0 �̂Ƃ��́C����(k+1��)����ƂȂ�D
		// �A�t�H�[�J���̂Ƃ��́Cdk()��0�ł��邽�߁C
		// ���ʂ͍ŏI��(k��)�̕ΐS��̍��W�ƈ�v����D
		// ���������āCdecentertype(k)=1 �̂Ƃ��� =2
		// �̂Ƃ��ł͊�ʂ��قȂ�C���Ȃ킿���ʂ��قȂ邱�Ƃɒ��ӁD
		//
		// �܂��C���ʂ�z���ł͂Ȃ������ɉ������l�ł���D
		A.a[1][1]=ha.x; A.a[1][2]=ha.y; A.a[1][3]=0;
		A.a[2][1]=0;    A.a[2][2]=ha.x; A.a[2][3]=ha.y;
		A.a[3][1]=hb.x; A.a[3][2]=hb.y; A.a[3][3]=0;
		B=vector<double>(dQa.x,dQa.y,dQb.x);
		X=inv(A)*B;
		a=X.x; b=X.y; c=X.z;
		cyl=-sqrt( (a-c)*(a-c)+4*b*b );
		sph=-(a+c+cyl)/2;
		axis= cyl==0 ? 0 : 0.5*atan2(-2*b/cyl,-(a-c)/cyl);
		
		Sph=sph*1000-defocus;
		Cyl=cyl*1000;
		Axis=axis*180/PI;
	}
	return 1;
}

double cLens1::Sph(double yObj,double xObj,double defocus,int findpupil/*=0*/,int i/*=0*/) {
	double s,c,a;
	point pupil;
	if(FindRay(pupil,yObj,xObj,"principal",findpupil,1)) {
		SCA(s,c,a,yObj,xObj,pupil.y,pupil.x,defocus,i,1);
		return s;
	}
	else return 0;
}

double cLens1::Sph(){
	return Sph(yObjectMax,xObjectMax,0);
}

double cLens1::Cyl(double yObj,double xObj,int findpupil/*=0*/,int i/*=0*/) {
	double s,c,a;
	point pupil;
	if(FindRay(pupil,yObj,xObj,"principal",findpupil,1)) {
		SCA(s,c,a,yObj,xObj,pupil.y,pupil.x,0,i,1);
		return c;
	}
	else return 0;
}

double cLens1::Cyl(){
	return Cyl(yObjectMax,xObjectMax);
}
double cLens1::Cyl90(){
	return Cyl(yObjectMax*0.9,xObjectMax*0.9);
}
double cLens1::Cyl70(){
	return Cyl(yObjectMax*0.7,xObjectMax*0.7);
}
double cLens1::Cyl50(){
	return Cyl(yObjectMax*0.5,xObjectMax*0.5);
}

double cLens1::CylSgn(double yObj,double xObj,int findpupil/*=0*/,int i/*=0*/){
	// Cyl()��0�t�߂ŕs�A���̂��ߍœK���̕]���֐��Ƃ��Ďg���Ȃ����C
	// �{�֐��͕����t�ł���D
	double axis;

	axis=Axis(yObj,xObj,findpupil,i);
	if(0<=axis && axis<90){  // Axis�̖߂�l�� -90�`+90deg
		return Cyl(yObj,xObj,i);
	}
	else{
		return -Cyl(yObj,xObj,i);
	}
}

double cLens1::Cyl0Deg(double yObj,double xObj,int findpupil/*=0*/,int j/*=0*/){
	// �c�������̃N���X�V�����_�[������Ԃ�
	double cyl,axis;

	cyl=Cyl(yObj,xObj,findpupil,j);
	axis=Axis(yObj,xObj,findpupil,j);
	return cyl*cos(2*axis*PI/180);
}
double cLens1::Cyl45Deg(double yObj,double xObj,int findpupil/*=0*/,int j/*=0*/){
	// 45�������̃N���X�V�����_�[������Ԃ�
	double cyl,axis;

	cyl=Cyl(yObj,xObj,findpupil,j);
	axis=Axis(yObj,xObj,findpupil,j);
	return cyl*sin(2*axis*PI/180);
}

double cLens1::Axis(double yObj,double xObj,int findpupil/*=0*/,int i/*=0*/) {
	double s,c,a;
	point pupil;
	if(FindRay(pupil,yObj,xObj,"principal",findpupil,1)) {
		SCA(s,c,a,yObj,xObj,pupil.y,pupil.x,0,i,1);
		return a;
	}
	else return 0;
}

double cLens1::Axis(){
	return Axis(yObjectMax,xObjectMax);
}

double cLens1::SphEq(double yObj,double xObj,double defocus,int findpupil/*=0*/,int i/*=0*/) {
	// �������ʓx�� Spherical Equivalent
	return Sph(yObj,xObj,defocus,findpupil,i)+Cyl(yObj,xObj,findpupil,i)/2;
}

double cLens1::SphEq(){
	return SphEq(yObjectMax,xObjectMax,0);
}

double cLens1::SphEq70(){
	return SphEq(yObjectMax*0.7,xObjectMax*0.7,0);
}

double cLens1::SphEq50(){
	return SphEq(yObjectMax*0.5,xObjectMax*0.5,0);
}

double cLens1::Mskew(char xy,double yObj,double xObj,double yPupil,double xPupil,int i1,int i2,int j,
                     int ObjImgTilt/*=1*/){
	// x(y)�����ɋ߂����������Ȕ��������̒��� l, ���̑��̒��� l' �̂Ƃ� Mskew=l'/l
	// i1,i2����1�ʁC�ŏI�ʂ̂Ƃ��́Cl,l'�̕��̖ʁC���ʂւ̓��e��L,L'�Ƃ��CMskew=L'/L�Ƃ���D
	// �������C�������������Ȃ��Ƃ�
	//�i�ߖT�����������Ŏ�o�������ɂȂ������ƌ����Ȃ��Ƃ��C
	//  �Ⴆ�Ύ�o��������xy�����łȂ��~�������Y)
	// �́C�s���m�Ȍ��ʂƂȂ�i�{���{�����`�ł��Ȃ��j�D
	// ObjImgTilt���^�̂Ƃ��́C�����ƕ��̖ʂ܂��͑��ʂ������łȂ����Ƃ̕␳���s���iScheimpflug���Ȃǁj�D
	
	double m;
	vector<double> dQ_i1,dQ1_i2,h1_i2,h_i1, dQ1_o,dQ_i;

	if( 1<=i1 && i1<=i2 && i2<=k ){
		RayTrace(yObj,xObj,yPupil,xPupil,0,j,0,0,0,1);

		if(xy=='x'){
			dQ_i1  =dQb[i1];
			dQ1_i2 =dQb1[i2];
			h_i1   =hb[i1];
			h1_i2  =hb1[i2];
			dQ1_o  =dQb1[0];
			dQ_i   =dQb[k+1];
		}
		else if(xy=='y'){
			dQ_i1  =dQa[i1];
			dQ1_i2 =dQa1[i2];
			h_i1   =ha[i1];
			h1_i2  =ha1[i2];
			dQ1_o  =dQa1[0];
			dQ_i   =dQa[k+1];
		}
		else{
			return 0;
		}

		if( abs(dQ_i1)==0 && abs(dQ1_i2)==0 ){
			m=abs(h1_i2)/abs(h_i1);
		}
		else{
			m=fabs(N(i1-1,j)/N(i2,j))*abs(dQ_i1)/abs(dQ1_i2);
		}

		if(ObjImgTilt!=0 && fabs(s)<LN && i1==1 && i2==k){
			// ���̖ʁC���ʂƌ����������łȂ��ꍇ�̕␳�iSheimpflug�J�����Ȃǁj
			double cs,cs1;
			vector<double> v,n, v1,n1;

			v=vector<double>(X1[0],Y1[0],Z1[0]);
			n=surface_normal(0,y[0],x[0],1);
			v1=dQ1_o-sProduct(dQ1_o,n)*n;     // dQ1�̖ʖ@���ɐ�������(�ʓ�)�̐���
			n1=dQ1_o-sProduct(dQ1_o,v)*v;     // dQ1�̌����ɐ��������̐���
			cs=cosine(v1,n1); 

			v=vector<double>(X[k+1],Y[k+1],Z[k+1]);
			n=surface_normal(k+1,y[k+1],x[k+1],1);
			v1=dQ_i-sProduct(dQ_i,n)*n;
			n1=dQ_i-sProduct(dQ_i,v)*v;
			cs1=cosine(v1,n1);

			m*=cs/cs1;			
		}

		return m;
	}
	else return 0;
}

double cLens1::My(double yObj,double xObj,double yPupil,double xPupil,int i1,int i2,int j){
	return Mskew('y',yObj,xObj,yPupil,xPupil,i1,i2,j);
}

double cLens1::My(int i1,int i2){
	point ymax,ymin,xmax,xmin, pupil;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax,0,1,1)){
		pupil=(ymax+ymin)/2;
		return My(yObjectMax,0,pupil.y,pupil.x,i1,i2,1);
	}
	else return 0;
}

double cLens1::My(){
	return My(1,k);
}

double cLens1::Mx(double yObj,double xObj,double yPupil,double xPupil,int i1,int i2,int j){
	return Mskew('x',yObj,xObj,yPupil,xPupil,i1,i2,j);
}

double cLens1::Mx(int i1,int i2){
	point ymax,ymin,xmax,xmin, pupil;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax,0,1,1)){
		pupil=(ymax+ymin)/2;
		return Mx(yObjectMax,0,pupil.y,pupil.x,i1,i2,1);
	}
	else return 0;
}

double cLens1::Mx(){
	return Mx(1,k);
}

double cLens1::Mpskew(char xy,int i1,int i2,int j){
	cLens1 X=*this;

	Swap(X.s,X.t);
	return X.Mskew(xy,0,0,0,0,i1,i2,j,0);
}

double cLens1::Mpy(int i1,int i2){
	return Mpskew('y',i1,i2,1);
}
double cLens1::Mpy(){
	return Mpy(1,k);
}

double cLens1::Mpx(int i1,int i2){
	return Mpskew('x',i1,i2,1);
}
double cLens1::Mpx(){
	return Mpx(1,k);
}

double cLens1::fskew(char xy,int j){
	// �����ɉ�������_�����ǐՂɂ��Cx,y�����̏œ_������Ԃ�	
	cLens1 X=*this;

	X.s=LN*10;
	X.RayTrace(0,0,0,0,0,j,0,0,0,1);

	if(xy=='x')      return 1/( fabs(X.N(k,j)) * abs(X.dQb1[k]) );
	else if(xy=='y') return 1/( fabs(X.N(k,j)) * abs(X.dQa1[k]) );
	else             return 0;
}

double cLens1::fx(){
	return fskew('x',1);
}
double cLens1::fy(){
	return fskew('y',1);
}

double cLens1::S1skew(char xy,double yObj,double xObj,double yPupil,double xPupil,int i,int j){
	// x(y)�����̋ߖT�����̑�i�ʌ�̌����ʒu(�����ɉ�������i�ʂ���̋����Ƃ���)�����߂�D
	// �ߖT���������̌����ƌ����Ȃ��Ƃ�(��]��Ώ̂Ȍn�Ȃ�)�͌��̌����ƋߖT�������ł��ߐڂ���ʒu�Ƃ���D
	vector<double> Q1,dQ1_i,h1_i, s1;

	if( 1<=i && i<=k+1 ){
		RayTrace(yObj,xObj,yPupil,xPupil,0,j,0,0,0,1);

		if(xy=='x'){
			dQ1_i =dQb1[i];
			h1_i  =hb1[i];
		}
		else if(xy=='y'){
			dQ1_i =dQa1[i];
			h1_i  =ha1[i];
		}
		else{
			return 0;
		}
		Q1=vector<double>(X1[i],Y1[i],Z1[i]);
		s1=NearestPointLineLine(h1_i,Q1+dQ1_i, vector<double>(0,0,0),Q1);
		// �� ���ӁF h_i, dQ1_i, Q1 �͒��_���W�ɂ���Ă���.
		return sProduct(s1, Q1/abs(Q1));   // �����ɉ���������
	}
	else return 0;
}

double cLens1::S1y(double yObj,double xObj,double yPupil,double xPupil,int i,int j){
	return S1skew('y',yObj,xObj,yPupil,xPupil,i,j);
}

double cLens1::S1y(int i){
	point ymax,ymin,xmax,xmin, pupil;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax,0,1,1)){
		pupil=(ymax+ymin)/2;
		return S1y(yObjectMax,0,pupil.y,pupil.x,i,1);
	}
	else return 0;
}

double cLens1::S1y(){
	return S1y(k);
}

double cLens1::S1x(double yObj,double xObj,double yPupil,double xPupil,int i,int j){
	return S1skew('x',yObj,xObj,yPupil,xPupil,i,j);
}

double cLens1::S1x(int i){
	point ymax,ymin,xmax,xmin, pupil;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax,0,1,1)){
		pupil=(ymax+ymin)/2;
		return S1x(yObjectMax,0,pupil.y,pupil.x,i,1);
	}
	else return 0;
}

double cLens1::S1x(){
	return S1x(k);
}

double cLens1::Dist(double yObj,double yPupil,double xPupil,int j){
	double ybar;
	double now_s1fix=s1fix; int now_Afocal=Afocal; // �ߎ����ʂłȂ��Ɨ��z������ݒ�ł��Ȃ��D
	s1fix=Afocal=0;
	RayTrace(yObj,0, yPupil,xPupil, 0,j,0,0,0,0);
	s1fix=now_s1fix; Afocal=now_Afocal;
	if(abs(s)>=LN) ybar=yObj*f(j); else ybar=yObj*M(j);
	if(ybar!=0) return (y[k+1]-ybar)/ybar; else return 0;
}
double cLens1::Dist(double yObjNormalized,int findpupil){
	point ymax,ymin,xmax,xmin, pupil;

	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax*yObjNormalized,0,1,findpupil,1)) {
		pupil=(ymax+ymin)/2;
		return Dist(yObjectMax*yObjNormalized,pupil.y,pupil.x,1);
	}
	else return 0;
}
double cLens1::Dist(){
	point ymax,ymin,xmax,xmin, pupil;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax,0,1,1)) {
		pupil=(ymax+ymin)/2;
		return Dist(yObjectMax,pupil.y,pupil.x,1);
	}
	else return 0;
}
double cLens1::Dist70(){
	point ymax,ymin,xmax,xmin, pupil;
	if(FindPupilEdge(ymax,ymin,xmax,xmin,yObjectMax,0,1,1)) {
		pupil=(ymax+ymin)/2;
		return Dist(yObjectMax*0.7,pupil.y,pupil.x,1);
	}
	else return 0;
}

double cLens1::DeltaH(int DYIs1DXIs2,
					 double yObj,double xObj,double yPupil,double xPupil,
	                 double yPupilPrincipal,double xPupilPrincipal,double defocus,int j,int j0,
					 int AbeTheory/*=0*/){
	double yhp,xhp,zhp,yh,xh,zh;

	if(AbeTheory==0){
		make_coordinate(defocus);
		// ��j0�g���̎����(xPupilPrincipal,yPupilPrincipal)�̑�����(xhp,yhp,zhp)����Ƃ���
		ImageHeight(yhp,xhp,zhp,yObj,xObj,"",yPupilPrincipal,xPupilPrincipal,defocus,j0,0,0,0,1,0);
		// ��j�g���̌���(xPupil,yPupil)�̑�����(xh,yh,zh)
		ImageHeight(yh,xh,zh,yObj,xObj,"",yPupil,xPupil,defocus,j,0,0,0,1,0);
	}
	else{
		ImageHeightAbe(yhp,xhp,yObj,xObj,"",yPupilPrincipal,xPupilPrincipal,defocus,j0,AbeTheory);
		ImageHeightAbe(yh,xh,yObj,xObj,"",yPupil,xPupil,defocus,j,AbeTheory);
		zhp=surface_sag(k+1,yhp,xhp,1);
		zh=surface_sag(k+1,yh,xh,1);
	}
	
	if(zh!=zhp){
		// ���ʁi���ȂǋȖʂ����蓾��j�̐ڕ��ʂ�x,y�ʂƂ���
		matrix<double> T,Rp(3,1),R(3,1);

		Rp[1][1]=xhp; Rp[2][1]=yhp; Rp[3][1]=zhp;
		R[1][1]=xh;   R[2][1]=yh;   R[3][1]=zh;
	
		T=Tmatrix(surface_normal(k+1,yhp,xhp,1));
		Rp=T*Rp;
		R= T*R;

		xhp=Rp[1][1]; yhp=Rp[2][1];
		xh=R[1][1];   yh=R[2][1];
	}

	switch(DYIs1DXIs2){
		case 1:  return yh-yhp;
		case 2:  return xh-xhp;
		default: return 0;
	}
}

double cLens1::DeltaH(int DYIs1DXIs2,
					 double yObj,double xObj,double yPupil,double xPupil,
	                 double yPupilPrincipal,double xPupilPrincipal,double defocus,int j,int j0,
					 int i1,int i2,int AbeTheory/*=0*/){
	// i1�ʂ���i2�ʂ̕����̉�������Ԃ��D
	//     �E�ߎ����{���ɂ��C����ԂɊ��Z�����l��Ԃ��D
	//     �Ei1==0����i2==0�̂Ƃ��́C�S�n�̉�������Ԃ��D
	double dh;

	if(1<=i1 && i1<=i2 && i2<=k){
		if(i2==k){
			dh=DeltaH(DYIs1DXIs2,yObj,xObj,yPupil,xPupil,yPupilPrincipal,xPupilPrincipal,defocus,j,j0,AbeTheory);
		}
		else{
			dh=Trimed(1,i2)
				.DeltaH(DYIs1DXIs2,yObj,xObj,yPupil,xPupil,yPupilPrincipal,xPupilPrincipal,defocus,j,j0,AbeTheory)
				*M(i2+1,k,1);
		}

		if(i1>1){
			dh-=Trimed(1,i1-1)
				.DeltaH(DYIs1DXIs2,yObj,xObj,yPupil,xPupil,yPupilPrincipal,xPupilPrincipal,defocus,j,j0,AbeTheory)
				*M(i1,k,1);
		}
	}
	else if(i1==0 && i2==0){
		dh=DeltaH(DYIs1DXIs2,yObj,xObj,yPupil,xPupil,yPupilPrincipal,xPupilPrincipal,defocus,j,j0,AbeTheory);
	}
	else{
		dh=0;
	}

	return dh;
}

double cLens1::DeltaH(int DYIs1DXIs2,double yObj,double xObj,
                      int FindPupil,double ypNormalized,double xpNormalized,double defocus,int j,int j0){
	// (ypNormalized,xpNormalized) : ���̔��a��1�Ƃ��������W
	point ymax,ymin,xmax,xmin, p,p0;

	make_coordinate(defocus);
	FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,j,FindPupil,0);
	p0=(ymax+ymin)/2;
	p=p0 +(ymax-ymin)/2*ypNormalized +(xmax-xmin)/2*xpNormalized;
	return DeltaH(DYIs1DXIs2,yObj,xObj,p.y,p.x,p0.y,p0.x,defocus,j,j0);
}

std::string cLens1::SPO(double yObj,double xObj,int findpupil,int colors,
						int IgnoreTC/*=0*/,double weight/*=1*/,int yfan/*=0*/){
	// optimize()��rms�X�|�b�g�a�ł͂Ȃ��e�����̃�x,��y��]���֐��Ƃ���ƁC
	// optimize()�̈���targets������������ƂȂ�C���͂��ʓ|�ł���D
	// �{�֐��ɂ�肱�̒���������������������D
	//   �i �e�����̃�x,��y��]���֐��Ƃ���̂́C����O�̒�Ă���g�X�|�b�g�����h�ŁC 
	//      �Ⴆ�� rms�X�|�b�g�a���^�[�Q�b�g�Ƃ���ƍœK���͂����Ɏ~�܂��Ă��܂����C
	//      �X�|�b�g�����ł���Ό��w�n�͑傫�������ƌ����Ă���D
	//      rms�X�|�b�g�a���ƍœK�l�t�߂Ŕ����W����0�ł��邪�C��x,��y�ł���ΐ��������l�ł��邽��
	//      ��x,��y=0 �̕t�߂Ŕ����W����0�Ƃ͂Ȃ�Ȃ����Ƃ����R���D)
	//
	// 2014.04.25
	//  �E�{���F�����ɔz�����Ď�����̃E�G�C�g�� 0.1 -> 1 �ɕύX
	//  �E�Ώ̎���(���ʘp�ȂȂ�)�C���Ȃ킿���̗��[�ł̍���ǉ��D���ꂪ�Ȃ��Ɣ�Ώ̎���(�R�})
	//    �ɔ�ׂāC�Ώ̎������}������ɂ����D�Ⴆ�΁Cdy=-1�`+1�̑Ώ̎�����dy=0�`+1�̔�Ώ�
	//    �����̔{�̃X�|�b�g�T�C�Y�����C�����b�g�֐��͓������Ȃ��Ă��܂��Ă����D
	// 2016.06.19
	//   ������ IgnoreTC (�{���F�����𖳎��j��ǉ��DSLO�Ń\�t�g�E�G�A�ɂ��{���F�����␳���\�ȂƂ��Ɏg����D
	// 2017.05.09
	//    ������ yfan ��ǉ��D�^�̂Ƃ���yz�ʓ��̂�(y�t�@���̂�)�]������D
	//    ���R�Ȗʃ~���[�n�ł܂��q�ߖʓ��̂ݐ݌v����Ƃ��Ɏg����D
	int j,j0;
	std::string s;
	char buf[1000];
	double w;

	if(colors>cn) colors=cn;
	for(j=1; j<=colors; j++){
		j0= IgnoreTC ? j : 1;
		w=colorweight[j]*weight;
		// �ȉ���y
		sprintf(buf,"(DeltaH 1 %g %g %d 1    0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w     ); s+=buf;
		sprintf(buf,"(DeltaH 1 %g %g %d .95  0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.95); s+=buf;
		sprintf(buf,"(DeltaH 1 %g %g %d .9   0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.9 ); s+=buf;
		sprintf(buf,"(DeltaH 1 %g %g %d .8   0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.8 ); s+=buf;
		sprintf(buf,"(DeltaH 1 %g %g %d .7   0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.7 ); s+=buf;
		sprintf(buf,"(DeltaH 1 %g %g %d .6   0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.6 ); s+=buf;
		sprintf(buf,"(DeltaH 1 %g %g %d .5   0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.5 ); s+=buf;
		sprintf(buf,"(DeltaH 1 %g %g %d .4   0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.4 ); s+=buf;
		sprintf(buf,"(DeltaH 1 %g %g %d .3   0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.3 ); s+=buf;
		if(!IsRotationallySymmetric() || yObj!=0){  // ��]�Ώ̌n����yObj=0�̂Ƃ��͈ȏ�ŏ\��
			sprintf(buf,"(DeltaH 1 %g %g %d 0    0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w     ); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -.3  0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.3 ); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -.4  0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.4 ); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -.5  0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.5 ); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -.6  0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.6 ); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -.7  0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.7 ); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -.8  0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.8 ); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -.9  0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.9 ); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -.95 0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.95); s+=buf;
			sprintf(buf,"(DeltaH 1 %g %g %d -1   0 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w     ); s+=buf;
			//    ��y�̑Ώ̎���
			sprintf(buf,"{- (DeltaH 1 %g %g %d 1 0 0 %d %d) (DeltaH 1 %g %g %d -1 0 0 %d %d)} 0 %g;\n", 
					yObj,xObj,findpupil,j,j0, yObj,xObj,findpupil,j,j0, w); s+=buf;

			// �ȉ���x
			if(yfan==0){ // yfan���^�ɐݒ肳��Ă��Ȃ��Ƃ�
				sprintf(buf,"(DeltaH 2 %g %g %d 0 1    0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w     ); s+=buf;
				sprintf(buf,"(DeltaH 2 %g %g %d 0 .95  0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.95); s+=buf;
				sprintf(buf,"(DeltaH 2 %g %g %d 0 .9   0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.9 ); s+=buf;
				sprintf(buf,"(DeltaH 2 %g %g %d 0 .8   0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.8 ); s+=buf;
				sprintf(buf,"(DeltaH 2 %g %g %d 0 .7   0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.7 ); s+=buf;
				sprintf(buf,"(DeltaH 2 %g %g %d 0 .6   0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.6 ); s+=buf;
				sprintf(buf,"(DeltaH 2 %g %g %d 0 .5   0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.5 ); s+=buf;
				sprintf(buf,"(DeltaH 2 %g %g %d 0 .4   0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.4 ); s+=buf;
				sprintf(buf,"(DeltaH 2 %g %g %d 0 .3   0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.3 ); s+=buf;
				if(!IsXAxisSymmetric() || xObj!=0){ // YZ�ʂɂ��đΏ̂���xObj=0 �̂Ƃ��͈ȏ�ŏ\��
					sprintf(buf,"(DeltaH 2 %g %g %d 0 0    0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w     ); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -.3  0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.3 ); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -.4  0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.4 ); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -.5  0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.5 ); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -.6  0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.6 ); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -.7  0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.7 ); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -.8  0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.8 ); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -.9  0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.9 ); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -.95 0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w*0.95); s+=buf;
					sprintf(buf,"(DeltaH 2 %g %g %d 0 -1   0 %d %d) 0 %g;\n", yObj,xObj,findpupil,j,j0,w     ); s+=buf;
					//    ��x�̑Ώ̎���
					sprintf(buf,"{- (DeltaH 2 %g %g %d 0 1 0 %d %d) (DeltaH 2 %g %g %d 0 -1 0 %d %d)} 0 %g;\n", 
							yObj,xObj,findpupil,j,j0, yObj,xObj,findpupil,j,j0, w); s+=buf;
				}
			}
		}
	}
	return s;
}

std::string cLens1::SPO2(double yObj,double xObj,int findpupil,int colors,int IgnoreTC/*=0*/,int n/*=0*/,double weight/*=1*/){
	// SPO()�Ƃقړ����@�\�����C������ł͓����W��XY����݂̂łȂ����S�̂Ɋi�q��ɔz�u����D
	// n : �i�q�� nxn. �������C���ۂ̊i�q�_����nxn��菭�Ȃ��D
	int i,j,j0;
	std::string s;
	char buf[1000];
	double w;

	if(n==0) n=10;   // �f�t�H���g�� n=10
	if(colors>cn) colors=cn;
	MakePupilGrid(yObj,xObj,findpupil,n,1,1,0,1,1);

	for(j=1; j<=colors; ++j){
		j0= IgnoreTC ? j : 1;
		w=colorweight[j]*weight;

		for(i=1; i<=this->gpoints.GetSize(); ++i){
			sprintf(buf,"(DeltaH 1 %g %g %d %g %g 0 %d %d) 0 %g;\n",
			        yObj,xObj,findpupil,gpoints[i].y,gpoints[i].x,j,j0,w);
			s+=buf;
			sprintf(buf,"(DeltaH 2 %g %g %d %g %g 0 %d %d) 0 %g;\n",
			        yObj,xObj,findpupil,gpoints[i].y,gpoints[i].x,j,j0,w);
			s+=buf;
		}
	}
	return s;
}

double cLens1::DYpr(double yObj,int j,int FindPupil){
	// ������̉��F����
	// = ��j�g����������� - ��g�����������
	point ymax,ymin,xmax,xmin, pr;

	if(FindPupilEdge(ymax,ymin,xmax,xmin, yObj,0, 1,FindPupil,1)) {
		pr=(ymax+ymin)/2;
		return DeltaH(1,yObj,0, pr.y,pr.x, pr.y,pr.x, 0,j,1);	
		// DeltaH()�ł͎�g�������������0�Ƃ��Ă���D
	}
	else return 0;
}
double cLens1::DYpr(int j,int FindPupil){
	return DYpr(yObjectMax,j,FindPupil);
}
double cLens1::DYpr75(int j,int FindPupil){
	return DYpr(yObjectMax*0.75,j,FindPupil);
}
double cLens1::DYpr50(int j,int FindPupil){
	return DYpr(yObjectMax*0.5,j,FindPupil);
}

double cLens1::DYymax(double yObj,int FindPupil/*=1*/){
	point ymax,ymin,xmax,xmin, pr;
	if(FindPupilEdge(ymax,ymin,xmax,xmin, yObj,0, 1,FindPupil,1)) {
		pr=(ymax+ymin)/2;
		return DeltaH(1,yObj,0, ymax.y,ymax.x, pr.y,pr.x, 0,1,1);
	}
	else return 0;
}
double cLens1::DYymax(){
	return DYymax(this->yObjectMax);
}

double cLens1::DYymin(double yObj,int FindPupil/*=1*/){
	point ymax,ymin,xmax,xmin, pr;
	if(FindPupilEdge(ymax,ymin,xmax,xmin, yObj,0, 1,FindPupil,1)) {
		pr=(ymax+ymin)/2;
		return DeltaH(1,yObj,0, ymin.y,ymin.x, pr.y,pr.x, 0,1,1);
	}
	else return 0;
}
double cLens1::DYymin(){
	return DYymin(this->yObjectMax);
}

double cLens1::DXxmax(double yObj){
	point ymax,ymin,xmax,xmin, pr;
	if(FindPupilEdge(ymax,ymin,xmax,xmin, yObj,0, 1,1)) {
		pr=(ymax+ymin)/2;
		return DeltaH(2,yObj,0, xmax.y,xmax.x, pr.y,pr.x, 0,1,1);
	}
	else return 0;
}
double cLens1::DXxmax(){
	return DXxmax(this->yObjectMax);
}

double cLens1::DYunsymmetric(double yObj,int FindPupil){
	// ��Ώ̎���������Ԃ�
	return ( DYymax(yObj,FindPupil)+DYymin(yObj,FindPupil) )/2 - DYpr(yObj,1,FindPupil);
}

int cLens1::ImageHeight(double& yh,double& xh,double& zh,double yObj,double xObj,
                        std::string SetRay,double yPupil,double xPupil,
					    double defocus,int j,int EA_enable,int mask_enable,int IsLambert,int InAngle,
						int MakeCoordinate) 
{
	// EA_enable,mask_enable,IsLambert : �X�|�b�g�_�C�A�O�����쐬�����̂��߂ɕK�v
	//
	point p;
	p=point(xPupil,yPupil);
	if( !FindRay(p,yObj,xObj,SetRay,1,j) ) return 0;
	if(RayTrace(yObj,xObj,p.y,p.x,defocus,j,EA_enable,mask_enable,IsLambert,0,0,vector<complex>(),MakeCoordinate)){
		return 0;
	}
	else{
		// ���̖ʖ������̂Ƃ��C
		//     InAngle=False : �����N * Y(�܂���X)/Z                    (���Ȃ킿Ntan��. ����͗��z�����Y�̑����ɑΉ�����)
		//     InAngle=True  : �����N * Y(�܂���X)/Z/AfocalImageUnit()  (���w�I�ɂ͌����Ɋp�x�ł͂Ȃ�)
		// ��Ԃ��D(X,Y,Z)�͌����̕����D
		//
		// ���z�����Y��̑����́C�g�􉽌��w(�O��)�h���C
		//      y' = -f(y/z) (=-f*tan��) (3.6)
		//      n'/f' = -n/f (=fideal)   (3.17)
		// ���������āC
		//      y' = fideal*n*tan��
		// �ł���D�܂��C���̎���
		//    ���ܗ������̋�Ԃł͌�����+z�����֐i�s
		//    ���ܗ������̋�Ԃł͌�����-z�����֐i�s
		// ���O��ƂȂ�D
		if(this->Afocal){
			// ����(���ʂ����̐i�s��)�̂Ƃ�defocus>0�Ƃ���ƁC
			// defocus�������ʂ܂ł̋��� = 1000/defocus*sgn(Z[k+1])
			// (���x�␳�@�\���t�������w�@��͎��x�␳�l>0�̂Ƃ���ɑ΂��Ď�����������D)
			yh=N(k,j)*( Y[k+1]/Z[k+1]+y[k+1]*defocus/1000*sgn(Z[k+1]) )/( InAngle ? AfocalImageUnit() : 1 );
			xh=N(k,j)*( X[k+1]/Z[k+1]+x[k+1]*defocus/1000*sgn(Z[k+1]) )/( InAngle ? AfocalImageUnit() : 1 );
			zh=0;
		}
		else{
			yh=y[k+1];
			xh=x[k+1];
			zh=z[k+1];
		}
		return 1;
	}
}

int cLens1::ImageHeightAbe(double& yh,double& xh,double yObj,double xObj,
                           std::string SetRay,double yPupil,double xPupil,double defocus,int j,int order/*=3*/) 
{
	// �����ǐՂł͂Ȃ��C�����W���ɂ����������v�Z����
	std::string normalizeunit0;
	int normalizetype0;
	point pupil,img;
	double g,g_hat, th,Y, y0,al1;
	double tw,R,phi;   // ���ꂼ�� "�����Y�݌v�@" (4.2)(4.4)���� Ntan��,R,��

	pupil=point(xPupil,yPupil);
	if( !FindRay(pupil,yObj,xObj,SetRay,1,1) ) return 0;
	
	normalizeunit0=NormalizeUnit; normalizetype0=NormalizeType;  // ���̒l��ۑ��ipush,pop���g���ƒx���Ȃ�)
	NormalizeUnit="1";   // �X�P�[���̋K�i���͂��Ȃ�
	NormalizeType=1;     // �g���K��1�h�Ƃ��C�g�����Y�݌v�@�h�� (4.2)(4.4) ��p����D

	CalcCoefficients(j); // �F�����͑�j�g�� vs ��1�g��

	if(fabs(s)<LN){
		g=s-t;
		g_hat=this->g_hat();
		pupil.y*=g_hat/g; pupil.x*=g_hat/g;  // �����W�͎啽�ʏ�Ƃ���
	}
		
	th=PI/2-atan2(yObj,xObj);  // ���_�̕��ʊp�Ƃ��������Ƃ̍�th

	Y=sqrt(yObj*yObj+xObj*xObj);
	if(fabs(s)<LN){
		tw=N(0,1)*Y/g_hat;     // ���̂�th�����񂵂�y����� y>0 �Ŕz�u����
	}
	else{
		tw=N(0,1)*Y;
	}

	pupil.rotate(th*180/PI);          // �����W��th������
	R=abs(pupil);
	phi=PI/2-atan2(pupil.y,pupil.x);  // �g�����Y�݌v�@�h �}4.1���Q��
	
	if(fabs(s)<LN) y0=Y*M(); else y0=Y*f();  // �ߎ������i���Ȃ݂�x�ߎ�������0)

	// �g�����Y�݌v�@�h(4.2) �ɂ�葜�ʒu��3�������W������v�Z����
	al1=N(k,1)*u[k+1];
	yh = y0 -(0.5/al1)*( SAt*R*R*R*cos(phi) +CMt*tw*R*R*(2+cos(2*phi)) 
		                +(3*ASt+PTt)*tw*tw*R*cos(phi) +DSt*tw*tw*tw 
						+LCt*R*cos(phi) +TCt*tw );
	xh = -(0.5/al1)*( SAt*R*R*R*sin(phi) +CMt*tw*R*R*sin(2*phi) +(ASt+PTt)*tw*tw*R*sin(phi)
		             +LCt*R*sin(phi) );

	if(order==5){
		// 5��������ǉ�����D�g�����Y�݌v�@�h(4.4)��� �D
		yh += -(0.5/al1)*( 0.25*SA5t*R*R*R*R*R*cos(phi)
			              +0.25*(CM41t*(3+2*cos(2*phi)+CM41Zt))*tw*R*R*R*R
						  +0.5*(SA32t+SA32Ft*(3+cos(2*phi))+2*SA32Zt)*tw*tw*R*R*R*cos(phi)
						  +0.5*(CM23t*(2+cos(2*phi))+CM23Pt*(1+cos(2*phi))+CM23Zt)*tw*tw*tw*tw*R*R
						  +0.25*(4*AS5t+SG5t)*tw*tw*tw*tw*R*cos(phi)
						  +0.25*DS5t*tw*tw*tw*tw*tw );
		xh += -(0.5/al1)*( 0.25*SA5t*R*R*R*R*R*sin(phi)
			              +0.5*CM41t*tw*R*R*R*R*sin(2*phi)
						  +0.5*(SA32t+SA32Ft*(1+cos(2*phi)))*tw*tw*R*R*R*sin(phi)
						  +0.5*CM23t*tw*tw*tw*R*R*sin(2*phi)
						  +0.25*SG5t*tw*tw*tw*tw*R*sin(phi) );
	}

	img=point(xh,yh); img.rotate(-th*180/PI); yh=img.y; xh=img.x; pupil.rotate(-th*180/PI); // th��]�����ɖ߂�

	if(s1fix!=0 || defocus!=0){ // �ߎ����ʂ���̃f�t�H�[�J�X������Ƃ�
		double g1_hat,dz;

		g1_hat=this->g1_hat();   // �啽�ʂ���ߎ����ʂ܂�
		dz=dk(defocus)-s1();
		yh=pupil.y + (g1_hat+dz)/g1_hat*(yh-pupil.y);
		xh=pupil.x + (g1_hat+dz)/g1_hat*(xh-pupil.x);
		// ���̎��́g�����_�h�� �}3.2.3 �̂悤�Ȑ}��`���Ɨ����ł���
	}

	NormalizeUnit=normalizeunit0; NormalizeType=normalizetype0;  // ���ɖ߂�
	return 1;
}

int cLens1::ImageHeight(point &height,const point& obj,const point& pupil,
                        double defocus,int j,int EA_enable,int mask_enable,int InAngle)
{
	double dummy;

	return ImageHeight(height.y,height.x,dummy,
	                   obj.y,obj.x,"",pupil.y,pupil.x,defocus,j,EA_enable,mask_enable,0,InAngle);
}

double cLens1::ImageHeightX(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
                            double defocus,int j,int InAngle)
{
	double yh,xh,zh;
	ImageHeight(yh,xh,zh,yObj,xObj,SetRay,yPupil,xPupil,defocus,j,0,0,0,InAngle);
	return xh;
}

double cLens1::ImageHeightY(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
                            double defocus,int j,int InAngle)
{
	double yh,xh,zh;
	ImageHeight(yh,xh,zh,yObj,xObj,SetRay,yPupil,xPupil,defocus,j,0,0,0,InAngle);
	return yh;
}

double cLens1::ImageHeightTh(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
                             double defocus,int j)
{
	// ���ʒu�̓��a������ rad �ŕԂ��DOTF�v�Z��S,M�����𓾂�̂Ɏg����D
	double x,y,z;

	ImageHeight(y,x,z,yObj,xObj,SetRay,yPupil,xPupil,defocus,j,0,0,0,0);
	if(x==0 && y==0){
		return 0;
	}
	else{
		return arg(complex(x,y));
	}
}

double cLens1::DImageHeight(double yObj,double xObj,
	                    std::string SetRay,double yPupil,double xPupil,
		                double defocus,int j,int EA_enable,int mask_enable,int IsLambert,int InAngle,
	                    double DyObj,double DxObj)
{
	point p0,p;
	double dummy;

	if(ImageHeight(p0.y,p0.x,dummy,yObj,xObj,SetRay,yPupil,xPupil,
	               defocus,j,EA_enable,mask_enable,IsLambert,InAngle)==0){
		return 0;
	}
	if(ImageHeight(p.y,p.x,dummy,yObj+DyObj,xObj+DxObj,SetRay,yPupil,xPupil,
	               defocus,j,EA_enable,mask_enable,IsLambert,InAngle)==0){
		return 0;
	}
	return distance(p0,p);
}

double cLens1::Get_x(int i){
	if(0<=i && i<=k+1) return x[i]; else return 0;
}

double cLens1::Get_y(int i){
	if(0<=i && i<=k+1) return y[i]; else return 0;
}

double cLens1::Get_z(int i){
	if(0<=i && i<=k+1) return z[i]; else return 0;
}

point cLens1::dObject(){
	double ry,rx;
	double dy,dx;

	if( SourcePhiY==0 ) return point(0,0);
	switch( SourceAreaType ){
	case 0:
		// ellipse
		ry= SourcePhiY/2;
		rx= SourcePhiX==0 ? SourcePhiY/2 : SourcePhiX/2;
		do{
			dy=Random(ry,0);
			dx=Random(rx,0);
		} while( dy*dy/ry/ry + dx*dx/rx/rx > 1 );
		break;
	case 1:
		// rectangle
		ry= SourcePhiY/2;
		rx= SourcePhiX==0 ? SourcePhiY/2 : SourcePhiX/2;
		dy=Random(ry,0);
		dx=Random(rx,0);
		break;
	case 2:
		// �K�E�X���z
		ry=SourcePhiY/2;
		rx=SourcePhiX==0 ? SourcePhiY/2 : SourcePhiX/2;
		dy=RandomGauss(ry);
		dx=RandomGauss(rx);
		break;
	default:
		dy=dx=0;
		break;
	}
	return point(dx,dy).rotate(SourceAxis);
}


int cLens1::MakeSpot(cSpot& spot,
					double yObj,double xObj,int FindPupil,double defocus,int footprint,
					int ColorStart,int ColorEnd,
					int IsAreaSource,int IsLambert,int OriginIsGravityCenter,int Add,int Randomize,
					int AutoInTangent/*=1*/) {
	
	// yObj(xObj)>=LN �̂Ƃ��� yObj(xObj)=yObjectMax(xObjectMax) �Ƃ���D

	// FindPupil % 10 ==0 : �����̋��E��EPphi�Ŏw�肷����˓��Ƃ���D
	// FindPupil % 10 !=0 : �����̋��E�͑S�Ă̗L���a��ʂ邱�Ƃ̂ł���ő�̌����̂���ł���D
	// FindPupil / 10 ==0 : �����S��
	// FindPupil / 10 !=0 : ���������̂�
	//    ��F FindPupil== 0 : EPphi�Ō��܂�����D�����S�́D
	//         FindPupil== 1 : �����Y�a�Ō��܂�����D�����S�́D
	//         FindPupil==10 : EPphi�Ō��܂�����D���������̂݁D
	//         FindPupil==11 : �����Y�a�Ō��܂�����D���������̂݁D

	// Randomize !=0 �Ƃ��邱�Ƃɂ��C�œ_���牓���f�t�H�[�J�X�ʒu�ł�OTF�v�Z�Ŋi�q�s�b�`�̉e����h���D
	// footprint >0 �̂Ƃ� �S�n��ʂ��������ő�footprint�ʂ̃t�b�g�v�����g���쐬����D
	//           <0 �̂Ƃ� ��|footprint|�ʂɒB�����S�Ă̌����Ńt�b�g�v�����g���쐬����D
	// AutoInTangent : �X�|�b�g�̑��ʐڕ��ʂ֓��e�������ōs�����ǂ���
	//
	// �߂�l : ���ʂ̃X�|�b�g�̌�����
	//
	// �y���z ���ʂ����ʂłȂ��Ƃ��́C�X�|�b�g�_�C�A�O��������ˈʒu�̐ڕ��ʂɓ��e����D
	
	if( ColorStart>cn || ColorEnd>cn ) return 0;

	int iy,ix,j, a, errcode, makecoordinate;
	double yPupil,xPupil, hy,hx,hz,hy0,hx0,hz0,hxp,hyp,hzp, th, yp0,xp0,dyp,dxp, ra,rb;
	point dr;
	int n=0;
	cPoint p;
	point ymax,ymin,xmax,xmin, yw,xw, principal;
	int all_ray;
	matrix<double> T;
	bool InTangent;

	if(!Add) spot.RemoveAll();

	all_ray= footprint<=-1 ? 1 : 0;
	footprint= footprint>=0 ? footprint : -footprint;
	a= (1<=footprint && footprint<=k) ? footprint : k+1;

	if(fabs(yObj)>=LN) yObj=yObjectMax;
	if(fabs(xObj)>=LN) xObj=xObjectMax;
	
	if(FindPupil%10 !=0){
		if( !FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,ColorStart,1) ) return 0;
		yw=ymax-ymin;
		ymax+=yw*PupilMargin/2;
		ymin-=yw*PupilMargin/2;
		xw=xmax-xmin;
		xmax+=xw*PupilMargin/2;
		xmin-=xw*PupilMargin/2;
	}
	else{
		ymax=point(0,EPD/2); ymin=point(0,-EPD/2);
		xmax=point((EPDx==0 ? EPD:EPDx)/2,0); xmin=point((EPDx==0 ? -EPD:-EPDx)/2,0);
	}

	make_coordinate(defocus);
	makecoordinate=0;            // �������̂��߁D���ۂɌ��ʂ�����D
	principal=(ymax+ymin)/2;

	if(AutoInTangent==0){
		InTangent=false;
	}
	else{
		InTangent=false;
		if(a==k+1 && is_plane(k+1)==0){
			if(RayTrace(yObj,xObj,principal.y,principal.x,defocus,ColorStart,0,0,0,0)==0){
				hxp=this->x[k+1]; hyp=this->y[k+1]; hzp=this->z[k+1];   //  ������Ƒ��ʂ̌�_
				T=Tmatrix(surface_normal(k+1,hyp,hxp,1));
				InTangent=true;
			}
		}
	}

	yp0=ymin.y; dyp=(ymax-ymin).y/this->nSpot;
	xp0=xmin.x; dxp=(xmax-xmin).x/this->nSpot;

	for(j=ColorEnd; j>=ColorStart; --j) {  // �`���ɑ�1�g�������g���ɖ�����Ȃ��悤, ColorEnd����J�n�D
		                                   // (�`���spot�̐擪�v�f����s�Ȃ���D)
		for(iy=1; iy<=this->nSpot; iy++) for(ix=1; ix<=nSpot; ix++){
			
			if(FindPupil/10 !=0){ // ���������̂�
				th=360.0/this->nSpot/this->nSpot* ((this->nSpot*(iy-1))+ix);
				if(FindPupil%10 !=0){
					FindMarginalRay(yPupil,xPupil,th,yObj,xObj,principal.y,principal.x,j,makecoordinate);
				}
				else{
					if(EPDx==0){
						yPupil=EPD/2*sin(th*PI/180);
						xPupil=EPD/2*cos(th*PI/180);
					}
					else{
						// �ȉ~�̕�����  x^2/a^2 + y^2/b^2 = (rcos��)^2 / a^2 + (rsin��)^2 / b^2 = 1
						// ��� r �����܂�̂ŁC x=rcos��, y=rsin�� �ɂ�� x,y �����܂�D
						double a,b,r, cs,sn;

						a=EPDx/2;
						b=EPD/2;
						cs=cos(th*PI/180);
						sn=sin(th*PI/180);
						r=sqrt( 1/(cs*cs/a/a+sn*sn/b/b) );
						xPupil=r*cs;
						yPupil=r*sn;			
					}
				}
			}
			else{
				yPupil=yp0+dyp*(iy-0.5+(Randomize ? Random(0.5,0):0));
				xPupil=xp0+dxp*(ix-0.5+(Randomize ? Random(0.5,0):0));
			}
			
			ra= EPDx==0 ? EPD/2 : EPDx/2;
			rb= EPD/2;

			if( FindPupil || yPupil*yPupil/rb/rb+xPupil*xPupil/ra/ra<=1 ){
				dr= IsAreaSource ? dObject() : point(0,0);
				// ���ӁF
				//  �ʌ������g���Ƃ��C���ʒut�����ۂƗ��ꂷ���Ă���(���̈ʒu�ɂ���Č����̈ʒu���قȂ�)�C
				//  ���� FindPupil=true �̂Ƃ��C�ʌ������ӂ���̒ʂ�ׂ��������ʂ�Ȃ����Ƃ�����D
				if( a==k+1 ){
					if(ImageHeight(hy,hx,hz, yObj+dr.y,xObj+dr.x,"",yPupil,xPupil,defocus,j,1,1,IsLambert,1,
						           makecoordinate)) {
						if(InTangent){
							// �ڕ��ʂɓ��e����D�Ⴆ�΁C�L�pSLO�̊����ӂł̕]���ɗL���D
							hx0=hx-hxp; hy0=hy-hyp; hz0=hz-hzp;      // ������Ƃ̌�_�𒆐S�ɉ�]
							hx=T[1][1]*hx0+T[1][2]*hy0+T[1][3]*hz0;
							hy=T[2][1]*hx0+T[2][2]*hy0+T[2][3]*hz0;
							hx+=hxp;
							hy+=hyp;
						}
						n++;
						p.p.y=hy; p.p.x=hx;
						p.wl=wl[j];
						p.weight=colorweight[j]*ApodizationIntensity();
						p.opl=0;
						spot.AddTail(p);
					}
					spot.Unit=AfocalImageUnitStr();
				}
				else{
					// �t�b�g�v�����g
					errcode=RayTrace(yObj+dr.y,xObj+dr.x,yPupil,xPupil,defocus,j,1,1,IsLambert,0,
					                 0,vector<complex>(),makecoordinate);
					if( errcode==0 || ((NotThruSurf(errcode)>=a) && (all_ray!=0)) ){
						n++;
						p.p.y=y[a]; p.p.x=x[a];
						p.wl=wl[j];
						p.weight=colorweight[j]*ApodizationIntensity();
						p.opl=optical_path[a];
						spot.AddTail(p);
					}
					spot.Unit="";
				}
			}
		}
	}

	if(this->IgnoreTC) spot.RemoveTC();

	if(OriginIsGravityCenter){
		spot.OriginToGravityCenter();
	}
	else{
		spot.OriginToNewCenter(0,0);
	}

	// �����f�ʐς̌v�Z
	if(FindPupil/10 !=0){ // ���������݂̂̏ꍇ�����v�Z����D
		int i,n, counter;
		double *x,*y;
		cPoint p;
		
		n=spot.GetSize();
		x=new double[n+1];
		y=new double[n+1];
		counter=0;
		for(i=1; i<=n; ++i){
			spot.GetData(p,i);
			if(p.wl==Wl(ColorStart)){  // ColorStart�̔g���ɂ��Čv�Z����
				counter++;
				x[counter]=p.p.x;
				y[counter]=p.p.y;
			}
		}

		spot_area=PolygonArea(counter,x,y);
		delete [] x;
		delete [] y;
	}
	else{
		spot_area=0;
	}
	
	return spot.GetSize();
}

int cLens1::MakeSpot(double yObj,double xObj,int FindPupil,double defocus,int footprint,
					int ColorStart,int ColorEnd,
					int IsAreaSource,int IsLambert,int OriginIsGravityCenter,int Add,
					int AutoInTangent/*=1*/) {
	return MakeSpot(this->spot,yObj,xObj,FindPupil,defocus,footprint,ColorStart,ColorEnd,
	                IsAreaSource,IsLambert,OriginIsGravityCenter,Add,0,AutoInTangent);
}

double cLens1::SpotTotalIntensity() {
	return spot.TotalIntensity();
}

double cLens1::SpotIntensity(double y,double x,double SensorPhi) {
	return spot.Intensity(y,x,SensorPhi);
}

double cLens1::SpotXIntensity(double x,double SensorFW) {
	return spot.XIntensity(x,SensorFW);
}

double cLens1::SpotYIntensity(double y,double SensorFW) {
	return spot.YIntensity(y,SensorFW);
}

double cLens1::SpotMTF(double nu_y,double nu_x) {
	return spot.MTF(nu_y,nu_x);
}

double cLens1::SpotRmsPhi() {
	return spot.RmsPhi();
}
	
double cLens1::SpotXMax()   { return spot.XMax(); }
double cLens1::SpotYMax()   { return spot.YMax(); }
double cLens1::SpotXMin()   { return spot.XMin(); }
double cLens1::SpotYMin()   { return spot.YMin(); }
double cLens1::SpotXWidth() { return spot.XWidth(); }
double cLens1::SpotYWidth() { return spot.YWidth(); }
double cLens1::SpotXGravityCenter() { return spot.XGravityCenter(); }
double cLens1::SpotYGravityCenter() { return spot.YGravityCenter(); }
double cLens1::SpotInertiaPrincipalAxis() { return spot.InertiaPrincipalAxis(); }
double cLens1::SpotPrincipalWx() { return spot.PrincipalWx(); }
double cLens1::SpotPrincipalWy() { return spot.PrincipalWy(); }
double cLens1::SpotArea() { return spot_area; } // ��������������`���ꍇ�C���̗֊s�ň͂܂��ʐς�Ԃ�

int cLens1::AddRayToSpot(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
                        double defocus,int footprint,int color)
{
	point p;
	if( FindRay(p,yObj,xObj,SetRay,1,color) ){
		yPupil=p.y;
		xPupil=p.x;
	}
	else{
		return 0;
	}
	int a= (1<=footprint && footprint<=k) ? footprint : k+1;
	if( RayTrace(yObj,xObj,yPupil,xPupil,defocus,color,1,1,0,0)==0 ) {
		spot.AddTail( cPoint(x[a],y[a],wl[color]) );
		return 1;
	}
	else return 0;
}

void cLens1::ClearSpot() { spot.RemoveAll(); }

double cLens1::rmsphi(double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,int FindPupil,
					  int OptimizeDefocus)
{
	cSpot spot;
	
	if(OptimizeDefocus){
		defocus=OptimizedDefocus(FindPupil); // ��1�g���C���㕨�_�ɑ΂���œK�f�t�H�[�J�X��
	}
	if( MakeSpot(spot,yObj,xObj,FindPupil,defocus,0,ColorStart,ColorEnd,0,0,1,0,0) ){
		return spot.RmsPhi();
	}
	else{
		return 0;
	}
}

double cLens1::rmsphiOn(int OptimizeDefocus){
	return rmsphi(0,0,0,1,cn,1,OptimizeDefocus);
}
double cLens1::rmsphiOff(int OptimizeDefocus){
	return rmsphi(yObjectMax,xObjectMax,0,1,cn,1,OptimizeDefocus);
}

double cLens1::EncircledEnergy(double SensorPhi,double yObj,double xObj,double defocus,
							   int ColorStart,int ColorEnd,int FindPupil,int OptimizeDefocus)
{
	cSpot spot;
	
	if(OptimizeDefocus){
		defocus=OptimizedDefocus(FindPupil); // ��1�g���C���㕨�_�ɑ΂���œK�f�t�H�[�J�X��
	}
	if( MakeSpot(spot,yObj,xObj,FindPupil,defocus,0,ColorStart,ColorEnd,0,0,1,0,0) ){
		return spot.EncircledEnergy(SensorPhi);
	}
	else{
		return 0;
	}
}

int cLens1::MakeRmsMap(double yObjMax,double xObjMax,int IsRectObj,double defocus,
                       int ColorStart,int ColorEnd,int FindPupil){
	const int N=17;    // �Б�8 + ����1
	int iy,ix;
	double yo,xo;

	if(yObjMax==0 && xObjMax==0) { yObjMax=yObjectMax; xObjMax=xObjectMax; }
	if(yObjMax==0 && xObjMax==0) return 0;
	if(yObjMax==0) { yObjMax=xObjMax; }
	if(xObjMax==0) { xObjMax=yObjMax; }

	this->rms_map.redim(N,N);

	for(iy=1; iy<=N; ++iy) for(ix=1; ix<=N; ++ix) {
		yo= yObjMax-yObjMax*2*(iy-1)/(N-1);      //  i=1��y�ő�ɑΉ�������
		xo=-xObjMax+xObjMax*2*(ix-1)/(N-1);      //  j=1��x�ŏ��ɑΉ�������
		if( IsRectObj==0 && (yo/yObjMax)*(yo/yObjMax)+(xo/xObjMax)*(xo/xObjMax)>1 ){
			rms_map[iy][ix]=0;
		}
		else{
			rms_map[iy][ix]=rmsphi(yo,xo,defocus,ColorStart,ColorEnd,FindPupil,0);
		}
	}

	return N;
}
double cLens1::RmsMap(int i,int j){
	if( 1<=i && i<=rms_map.rows() && 1<=j && j<=rms_map.columns() ){
		return this->rms_map[i][j];
	}
	else{
		return 0;
	}
}

complex cLens1::OTF(double nu_y,double nu_x,double yObj,double xObj,double defocus, 
				   int Geometrical, int AberrationFree, int ColorStart, int ColorEnd, int FindPupil){
	// Apodization�ɖ��Ή��i06.09.14�j
	if( ColorStart>cn || ColorEnd>cn || ColorStart>ColorEnd ) return 0;
	if( nu_y==0 && nu_x==0 ) return 1;
	
	if(Geometrical) {
		// �􉽌��wOTF
		cSpot spot;
		MakeSpot(spot,yObj,xObj,FindPupil,defocus,0,ColorStart,ColorEnd,0,0,0,0,1);
		return spot.OTF(nu_y,nu_x);
	}
	else {
		// �g�����wOTF
		int iy,ix, j;
		double yPupil,xPupil;
		double n,Rc,Rs, w,RWc,RWs;
		point principal,ymax,ymin,xmax,xmin;
		double epphi,epyo,epxo;
		double fno_y,fno_x;
		double dyPupil,dxPupil;
		double op,op0;
		vector<complex> dummy;

		w=RWc=RWs=0;

		if(FindPupil){
			if( !(FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,ColorStart,1)) ) return 0;
			epphi=Max( distance(ymax,ymin),distance(xmax,xmin) ) *(1+PupilMargin);
			epyo=((ymax+ymin)/2).y; epxo=((ymax+ymin)/2).x;
		}
		else{
			epphi= EPDx==0 ? EPD : Max(EPD,EPDx);
			epyo=epxo=0;
		}

		for(j=ColorStart; j<=ColorEnd; ++j) {
			if( !Afocal ){
				/*
				if(t<LN && s<LN){ 
					dyPupil=nu_y*wl[j]/1000000*(dk(defocus)-delta1())*(s-t)/(s-delta());
					dxPupil=nu_x*wl[j]/1000000*(dk(defocus)-delta1())*(s-t)/(s-delta());
				}
				else if(t<LN && s>=LN){
					dyPupil=nu_y*wl[j]/1000000*(dk(defocus)-delta1());
					dxPupil=nu_x*wl[j]/1000000*(dk(defocus)-delta1());
				}
				else if(t>=LN && s<LN){
					dyPupil=nu_y*wl[j]/1000000*(dk(defocus)-delta1())/(s-delta());
					dxPupil=nu_x*wl[j]/1000000*(dk(defocus)-delta1())/(s-delta());
				}
				else{
					return 0;
				}
				*/
				if(fno_for_unit_pupil(fno_y,fno_x,yObj,xObj,defocus,j,0)){
					dyPupil=nu_y*wl[j]/1000000*fno_y;
					dxPupil=nu_x*wl[j]/1000000*fno_x;
				}
				else{
					return 0;
				}
			}
			else{
				double A=1/AfocalImageUnit();
				/*
				if(t<LN && s<LN){ 
					dyPupil=nu_y*A*wl[j]/1000000*(s-t)/(s-delta());
					dxPupil=nu_x*A*wl[j]/1000000*(s-t)/(s-delta());
				}
				else if(t<LN && s>=LN){
					dyPupil=nu_y*A*wl[j]/1000000/gamma();
					dxPupil=nu_x*A*wl[j]/1000000/gamma();
				}
				else if(t>=LN && s<LN){
					dyPupil=nu_y*A*wl[j]/1000000/(s-delta());
					dxPupil=nu_x*A*wl[j]/1000000/(s-delta());
				}
				else{
					return 0;
				}
				*/
				if(fno_for_unit_pupil(fno_y,fno_x,yObj,xObj,defocus,j,0)){
					dyPupil=nu_y*A*wl[j]/1000000*fno_y;
					dxPupil=nu_x*A*wl[j]/1000000*fno_x;
				}
				else{
					return 0;
				}
			}			

			if( !make_ref_sphere(yObj,xObj,epyo,epxo,defocus, this->IgnoreTC ? j : 1) ) return 0;
			n=Rc=Rs=0;
			ref_sphere->make_coordinate(0);

			for(iy=1; iy<=this->nSpot; ++iy) for(ix=1; ix<=this->nSpot; ++ix){
				yPupil=-epphi/2+epphi/this->nSpot*(iy-0.5 +Random(0.5,0))+epyo;
				xPupil=-epphi/2+epphi/this->nSpot*(ix-0.5 +Random(0.5,0))+epxo;
				if( FindPupil || (yPupil-epyo)*(yPupil-epyo)+(xPupil-epxo)*(xPupil-epxo)<=epphi*epphi/4 ){
					if( ref_sphere->RayTrace(yObj,xObj,yPupil,xPupil,0,j,1,1,0,0,
						                     0,dummy,0)==0 ){
						op0=ref_sphere->optical_path[(this->k)+2]; n+=1;
						if( ref_sphere->RayTrace(yObj,xObj,yPupil+dyPupil,xPupil+dxPupil,0,j,1,1,0,0,
							                     0,dummy,0)==0 ){
							op=ref_sphere->optical_path[(this->k)+2];
							if(AberrationFree) op0=op=0;
							Rc+=cos( 2*PI/(wl[j]/1000000)*(op-op0) ); 
							Rs+=sin( 2*PI/(wl[j]/1000000)*(op-op0) );
						}
					}
				}
			}

			RWc+=Rc/n*colorweight[j];
			RWs+=Rs/n*colorweight[j];
			w+=colorweight[j];
		}
		if(n!=0) return complex(RWc/w,RWs/w); else return 0;
	}
}

complex cLens1::OTFs(double nu,double yObj,double xObj,double defocus, 
				    int Geometrical, int AberrationFree, int ColorStart, int ColorEnd, int FindPupil){
	double th, nux,nuy;
	std::string setray;

	setray= FindPupil ? "principal" : "";
	th=ImageHeightTh(yObj,xObj,setray,0,0,defocus,1);
	nux=-nu*sin(th);
	nuy= nu*cos(th);

	return OTF(nuy,nux,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil);
}

complex cLens1::OTFm(double nu,double yObj,double xObj,double defocus, 
				    int Geometrical, int AberrationFree, int ColorStart, int ColorEnd, int FindPupil){
	double th, nux,nuy;
	std::string setray;

	setray= FindPupil ? "principal" : "";
	th=ImageHeightTh(yObj,xObj,setray,0,0,defocus,1);
	nux=nu*cos(th);
	nuy=nu*sin(th);

	return OTF(nuy,nux,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil);
}

double cLens1::MTFy(double nu,double yObj,double xObj,double defocus,
				   int Geometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil){
	return abs( OTF(nu,0,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil) );
}

double cLens1::MTFy(double nu,double yObj,double xObj){
	return MTFy(nu,yObj,xObj,OptimizedDefocus(0,0,1,1),0,0,1,cn,1);
}

double cLens1::MTFx(double nu,double yObj,double xObj,double defocus,
				   int Geometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil){
	return abs( OTF(0,nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil) );
}

double cLens1::MTFx(double nu,double yObj,double xObj){
	return MTFx(nu,yObj,xObj,OptimizedDefocus(0,0,1,1),0,0,1,cn,1);
}

double cLens1::MTFxyave(double nu,double yObj,double xObj,double defocus,
				        int Geometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil){
	// ���ӁF�𑜗͂̕��ς�MTF�̕��ς��瓱�o�ł��Ȃ��D
	return (MTFy(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil)
		   +MTFx(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil)) /2;
}

double cLens1::MTFxyave(double nu,double yObj,double xObj){
	// ���ӁF�𑜗͂̕��ς�MTF�̕��ς��瓱�o�ł��Ȃ��D
	return (MTFy(nu,yObj,xObj)+MTFx(nu,yObj,xObj))/2;
}

double cLens1::MTFs(double nu,double yObj,double xObj,double defocus,
				   int Geometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil){
	return abs( OTFs(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil) );
}

double cLens1::MTFs(double nu,double yObj,double xObj){
	return MTFs(nu,yObj,xObj,OptimizedDefocus(0,0,1,1),0,0,1,cn,1);
}

double cLens1::MTFm(double nu,double yObj,double xObj,double defocus,
				   int Geometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil){
	return abs( OTFm(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil) );
}

double cLens1::MTFm(double nu,double yObj,double xObj){
	return MTFm(nu,yObj,xObj,OptimizedDefocus(0,0,1,1),0,0,1,cn,1);
}

double cLens1::MTFsmave(double nu,double yObj,double xObj,double defocus,
				        int Geometrical,int AberrationFree, int ColorStart,int ColorEnd,int FindPupil){
	// ���ӁF�𑜗͂̕��ς�MTF�̕��ς��瓱�o�ł��Ȃ��D
	return (MTFs(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil)
		   +MTFm(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil)) /2;
}

double cLens1::MTFsmave(double nu,double yObj,double xObj){
	// ���ӁF�𑜗͂̕��ς�MTF�̕��ς��瓱�o�ł��Ȃ��D
	return (MTFs(nu,yObj,xObj)+MTFm(nu,yObj,xObj))/2;
}

double cLens1::ResPower(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                    int ColorStart,int ColorEnd,int FindPupil,double NuStepY,double NuStepX){
	const double C0=0.1;         // �R���g���X�gC0���𑜂̌��E�Ƃ���
	const int MAXTIMES=200;
	int i;
	double nuy,nux, mtf,nu;
	cXYList li;

	nuy=NuStepY;
	nux=NuStepX;
	for(i=1; i<=MAXTIMES; i++){
		nu=sqrt(nuy*nuy+nux*nux);
		mtf=abs(OTF(nuy,nux,yObj,xObj,defocus,IsGeometrical,AberrationFree,ColorStart,ColorEnd,FindPupil));
		li.AddData(nu,mtf);
		nuy+=NuStepY;
		nux+=NuStepX;
		if(mtf<C0) break;
	}
	return li.TransitionPoint(C0);
}

double cLens1::ResPowerY(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                     int ColorStart,int ColorEnd,int FindPupil,double NuStep){
	return ResPower(yObj,xObj,defocus,IsGeometrical,AberrationFree,ColorStart,ColorEnd,FindPupil,NuStep,0);
}

double cLens1::ResPowerX(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                     int ColorStart,int ColorEnd,int FindPupil,double NuStep){
	return ResPower(yObj,xObj,defocus,IsGeometrical,AberrationFree,ColorStart,ColorEnd,FindPupil,0,NuStep);
}

double cLens1::ResPowerS(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                     int ColorStart,int ColorEnd,int FindPupil,double NuStep){
	double th, nux,nuy;
	std::string setray;

	setray= FindPupil ? "principal" : "";
	th=ImageHeightTh(yObj,xObj,setray,0,0,defocus,1);
	nux=-NuStep*sin(th);
	nuy= NuStep*cos(th);
	return ResPower(yObj,xObj,defocus,IsGeometrical,AberrationFree,ColorStart,ColorEnd,FindPupil,nuy,nux);
}

double cLens1::ResPowerM(double yObj,double xObj,double defocus,int IsGeometrical,int AberrationFree,
	                     int ColorStart,int ColorEnd,int FindPupil,double NuStep){
	double th, nux,nuy;
	std::string setray;

	setray= FindPupil ? "principal" : "";
	th=ImageHeightTh(yObj,xObj,setray,0,0,defocus,1);
	nux=NuStep*cos(th);
	nuy=NuStep*sin(th);
	return ResPower(yObj,xObj,defocus,IsGeometrical,AberrationFree,ColorStart,ColorEnd,FindPupil,nuy,nux);
}

double cLens1::PTFy(double nu,double yObj,double xObj,double defocus,
				   int IsGeometrical,int ColorStart,int ColorEnd,int FindPupil){
	return arg( OTF(nu,0,yObj,xObj,defocus,IsGeometrical,0,ColorStart,ColorEnd,FindPupil) );
}

double cLens1::PTFx(double nu,double yObj,double xObj,double defocus,
				   int IsGeometrical,int ColorStart,int ColorEnd,int FindPupil){
	return arg( OTF(0,nu,yObj,xObj,defocus,IsGeometrical,0,ColorStart,ColorEnd,FindPupil) );
}

void cLens1::ObjectListAdd(double y, double x) { object_list.AddTail(point(x,y)); }
void cLens1::ObjectListClear() { object_list.RemoveAll(); }
int cLens1::ObjectListSize() { return object_list.GetSize(); }
double cLens1::ObjectListY(int i) { 
	point p;
	object_list.GetData(p,i);
	return p.y;
}
double cLens1::ObjectListX(int i) { 
	point p;
	object_list.GetData(p,i);
	return p.x;
}
void cLens1::SetDefaultObjects(double h1,double h2,double h3) {
	ObjectListClear();
	double h[4]; h[1]=h1; h[2]=h2; h[3]=h3;
	int i;
	ObjectListAdd(0,0);
	if(IsRotationallySymmetric()){
		for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(yObjectMax*h[i],xObjectMax*h[i]);
	}
	else if( IsYAxisSymmetric() ){
		if(yObjectMax!=0 && xObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(-yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,xObjectMax*h[i]);
		}
		else if(yObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(-yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,yObjectMax*h[i]);
		}
		else if(xObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,xObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(xObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(-xObjectMax*h[i],0);
		}
	}
	else if(IsXAxisSymmetric()){
		if(yObjectMax!=0 && xObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,xObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,-xObjectMax*h[i]);
		}
		else if(yObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,yObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,-yObjectMax*h[i]);
		}
		else if(xObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,xObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,-xObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(xObjectMax*h[i],0);
		}
	}
	else{
		if(yObjectMax!=0 && xObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(-yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,xObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,xObjectMax*h[i]);
		}
		else if(yObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(-yObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,yObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,-yObjectMax*h[i]);
		}
		else if(xObjectMax!=0){
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,xObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(0,-xObjectMax*h[i]);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(xObjectMax*h[i],0);
			for(i=1; i<=3; ++i) if(h[i]!=0) ObjectListAdd(-xObjectMax*h[i],0);
		}
	}
}

double cLens1::MinMTFAmongObjects(double nu,double defocus,int IsGeometrical){
	list<point> object_list=this->object_list;
	object_list.EraseDuplicates();
	if( object_list.GetSize()==0 ) return 0;
	int i;
	list<double> R;
	point object;
	for(i=1; i<=object_list.GetSize(); ++i){
		object_list.GetData(object,i);
		R.AddTail( MTFy(nu,object.y,object.x,defocus,IsGeometrical,0,1,cn,0) ); // FindPupil=0
		R.AddTail( MTFx(nu,object.y,object.x,defocus,IsGeometrical,0,1,cn,0) );
	}
	return R.Min();
}
double cLens1::MinMTF(double nu)  { return MinMTFAmongObjects(nu,0,0); }

double cLens1::MTFArea(double defocus,int IsGeometrical) {
	if( nu_list.GetSize()==0 || object_list.GetSize()==0 ) return 0;
	int i;
	double nu;
	double area=0;
	for(i=1; i<=nu_list.GetSize(); ++i){
		nu_list.GetData(nu,i);
		area+=MinMTFAmongObjects(nu,defocus,IsGeometrical);
	}
	return area;
}

double cLens1::MaximizeMTFArea(double defocus_step,int IsGeometrical){
	double max=MTFArea(0,IsGeometrical);
	double area, defocus;
	if( ( area=MTFArea(defocus_step,IsGeometrical) ) > max ) {
		max=area;
		defocus=defocus_step*2;
	}
	else {
		defocus_step*=(-1);
		defocus=defocus_step;
	}
	while(true){
		area=MTFArea(defocus,IsGeometrical);
		if(area>max){
			max=area;
			defocus+=defocus_step;
		}
		else{
			defocus-=defocus_step;
			break;
		}
	}
	return defocus;
}

double cLens1::OPL(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
                   double defocus,int j,int i1/*=0*/,int i2/*=0*/,int makecoordinate/*=1*/){
	// ��j�g���P�����ɂ��đ��ʂ܂ł̌��H�����v�Z����D
	// �Ⴆ�� OCT B�X�L�����摜�̘p�Ȃ���������̂Ɏg���i15.05.15)�D
	// i1==0 && i2==0 �łȂ��Ƃ��́C������i1�ʌ�_����i2�ʌ�_�܂ł̌��H����Ԃ��D
	point p;
	p=point(xPupil,yPupil);
	if( !FindRay(p,yObj,xObj,SetRay,1,j) ) return 0;
	if(RayTrace(yObj,xObj,p.y,p.x,defocus,j,0,0,0,0,0,vector<complex>(),makecoordinate)){
		return 0;
	}
	else{
		if(i1==0 && i2==0){
			return optical_path[k+1];
		}
		else{
			if(0<=i1 && i1<i2 && i2<=k+1){
				return this->optical_path[i2]-optical_path[i1];
			}
			else{
				return 0;
			}
		}
	}
}

int cLens1::OPD_calc(double& opd, int ref_sphere_make,double yObj,double xObj,double yPupil,double xPupil,
                    double yPupil_principal,double xPupil_principal,double defocus,int j0,int j){
	// j0�g���̎��������ɍ쐬�����Q�Ƌ����g���āC
	// j�g���̔g�ʎ������v�Z����D
	// �Q�ƌ����C�Ώی����̏��ɒǐՂ��C��őΏی����̃f�[�^(������)�Ȃǂ��Q�Ƃł���D
	double op0,op;

	if(ref_sphere_make){ // ����ref_sphere���쐬�ς̂Ƃ��͂Ƃ΂��Čv�Z���Ԃ�Z�k����DMakeOPDMap(...)�ȂǂɗL���D
		                 // ref_sphere_make=true  : make_ref_spher�����s
		                 //                =false : make_ref_sphere�����s���Ȃ�
		if( !make_ref_sphere(yObj,xObj,yPupil_principal,xPupil_principal,defocus,j0) ) return 0;
	}
	
	ref_sphere->make_coordinate(0);

	if( ref_sphere->RayTrace(yObj,xObj,yPupil_principal,xPupil_principal,0,j,1,0,0,0,
                             0,vector<complex>(),0)==0 ) {
		// ������ǐՂł�mask�𖳌��ɂ��Ă����Ȃ��ƒ����ɎՕ��̂�����w�n���v�Z�ł��Ȃ��D
		op0=ref_sphere->optical_path[(this->k)+2];
	}
	else{
		return 0;
	}
	
	if( ref_sphere->RayTrace(yObj,xObj,yPupil,xPupil,0,j,1,1,0,0,
		                     0,vector<complex>(),0)==0 ){
		op =ref_sphere->optical_path[(this->k)+2];
	}
	else{
		return 0;
	}

	opd=op-op0;    // �g�ʂ��Q�Ƌ����x���Ƃ��𐳂Ƃ���(�����C�����Ɠ����D�O��Ƌt)�D
	return 1;
}

int cLens1::OPD_calc(double& opd, double yObj,double xObj,double yPupil,double xPupil,
                    double yPupil_principal,double xPupil_principal,double defocus,int j0,int j){
	return OPD_calc(opd,0,yObj,xObj,yPupil,xPupil,yPupil_principal,xPupil_principal,defocus,j0,j);
}

double cLens1::OPD(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
                   double yPupil_principal,double xPupil_principal,double defocus,int j0,int j){
	// �P������OPD���v�Z����D
	point p;
	double opd;
	p=point(xPupil,yPupil);
	FindRay(p,yObj,xObj,SetRay,1,j);
	if( OPD_calc(opd,1,yObj,xObj,p.y,p.x,yPupil_principal,xPupil_principal,defocus,j0,j) != 0 ){
		return opd;
	}
	else{
		return 0;
	}
}

double cLens1::OPD(double yObj,double xObj,
                   int FindPupil,double ypNormalized,double xpNormalized,double defocus,int j0,int j){
	// (ypNormalized,xpNormalized) : ���̔��a��1�Ƃ��������W
	point ymax,ymin,xmax,xmin, p,p0;

	make_coordinate(defocus);
	FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,j,FindPupil,0);
	p0=(ymax+ymin)/2;
	p=p0 +(ymax-ymin)/2*ypNormalized +(xmax-xmin)/2*xpNormalized;
	return OPD(yObj,xObj,"",p.y,p.x,p0.y,p0.x,defocus,j0,j);
}

std::string cLens1::OPDs(double yObj,double xObj,int findpupil,int colors,int IgnoreTC/*=0*/,int n/*=0*/,double weight/*=1*/){
	// n : �i�q�� nxn. �������C���ۂ̊i�q�_����nxn��菭�Ȃ��D
	int i,j,j0;
	std::string s;
	char buf[1000];
	double w;

	if(n==0) n=10;   // �f�t�H���g�� n=10
	if(colors>cn) colors=cn;
	MakePupilGrid(yObj,xObj,findpupil,n,1,1,0,1);

	for(j=1; j<=colors; ++j){
		j0= IgnoreTC ? j : 1;
		w=colorweight[j]*weight;

		for(i=1; i<=this->gpoints.GetSize(); ++i){
			sprintf(buf,"(OPD %g %g %d %g %g 0 %d %d) 0 %g;\n",
			        yObj,xObj,findpupil,gpoints[i].y,gpoints[i].x,j0,j,w);
			s+=buf;
		}
	}
	return s;
}

double cLens1::OPD2(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
                    double yPupil_principal,double xPupil_principal,double defocus,
					double wl_nm0,double wl_nm){
	// �P������OPD���v�Z����D�F�͔g��(nm)�Ŏw�肷��D
	cLens1 X=*this;

	X.IndexDigits=10; // ���ܗ��������_�ȉ�10���Ƃ���D
	                  // �Ⴆ�΁CSD-OCT�̕����v�̃Z���T��846����913nm��1024��������D
	                  // ���̂Ƃ������_5���ł͋��ܗ���1024�i�K�̊��炩�ȕω��ŕ\���Ȃ��D
	X.Set_cn(2);
	X.Set_color(1,str(wl_nm0));
	X.Set_color(2,str(wl_nm));
	return X.OPD(yObj,xObj,SetRay,yPupil,xPupil,yPupil_principal,xPupil_principal,defocus,1,2);
}

std::string cLens1::DispersionTable(double yObj,double xObj,std::string SetRay,double yPupil,double xPupil,
	                                double yPupil_principal,double xPupil_principal,
								    double wl_start,double wl_end,int wl_points){
	// ���H���̕��U�̕\���쐬����D
	// (yObj,xObj, yPupil_principal,xPupil_principal)�Ō��܂�����ɉ��������H���ɑ΂���C
	// (yObj,xObj, yPupil,xPupil)�Ō��܂�����ɉ��������H�����̕\���쐬����D
	// �\�͓��Ԋu�̔g���i1/wl_end�`1/wl_start)�ɑ΂��č쐬����D���H����deg�ŕ\���D
	// �܂��C�萔�����C�g���ɑ΂���1���������͏����������̂𕹂��ďo�͂���D
	// ����ɂ��COCT�ɉe�����鎋�쒆�S�ȊO�ł̕��U��]���ł���D
	cFitting X;
	std::string s;
	int i;
	double k0,dk,opd;
	double *wl,*p0,*p;
	char buf[1000];

	wl=new double [wl_points+1];
	p0=new double [wl_points+1];  // ���U
	p =new double [wl_points+1];  // �ꎟ�������������������U


	k0=1/wl_end;
	dk=((1/wl_start)-(1/wl_end))/(wl_points-1);

	for(i=1; i<=wl_points; i++){
		wl[i]=1/(k0+dk*(i-1));
		opd=OPD2(yObj,xObj,SetRay,yPupil,xPupil,yPupil_principal,xPupil_principal,0,wl_start,wl[i]);
		p0[i]=(opd/(wl[i]/1000000))*360;
		X.AddData(wl[i],0,p0[i]);
	}

	X.SetNumberOfTerms(2);
	X.dimensionSet(1,0,0);
	X.dimensionSet(2,1,0);
	X.RemoveTerms();

	for(i=1; i<=wl_points; i++){
		p[i]=X.ydataGet(i);
	}

	sprintf(buf,"  i\t   wl   \t    k   \t    p0(deg) \t p(deg) \n");
	s+=buf;
	for(i=1; i<=wl_points; i++){
		sprintf(buf,"%3d\t%8.3f\t%8.2f\t%12.5f\t%8.5f\n",
		        i,wl[i],1000000/wl[i]*2*PI,p0[i],p[i]);
		s+=buf;
	}

	delete [] wl; delete [] p0; delete [] p;
	return s;
}

int cLens1::MakeOPDMap(matrix<double>& X,matrix<double>& Y,
                      matrix<double>& W,matrix<double>& P,matrix<double>& A,matrix<int>& E,
					  cZernike& zernike,
                      double yObj,double xObj,double defocus,int j0,int j,int FindPupil,int InLambda,
					  int ZernikeOrder,int AdjustSph,int ZernikePupil,double ShearDx,double ShearDy,
					  std::string ZernikeFix)
{
	int iy,ix, n;
	point ymax,ymin,xmax,xmin, ymax1,ymin1,xmax1,xmin1;
	double epphi,epphi_z;
	point epo,epo_z;
	double yPupil,xPupil, a,b;
	double opd;
	double unit,wl;
	matrix<double> X1,Y1;
	double zo;
	vector<double> eo,qo;

	n=this->nSpot;
	X.redim(n,n);  Y.redim(n,n);   // ���˓����W(���˓��ƌ����̌�_)
	X1.redim(n,n); Y1.redim(n,n);  // Zernike�W�J�p�������W
	W.redim(n,n);   // �g�ʎ���
	P.redim(n,n);   // ���֐��ʑ�
	A.redim(n,n);   // ���֐��U��
	E.redim(n,n);   // �������ʉ߂����ꍇ1, ����ȊO��0

	if(FindPupil || EPD==0){
		if( !(FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,j,1)) ) return 0;
		epphi=Max( distance(ymax,ymin),distance(xmax,xmin) )*(1+PupilMargin);
		epo=((ymax+ymin)/2);
	}
	else{
		ymax=point(0, EPD/2); ymin=point(0,-EPD/2);
		xmax=point((EPDx==0 ? EPD:EPDx)/2,0); xmin=point((EPDx==0 ? -EPD:-EPDx)/2,0);
		epphi= EPDx==0 ? EPD:Max(EPD,EPDx);
		epo=point(0,0);
	}
	if( !make_ref_sphere(yObj,xObj,epo.y,epo.x,defocus,j0) ) return 0;

	// ����ZernikePupil : ��ZernikePupil�ʂƌ����̌�_(x,y)��Zernike�W�J����D
	//   �t�b�g�v�����g�ɑ傫�Șc�݂��Ȃ����r0���ς邾���ŌW���͕ς�Ȃ��D
	//   ���������ĒP�Ɍ����S�̂̎��������邾���Ȃ��{�I�ɂ͕s�v�ł���D
	//   �������C�Ⴆ��SetOPDMapZernikeR0Fix()��R0���w�肷��Ƃ�
	//   �w�肵�₷����(�����a���������Ă���ʁj��I�Ԃ��Ƃ��ł��ĕ֗��ł���D

	if(ZernikePupil==-1){
		// ���˓��ƌ����̌�_(x,y)�œW�J����
		ymax1=ymax;
		ymin1=ymin;
		xmax1=xmax;
		xmin1=xmin;
	}
	else if(ZernikePupil==0){
		// �Q�Ƌ��ƌ����̌�_(x,y)�œW�J����i���̒ǋL���Q�Ɓj
		// <note>
		//   �Ⴆ�Γ��˓����(x,y)��Zernike�W�J������ƁC
		//   �ˏo�������˓��ɑ΂��Ęc�ނ悤�Ȍ��w�n�ł�
		//   �W�J���ʂɑ��݂��Ȃ�cyl�������t������邱�ƂɂȂ�C
		//   AdjustSph��True�ŋ��ʐ����ȉ�����������Ƃ���cyl�������c��s����������邱�Ƃ�����D
		//   ���̏ꍇ�͎ˏo��(�Q�Ƌ�)���(x,y)�œW�J�Ƃ悢�D
		//
		// �ǋL : �Q�Ƌ��ƌ����̌�_(x,y)�́gz����������ƈ�v������W�ւ̓��e�h�ł̓W�J�ɕύX�D
		//        �������Ȃ��ƁC���O�i�����ȉ~�ɂȂ�j�Ŕ�_����������������Ă��܂��D
		//        (ref_sphere->k+1 �ʂ͌X���ĂȂ����C�ȉ���X1�CY1�ւ̑�����ɍ��W�ϊ����Ă���j
		ymax1=ref_sphere->RayHeight(ref_sphere->k+1,yObj,xObj,"",ymax.y,ymax.x,j);
		ymin1=ref_sphere->RayHeight(ref_sphere->k+1,yObj,xObj,"",ymin.y,ymin.x,j);
		xmax1=ref_sphere->RayHeight(ref_sphere->k+1,yObj,xObj,"",xmax.y,xmax.x,j);
		xmin1=ref_sphere->RayHeight(ref_sphere->k+1,yObj,xObj,"",xmin.y,xmin.x,j);
	}
	else if(1<=ZernikePupil && ZernikePupil<=this->k){
		// ��ZernikePupil��(�ڕ��ʂł͂Ȃ�)�ƌ����̌�_(x,y)�œW�J����
		// �y���Ӂz
		//   ������Ƒ�ZernikePupil�ʂ̖@�����傫���p�x���Ȃ��ꍇ�C
		//   ����AdjustSph�ɂ��f�t�H�[�J�X�������������悤�Ƃ���ƁC
		//   ��_�����������������Ă��܂����ƂɂȂ�D
		//   �Ⴆ�Δ�_�������x�z�I�ȏꍇ�C���ۂ����������ǂ������Ă��܂����Ƃ�����D
		ymax1=ref_sphere->RayHeight(ZernikePupil,yObj,xObj,"",ymax.y,ymax.x,j);
		ymin1=ref_sphere->RayHeight(ZernikePupil,yObj,xObj,"",ymin.y,ymin.x,j);
		xmax1=ref_sphere->RayHeight(ZernikePupil,yObj,xObj,"",xmax.y,xmax.x,j);
		xmin1=ref_sphere->RayHeight(ZernikePupil,yObj,xObj,"",xmin.y,xmin.x,j);
	}
	else{
		return 0;
	}

	// epo_z=Zernike�W�J�̌��_
	// epphi_z=Zernike�W�J����a�DcZernike::CalcR0()���g���ƌ������ɂ���ĕς���Ă��܂��̂ł������������悢�D
	epo_z=(ymax1+ymin1+xmax1+xmin1)/4;
	epphi_z=Max( distance(ymax1,epo_z),distance(ymin1,epo_z),distance(xmax1,epo_z),distance(xmin1,epo_z) )*2;

	if(ZernikePupil==0){
		// iy,ix���[�v�ɓ���O�ɏ�������D
		zo=ref_sphere->surface_sag(ref_sphere->k+1,epo_z.y,epo_z.x,0); // �Q�Ƌ��ʋǏ����W�ŁCepo_z�ɂ�����T�O
		eo=ref_sphere->surface_normal(ref_sphere->k+1,epo_z.y,epo_z.x,0); 
		eo=eo/abs(eo);  // �Q�Ƌ��ʍ��W�ŁCepo_z�ɂ�����@�������P�ʃx�N�g�������������
	}

	zernike.DataClear();
	zernike.SetMaxOrder( ZernikeOrder>4 ? ZernikeOrder : 4 ); // ��Ńs�X�g�������̂���Zernike�W�J�K�{
	// �ȑO��ZernikeOrder�̍Œ�l��4�ɂ��Ă������C������w�n��3�ɂ��Ȃ��Ɣg�ʎ�����0�ɂȂ��Ă��܂�
	// (�ŏ����@�ŋt�s�񂪔񐳑��Ōv�Z�ł��Ȃ�?)���߁C�Œ�l������2�ɐݒ肵���D(080910)
	// -> �������CAdjustSph=1 �̂Ƃ��CZernikeOrder=2�ł͍��������������Ă��܂�����,
	//    �g�ʎ�����0�ɂȂ��Ă��܂��D�ĂэŒ�l4�ɐݒ�D�񐳑��̌��ɂ��Ă� matrix.h ��
	//    �t�s��v�Z�ɂ�����Â��������ύX�����̂ŁC�܂��l�q������D (090819)
	zernike.SetR0(epphi_z/2);
	wl  = Wl(j)/1000000;
	// InLambda==0 �̂Ƃ��͔g�ʎ�������m�P�ʂŕ\��. ����ƁC
	// InLambda!=0 �̂Ƃ��Ɠ����ʂ̃I�[�_�[�ƂȂ�C���^�[�t�F�[�X�쐬��s�����悢�D
	unit= InLambda ? wl : 0.001;
	a=EPDx==0 ? EPD/2 : EPDx/2;
	b=EPD/2;
	for(iy=1; iy<=n; iy++) for(ix=1; ix<=n; ix++){
		yPupil=-epphi/2+epphi/n*(iy-0.5)+epo.y;
		xPupil=-epphi/2+epphi/n*(ix-0.5)+epo.x;
		X[n-iy+1][ix]=xPupil-epo.x;
		Y[n-iy+1][ix]=yPupil-epo.y;
		X1[n-iy+1][ix]=0;
		Y1[n-iy+1][ix]=0;
		W[n-iy+1][ix]=0;
		P[n-iy+1][ix]=0;
		A[n-iy+1][ix]=0;
		E[n-iy+1][ix]=0;
		if( FindPupil || EPD==0 || yPupil*yPupil/b/b+xPupil*xPupil/a/a<=1 ){
			if(OPD_calc(opd,0,yObj,xObj,yPupil,xPupil,epo.y,epo.x,defocus,j0,j)){
				// ���ӁFOPD_calc()��ref_sphere��RayTrace()�����s����D
				if(ZernikePupil==-1){
					X1[n-iy+1][ix]=X[n-iy+1][ix];
					Y1[n-iy+1][ix]=Y[n-iy+1][ix];
				}
				else if(ZernikePupil==0){
					qo.x=ref_sphere->x[ref_sphere->k+1]-epo_z.x; // Zernike�W�J���_�ɑ΂�������ʉߓ_�ʒu
					qo.y=ref_sphere->y[ref_sphere->k+1]-epo_z.y;
					qo.z=ref_sphere->z[ref_sphere->k+1]-zo;
					qo=Tmatrix(eo)*qo;    // �@��������z���Ƃ�����W�ɕϊ�
					X1[n-iy+1][ix]=qo.x;
					Y1[n-iy+1][ix]=qo.y;
				}
				else if(1<=ZernikePupil && ZernikePupil<=this->k){
					X1[n-iy+1][ix]=ref_sphere->x[ZernikePupil]-epo_z.x;
					Y1[n-iy+1][ix]=ref_sphere->y[ZernikePupil]-epo_z.y;
				}
				W[n-iy+1][ix]=opd/unit;
				A[n-iy+1][ix]=ref_sphere->ApodizationAmplitude();
				E[n-iy+1][ix]=1;
				zernike.SetData(X1[n-iy+1][ix],Y1[n-iy+1][ix],W[n-iy+1][ix],n-iy+1,ix);
			}
		}
	}
	
	{
		int i;
		
		switch(AdjustSph){
		case -1:
			break;                       // �������Ȃ��Dpiston,tilt,defocus�W�������̂܂ܓ�����
		case 0:
			zernike.AdjustTerms(1,1,0);  // piston,tilt�𒲐����Ĕg�ʎ���RMS���ŏ��ɂ���
			break;
		default:
			zernike.AdjustTerms(1,1,1);  // piston,tilt,defocus�𒲐����Ĕg�ʎ���RMS���ŏ��ɂ���
			break;
		}

		for(i=1; i<=zernike.NumberOfData(); i++){
			W[zernike.GetIData(i)][zernike.GetJData(i)]=zernike.GetZData(i);
		}
	}

	if( words(ZernikeFix) >= 7 ){
		// Zernike�W����S��0�Ƃ�����C
		// �w�肷��l�ɌŒ肷��D�g�ʎ���������ɉ����ČŒ肳���D
		// �����F Inlambda r0 Normalize IsFringeOrder jBase j1 val1 j2 val2 ....
		// ��  �F opd_map_zernike_fix="0 3 1 0 0 4 1 12 1" �̂Ƃ��C
		//           InLambda=0
		//           zernike.r0=3
		//           zernike.Normalize=1
		//           zernike.IsFringeOrder=0
		//           zernike.jBase=0
		//           �f�t�H�[�J�X��(��4��)�Ƌ��ʎ�����(��12��)��1
		//        �Ɏw�肷��D
		//
		// MakePSFMap()�Ȃǂ�����Ă΂�邽�߁C�����������ʓ|�����w�肵�Ȃ���΂Ȃ�Ȃ��D

		int i, n;
		std::string *s;
		int in_lambda;
		double unit1;
		int order,order_max;
		cZernike a;

		n=words(ZernikeFix);
		s=new std::string [n+1];
		for(i=1; i<=n; i++){
			s[i]=word(ZernikeFix,i);
		}

		in_lambda=atoi(s[1].c_str());
		zernike.SetR0(atof(s[2].c_str()));
		zernike.SetNormalize(atoi(s[3].c_str()));
		zernike.SetIsFringeOrder(atoi(s[4].c_str()));
		zernike.SetJBase(atoi(s[5].c_str()));

		order_max=0;
		for(i=6; i<=n-1; i+=2){
			order=zernike.Order(atoi(s[i].c_str()));
			if(order>order_max) order_max=order;
		}
		zernike.SetMaxOrder(order_max);
		
		for(i=6; i<=n-1; i+=2){
			zernike.SetC( atoi(s[i].c_str()),atof(s[i+1].c_str()));
		}
		zernike.FitZToCoefficient();
		unit1= in_lambda ? wl : 0.001;
		
		a=zernike;
		a.DataClear();
		for(i=1; i<=zernike.NumberOfData(); i++){	
			a.SetData(zernike.GetXData(i),zernike.GetYData(i),zernike.GetZData(i)*unit1/unit,
			          zernike.GetIData(i),zernike.GetJData(i));
		}
		zernike=a;
		for(i=1; i<=zernike.NumberOfData(); i++){	
			W[zernike.GetIData(i)][zernike.GetJData(i)]=zernike.GetZData(i);
		}

		delete [] s;
	}

	if(ShearDx!=0 || ShearDy!=0){
		// zernike�W�J��x,y���W��(ShearDx,SherDy)�������炵���g�ʂƂ̔g�ʎ����̍����Ƃ�
		int i;

		zernike.Shear(ShearDx,ShearDy,ZernikeOrder);	
		for(i=1; i<=zernike.NumberOfData(); i++){
			W[zernike.GetIData(i)][zernike.GetJData(i)]=zernike.GetZData(i);
		}
	}	

	for(iy=1; iy<=n; iy++) for(ix=1; ix<=n; ix++){
		if(E[iy][ix]){
			opd=W[iy][ix]*unit;
			P[iy][ix]=2*PI*opd/wl;
		}
	}
	
	return n; // �摜�̑傫����Ԃ�
}

int cLens1::MakeOPDMap(std::string filename,double yObj,double xObj,double defocus,int j0,int j,
					   int FindPupil,int InLambda,int AdjustSph,int ZernikePupil){
	matrix<double> X,Y,W,P,A;
	matrix<int> E;
	cZernike zernike;
	list<double> li;
	double rms,pv;
	std::string s;
	if( MakeOPDMap(X,Y,W,P,A,E,zernike,yObj,xObj,defocus,j0,j,FindPupil,InLambda,4,
		           AdjustSph,ZernikePupil,0,0,"") )
	{
		rms=OPDMapRMS();
		pv=OPDMapPV();
	}
	else{
		return 0;
	}
	s=remove_path(this->filename)+" y="+str(yObj)+" x="+str(xObj)+" defocus="+str(defocus)+" j="+str(j);
	if(filename==""){
		filename=s;
	}
	std::ofstream to(filename.c_str());
	if(to){
		to << s << '\n';
		to << "rms=" << Round(rms,-3) << " wave" << '\n';
		to << "p-v=" << Round(pv,-3) << " wave" << '\n';
		to << W;
		return 1;
	}
	else{
		return 0;
	}
}

int cLens1::MakeOPDMap(double yObj,double xObj,double defocus,int j0,int j,int FindPupil,int InLambda,
					  int ZernikeOrder,int AdjustSph,int ZernikePupil,double ShearDx,double ShearDy,
					  std::string ZernikeFix)
{
	return MakeOPDMap(this->opd_map_x,this->opd_map_y,
		              this->opd_map_w,this->opd_map_p,this->opd_map_a,this->opd_map_e,
					  this->opd_map_zernike,
	                  yObj,xObj,defocus,j0,j,FindPupil,InLambda,ZernikeOrder,AdjustSph,ZernikePupil,
					  ShearDx,ShearDy,ZernikeFix);
	// �y���Ӂzthis->opd_map_zernike �ւ̃A�N�Z�X�͂��������Ƃ��邱��
	// 
	//      MakePSFMap(double yObj,double xObj, ...),MakeOTFMap(double yObj,double xObj, ...)�� 
	//      this->opd_map_zernike �ɂ̓A�N�Z�X���Ȃ��D
	//      ���������āC�֐�����MakeOPDMap(matrix<double>& X,matrix<double>& Y, ...)
	//      ���قȂ�����Ŏ��s���Ă��邪 this->opd_map_zernike ��Zernike�W���̏㏑���͂Ȃ��D
}

double cLens1::OPDMap(int i,int j){
	if( 1<=i && i<=opd_map_w.rows() && 1<=j && j<=opd_map_w.columns() ){
		return opd_map_w[i][j];
	}
	else{
		return 0;
	}
}

double cLens1::OPDMapRMS(){
	// �grms�g�ʎ����h�ƌĂ΂�邪�C�v����ɕW���΍��ł���D
	// �W���΍�^2 = < (x-<x>)^2 >
	//            = < x^2 -2x<x> +<x>^2 >
	//            = <x^2> -<2x<x>> +<<x>^2>
	//            = <x^2> -2<x><x> +<x>^2
	//            = <x^2> -<x>^2
	int i,j;
	list<double> li;
	for(i=1; i<=opd_map_w.rows(); i++) for(j=1; j<=opd_map_w.columns(); j++){
		if( opd_map_e[i][j] ) li.AddTail(opd_map_w[i][j]);
	}
	return li.Stdev();
}

double cLens1::OPDMapPV(){
	int i,j;
	list<double> li;
	for(i=1; i<=opd_map_w.rows(); i++) for(j=1; j<=opd_map_w.columns(); j++){
		if( opd_map_e[i][j] ) li.AddTail(opd_map_w[i][j]);
	}
	return li.PV();
}

double cLens1::OPDMapStrehlDef(){
	const complex I(0,1);
	cZernike z;
	int i,j;
	list<complex> li;
	complex ave;
	// �g�ʎ����̃v���Y����␳����D
	for(i=1; i<=opd_map_w.rows(); i++) for(j=1; j<=opd_map_w.columns(); j++){
		if( opd_map_e[i][j] ) z.SetData(opd_map_x[i][j],opd_map_y[i][j],opd_map_p[i][j]);
	}
	z.SetMaxOrder(4);
	z.CalcR0();
	z.RemoveTerms(1,1,0);
	// Strehl's definition = |<f>|^2 ; f=���֐��C<>=����
	// �Q�l���g���@��̌��w�U(����)�hp.509
	for(i=1; i<=z.NumberOfData(); i++) {
		li.AddTail(exp(I*z.GetZData(i)));
	}
	ave=li.Ave();
	return Re(ave*conj(ave));
}

double cLens1::OPDMapZernikeR0(){
	return opd_map_zernike.GetR0();
}

double cLens1::GetOPDMapZernikeR0Fix(){
	return opd_map_zernike.GetR0Fix();
}
void cLens1::SetOPDMapZernikeR0Fix(double value){
	opd_map_zernike.SetR0Fix(value);
}

int  cLens1::GetOPDMapZernikeNormalize(){
	return opd_map_zernike.GetNormalize();
}
void cLens1::SetOPDMapZernikeNormalize(int value){
	opd_map_zernike.SetNormalize(value);
}

int  cLens1::GetOPDMapZernikeIsFringeOrder(){
	return opd_map_zernike.GetIsFringeOrder();
}
void cLens1::SetOPDMapZernikeIsFringeOrder(int value){
	opd_map_zernike.SetIsFringeOrder(value);
}

int  cLens1::GetOPDMapZernikeJBase(){
	return opd_map_zernike.GetJBase();
}
void cLens1::SetOPDMapZernikeJBase(int value){
	opd_map_zernike.SetJBase(value);
}

int cLens1::OPDMapZernikeTotalTerms(){
	return opd_map_zernike.GetNumberOfTerms();
}

double cLens1::OPDMapZernikeCoefficient(int term_no){
	opd_map_zernike.Fit();
	return opd_map_zernike.GetC(term_no);
}

double cLens1::OPDMapZernikeRMSError(){
	opd_map_zernike.Fit();
	return opd_map_zernike.RMSError();
}

double cLens1::OPDMapZernikePVError(){
	opd_map_zernike.Fit();
	return opd_map_zernike.PVError();
}

double cLens1::OPDMapZernikeRMS(int tilt,int sph,int cyl,int high_order){
	opd_map_zernike.Fit();
	return opd_map_zernike.RMS(tilt,sph,cyl,high_order);
}

int cLens1::CopyOPDMapToImage(cImage& image,int PixelWidth){
	// OPDMap�����ȂƂ��ĉ摜����image�ɕۑ�����
	// PixelWidth>0  �̂Ƃ��͉摜�̃s�N�Z������PixelWidth�ɕς���D
	//           <=0 �̂Ƃ���OPDMap�Ɠ����ɂȂ�D
	int i,j,n, i0,j0,n0;

	image.Clear();
	n0=opd_map_p.rows();
	n= PixelWidth<=0 ? n0 : PixelWidth;
	image.SetXpixels(n);
	image.SetYpixels(n);
	for(i=1; i<=n; i++) for(j=1; j<=n; j++) {
		// OPDMap�̎ȋ��x��image�̒����ɃR�s�[����
		i0=i-(n-n0)/2;
		j0=j-(n-n0)/2;
		if( 1<=i0 && i0<=n0 && 1<=j0 && j0<=n0 ){
			if(opd_map_e[i0][j0]!=0){
				image.SetIntensity(i,j,1+cos(opd_map_p[i0][j0]));
			}
			else{
				image.SetIntensity(i,j,1);
			}
		}
		else{
			image.SetIntensity(i,j,1);
		}
	}
	// �摜�̎��������킹��
	image.Zoom((double)n/(double)n0);
	return n;
}

int cLens1::SaveOPDMapAsBmp(std::string filename,int PixelWidth){
	cImage image;
	CopyOPDMapToImage(image,PixelWidth);
	image.Normalize(255);
	return image.SaveAsBmp(filename,0);
}

double cLens1::OPDRMS(double yObj,double xObj,double defocus,int j,int FindPupil,int InLambda,int AdjustSph,
                      int OptimizeDefocusOnAxis,std::string ZernikeFix){
	if(OptimizeDefocusOnAxis){
		defocus=OptimizedDefocus(FindPupil); // ���㕨�_�ɑ΂���œK�f�t�H�[�J�X��
	}
	MakeOPDMap(yObj,xObj,defocus,j,j,FindPupil,InLambda,AdjustSph!=0 ? 4:-1,AdjustSph,0,0,0,ZernikeFix);
	return OPDMapRMS();
}

double cLens1::OPDRMS0(double yObj,double xObj,int FindPupil){
	// �f�t�H�[�J�X�����̏����C�������̒������s��Ȃ��D
	// �œK���ŏœ_���킹������Ƃ��̖ڕW�l�ȂǂɎg����D
	return OPDRMS(yObj,xObj,0,1,FindPupil,1,0,0,"");
}

double cLens1::OPDPV(double yObj,double xObj,double defocus,int j,int FindPupil,int InLambda,int AdjustSph,
					 int OptimizeDefocusOnAxis,std::string ZernikeFix){
	if(OptimizeDefocusOnAxis){
		defocus=OptimizedDefocus(FindPupil); // ���㕨�_�ɑ΂���œK�f�t�H�[�J�X��
	}
	MakeOPDMap(yObj,xObj,defocus,j,j,FindPupil,InLambda,AdjustSph!=0 ? 4:-1,AdjustSph,0,0,0,ZernikeFix);
	return OPDMapPV();
}

double cLens1::StrehlDef(double yObj,double xObj,double defocus,int j,int FindPupil,int AdjustSph,
						 int OptimizeDefocusOnAxis,std::string ZernikeFix){
	if(OptimizeDefocusOnAxis){
		defocus=OptimizedDefocus(FindPupil); // ���㕨�_�ɑ΂���œK�f�t�H�[�J�X��
	}
	MakeOPDMap(yObj,xObj,defocus,j,j,FindPupil,0,AdjustSph!=0 ? 4:-1,AdjustSph,0,0,0,ZernikeFix);
	return OPDMapStrehlDef();
}
double cLens1::StrehlDef(double yObj,double xObj){
	// ����ő��������œK���������ʂŃX�g���[�����]������D
	// (�f�t�H�[�J�X���̏����͍s��Ȃ��D�������ɑ������𒲐����邱�Ƃ͂��Ȃ��D�j
	return StrehlDef(yObj,xObj,0,1,1,0,1,"");
}


double cLens1::ZernikeC(int term_no,double yObj,double xObj,double defocus,int j,int FindPupil,
                        int InLambda,int AdjustSph,int OptimizeDefocus){
	int order;
	matrix<double> X,Y,W,P,A;
	matrix<int> E;
	cZernike z;

	z.SetR0Fix(0);
	z.SetNormalize(1);
	z.SetIsFringeOrder(0);
	z.SetJBase(0);
	if(OptimizeDefocus){
		defocus=OptimizedDefocus(FindPupil); // ���㕨�_�ɑ΂���œK�f�t�H�[�J�X��
	}
	order=cZernike::nNumber(term_no,0,0);
	MakeOPDMap(X,Y,W,P,A,E,z,yObj,xObj,defocus,j,j,FindPupil,InLambda,order,AdjustSph,0,0,0,"");
	return z.GetC(term_no);
}

double cLens1::ZernikeCHigh(double yObj,double xObj,double defocus,int j,int FindPupil,
                            int InLambda,int AdjustSph,int OptimizeDefocus){
	// ���������W���i3���ȏ�̌W���j�̓��a��������Ԃ��D
	const int MAX_ORDER=8;  // 8���܂łƂ���
	int i,n;
	double c,sum;
	matrix<double> X,Y,W,P,A;
	matrix<int> E;
	cZernike z;

	z.SetR0Fix(0);
	z.SetNormalize(1);
	z.SetIsFringeOrder(0);
	z.SetJBase(0);
	if(OptimizeDefocus){
		defocus=OptimizedDefocus(FindPupil); // ���㕨�_�ɑ΂���œK�f�t�H�[�J�X��
	}
	MakeOPDMap(X,Y,W,P,A,E,z,yObj,xObj,defocus,j,j,FindPupil,InLambda,MAX_ORDER,AdjustSph,0,0,0,"");

	n=cZernike::TotalTerms(MAX_ORDER,0);
	
	sum=0;
	for(i=0; i<=n-1; ++i){
		c=z.GetC(i);
		sum+=c*c;
	}
	return sqrt(sum);
}

int cLens1::MakePSFMap(cImage& B,double& wx,double& wy,double& phi_x,double& phi_y,double phi_threshold,
                      double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,
                      int FindPupil,int StrehlNormalize,double Zoom,int AberrationFree,int AdjustSph,
					  std::string ZernikeFix)
{
	// �y���Ӂz ���F�̏ꍇ�� AdjustSph �͖����i���L�Q�Ɓj

	// StrehlNormalize = true  : �ő�l��Strehl��ƂȂ�悤�ɋK�i������
	//                 = false : �S���ʂ�1�ƂȂ�悤�ȏƓx�Ƃ���D
	//                           (�Ⴆ�΁C��Ŏ��Ƃɂē��˓��ʐς��|����Γ��˓��P�x��1�̂Ƃ��̏Ɠx�ɕϊ��ł���)

	matrix<double> X,Y,W,P,A;
	matrix<int> E;
	cZernike zernike;
	int ii,jj,n, j,j0;
	double max0,max, w;
	cImage *pX;
	double *pwx,*pwy;

	if(ColorStart>ColorEnd) ColorEnd=ColorStart;
	pX=new cImage[ColorEnd+1];
	pwx=new double[ColorEnd+1];
	pwy=new double[ColorEnd+1];

	if(ColorStart<ColorEnd){
		// ���F�̏ꍇ�͉��Əc�̐F�����𖳎��ł��Ȃ��̂� AdjustSph=true �͖���.
		// ���Ȃ킿 tilt,defocus ���l������D
		// �K�v�ɉ����āC����defocus�Œ������邱�Ƃ͂ł���D
		AdjustSph=-1;
	}

	// n���v�Z����D�����Strehl��v�Z�p�ɖ������C��ColorStart�g���P�F�ł̉�ܑ��̍ő勭�xmax0���v�Z����D
	n=MakeOPDMap(X,Y,W,P,A,E,zernike,yObj,xObj,defocus,ColorStart,ColorStart,FindPupil,0,0,0,0,0,0,ZernikeFix);
	if(StrehlNormalize){
		pX[0].SetPixels(n);
		for(ii=1; ii<=n; ii++) for(jj=1; jj<=n; jj++){
			pX[0].SetAmplitude(ii,jj,A[ii][jj]);
			pX[0].SetPhase(ii,jj,0);
		}
		pX[0].DFTAmplitude(1,1);
		pX[0]=(1/wl[ColorStart])*pX[0];
		max0=0;
		for(ii=1; ii<=n; ii++) for(jj=1; jj<=n; jj++){
			if( pX[0].GetIntensity(ii,jj)>max0 ) max0=pX[0].GetIntensity(ii,jj);
		}

		{
			// ���F�̂Ƃ��́C�ő勭�x���C
			//     ColorStart�g���P�F�ł̍ő勭�x * colorweight[j] / colorweight[ColorStart] 
			// �̘a�ɂȂ���̂ƍl����D
			// �e�F���ɍő勭�x����܌v�Z���a���Ƃ�̂��{�������ʓ|�Ȃ̂ł�������D
			double max00;
			max00=max0;
			for(j=ColorStart+1; j<=ColorEnd; j++){
				max0+=max00*colorweight[j]/colorweight[ColorStart];
			}
		}
	}

	Zoom= Zoom==0 ? 1 : Zoom;

	w=0;
	for(j=ColorStart; j<=ColorEnd; j++){
		w+=colorweight[j];
	}

	for(j=ColorStart; j<=ColorEnd; j++){
		j0= this->IgnoreTC ? j : ColorStart;   // IgnoreTC���^�̂Ƃ��͊e�F���ɎQ�Ƌ����쐬
		MakeOPDMap(X,Y,W,P,A,E,zernike,yObj,xObj,defocus,j0,j,FindPupil,0,4,AdjustSph,0,0,0,ZernikeFix);
		pX[j].SetPixels(n);

		for(ii=1; ii<=n; ii++) for(jj=1; jj<=n; jj++){
			pX[j].SetAmplitude(ii,jj,A[ii][jj]);
			pX[j].SetPhase(ii,jj,AberrationFree ? 0 : P[ii][jj]);
		}
		pX[j].DFTAmplitude(1,Zoom);
		pX[j]=(sqrt(colorweight[j]/w)/wl[j])*pX[j];
		{
			double wl0;
			double fno_x,fno_y;
			if(fno_for_unit_pupil(fno_y,fno_x,yObj,xObj,defocus,j,0)){
				fno_x/=fabs(X[n/2][1]-X[n/2][n]);
				fno_y/=fabs(Y[1][n/2]-Y[n][n/2]);
				wl0=wl[j]/1000000;
				// wx,wy=A�̈�ӂɑΉ����钷��
				//   �g�����Y�݌v�̂��߂̔g�ʌ��w(����)�h(6-43)�����C
				//    ��x = ��R/(n����)
				//        = ��R(n-1)/{n����(n-1)}     ... n�����ł͂Ȃ�n-1�����ł��邱�Ƃɒ���
				//        = ��0fno(n-1)/n, ��0=�^�󒆔g��
				//    �� w = (n-1)��x = ��0fno(n-1)^2/n
				pwx[j]= !Afocal ? wl0*fno_x*(n-1)*(n-1)/n/Zoom : wl0*fno_x*(n-1)*(n-1)/n/Zoom/AfocalImageUnit();
				pwy[j]= !Afocal ? wl0*fno_y*(n-1)*(n-1)/n/Zoom : wl0*fno_y*(n-1)*(n-1)/n/Zoom/AfocalImageUnit();
			}
			else{
				pwx[j]=pwy[j]=0;
			}
		}
	}

	// �摜�̕������킹��
	wx=wy=1E30;
	for(j=ColorStart; j<=ColorEnd; j++){
		wx= pwx[j]<wx ? pwx[j] : wx;
		wy= pwy[j]<wy ? pwy[j] : wy;
	}
	for(j=ColorStart; j<=ColorEnd; j++){
		pX[j].Zoom(pwy[j]/wy,pwx[j]/wx,n/2+1,n/2+1);
		// �����ɂ͊m�F���ĂȂ����C�摜�̒�����(n/2+1,n/2+1)�̗l�q(07.11.05)
	}

	// ���x�a��葽�FPSF���v�Z����
	B.SetPixels(n);
	for(j=ColorStart; j<=ColorEnd; j++){
		if(j==ColorStart){
			B=pX[j]; // �P�F�̏ꍇ�͈ʑ�����B�ɕۑ������
		}
		else{
			B=SumOnIntensity(B,pX[j]);  // 2�F�ڈȍ~�͋��x�a
		}
	}

	// �ő勭�xmax�����߂�
	max=0;
	for(ii=1; ii<=n; ii++) for(jj=1; jj<=n; jj++){
		if( B.GetIntensity(ii,jj)>max ) max=B.GetIntensity(ii,jj);
	}
	// ���K��
	if(StrehlNormalize){
		B.Normalize(max/max0);
	}
	else{
		// ���ˑ�(�S����)��1�ƂȂ�悤�ɋK�i������D
		// ���ˑ� = 1 = �� I��x��y,    I=�e��f�̋��x�C ��x,��y = ��f�̕�
		// �� ��I = 1 / (��x��y)
		B.NormalizeTotalIntensity(1/(wx/n)/(wy/n));
	}
	// PSF�̌aphi_x,phi_y���v�Z����
	phi_x= phi_threshold>0 ? B.FWxIntensity(phi_threshold,wx) : 0;
	phi_y= phi_threshold>0 ? B.FWyIntensity(phi_threshold,wy) : 0;

	delete [] pX; delete [] pwx; delete [] pwy;

	// �摜�̑傫����Ԃ�
	return n;
}

int cLens1::MakePSFMap(double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,
                      int FindPupil,int StrehlNormalize,double Zoom,int AberrationFree,int AdjustSph,
					  std::string ZernikeFix)
{
	int n;

	n=MakePSFMap(psf_map,psf_map_w_x,psf_map_w_y,psf_map_phi_x,psf_map_phi_y,psf_map_phi_threshold,
	             yObj,xObj,defocus,ColorStart,ColorEnd,
                 FindPupil,StrehlNormalize,Zoom,AberrationFree,AdjustSph,ZernikeFix);
	psf_map.SetHeight(psf_map_w_y);
	psf_map.SetWidth(psf_map_w_x);
	return n;
}

double cLens1::PSFMap(int i,int j){
	return psf_map.GetIntensity(i,j);
}

double cLens1::PSFMapWX(){
	return psf_map_w_x;
}
double cLens1::PSFMapWY(){
	return psf_map_w_y;
}
double cLens1::PSFMapPhiX(){
	return psf_map_phi_x;
}
double cLens1::PSFMapPhiY(){
	return psf_map_phi_y;
}
double& cLens1::PSFMapPhiThreshold(){
	return psf_map_phi_threshold;
}
int cLens1::CopyPSFMapToImage(cImage& image,int PixelWidth,double Width) {
	// PixelWidth>0  �̂Ƃ��͉摜�̃s�N�Z������PixelWidth�ɕς���D
	//           <=0 �̂Ƃ���PSFMap�Ɠ����ɂȂ�D
	// Width>0  �̂Ƃ��͉摜�̎�����Width�ɕς���D
	//      <=0 �̂Ƃ��͉摜�̎�����PSFMap�Ɠ����ɂȂ�D
	// ������PSF�摜�Ƒ��̉摜�Ƃ̃R���{�����[�V�����̏����Ɏg�����Ƃ��ł���D
	int i,j,n,n0;
	double w0;

	image.Clear();
	n0=psf_map.GetPixels();
	n= PixelWidth<=0 ? n0 : PixelWidth;
	image.SetXpixels(n);
	image.SetYpixels(n);
	w0=(PSFMapWX()+PSFMapWY())/2;  // ���ς��Ƃ�
	image.SetPitch(w0/n0);
	for(i=1; i<=n; i++) for(j=1; j<=n; j++) {
		// PSFMap��image�̒����ɃR�s�[����
		image.SetIntensity(i,j,PSFMap(i-(n-n0)/2,j-(n-n0)/2));
	}
	// �摜�̎�����Width�ɍ��킹��
	if(Width<=0) Width=w0;
	image.Zoom(image.GetPitch()*n/Width);
	return n;
}
int cLens1::CopyPSFMapToImage(int PixelWidth,double Width) {
	return CopyPSFMapToImage(this->image,PixelWidth,Width);
}
int cLens1::SavePSFMapAsBmp(std::string filename,int PixelWidth,double Width,double gamma){
	cImage image;
	CopyPSFMapToImage(image,PixelWidth,Width);
	image.Normalize(255);
	return image.SaveAsBmp(filename,gamma);
}

double cLens1::PSFMapFiberEfficiency(double MFD,int DoublePass,int Incoherent){
	// PSF�����[�h�t�B�[���h���aMFD�̃V���O�����[�h�t�@�C�o�̒[�ʂ�
	// �Ǝ˂����Ƃ��̌������������߂�D
	// this->psf_map�ɑ΂��Čv�Z����̂ŁC
	// MakePSFMap()�̎��s����̂ݗL���D
	// DoublePass=true�̂Ƃ���PSF�̎��ȑ��ւɑ΂��Čv�Z����D
	//   �Ⴆ�΁COCT���w�n�̊��PSF��m��Ƃ��C��ꔽ�ˌ������ʔ��˂ł���Ƃ�
	//   ��ꔽ�ˌ��̃t�@�C�o�[�ł�PSF�͊��PSF�̎��ȑ��ւƍl������D
	//   ���̂Ƃ��CMFD�̓t�@�C�o�̎�MFD�Ɍ����{�����悶�����̂��w�肷��D
	// Incoherent��DoublePass=true�̂Ƃ��̂ݗL���ŁC
	// OCT�̗�ł͊�ꔽ�˂��C���R�q�[�����g�̏ꍇ�ɑΉ�����D
	cImage G,A,X;
	double max;
	cImage B=this->psf_map;
	const int DEBUG=0;

	if(DoublePass && Incoherent){ // �_�u���p�X�Ő܂�Ԃ��ʂ��U���ʂł���ꍇ
		B.GCenterToCenter();
		B.NormalizeTotalIntensity();
		if(DEBUG){
			X=B; X.Normalize(255); X.SaveAsBmp("1 psf.bmp");
		}

		G=B;  // ��f���C��f�s�b�`��������D
		G.SetGauss(MFD);
		G.NormalizeTotalIntensity();
		if(DEBUG){
			X=B; X.Normalize(255); X.SaveAsBmp("2 gauss mode.bmp");
		}

		max=sqabs((B*G).Total());  // ���˖ʂ̏d�S����̌������Ɋւ���t�@�C�o���������D
		                           // "���f�o�C�X�̂��߂̌������n�̊�b�Ɖ��p (�͖쌒��)" (3.1-7��)���D
		                           // ���ƂŋK�i���ɗp����D
		
		// A�ipsf��gauss�̏�ݍ��݁j���t�[���G�ϊ��ɂ��v�Z����D
		B=B.DFTAmplitude(1,1);
		if(DEBUG){
			X=B; X.Normalize(255); X.SaveAsBmp("3 F(psf).bmp");
		}
		G=G.DFTAmplitude(1,1);
		if(DEBUG){
			X=G; X.Normalize(255); X.SaveAsBmp("4 F(gauss).bmp");
		}
		A=B*G;
		if(DEBUG){
			X=A; X.Normalize(255); X.SaveAsBmp("5 F(psf)xF(gauss).bmp");
		}
		A=A.InvDFTAmplitude(1,1);
		if(DEBUG){
			X=A; X.Normalize(255); X.SaveAsBmp("6 invF(F(B)xF(G)).bmp");
		}
		
		// �ȏ��A�� "���f�o�C�X�̂��߂̌������n�̊�b�Ɖ��p (�͖쌒��)" (3.1-7��) �̕��q��
		// �ϕ��ł���D��Βl�̓����Ƃ邱�ƂŁCA[i][j]�͎U���ʏ�̊e�_����̌��ɑ΂��錋�������ƂȂ�D
		A.ToSquare();
		A.Normalize(max); // ��Ōv�Z�������˖ʏd�S����̌��ɑ΂��錋�������ɍ��킹��D

		B=this->psf_map;
		B.GCenterToCenter();
		B.NormalizeTotalIntensity();
		B.ToSquare();

		A=A*B;  // A���U���ʏ�̋��x�ŏd�ݕt������D
		return abs(A.TotalAmplitude());  // �C���R�q�[�����g�ɉ��Z����D
	}
	else{
		// ����ȊO�̏ꍇ  (���� Coherent �͖����j
		if(DoublePass) B=Convolution(B,B);
		B.GCenterToCenter();  // ���ꂪ�Ȃ���B�̏d�S�����ꂽ��ԂƂȂ�\��������C
							  // �Ƃ���B,G�̊g���肪�������Ƃ��d�S�̂���͑��ΓI��
							  // �傫�����̂ƂȂ邩��덷�̌����ƂȂ�D
		B.NormalizeTotalIntensity();

		G=B;  // ��f���C��f�s�b�`��������D
		G.SetGauss(MFD);
		G.NormalizeTotalIntensity();

		A=B*G;
		return sqabs(A.Total());  // "���f�o�C�X�̂��߂̌������n�̊�b�Ɖ��p (�͖쌒��)" (3.1-7��)���D
	}
}

int cLens1::PSFMapToShadow(double MaskArea){
	// PSFMap�𑊕�I�}�X�N�ł̕��z�ɕϊ�����D
	//   �E�}�X�N���Ȃ��Ƃ��X�N���[���͐U��1(���x1)�̋ψ�Ȍ����z�Ƃ���D
	//   �E�Ⴆ�΃�0.3�o�̊J���̉�ܑ��ɓK�p����΁C
	//     ���s���ɂ���0.3�̎Ռ���(FT�̃v���Y���̖A�Ȃ�)�̉e���v�Z�ł���D
	int i,j,n;
	double p,dx,dy;
	
	if(MaskArea<=0) return 0;
	p=MaskArea;   // �}�X�N���Ȃ��Ƃ��̑��ʏ�̐U�����ψ��1�Ƃ��C�Ռ������p���[ = MaskArea �Ƃ���D
	              //   ���Ȃ킿�C���ʂƌ����a���قȂ�Ƃ���Ƀ}�X�N��ݒu����Ƃ��́C
	              //   ���ʂɓ��e�����ʐς�ݒ肷�邱�ƁD
	              //   ��FFT(flatnesstester�ŃR�����[�e�B���O�����Y���㗬�ɐݒu����Ƃ�)
	n=psf_map.GetPixels();
	dx=this->PSFMapWX()/(n-1);
	dy=this->PSFMapWY()/(n-1);
	psf_map.NormalizeTotalIntensity(p/dx/dy);   // �p���[��p�ɂȂ�悤��PSF�}�b�v�̋K�i������
	n=psf_map.GetPixels();
	for(i=1; i<=n; ++i) for(j=1; j<=n; ++j){
		psf_map.SetComplexAmplitude(i,j, 1-psf_map.GetComplexAmplitude(i,j));  // Babinet�藝
		// "1-psf_map.Get... " �� "1" �� exp(i��) �łȂ��Ă悢�̂��H
		//    -> PSF��OPD�}�b�v����v�Z����邪�COPD�͎�����̌��H����0�i���Ȃ킿�ʑ���0�j�Ƃ��Ă���D
		//       �w�i�̈ʑ��͎�����̈ʑ��Ɠ����ƍl������̂ŁC"1" �ł悢�D
	}
	return n;
}

int cLens1::PSFMapAddAmplitude(double amp,double phase_deg){
	// PSFMap�̑S�̂ɋψ�ȕ��f�U����������D
	// ��F PSFMapToShadow��FT(flatness tester)�̃X�N���[����̖A�̉e���v�Z����D
	//      ����ɖ{�֐�����p�����Ċ��ȓ��ł̖A�̉e���v�Z����D
	int i,j,n;
	complex z;

	z=complex( amp*cos(phase_deg*PI/180), amp*sin(phase_deg*PI/180) );
	n=psf_map.GetPixels();
	for(i=1; i<=n; ++i) for(j=1; j<=n; ++j){
		psf_map.SetComplexAmplitude(i,j, psf_map.GetComplexAmplitude(i,j)+z); 
	}
	return n;
}

int cLens1::MakeOTFMap(cImage& B,double& wx,double& wy,double& rayleigh_nu_x,double& rayleigh_nu_y,
	               double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,
				   int FindPupil,int AberrationFree,int AdjustSph,std::string ZernikeFix)
{
	double dummy, temp;
	int n, size;
	cImage X;
	temp=PupilMargin;
	PupilMargin=1;     // OTF�ɐ܂肽���݂��������Ȃ��悤�ɂ���D
	n=MakePSFMap(X,wx,wy,dummy,dummy,0,yObj,xObj,defocus,ColorStart,ColorEnd,
	             FindPupil,0,0,AberrationFree,AdjustSph,ZernikeFix);
	PupilMargin=temp;  // cLens1�̃����oPupilMargin�����ɖ߂�
	if(n>0){
		X.DFTIntensity(0,0);
		X.Normalize(1);
		B=X;
		size=n/2;
		// n�_���U�t�[���G�ϊ��̔g���́C�ϊ�����i=0�`n-1�ɑΉ����钷����L�Ƃ���ƁC
		//   k=0   -> ��
		//   k=1   -> L*n/(n-1)
		//   k=2   -> L*n/(n-1)/2
		//   ...
		//   k=n-1 -> L*n/(n-1)/(n-1)
		//   �� ���g�� = 0 �` (n-1)^2/n/L
		// �ƂȂ�D
		wx=(double)(n-1)*(double)(n-1)/(double)n/wx/((double)(n-1)/(double)(size-1));
		wy=(double)(n-1)*(double)(n-1)/(double)n/wy/((double)(n-1)/(double)(size-1));
		
		{
			// Rayleigh�̕���\�ɑΉ�������g�������߂�D
			const double MTFC=0.08940; // 0.61��/NA�ɑΉ����闝�z�����Y��MTF
			int i;
			cXYList li;
			for(i=1; i<=size; i++){
				li.AddData(i,B.GetAmplitude(1,i));
				rayleigh_nu_x=wx*(li.TransitionPoint(MTFC)-1)/((double)(size-1));
			}
			li.RemoveAll();
			for(i=1; i<=size; i++){
				li.AddData(i,B.GetAmplitude(i,1));
				rayleigh_nu_y=wy*(li.TransitionPoint(MTFC)-1)/((double)(size-1));
			}
		}

		return size;
	}
	else{
		wx=wy=0;
		return 0;
	}
}

int cLens1::MakeOTFMap(double yObj,double xObj,double defocus,int ColorStart,int ColorEnd,
					   int FindPupil,int AberrationFree,int AdjustSph,std::string ZernikeFix){
	return MakeOTFMap(otf_map,otf_map_w_x,otf_map_w_y,otf_map_rayleigh_nu_x,otf_map_rayleigh_nu_y,
	                  yObj,xObj,defocus,ColorStart,ColorEnd,
					  FindPupil,AberrationFree,AdjustSph,ZernikeFix);
}

double cLens1::OTFMapMTF(int i,int j){
	return otf_map.GetAmplitude(i,j);
}

double cLens1::OTFMapWX(){
	return otf_map_w_x;
}
double cLens1::OTFMapWY(){
	return otf_map_w_y;
}

double cLens1::OTFMapRayleighNuX(){
	return otf_map_rayleigh_nu_x;
}
double cLens1::OTFMapRayleighNuY(){
	return otf_map_rayleigh_nu_y;
}
/*
double cLens1::yVignetting(double yObj,int findpupil){
	// �y���̎����͌��z �����̏㉺�[���ʁX�̖ʂł�����ꍇ��z�肵�Ă��Ȃ��D

	//   findpupil==0 �̂Ƃ� :  �����f�ʍ��� / ���˓����a
	//   findpupil==1 �̂Ƃ� :  �����f�ʍ��� / stop�ʂ̗L���a�Ō��܂�����f�ʍ���
	// ��Ԃ�
	//
	// FindAThruRay()�Ŏ��Ԃ��Ƃ��Ȃ��悤�ɁC�ꎞ�I��StopDominate=1�Ƃ��Čv�Z����D
	// StopDominate=0�̂܂܂Ōv�Z������@�ł́C�L���a�ɗ]�T���Ȃ��C���������ׂ��Ƃ�
	// (�L�p�̐ڊ჌���Y�Ȃ�)�C�Ŏ��O�̎������œK������ꍇ�C���΂��Ό������L���a�Ŋ��S�ɎՕ�����邽�߁C
	// FindAThruRay()�Ŏ��Ԃ��Ƃ�ꎩ���݌v���x���Ȃ�D
	point ymax,ymin,xmax,xmin;
	int i,j,temp;
	double y1,y2, a,min,result;

	temp=StopDominate;
	StopDominate=1;

	if(findpupil!=0 && stop==0){
		result=0;
	}
	else{
		min=1e30;

		for(j=1; j<=cn; j++){
			FindPupilEdge(ymax,ymin,xmax,xmin,yObj,0,j,findpupil,1);
			for(i=1; i<=k; i++){
				if(EAy(i)>0){
					y1=RayPosY(i,yObj,0,"",ymax.y,ymax.x,j);
					y2=RayPosY(i,yObj,0,"",ymin.y,ymin.x,j);

					if( (y1-EAy(i)/2)*(y2-EAy(i)/2)<0 && (y1+EAy(i)/2)*(y2+EAy(i)/2)<0 ){
						// y1,y2 ��EAy(i)/2�C-EAy(i)/2 �̗������͂���
						a=fabs( EAy(i)/(y1-y2) );
					}
					if( (y1-EAy(i)/2)*(y2-EAy(i)/2)<0 ){
						// y1,y2�� EAy(i)/2���͂���
						a= y1>y2 ? (EAy(i)/2-y2)/(y1-y2) : (EAy(i)/2-y1)/(y2-y1);
					}
					else if( (y1+EAy(i)/2)*(y2+EAy(i)/2)<0  ){
						// y1,y2�� -EAy(i)/2���͂���
						a= y1>y2 ? (y1+EAy(i)/2)/(y1-y2) : (y2+EAy(i)/2)/(y2-y1);
					}
					else if( fabs(y1)<=EAy(i)/2 && fabs(y2)<=EAy(i)/2 ){
						// y1,y2�Ƃ�EAy(i)�̓�
						a=1+( EAy(i)/2-Max(fabs(y1),fabs(y2)) )/abs(y1-y2);
						// �����݌v�Ő���ł���悤�C�ꗥa=1�Ƃ��Ȃ��ŗ]�T�x���܂񂾒l��Ԃ��D
						// abs(y1-y2)�Ŋ���̂� y1,y2���L����̉����͂��ޗ̈�Ɗ֐��̌��z�����킹�邽�߁D
					}
					else{
						// y1,y2��EAy(i)�̊O
						a=-(Min(fabs(y1),fabs(y2))-EAy(i)/2)/abs(y1-y2);
						// �����݌v�Ő���ł���悤�C�ꗥa=0�Ƃ��Ȃ��ŁC���̒l��Ԃ��D
						// abs(y1-y2)�Ŋ���̂� y1,y2���L����̉����͂��ޗ̈�Ɗ֐��̌��z�����킹�邽�߁D
					}

				if(a<min) min=a;
				}
			}
		}
		
		result=min;
	}

	StopDominate=temp;
	return result;
}
*/
double cLens1::yVignetting(double yObj,int findpupil){
	//   findpupil==0 �̂Ƃ� :  �����f�ʍ��� / ���˓����a
	//   findpupil==1 �̂Ƃ� :  �����f�ʍ��� / stop�ʂ̗L���a�Ō��܂�����f�ʍ���
	// ��Ԃ�
	//
	// FindAThruRay()�Ŏ��Ԃ��Ƃ��Ȃ��悤�ɁC�ꎞ�I��StopDominate=1�Ƃ��Čv�Z����D
	// StopDominate=0�̂܂܂Ōv�Z������@�ł́C�L���a�ɗ]�T���Ȃ��C���������ׂ��Ƃ�
	// (�L�p�̐ڊ჌���Y�Ȃ�)�C�Ŏ��O�̎������œK������ꍇ�C���΂��Ό������L���a�Ŋ��S�ɎՕ�����邽�߁C
	// FindAThruRay()�Ŏ��Ԃ��Ƃ�ꎩ���݌v���x���Ȃ�D
	point ymax,ymin,xmax,xmin;
	int i,j,temp;
	double y1,y2, edge, a1,a2, max1,max2, result;

	temp=StopDominate;
	StopDominate=1;

	if(findpupil!=0 && stop==0){
		result=0;
	}
	else{
		max1=max2=-1e30;

		for(j=1; j<=cn; j++){
			if(findpupil==0 || i!=stop){   // stop�ʂ�����ƁC�߂�l��1�ȏ�ɂȂ�Ȃ��Ȃ��Ă��܂��D
				FindPupilEdge(ymax,ymin,xmax,xmin,yObj,0,j,findpupil,1);
				for(i=1; i<=k; i++){
					if(EAy(i)>0){
						y1=RayPosY(i,yObj,0,"",ymax.y,ymax.x,j);  // �����̂������[
						y2=RayPosY(i,yObj,0,"",ymin.y,ymin.x,j);  //      �V
						
						edge = y1>y2 ? EAy(i)/2 : -EAy(i)/2;
						a1= (y1-edge)/(y1-y2);
						if(a1>max1) max1=a1;  // y1���̂����ʁD�]�T������Ƃ��͕��ɂȂ�D

						edge = y2>y1 ? EAy(i)/2 : -EAy(i)/2;
						a2= (y2-edge)/(y2-y1);
						if(a2>max2) max2=a2;  // y2���̂����ʁD�]�T������Ƃ��͕��ɂȂ�D
					}
				}
			}
		}

		if(max1>0 && max2>0){         // ���������Ă���
			result= 1-max1-max2;
		}
		else if(max1>0 && max2<=0){   // y1�����������Ă���
			result= 1-max1;
		}
		else if(max1<=0 && max2>0){   // y2�����������Ă���
			result= 1-max2;
		}
		else{                         // ���������Ă��Ȃ�(1�ȏ�̒l��Ԃ�)
			result= 1-Max(max1,max2);
		}
	}

	StopDominate=temp;
	return result;
}

double cLens1::yVignetting(int findpupil/*=1*/){
	return yVignetting(this->yObjectMax,findpupil);
}

double cLens1::AbyA0(){
	// ���A���S���Y���i�����̌�����ǐՂ��C�ʉߌ������^�ǐՌ����� ��Ԃ��j
	// ��菈�����x���x�����C�����݌v�� AbyA0 ���ǂ���������D
	double A0=0, A=0;
	int findpupil;
	cLens1 X=*this;   // �R�s�[���Ȃ��� *this �𒼐ڑ��삵�Ă��������x�͂���قǕς��Ȃ��D
	X.nSpot=5;        // �����l20���Ǝ��Ԃ������肷����D

	findpupil= ea_is_defined()==0 ? 0 : 1;	
	X.MakeSpot(0,0,11,0,1,1,1,0,0,0,0);          // ��������̑�1�ʂł̒f�ʗ֊s
	A0=X.SpotArea();
	X.MakeSpot(yObjectMax,0,11,0,1,1,1,0,0,0,0); // ���O�����̑�1�ʂł̒f�ʗ֊s
	A=X.SpotArea();
	if( A0!=0 ) return A/A0; else return 0;

	/*  �y���A���S���Y���z
	int iy,ix, makecoordinate;
	double yPupil,xPupil, a,b;
	double A0=0; double A=0;
	int n=this->nSpot;

	make_coordinate(0); makecoordinate=0;

	a= EPDx==0 ? EPD/2:EPDx/2;  // x���������a
	b= EPD/2;                   // y���������a

	for(iy=1; iy<=n; iy++) for(ix=1; ix<=n; ix++){
		yPupil=-b+b*2/n*(iy-0.5);
		xPupil=-a+a*2/n*(ix-0.5);
		if( yPupil*yPupil/b/b+xPupil*xPupil/a/a<=1 ){
			if( RayTrace(0,0,yPupil,xPupil,0,1,1,1,0,0,
				         0,vector<complex>(),makecoordinate) == 0 ){
				A0++;
			}
			if( RayTrace(yObjectMax,xObjectMax,yPupil,xPupil,0,1,1,1,0,0,
				         0,vector<complex>(),makecoordinate) == 0 ){
				A++;
			}
		}
	}
	if( A0!=0 ) return A/A0; else return 0;
	*/
}

double cLens1::Transmittance(double yObject,double xObject,int IsAreaSource,int IsLambert,
                             int MakeSurfaceList,list<int>& surfaces){
	// ���ߗ��i���˓��ɓ����������ɑ΂���n��ʉ߂ł����������̊����j���C
	// ��1�g���ɂ��Čv�Z����D

	int iy,ix, i;
	double yPupil,xPupil;
	point dr;
	double A0=0; double A=0;
	double t;
	int n=this->nSpot;
	vector<complex> dummy;
	double a,b,dx,dy;
	
	a= EPDx==0 ? EPD/2:EPDx/2;  // x���������a
	b= EPD/2;                   // y���������a
	dx=a*2/n;
	dy=b*2/n;
	make_coordinate(0);

	for(iy=1; iy<=n; ++iy) for(ix=1; ix<=n; ++ix){
		yPupil=-b+dy*(iy-0.5);
		xPupil=-a+dx*(ix-0.5);
		if( yPupil*yPupil/b/b+xPupil*xPupil/a/a<=1 ){
			dr= IsAreaSource ? dObject() : point(0,0);
			switch( i=RayTrace(yObject+dr.y,xObject+dr.x,yPupil,xPupil,0,1,1,1,IsLambert,0,
				               0,dummy,0) )  // ��1�g���Œǐ�
			{
				case 0: 
					t=ApodizationIntensity();
					A0+=t;
					A+= t;
					break;
				case INVALID: 
					break;
				default:
					A0+=ApodizationIntensity();
					if( MakeSurfaceList!=0 && NotThruSurf(i)!=0 ){
						surfaces.AddTail(NotThruSurf(i)); // �������ʉ߂��Ȃ��ʂ̃��X�g���쐬����
					}
					break;
			}
		}
	}
	if( MakeSurfaceList!=0 ){
		// ���ߗ��̒Ⴂ�n�ő����̌�����ǐՂ����ꍇ�C���X�g�̃T�C�Y���傫���Ȃ���Ɏ��Ԃ�������D
		surfaces.EraseDuplicates();
		surfaces.Sort();
	}
	if( A0!=0 ) return A/A0; else return 0;
}

double cLens1::Transmittance(double yObject,double xObject,int IsAreaSource,int IsLambert){
	list<int> dummy;
	return Transmittance(yObject,xObject,IsAreaSource,IsLambert,0,dummy);
}

double cLens1::Ton(){
	return Transmittance(0,0,1,0);
}

double cLens1::Toff(){
	return Transmittance(yObjectMax,xObjectMax,1,0);
}

std::string cLens1::NotThruSurfaces(double yObject,double xObject,int IsAreaSource){
	list<int> li;
	std::string s;
	int i,data;
	char buf[100];
	Transmittance(yObject,xObject,IsAreaSource,0,1,li);
	for(i=1; i<=li.GetSize(); i++){
		li.GetData(data,i);
		sprintf(buf,"%d\n", data); s+=buf;
	}
	return s;
}

double cLens1::tLack(){
	// d�̕s��
	int i; double lack_c,lack_e;  double lack=0;
	for(i=1; i<=k-1; i++){
		if(N(i,1)==1 || N(i,1)==-1){
			lack_c=CenterThickness(i)-TcMinAir;
			lack_e=koba(i)-TeMinAir;
		}
		else{
			lack_c=CenterThickness(i)-TcMinGlass;
			lack_e=koba(i)-TeMinGlass;
		}
		if(lack_c<0 && dVariable(i)!=0) lack+=lack_c;
		if(lack_e<0 && (rVariable(i)!=0 || rVariable(i+1)!=0 || dVariable(i)!=0)) lack+=lack_e;
	}
    return lack;
}

cLens1 cLens1::perturbed(int IsEndNotUni/*=0*/){
	int i;
//	double AsLimit;
	double ea;
	cLens1 x=*this;
/*	for(i=1; i<=k; i++) {
		if(x.asph_type(i)==SPH) x.asph_type(i)=CONIC;
		x.Newton(i)=Random(x.NewtonTol(i),IsEndNotUni);
		AsLimit= fabs(x.dNewtonTol(i))<fabs(x.Newton(i)*2) ? 
		         fabs(x.dNewtonTol(i)) : fabs(x.Newton(i)*2);
		x.dNewton(i)=Random(0,AsLimit,IsEndNotUni);
		x.Axis(i)=Random(0,360,0);
	}
*/	
	for(i=1; i<=k; i++) {
		ea=ea_max(i);
		x.r(i)=NewtonToR( r(i), ea, FRW, Random(NewtonTol(i),IsEndNotUni) );
	}
	for(i=1; i<=k-1; i++) {
		x.d(i)+=Random(delta_d(i),IsEndNotUni);
	}
	for(i=0; i<=k+1; ++i) {
		x.dN(i)=Random(x.dN(i),IsEndNotUni);
	}
	return x;
}

void cLens1::perturbe_this(int IsEndNotUni/*=0*/){
	*this=perturbed(IsEndNotUni);
}

double cLens1::PerturbDxDy(int i,double dr/*=0*/){
	// r��0�`dr�̈�l���z�C�ӂ�0�`2�΂̈�l���z�̂Ƃ���
	// (dx,dy)=(rcos��,rsin��)�̃T���v����Ԃ��D
	double r,phi;

	Srand();
	r=Random(0,dr,0);
	phi=Random(0,2*PI,0);
	dx(i)+=r*cos(phi);
	dy(i)+=r*sin(phi);
	// �Ⴆ�΃����e�J�����V�~�����[�V�����Ŏ��̌����ƌ��̌�����ʁX�ɐݒ肷�邽�߂ɁC
	// dx(i)=.. �ł͂Ȃ� dx(i)+=.. �Ƃ���D
	return 0;
	// ������͓��Ńv���p�e�B�Ɠ��������Ƃ��邽�߂�0��Ԃ��D
	// ������̗͂���F
	//   (1)�����l���擾 (����dr���ȗ��\�ł��邱�Ƃ��K�v�j
	//   (2)�����l+�ϕ� �ɐݒ�
	//   (3)�]���֐����v�Z
	//  �X�e�b�v(1)�Ŗ߂�l���K�v�D
}

double cLens1::PerturbXYObjectMax(double dr/*=0*/){
	// PerturbDxDy �Ɠ��l
	double r,phi;

	Srand();
	r=Random(0,dr,0);
	phi=Random(0,2*PI,0);
	xObjectMax+=r*cos(phi);
	yObjectMax+=r*sin(phi);
	return 0;
}



void cLens1::AddGroup(int i1,int i2){
	if( 0<=i1 && i1<i2 && i2<=k+1 ){
		groups.AddTail( group(i1,i2) );
	}
}

void cLens1::GroupsClear(){
	groups.RemoveAll();
}

int cLens1::RowsNumber(){
	int n;

	n=k;
	groups.EraseDuplicates(); n+=groups.GetSize();
	if(rObj()!=0)   n+=1;
	if(rImage()!=0) n+=1;
	return n;
}

int cLens1::i1Row(int i){
	int ii; group g;
	list<group> rows;

	if(rObj()!=0)   rows.AddTail( group(0,0) );
	if(rImage()!=0) rows.AddTail( group(k+1,0) );
	for(ii=1; ii<=k; ++ii) rows.AddTail( group(ii,0) );
	for(ii=1; ii<=groups.GetSize(); ++ii) {
		groups.GetData(g,ii);
		rows.AddTail(g);
	}
	rows.EraseDuplicates();
	rows.Sort();
	rows.GetData(g,i);
	return g.i1;
}

int cLens1::i2Row(int i){
	int ii; group g;
	list<group> rows;

	if(rObj()!=0)   rows.AddTail( group(0,0) );
	if(rImage()!=0) rows.AddTail( group(k+1,0) );
	for(ii=1; ii<=k; ++ii) rows.AddTail( group(ii,0) );
	for(ii=1; ii<=groups.GetSize(); ++ii) {
		groups.GetData(g,ii);
		rows.AddTail(g);
	}
	rows.EraseDuplicates();
	rows.Sort();
	rows.GetData(g,i);
	return g.i2;
}

std::string cLens1::strRow(int i){
	int ii; group g;
	list<group> rows;

	if(rObj()!=0)   rows.AddTail( group(0,0) );
	if(rImage()!=0) rows.AddTail( group(k+1,0) );
	for(ii=1; ii<=k; ++ii) rows.AddTail( group(ii,0) );
	for(ii=1; ii<=groups.GetSize(); ++ii) {
		groups.GetData(g,ii);
		rows.AddTail(g);
	}
	rows.EraseDuplicates();
	rows.Sort();
	rows.GetData(g,i);
	return g.str(this->k);
}

void cLens1::MakeThinDoubletOnNud2(double beta,
	                               std::string color1,std::string color2,std::string color3, 
	                               std::string glass1,
	                               double L,double B,double A, double nud2, double f)
{
	*this=ThinDoubletOnNud2(beta, color1,color2,color3, glass1, L,B,A, nud2);
	AdjustFocalLength(f,0);
	SetM(beta);
}
void cLens1::MakeThinDoublet(double beta,
	                        std::string color1,std::string color2,std::string color3, 
							std::string glass1,std::string glass2,
	                        double L,double B, double f)
{
	*this=ThinDoublet(beta, color1,color2,color3, glass1,glass2, L,B);
	AdjustFocalLength(f,0);
	SetM(beta);
}
void cLens1::MakeThinSeparateDoublet(double beta,
	                                 std::string color1,std::string color2,std::string color3, 
							         std::string glass1,std::string glass2,
	                                 double L,double B,double A, double f)
{
	*this=ThinSeparateDoublet(beta, color1,color2,color3, glass1,glass2, L,B,A);
	AdjustFocalLength(f,0);
	SetM(beta);
}
void cLens1::MakeThinTriplet(int FrontIsCemented,double beta,
	                        std::string color1,std::string color2,std::string color3, 
	                        std::string glass1,std::string glass2,std::string glass3,
	                        double Phi3Phi1Ratio, double L,double B,double A, double f)
{
	*this=ThinTriplet(FrontIsCemented, beta, 
	                  color1,color2,color3, glass1,glass2,glass3, Phi3Phi1Ratio,L,B,A);
	AdjustFocalLength(f,0);
	SetM(beta);
}
void cLens1::MakeDoubleThinDoublet(double beta,double t,double e,
	                              std::string color1,std::string color2,std::string color3,
	                              std::string glass1,std::string glass2,
				                  std::string glass3,std::string glass4,
			                      double Length, double L,double T,double CM,double SA,
						          int SolutionNo, double f,double d1,double d2,double d4,double d5)
{
	// t,e : values with f=1
	// Length : value with s=��
	*this=DoubleThinDoublet(beta,t,e, color1,color2,color3, glass1,glass2,glass3,glass4,
	                        Length,L,T,CM,SA,SolutionNo);
	double fA,fB;
	AdjustFocalLength(f,0);
	fA=this->f(1,3); fB=this->f(4,6);
	d(1)=d1; d(2)=d2; d(4)=d4; d(5)=d5;
	AdjustFocalLength(fA,1,3,1,0,0); AdjustFocalLength(fB,4,6,1,0,0);
	d(3)+=Trimed(1,3).delta1();
	d(3)-=Trimed(4,6).delta();
	SetM(beta);
	if(this->t<LN) this->t=Trimed(1,3).delta()+this->t;
}
void cLens1::MakeTriplet(double beta, double SigmaPhi,double Length,double Kappa,double Epsilon,
	                    std::string color1,std::string color2,std::string color3,
	                    std::string glass1,std::string glass2,std::string glass3,
                        double CM,double AS,double DS, double f,double d1,double d3,double d5)
{
	*this=Triplet(beta, SigmaPhi,Length,Kappa,Epsilon, color1,color2,color3, glass1,glass2,glass3, 
	              CM,AS,DS, f,d1,d3,d5);
}


double cLens1::BeamFullDiv(double waist_dia,double wl_nm,double N,double M2){
	return waist_dia==0 ? 0 : atan( (2/PI)*(wl_nm*M2/1000000/N)/waist_dia ) *2;
	// �g�g���ɂ�M2*�ɂƂ���Ώ]����TEM00�K�E�V�A���r�[���ɑ΂��č��ꂽ
	//   �V�~�����[�^�����̂܂܊��p�ł���h
	//  (���[�U�r�[���i������̊�b ������� ���[�U�[�r�[���i������̊�b ��26����10���j
}

double cLens1::WaistDia(double beam_full_div,double wl_nm,double N,double M2){
	return beam_full_div==0 ? 0 : (2/PI)*(wl_nm*M2/1000000/N)/tan(beam_full_div/2);
}

double cLens1::RaylaySemiDOF(double waist_dia,double wl_nm,double N,double M2){
	return (PI/4)*waist_dia*waist_dia/(wl_nm*M2/1000000/N);
}

double cLens1::GetM2(){
	return M2;
}
void cLens1::SetM2(double value){
	if(value>=1) M2=value;
}

double cLens1::GetWaistDiaIn(){
	if(waist_dia_in==0){
		waist_dia_in=WaistDia(beam_full_div_in,wl[1],N(0,1),M2);
	}
	return waist_dia_in;
}
void cLens1::SetWaistDiaIn(double value){
	waist_dia_in=value;
}

double cLens1::GetBeamFullDivIn(){
	if(waist_dia_in==0){
		waist_dia_in=WaistDia(beam_full_div_in,wl[1],N(0,1),M2);
		return beam_full_div_in;
	}
	else{
		beam_full_div_in=BeamFullDiv(waist_dia_in,wl[1],N(0,1),M2);
		return beam_full_div_in;
	}
}
void cLens1::SetBeamFullDivIn(double value){
	beam_full_div_in=value;
}

double cLens1::RaylaySemiDOFIn(){
	return RaylaySemiDOF(waist_dia_in,wl[1],N(0,1),M2);
}

complex cLens1::q0(){
	// �Q�l�� �g���@��̌��w�U(����)�h
	complex q; 
	const complex i=complex(0,1);
	q=i*(waist_dia_in*waist_dia_in/4)*PI/(wl[1]*M2/1000000/N(0,1));
	q+=WaistPosIn;
	return q;
}

complex cLens1::q1(){
	complex q=q0()/N(0,1); 
	const complex i=complex(0,1);
	AMatrixCalc(1,k,1);
	q=(A2[2]+A2[1]*q)/(A2[4]+A2[3]*q);  // (7.1.27)
	q=q*N(k,1);
	return q;
}

double cLens1::WaistPosOut(){
	// (7.1.15)���r�[���E�G�X�g�ł�q�͏������ł��邱�Ƃ� dq/dz=-1 (q'=q-z) ���C
	return Re(q1());
}

double cLens1::WaistDiaOut(){
	// (7.1.15)���r�[���E�G�X�g�ł�q�͏������ł��邱�Ƃ� dq/dz=-1 (q'=q-z) ���C
	// �r�[���E�G�X�g�ł� q0=Im(q) �ł���D
	// ���������� (7.1.16)���C
	return 2*sqrt( Im(q1())*(wl[1]*M2/1000000/N(k,1))/PI );
}

double cLens1::BeamDiaIn(){
	const complex i=complex(0,1);
	complex q=q0()-( fabs(s)>=LN ? 0 : s );
	return 2*sqrt( -(wl[1]*M2/1000000/N(k,1))/PI/Im(1/q) );  // (7.1.15)
}

double cLens1::BeamDiaOut(double defocus){
	const complex i=complex(0,1);
	complex q=q1()-dk(defocus);
	return 2*sqrt( -(wl[1]*M2/1000000/N(k,1))/PI/Im(1/q) );  // (7.1.15)
}

double cLens1::BeamFullDivOut(){
	return BeamFullDiv(WaistDiaOut(),wl[1],N(k,1),M2);
}

double cLens1::RaylaySemiDOFOut(){
	return RaylaySemiDOF(WaistDiaOut(),wl[1],N(k,1),M2);
}

int cLens1::PolarizationTrace(double yObj,double xObj,
							  std::string SetRay,int findpupil,double yPupil,double xPupil,
                              int j, double a,double b,double phi_deg)
{
	point pupil;
	double a1,a2,delta_deg;
	matrix<double> T;
	vector<complex> E0;

	pupil=point(xPupil,yPupil);
	if( !FindRay(pupil,yObj,xObj,SetRay,findpupil,j) ) return 0;
	cPolarization::Analyze(a1,a2,delta_deg,a,b,phi_deg);
	E0=vector<complex>( complex(a1,0), complex(a2*cos(delta_deg*PI/180),a2*sin(delta_deg*PI/180)), 0 ); 
	if( RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,1,0,0) ){
		return 0;
	}
	else{
		T=Tmatrix( vector<double>(X1[0],Y1[0],Z1[0]) );
		E0=inv(T)*E0;
		RayTrace(yObj,xObj,pupil.y,pupil.x,0,j,0,1,0,0, 1,E0);
		return 1;
	}
}

vector<complex> cLens1::GetE(int i){
	// i�ʂւ̓��˓d��i���O�̌����ǐՂ��K�v�j
	if(0<=i && i<=k+1){
		return Tmatrix(vector<double>(X[i],Y[i],Z[i]))*E[i];  // z���������̐i�s�����Ɍ�����
	}
	else{
		return vector<complex>(0,0,0);
	}
}

double cLens1::Ex(int i){
	return abs(GetE(i).x);  // x�����̓d�ꋭ�x
}

double cLens1::Ey(int i){
	return abs(GetE(i).y);  // y�����̓d�ꋭ�x
}

double cLens1::Edelta(int i){
	return (arg(GetE(i).y)-arg(GetE(i).x))*180/PI;  // Ex�̈ʑ���0�̂Ƃ���Ey�̈ʑ�(deg)
}

vector<complex> cLens1::GetE1(int i){
	// i�ʂ���̏o�˓d��i���O�̌����ǐՂ��K�v�j
	if(0<=i && i<=k+1){
		return Tmatrix(vector<double>(X[i],Y[i],Z[i]))*E1[i];  // z������ˌ����̐i�s�����Ɍ�����
	}
	else{
		return vector<complex>(0,0,0);
	}
}
	
void cLens1::Ellipse(int i,double &A,double &B,double &Phi_deg){
	// i�ʂ̓��ˑ��d��ɂ��Čv�Z����.
	// A=���������d�ꋭ�x�CB=�Z�������d�ꋭ�x�CPhe_deg=x���ɑ΂��钷���̊p�x
	matrix<double> T;
	vector<complex> E;
	double a1,a2,delta;

	if( i<0 || k+1<i ) i=k+1;
	T=Tmatrix( vector<double>(X[i],Y[i],Z[i]) );  // z������ˌ����̐i�s�����Ɍ�����
	E=T*this->E[i];
	a1=abs(E.x); a2=abs(E.y); delta=(arg(E.y)-arg(E.x))*180/PI;
	cPolarization::EllipseShape(A,B,Phi_deg, a1,a2,delta);
}

int cLens1::Ellipse(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
					int i,int j, 
					double a,double b,double phi_deg, double& A,double& B,double& Phi_deg){
	if(!PolarizationTrace(yObj,xObj,SetRay,findpupil,yPupil,xPupil,j,a,b,phi_deg)){
		return 0;
	}
	else{
		Ellipse(i,A,B,Phi_deg);
		return 1;	
	}
}

double cLens1::EllipseA(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
                        int i,int j, double a,double b,double phi_deg){
	// �ȉ~�Ό��̒�����
	double A,B,Phi;

	if(Ellipse(yObj,xObj,SetRay,findpupil,yPupil,xPupil,i,j,a,b,phi_deg,A,B,Phi)){
		return A;
	}
	else{
		return 0;
	}
}

double cLens1::EllipseB(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
                        int i,int j, double a,double b,double phi_deg){
	// �ȉ~�Ό��̒Z����
	double A,B,Phi;

	if(Ellipse(yObj,xObj,SetRay,findpupil,yPupil,xPupil,i,j,a,b,phi_deg,A,B,Phi)){
		return fabs(B);
	}
	else{
		return 0;
	}
}

double cLens1::EllipseRatio(int i){
	// �ȉ~�Ό��̒������ɑ΂���Z�����̓d��U����
	double A,B,Phi;
	Ellipse(i,A,B,Phi);
	return fabs(B)/A;
}

double cLens1::EllipseRatio(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
                            int i,int j, double a,double b,double phi_deg){
	double A,B,Phi;

	if(Ellipse(yObj,xObj,SetRay,findpupil,yPupil,xPupil,i,j,a,b,phi_deg,A,B,Phi)){
		return fabs(B)/A;
	}
	else{
		return 0;
	}
}

double cLens1::EllipseIntensityRatio(int i){
	// �ȉ~�Ό��̒������ɑ΂���Z�����̓d�ꋭ�x��
	double x;

	x=EllipseRatio(i);
	return x*x;
}

double cLens1::EllipsePhi(int i){
	// �ȉ~�Ό���x������Ƃ��钷������
	double A,B,Phi;
	Ellipse(i,A,B,Phi);
	return Phi;
}

double cLens1::EllipsePhi(double yObj,double xObj,std::string SetRay,int findpupil,double yPupil,double xPupil,
                          int i,int j, double a,double b,double phi_deg){
	double A,B,Phi;

	if(Ellipse(yObj,xObj,SetRay,findpupil,yPupil,xPupil,i,j,a,b,phi_deg,A,B,Phi)){
		return Phi;
	}
	else{
		return 0;
	}
}


cShapes* cLens1::CurrentShapes(int n){
	if(n<SHAPES_SIZE) current_shapes=&shapes[n];
	return current_shapes;
}
void cLens1::ShapesClear() { current_shapes->Clear(); }
void cLens1::ShapesViewTopLeft(double x,double y) { current_shapes->SetViewTopLeft(x,y); }
void cLens1::ShapesViewWidth(double xw,double yw) { current_shapes->SetViewWidth(xw,yw); }
int cLens1::ShapesSize() const { return current_shapes->Size(); }
int cLens1::ShapeKind(int i) { return current_shapes->Kind(i); }
double cLens1::ShapeX1(int i) { return current_shapes->X1(i); }
double cLens1::ShapeY1(int i) { return current_shapes->Y1(i); }
double cLens1::ShapeX2(int i) { return current_shapes->X2(i); }
double cLens1::ShapeY2(int i) { return current_shapes->Y2(i); }
int cLens1::ShapeLineStyle(int i) { return current_shapes->LineStyle(i); }
std::string cLens1::ShapeText(int i) { return current_shapes->Text(i); }
double cLens1::ShapeFontSize(int i) { return current_shapes->FontSize(i); }
long cLens1::ShapeColor(int i) { return current_shapes->Color(i); }
void cLens1::ShapesZoom(double m) { current_shapes->Zoom(m); }
void cLens1::ShapesTranslate(double dx,double dy) { current_shapes->Translate(dx,dy); }
void cLens1::ShapesAllView() { current_shapes->AllView(); }
int cLens1::ShapesSaveAsBmp(std::string filename,int yPixels){
	return current_shapes->SaveAsBmp(filename,yPixels);
}

int cLens1::SAGraph() { return SA_GRAPH; }

void cLens1::MakeSAGraph(int FindPupil,double yPupilMax,double FullScale,int PupilSA,int AbeTheory/*=0*/){
	int i,j;
	point p0;
	double **sa, *sc;
	double max=0;
	int xdiv;
	const int YDIV=20;
	const double W=80, H=60;
	double xPupil;
	cShapes* ps=CurrentShapes(SA_GRAPH);

	ps->Reset();
	sa=new double*[YDIV+1]; for(i=0; i<=YDIV; ++i) sa[i]=new double[cn+1]; 
	sc=new double[YDIV+1];

	if(PupilSA){
		push();
		SwapObjPupil();
	}

	if(yPupilMax==0) yPupilMax=this->EPD/2;
	
	if(FindPupil || yPupilMax==0) {
		point ymax,ymin,xmax,xmin;
		FindPupilEdge(ymax,ymin,xmax,xmin,0,0,1,1);
		yPupilMax=ymax.y; xPupil=ymax.x;
		p0=(ymax+ymin)/2;
	}
	else {
		xPupil=0;
		p0=point(0,0);
	}
	for(j=1; j<=cn; ++j) for(i=0; i<=YDIV; ++i) {
		sa[i][j]=LSA(yPupilMax*i/YDIV,xPupil, p0.y,p0.x, j,AbeTheory);
		max= fabs(sa[i][j])>max ? fabs(sa[i][j]) : max;
	}
	for(i=0; i<=YDIV; ++i) {
		sc[i]=SC(yPupilMax*i/YDIV,xPupil);
		max= fabs(sc[i])>max ? fabs(sc[i]) : max;
	}

	if(FullScale==0){
		max=fceil(max);
		xdiv=(int)fceil_mantissa(max);
	}
	else{
		max=FullScale;  // �ۂ߂Ȃ�
		if(FullScale==fceil(FullScale)){
			xdiv=(int)fceil_mantissa(FullScale);
		}
		else{
			xdiv=4;
		}
	}

	ps->SetWindow(-50,-20,50,80);
	ps->AddLine(-W/2,0, W/2,0);
	for(i=-xdiv; i<=xdiv; ++i) {
		ps->AddLine((W/2/xdiv)*i,0.7, (W/2/xdiv)*i,-0.7);
	}
	ps->AddLine(0,H, 0,-2);
	for(i=1; i<=10; ++i) {
		ps->AddLine(-0.7,(H/10)*i, 0.7,(H/10)*i);
	}
	for(j=cn; j>=1; --j) {
		// for(j=1; j<=cn; ++j) �ɂ���ƁC�O���t���d�Ȃ����Ƃ�
		// ��1�g���̐F����cn�g���̐F�ɉB����Ă��܂��D
		long col=cSpectrum::ApproximateColor(wl[j]);
		for(i=0; i<=YDIV-1; ++i) {
			ps->AddLine(sa[i][j]*(W/2/max),(H/YDIV)*i, sa[i+1][j]*(W/2/max),(H/YDIV)*(i+1),col);
		}
	}
	for(i=0; i<=YDIV-1; ++i) {
		ps->AddLine(sc[i]*(W/2/max),(H/YDIV)*i, sc[i+1]*(W/2/max),(H/YDIV)*(i+1),0,DASH);
	}
	ps->AddText(str(max)+( Afocal ? " Dptr":"" ), W/2*0.95,-1, 3);
	ps->AddText("SA", -3,H+8, 4);
	ps->AddText("Y*max="+str(Round(yPupilMax,-3)),4,H+13, 3);
	ps->AddText(" Fno=" + str(Round(FNumber(0,0,FindPupil),-3)) 
				+ "(PARAX" + str(Round(FNumberParaxial(FindPupil),-3)) + ")",
		        4,H+10, 3);
	ps->AddText(" NA=" + str(Round(NA(0,0,FindPupil),-4)) 
				+ "(PARAX" + str(Round(NAParaxial(FindPupil),-4)) + ")",
		        4,H+7, 3);
	for(j=1; j<=cn; ++j) ps->AddText(str(j), sa[YDIV][j]*(W/2/max)-1,H+4, 3);

	for(i=0; i<=YDIV; ++i) delete[] sa[i]; 
	delete[] sa;
	delete[] sc;

	if(PupilSA){
		pop();
	}
}
	
void cLens1::MakeSAGraph(int FindPupil,double yPupilMax,double FullScale) {
	int i,j;
	point p0;
	double **sa, *sc;
	double max=0;
	int xdiv;
	const int YDIV=20;
	const double W=80, H=60;
	double xPupil;
	cShapes* ps=CurrentShapes(SA_GRAPH);

	ps->Reset();
	sa=new double*[YDIV+1]; for(i=0; i<=YDIV; ++i) sa[i]=new double[cn+1]; 
	sc=new double[YDIV+1];

	if(yPupilMax==0) yPupilMax=this->EPD/2;
	
	if(FindPupil || yPupilMax==0) {
		point ymax,ymin,xmax,xmin;
		if(FindPupilEdge(ymax,ymin,xmax,xmin,0,0,1,1)){
			yPupilMax=ymax.y; xPupil=ymax.x;
			p0=(ymax+ymin)/2;
		}
		else return;
	}
	else {
		xPupil=0;
		p0=point(0,0);
	}
	for(j=1; j<=cn; ++j) for(i=0; i<=YDIV; ++i) {
		sa[i][j]=LSA(p0.y+(yPupilMax-p0.y)*i/YDIV,xPupil, p0.y,p0.x, j);
		max= fabs(sa[i][j])>max ? fabs(sa[i][j]) : max;
	}
	for(i=0; i<=YDIV; ++i) {
		sc[i]=SC(yPupilMax*i/YDIV,xPupil);
		max= fabs(sc[i])>max ? fabs(sc[i]) : max;
	}

	if(FullScale==0){
		max=fceil(max);
		xdiv=(int)fceil_mantissa(max);
	}
	else{
		max=FullScale;  // �ۂ߂Ȃ�
		if(FullScale==fceil(FullScale)){
			xdiv=(int)fceil_mantissa(FullScale);
		}
		else{
			xdiv=4;
		}
	}

	ps->SetWindow(-50,-20,50,80);
	ps->AddLine(-W/2,0, W/2,0);
	for(i=-xdiv; i<=xdiv; ++i) {
		ps->AddLine((W/2/xdiv)*i,0.7, (W/2/xdiv)*i,-0.7);
	}
	ps->AddLine(0,H, 0,-2);
	for(i=1; i<=10; ++i) {
		ps->AddLine(-0.7,(H/10)*i, 0.7,(H/10)*i);
	}
	for(j=cn; j>=1; --j) {
		// for(j=1; j<=cn; ++j) �ɂ���ƁC�O���t���d�Ȃ����Ƃ�
		// ��1�g���̐F����cn�g���̐F�ɉB����Ă��܂��D
		long col=cSpectrum::ApproximateColor(wl[j]);
		for(i=0; i<=YDIV-1; ++i) {
			ps->AddLine(sa[i][j]*(W/2/max),(H/YDIV)*i, sa[i+1][j]*(W/2/max),(H/YDIV)*(i+1),col);
		}
	}
	for(i=0; i<=YDIV-1; ++i) {
		ps->AddLine(sc[i]*(W/2/max),(H/YDIV)*i, sc[i+1]*(W/2/max),(H/YDIV)*(i+1),0,DASH);
	}
	ps->AddText(str(max)+( Afocal ? " Dptr":"" ), W/2*0.95,-1, 3);
	ps->AddText("SA", -3,H+8, 4);
	ps->AddText("Y*max="+str(Round(yPupilMax-p0.y,-3)),4,H+13, 3);
	ps->AddText(" Fno=" + str(Round(FNumber(0,0,FindPupil),-3)) 
				+ "(PARAX" + str(Round(FNumberParaxial(FindPupil),-3)) + ")",
		        4,H+10, 3);
	ps->AddText(" NA=" + str(Round(NA(0,0,FindPupil),-4)) 
				+ "(PARAX" + str(Round(NAParaxial(FindPupil),-4)) + ")",
		        4,H+7, 3);
	for(j=1; j<=cn; ++j) ps->AddText(str(j), sa[YDIV][j]*(W/2/max)-1,H+4, 3);

	for(i=0; i<=YDIV; ++i) delete[] sa[i]; 
	delete[] sa;
	delete[] sc;
}

int cLens1::ASGraph() { return AS_GRAPH; }

void cLens1::MakeASGraph(double yObjMax,int FindPupil,double FullScale,double defocus,int colors/*=1*/,int AbeTheory/*=0*/,
						 int Scan/*=0*/,double yScanThMax/*=0*/,double yScanThMin/*=0*/){
	int i,j,N,ii;
	double ***dm, ***ds;
	double yobjmax,max=0;
	int xdiv;
	long col;
	const int YDIV=40;
	const double W=80, H=60;
	cShapes* ps=CurrentShapes(AS_GRAPH);
	double yScanTh0;

	N= IsXAxisSymmetric() ? 1 : 2;     // y�����ɔ�Ώ̂ł���΁C-yObjMax �����̃O���t���\��

	ps->Reset();
	ps->SetWindow(-50, N==1 ? -20:-95, 50, 80);
	if(colors<=0) colors=1;
	Dim3(dm,N,colors,YDIV);
	Dim3(ds,N,colors,YDIV);
	if(yObjMax==0) yObjMax= yObjectMax;
	if(yObjMax==0 && Scan==0) return;
	if(Scan) yScanTh0=this->yScan();  // ��Ō��ɖ߂����߂ɏ����l��ۑ�����D

	for(j=1; j<=N; ++j){
		if(Scan==0){
			yobjmax= j==1 ? yObjMax : -yObjMax;
		}
		else{
			yobjmax=0;   // ��������ꍇ�C����yObjMax�ɂ�����炸���̍���0
		}
		
		for(i=0; i<=YDIV; ++i){
			double yobj,yp,xp;
			int found;

			yobj=yobjmax*i/YDIV;
			
			if(Scan){
				if(j==1) yScan(yScanThMax*i/YDIV);
				if(j==2) yScan( yScanThMin==0 ? -yScanThMax*i/YDIV : yScanThMin*i/YDIV );
			}

			if(FindPupil){
				point ymax,ymin,xmax,xmin;
				if( found=FindPupilEdge(ymax,ymin,xmax,xmin,yobj,0,1,1) ){ // ���ݒ�͑�1�g���ɂ��
					yp=((ymax+ymin)/2).y; xp=((ymax+ymin)/2).x;
				}
			}
			else{
				found=1; yp=xp=0;
			}

			for(ii=1; ii<=colors; ++ii){
				dm[j][ii][i]= found ? DeltaM(yobj,0, yp,xp, defocus, ii,AbeTheory) : 0;
				max= fabs(dm[j][ii][i])>max ? fabs(dm[j][ii][i]) : max;
				ds[j][ii][i]= found ? DeltaS(yobj,0, yp,xp, defocus, ii,AbeTheory) : 0;
				max= fabs(ds[j][ii][i])>max ? fabs(ds[j][ii][i]) : max;
			}
		}
	}

	if(FullScale==0){
		max=fceil(max);
		xdiv=(int)fceil_mantissa(max);
	}
	else{
		max=FullScale;  // �ۂ߂Ȃ�
		if(FullScale==fceil(FullScale)){
			xdiv=(int)fceil_mantissa(FullScale);
		}
		else{
			xdiv=4;
		}
	}

	for(j=1; j<=N; ++j){
		if(j==1){
			ps->OffsetY=0;
			yobjmax=yObjMax;
		}
		else if(j==2){
			ps->OffsetY=-75;
			yobjmax=-yObjMax;
		}
		ps->AddLine(-W/2,0, W/2,0);
		for(i=-xdiv; i<=xdiv; ++i) {
			ps->AddLine((W/2/xdiv)*i,0.7, (W/2/xdiv)*i,-0.7);
		}
		ps->AddLine(0,H, 0,-2);
		for(i=1; i<=10; ++i) {
			ps->AddLine(-0.7,(H/10)*i, 0.7,(H/10)*i);
		}
		ps->AddText( str(max)+(Afocal ? " Dptr":""), W/2*0.95,-1, 3);
		ps->AddText("AS", -3,H+8, 4);
		if(Scan==0){
			ps->AddText("Ymax="+str(yobjmax)+"("+str(Round(AngleOfView(yobjmax),-1))+"�)",4,H+7,3);
		}
		else{
			if(j==1) ps->AddText("YScanMax="+str(yScanThMax),4,H+7,3);
			if(j==2) ps->AddText("YScanMin="+str(yScanThMin==0 ? -yScanThMax:yScanThMin),4,H+7,3);
		}
		ps->AddText("DM", dm[j][1][YDIV]*(W/2/max)-3,H+4, 3);
		ps->AddText("DS", ds[j][1][YDIV]*(W/2/max)-3,H+4, 3);

		for(ii=1; ii<=colors; ++ii){
			for(i=0; i<=YDIV-1; ++i){
				col= colors>1 ? cSpectrum::ApproximateColor(wl[ii]) : rgb(230,0,100);
				ps->AddLine(dm[j][ii][i]*(W/2/max),(H/YDIV)*i,dm[j][ii][i+1]*(W/2/max),(H/YDIV)*(i+1),col,DASH);
				ps->AddLine(ds[j][ii][i]*(W/2/max),(H/YDIV)*i,ds[j][ii][i+1]*(W/2/max),(H/YDIV)*(i+1),col);
			}
		}
	}

	Erase3(dm,N,colors);
	Erase3(ds,N,colors);
	if(Scan) yScan(yScanTh0);
}

int cLens1::DistGraph() { return DIST_GRAPH; }

void cLens1::MakeDistGraph(double yObjMax,int FindPupil,double FullScale) {
	int i;
	double *dist;
	double max=0;
	int xdiv;
	const int YDIV=20;
	const double W=80, H=60;
	cShapes* ps=CurrentShapes(DIST_GRAPH);

	ps->Reset();
	ps->SetWindow(-50,-20,50,80);
	if(yObjMax==0) yObjMax=yObjectMax;
	if(yObjMax==0) return;

	dist=new double[YDIV+1];
	
	for(i=0; i<=YDIV; ++i) {
		double yObj,yPupil,xPupil;
		int found;
		yObj=yObjMax*i/YDIV;
		if(FindPupil){
			point ymax,ymin,xmax,xmin;
			if( found=FindPupilEdge(ymax,ymin,xmax,xmin,yObj,0,1,1) ){
				yPupil=((ymax+ymin)/2).y; xPupil=((ymax+ymin)/2).x;
			}
		}
		else{
			found=1; yPupil=xPupil=0;
		}
		dist[i]= found ? Dist(yObj,yPupil,xPupil,1)*100 : 0;
		max= fabs(dist[i])>max ? fabs(dist[i]) : max;
	}

	if(FullScale==0){
		max=fceil(max);
		xdiv=(int)fceil_mantissa(max);
	}
	else{
		max=FullScale;  // �ۂ߂Ȃ�
		if(FullScale==fceil(FullScale)){
			xdiv=(int)fceil_mantissa(FullScale);
		}
		else{
			xdiv=4;
		}
	}

	ps->AddLine(-W/2,0, W/2,0);
	for(i=-xdiv; i<=xdiv; ++i) {
		ps->AddLine((W/2/xdiv)*i,0.7, (W/2/xdiv)*i,-0.7);
	}
	ps->AddLine(0,H, 0,-2);
	for(i=1; i<=10; ++i) {
		ps->AddLine(-0.7,(H/10)*i, 0.7,(H/10)*i);
	}
	for(i=0; i<=YDIV-1; ++i) {
		ps->AddLine(dist[i]*(W/2/max),(H/YDIV)*i,dist[i+1]*(W/2/max),(H/YDIV)*(i+1),rgb(0,0,220));
	}
	ps->AddText(str(max)+"%", W/2*0.95,-1, 3);
	ps->AddText("DIST", -5,H+8, 4);
	ps->AddText("Ymax="+str(yObjMax)+"("+str(Round(AngleOfView(yObjMax),-1))+"�)",4,H+7,3);

	delete[] dist;
}

int cLens1::DistChart() { return DIST_CHART; }

std::string cLens1::MakeDistChart(double yObjMax,double xObjMax,int DIVy,int DIVx,int FindPupil,double defocus,
								  int Scan,
								  std::string yScanAxisXY,double yScanThMin,double yScanThMax,double yScanRatio,
							      std::string xScanAxisXY,double xScanThMin,double xScanThMax,double xScanRatio)
{
	int i,j,k;
	std::string setray;
	point **p,**p1,**p2, center;
	double max;
	const double L=100;
	const point NullPoint=point(999,999);
	list<double> p1_x,p1_y;
	point top_left,top_center,top_right;
	point center_left,center_right;
	point bottom_left,bottom_center,bottom_right;
	double dist_top,dist_left,dist_right,dist_bottom;
	double my_ratio_top_center;
	double mx_ratio_center_left,mx_ratio_center_right;
	double my_ratio_bottom_center;
	cShapes* ps=CurrentShapes(DIST_CHART);
	std::string s;
	char buf[1000];
	double dummy;
	double thx,thy;

	ps->Reset();
	ps->SetWindow(-L/2*1.2,-L/2*1.2,L/2*1.2,L/2*1.2);
	
	if(yObjMax==0 && xObjMax==0) { yObjMax=yObjectMax; xObjMax=xObjectMax; }
	if(yObjMax==0 && xObjMax==0 && Scan==0) return s;
	if(yObjMax==0) { yObjMax=xObjMax; }
	if(xObjMax==0) { xObjMax=yObjMax; }

	if(DIVx==0) DIVx=10;
	if(DIVy==0) DIVy=10;
	if(DIVx%2!=0) DIVx+=1;
	if(DIVy%2!=0) DIVy+=1;

	p =new point* [DIVx+2]; for(i=0; i<=DIVx+1; ++i) p[i] =new point[DIVy+2];
	p1=new point* [DIVx+2]; for(i=0; i<=DIVx+1; ++i) p1[i]=new point[DIVy+2];
	p2=new point* [DIVx+2]; for(i=0; i<=DIVx+1; ++i) p2[i]=new point[DIVy+2];
	
	setray= FindPupil==0 ? "" : "principal";

	push();   // �X�L�����n�ŕΐS�ɕύX��������̂ŕۑ�

	for(i=1; i<=DIVx+1; ++i) for(j=1; j<=DIVy+1; ++j) {
		if(Scan==0){
			p[i][j]=point(-xObjMax+(xObjMax)*2*(i-1)/DIVx, -yObjMax+(yObjMax)*2*(j-1)/DIVy);
			if(ImageHeight(p1[i][j].y,p1[i][j].x,dummy, p[i][j].y,p[i][j].x,setray,0,0,defocus,1,0,0,0,0)==0){
				p1[i][j]=NullPoint;
			}
		}
		else{
			// �X�L�����n�̏ꍇ
			thy=yScanThMin+(yScanThMax-yScanThMin)*(j-1)/DIVy;
			thx=xScanThMin+(xScanThMax-xScanThMin)*(i-1)/DIVx;
			p[i][j]=point(thx,thy);

			if(0<yScanRatio && yScanRatio<1){
				// ���U�~���[�ɂ��c�ȁi���Ԃɑ΂��Ċp�x�����`�łȂ����Ƃɂ��j
				//   xScanRatio = ���U�~���[�̐U���ɑ΂���g�p�p�x�͈͂̊����i���ʂ�0.8, 0.9�Ȃ�)
				//
				// yScanThMax=-yScanThMin�̂Ƃ��i�I�t�Z�b�g���Ȃ��Ƃ��j�C
				//      thy = (yScanThMax/yScanRatio)sin(x)
				//           x= -asin(yScanRaito) �` +asin(yScanRatio), x�͎��Ԃ�\����������
				// �ƂȂ�D
				double thmax,thave,th0max,th0;

				thmax=(yScanThMax-yScanThMin)/2;  // deg
				thave=(yScanThMax+yScanThMin)/2;  // deg
				th0max=asin(yScanRatio);          // rad
				th0=-th0max+2*th0max*(j-1)/DIVy;  // rad
				thy=(thmax/yScanRatio)*sin(th0);
			}

			if(0<xScanRatio && xScanRatio<1){
				double thmax,thave,th0max,th0;

				thmax=(xScanThMax-xScanThMin)/2;
				thave=(xScanThMax+xScanThMin)/2;
				th0max=asin(xScanRatio);
				th0=-th0max+2*th0max*(i-1)/DIVx;
				thx=(thmax/xScanRatio)*sin(th0);
			}

			yScan(yScanAxisXY,thy);
			xScan(xScanAxisXY,thx);
			if(ImageHeight(p1[i][j].y,p1[i][j].x,dummy, 0,0,setray,0,0,defocus,1,0,0,0,0)==0){
				p1[i][j]=NullPoint;
			}			
		}
	}

	pop();
	
	for(i=1; i<=DIVx+1; ++i) for(j=1; j<=DIVy+1; ++j){
		if(p1[i][j]!=NullPoint){
			p1_x.AddTail(p1[i][j].x);
			p1_y.AddTail(p1[i][j].y);
		}
	}
	center=point( (p1_x.Max()+p1_x.Min())/2, (p1_y.Max()+p1_y.Min())/2 );
	max=Max( p1_x.Max()-p1_x.Min(), p1_y.Max()-p1_y.Min() ); 
	for(i=1; i<=DIVx+1; ++i) for(j=1; j<=DIVy+1; ++j){
		if(p1[i][j]!=NullPoint){
			p2[i][j]=(p1[i][j]-center)*L/max;
		}
		else{
			p2[i][j]=NullPoint;
		}
	}
	
	for(j=1; j<=DIVy+1; ++j) for(i=1; i<=DIVx; ++i) {
		if( p2[i][j]!=NullPoint && p2[i+1][j]!=NullPoint ) {
			ps->AddLine(p2[i][j].x,p2[i][j].y, p2[i+1][j].x,p2[i+1][j].y);
		}
	}
	for(i=1; i<=DIVx+1; ++i) for(j=1; j<=DIVy; ++j) {
		if( p2[i][j]!=NullPoint && p2[i][j+1]!=NullPoint ) {
			ps->AddLine(p2[i][j].x,p2[i][j].y, p2[i][j+1].x,p2[i][j+1].y);
		}
	}
	
	if(p2[DIVx/2][DIVy/2].x<=p2[DIVx/2+1][DIVy/2].x && p2[DIVx/2][DIVy/2].y<=p2[DIVx/2][DIVy/2+1].y){
		top_left     =p2[1][DIVy+1]; 
		top_center   =p2[DIVx/2+1][DIVy+1]; 
		top_right    =p2[DIVx+1][DIVy+1];
		center_left  =p2[1][DIVy/2+1]; 
		center_right =p2[DIVx+1][DIVy/2+1];
		bottom_left  =p2[1][1]; 
		bottom_center=p2[DIVx/2+1][1]; 
		bottom_right =p2[DIVx+1][1];
	}
	else if(p2[DIVx/2][DIVy/2].x<=p2[DIVx/2+1][DIVy/2].x && p2[DIVx/2][DIVy/2].y>=p2[DIVx/2][DIVy/2+1].y){
		bottom_left  =p2[1][DIVy+1]; 
		bottom_center=p2[DIVx/2+1][DIVy+1]; 
		bottom_right =p2[DIVx+1][DIVy+1];
		center_left  =p2[1][DIVy/2+1]; 
		center_right =p2[DIVx+1][DIVy/2+1];
		top_left     =p2[1][1]; 
		top_center   =p2[DIVx/2+1][1]; 
		top_right    =p2[DIVx+1][1];
	}
	else if(p2[DIVx/2][DIVy/2].x>=p2[DIVx/2+1][DIVy/2].x && p2[DIVx/2][DIVy/2].y<=p2[DIVx/2][DIVy/2+1].y){
		top_right    =p2[1][DIVy+1];
		top_center   =p2[DIVx/2+1][DIVy+1];
		top_left     =p2[DIVx+1][DIVy+1];
		center_right =p2[1][DIVy/2+1];
		center_left  =p2[DIVx+1][DIVy/2+1];
		bottom_right =p2[1][1];
		bottom_center=p2[DIVx/2+1][1];
		bottom_left  =p2[DIVx+1][1];
	}
	else{
		bottom_right =p2[1][DIVy+1];
		bottom_center=p2[DIVx/2+1][DIVy+1];
		bottom_left  =p2[DIVx+1][DIVy+1];
		center_right =p2[1][DIVy/2+1];
		center_left  =p2[DIVx+1][DIVy/2+1];
		top_right    =p2[1][1];
		top_center   =p2[DIVx/2+1][1];
		top_left     =p2[DIVx+1][1];
	}
	if(top_left!=NullPoint && top_center!=NullPoint && top_right!=NullPoint ){
		dist_top=distance(top_center,line(top_left,top_right))/distance(top_left,top_right);
	}
	else dist_top=0;
	if(top_left!=NullPoint && center_left!=NullPoint && bottom_left!=NullPoint){
		dist_left=distance(center_left,line(top_left,bottom_left))/distance(top_left,bottom_left);
	}
	else dist_left=0;
	if(top_right!=NullPoint && center_right!=NullPoint && bottom_right!=NullPoint){
		dist_right=distance(center_right,line(top_right,bottom_right))/distance(top_right,bottom_right);
	}
	else dist_right=0;
	if(bottom_left!=NullPoint && bottom_center!=NullPoint && bottom_right!=NullPoint){
		dist_bottom=distance(bottom_center,line(bottom_left,bottom_right))/distance(bottom_left,bottom_right);
	}
	else dist_bottom=0;
	
	{
		// �Ǐ��{�� dy1/dy, dx1/dx �̎���[�̎���ɑ΂���䗦���v�Z����D

		// �ŏ��͍ŏ����@��p1[][]�ɂ��x1,y1��p[][]�ɂ��x,y�̑������W�J�ŕ\���C���̔������
		// dy1/dy, dx1/dx �����߂悤�Ƃ������C�Ȃ������x���o�Ȃ��̂Œf�O���� (090904)�D

		const double K=0.0001;
		point top_center;
		point center_left,center_right;
		point bottom_center;
		double dy,dx, dy1o,dy1,dx1o,dx1;

		if(p1[DIVx/2][DIVy/2].x<=p1[DIVx/2+1][DIVy/2].x && p1[DIVx/2][DIVy/2].y<=p1[DIVx/2][DIVy/2+1].y){
			top_center   =p[DIVx/2+1][DIVy+1]; 
			center_left  =p[1][DIVy/2+1]; 
			center_right =p[DIVx+1][DIVy/2+1];
			bottom_center=p[DIVx/2+1][1]; 
		}
		else if(p1[DIVx/2][DIVy/2].x<=p1[DIVx/2+1][DIVy/2].x && p1[DIVx/2][DIVy/2].y>=p1[DIVx/2][DIVy/2+1].y){
			bottom_center=p[DIVx/2+1][DIVy+1]; 
			center_left  =p[1][DIVy/2+1]; 
			center_right =p[DIVx+1][DIVy/2+1];
			top_center   =p[DIVx/2+1][1]; 
		}
		else if(p1[DIVx/2][DIVy/2].x>=p1[DIVx/2+1][DIVy/2].x && p1[DIVx/2][DIVy/2].y<=p1[DIVx/2][DIVy/2+1].y){
			top_center   =p[DIVx/2+1][DIVy+1];
			center_right =p[1][DIVy/2+1];
			center_left  =p[DIVx+1][DIVy/2+1];
			bottom_center=p[DIVx/2+1][1];
		}
		else{
			bottom_center=p[DIVx/2+1][DIVy+1];
			center_right =p[1][DIVy/2+1];
			center_left  =p[DIVx+1][DIVy/2+1];
			top_center   =p[DIVx/2+1][1];
		}

		dy=distance(bottom_center,top_center)*K;
		dx=distance(center_right,center_left)*K;

		dy1o=DImageHeight(0,0,setray,0,0,defocus,1,0,0,0,0,dy,0);
		dy1 =DImageHeight(top_center.y,top_center.x,setray,0,0,defocus,1,0,0,0,0,dy,0);
		my_ratio_top_center= dy1/dy1o;

		dx1o=DImageHeight(0,0,setray,0,0,defocus,1,0,0,0,0,0,dx);
		dx1 =DImageHeight(center_left.y,center_left.x,setray,0,0,defocus,1,0,0,0,0,0,dx);
		mx_ratio_center_left= dx1/dx1o;

		dx1o=DImageHeight(0,0,setray,0,0,defocus,1,0,0,0,0,0,dx);
		dx1 =DImageHeight(center_right.y,center_right.x,setray,0,0,defocus,1,0,0,0,0,0,dx);
		mx_ratio_center_right= dx1/dx1o;

		dy1o=DImageHeight(0,0,setray,0,0,defocus,1,0,0,0,0,dy,0);
		dy1 =DImageHeight(bottom_center.y,bottom_center.x,setray,0,0,defocus,1,0,0,0,0,dy,0);
		my_ratio_bottom_center= dy1/dy1o;
	}

	{
		const double H=L*0.07;
		const double P=H*1.4;

		ps->AddText("DistTop   ="+str(Round(dist_top*100   ,-3))+"%",L*0.53,L/2    ,H);
		ps->AddText("DistLeft  ="+str(Round(dist_left*100  ,-3))+"%",L*0.53,L/2-P  ,H);
		ps->AddText("DistRight ="+str(Round(dist_right*100 ,-3))+"%",L*0.53,L/2-P*2,H);
		ps->AddText("DistBottom="+str(Round(dist_bottom*100,-3))+"%",L*0.53,L/2-P*3,H);
		
		ps->AddText("my/my0 Top   ="+str(Round((my_ratio_top_center   -1)*100,-3))+"%",L*0.53,L/2-P*5,H);
		ps->AddText("mx/mx0 Left  ="+str(Round((mx_ratio_center_left  -1)*100,-3))+"%",L*0.53,L/2-P*6,H);
		ps->AddText("mx/mx0 Right ="+str(Round((mx_ratio_center_right -1)*100,-3))+"%",L*0.53,L/2-P*7,H);
		ps->AddText("my/my0 Bottom="+str(Round((my_ratio_bottom_center-1)*100,-3))+"%",L*0.53,L/2-P*8,H);
	}

	for(k=1; k<=2; k++){
		// �c�Ȃ�x,y�ׂ̂������œW�J����
		//   1���(k=1) : x,y�œW�J����
		//   2���(k=2) : x/xObjMax, y/yObjMax�œW�J����
		const int ORDER=4;
		int i,j, m;
		double dx,dy,dr,max;
		cFitting fx,fy;
		list<point> lp,lp1;

		for(i=1; i<=DIVx+1; ++i) for(j=1; j<=DIVy+1; ++j) {
			if( p1[i][j]!=NullPoint ) {
				switch(k){
				case 1: lp.AddTail(p[i][j]);                                    break;
				case 2: lp.AddTail(point(p[i][j].x/xObjMax,p[i][j].y/yObjMax)); break;
				}
				lp1.AddTail(p1[i][j]);
			}
		}

		fx.SetOrder(ORDER);
		fy.SetOrder(ORDER);
		
		m=lp.GetSize();
		fx.SetNumberOfData(m);
		fy.SetNumberOfData(m);

		for(i=1; i<=m; i++){
			fx.dataSet(i,lp[i].x,lp[i].y,lp1[i].x);
			fy.dataSet(i,lp[i].x,lp[i].y,lp1[i].y);			
		}

		fx.CalcCoefficients();
		switch(k){
		case 1: s+="x1 expanded by x,y power series\n";                 break;
		case 2: s+="x1 expanded by x/xObjMax,y/yObjMax power series\n"; break;
		}
		for(j=1; j<=fx.GetNumberOfTerms(); j++){
			sprintf(buf,"  %d %d %g\n", 
			        fx.x1dimensionGet(j),fx.x2dimensionGet(j),fx.coefficientGet(j,0));
			s+=buf;
		}
		sprintf(buf,"  magnification or f=%g\n",fx.coefficientGet(1,0,0)); s+=buf;
		sprintf(buf,"  pvError=%g\n\n",fx.pvError(0)); s+=buf;

		fy.CalcCoefficients();
		switch(k){
		case 1: s+="y1 expanded by x,y power series\n";                 break;
		case 2: s+="y1 expanded by x/xObjMax,y/yObjMax power series\n"; break;
		}
		for(j=1; j<=fy.GetNumberOfTerms(); j++){
			sprintf(buf,"  %d %d %g\n",
			        fy.x1dimensionGet(j),fy.x2dimensionGet(j),fy.coefficientGet(j,0));
			s+=buf;
		}
		sprintf(buf,"  magnification or f=%g\n",fy.coefficientGet(0,1,0)); s+=buf;
		sprintf(buf,"  pvError=%g\n\n",fy.pvError(0)); s+=buf;

		max=0;
		for(i=1; i<=m; i++){
			dx=fx.dyApproximate(fx.x1dataGet(i),fx.x2dataGet(i),0);
			dy=fy.dyApproximate(fy.x1dataGet(i),fy.x2dataGet(i),0);
			dr=sqrt(dx*dx+dy*dy);
			if(dr>max) max=dr;
		}
		sprintf(buf,"max of deviation from linear part=%g\n\n",max); s+=buf;
	}

	{
		// 1���̌W��(�ߎ��{��)��x,y�����ŋ��ʂƂ����W�J���s�Ȃ��D5���܂œW�J����Ƃ��C
		//   x1=beta*x +a00 +a20*x*x +a11*x*y +a02*y*y +a30*x*x*x +... +a05*y*y*y*y*y
		//   y1=beta*y +b00 +b20*x*x +b11*x*y +b02*y*y +b30*x*x*x +... +b05*y*y*y*y*y
		// �ŏ�2��@��AX=B������
		//     | x 1 x*x x*y y*y x*x*x .....  0 0      0         0         |
		//     | x 1 x*x x*y y*y x*x*x .....  0 0      0         0         |
		//   A=|     ....                      ....                        |
		//     | y 0 0   0   0   0     .....  1 x*x .. x*y*y*y*y y*y*y*y*y |
		//     | y 0 0   0   0   0     .....  1 x*x .. x*y*y*y*y y*y*y*y*y |
		//     |     ....                      ....                        |
		//
		//     | x1 |
		//     | x1 |
        //   B=| .. |
		//     | y1 |
		//     | y1 |
		//     | .. |
		//
		// ����ꂽX�̐����͓W�J�W����\��
		//     | beta |
		//     | a00  |
		//     | a20  |
		//     | a11  |
		//   X=| a02  |
		//     | a30  |
		//     | ..   |
		//     | b00  |
		//     | b20  |
		//     | ..   |
		//     | b14  |
		//     | b05  |

		const int ORDER=5;
		int m_half,n, i,j,jj,k;
		int a,b, *x_order,*y_order;
		double x,y,x1,y1;
		cLeneq ln;
		list<point> lp,lp1;

		for(i=1; i<=DIVx+1; ++i) for(j=1; j<=DIVy+1; ++j) {
			if( p1[i][j]!=NullPoint ) {
				lp.AddTail(p[i][j]);
				lp1.AddTail(p1[i][j]);
			}
		}

		m_half=lp.GetSize();
		n=(ORDER+1)*(ORDER+2)-3;

		x_order=new int[n+1];
		y_order=new int[n+1];

		ln.SetNumberOfEq(m_half*2);
		ln.SetNumberOfVar(n);

		for(i=1; i<=m_half; i++){
			x1=lp1[i].x; x=lp[i].x; y=lp[i].y;

			ln.SetA(i,1,x);
			ln.SetA(i,2,1); x_order[2]=y_order[2]=0;

			j=3;
			for(jj=2; jj<=ORDER; jj++){
				a=jj;
				b=0;
				for(k=1; k<=jj+1; k++){
					ln.SetA(i,j,pow(x,a)*pow(y,b));
					x_order[j]=a;
					y_order[j]=b;
					a--;
					b++;
					j++;
				}
			}

			ln.SetB(i,x1);			
		}

		for(i=1; i<=m_half; i++){
			y1=lp1[i].y; x=lp[i].x; y=lp[i].y;
			
			ln.SetA(m_half+i,1,y);
			
			j=(ORDER+1)*(ORDER+2)/2;

			ln.SetA(m_half+i,j,1); x_order[j]=y_order[j]=0; j++;

			for(jj=2; jj<=ORDER; jj++){
				a=jj;
				b=0;
				for(k=1; k<=jj+1; k++){
					ln.SetA(m_half+i,j,pow(x,a)*pow(y,b));
					x_order[j]=a;
					y_order[j]=b;
					a--;
					b++;
					j++;
				}
			}

			ln.SetB(m_half+i,y1);
		}

		s+="expand x1,y1 with same beta in x,y direction\n";
		s+=" x1=beta*x+a00+a20xx+a11xy+a02yy+a30xxx+...\n";
		s+=" y1=beta*y+b00+b20xx+b11xy+b02yy+b30xxx+...\n";
		sprintf(buf,"  beta= %g\n",ln.GetX(1)); s+=buf;
		for(j=2; j<=(ORDER+1)*(ORDER+2)/2-1; j++){
			sprintf(buf,"  a %d %d = %g\n",x_order[j],y_order[j],ln.GetX(j)); s+=buf;
		}
		for(j=(ORDER+1)*(ORDER+2)/2; j<=n; j++){
			sprintf(buf,"  b %d %d = %g\n",x_order[j],y_order[j],ln.GetX(j)); s+=buf;
		}
		sprintf(buf,"  pvError=%g\n\n",ln.pvError()); s+=buf;

		delete[] x_order;
		delete[] y_order;
	}

	{
		// �������W���o�͂���
		double x,y,x1,y1;

		//          " ########## ########## ########## ########## ############# #############"
		sprintf(buf,"        x          y    x/xObjMax  y/yObjMax          x1            y1 \n"); 
		s+=buf;
		
		for(i=1; i<=DIVx+1; ++i) for(j=1; j<=DIVy+1; ++j) {	
			x=p[i][j].x;
			y=p[i][j].y;
			x1=p1[i][j].x;
			y1=p1[i][j].y;
			
			if( p1[i][j]!=NullPoint ) {
				sprintf(buf," %10g %10g %10g %10g %13.8g %13.8g\n",x,y,x/xObjMax,y/yObjMax,x1,y1);
				s+=buf;
			}
			else{
				sprintf(buf," %10g %10g %10g %10g   no ray\n",x,y,x/xObjMax,y/yObjMax);
				s+=buf;
			}
		}
	}

	for(i=0; i<=DIVx+1; i++) delete[] p[i];  delete[] p;
	for(i=0; i<=DIVx+1; i++) delete[] p1[i]; delete[] p1;
	for(i=0; i<=DIVx+1; i++) delete[] p2[i]; delete[] p2;

	return s;
}

std::string cLens1::MakeDistChart(double yObjMax,double xObjMax,int DIVy,int DIVx,int FindPupil,double defocus){
	return MakeDistChart(yObjMax,xObjMax,DIVy,DIVx,FindPupil,defocus,0,"",0,0,0,"",0,0,0);
}

std::string cLens1::MakeDistChart(double yObjMax,double xObjMax,int FindPupil,double defocus) {
	return MakeDistChart(yObjMax,xObjMax,10,10,FindPupil,defocus);
}

std::string cLens1::MakeScanDistChart(std::string yScanAxisXY,double yScanThMin,double yScanThMax,double yScanRatio,
									  std::string xScanAxisXY,double xScanThMin,double xScanThMax,double xScanRatio,
									  int DIVy,int DIVx,int FindPupil,double defocus){
	return MakeDistChart(0,0,DIVy,DIVx,FindPupil,defocus,
	                     1,yScanAxisXY,yScanThMin,yScanThMax,yScanRatio,xScanAxisXY,xScanThMin,xScanThMax,xScanRatio);
}

std::string cLens1::ExpandDist(double yObjMax,double xObjMax,int FindPupil,double defocus){
	return MakeDistChart(yObjMax,xObjMax,FindPupil,defocus);
}

int cLens1::DeltaHGraph() { return DELTAH_GRAPH; }

void cLens1::MakeDeltaHGraph(double yObjMax,double xObjMax,int FindPupil,double PupilMax,
							 double defocus,int BothSide,double FullScale,std::string command,
							 int i1,int i2,int AbeTheory/*=0*/) 
{
	// ��i�����̌v�Z������O�Ɉ���command�̑�isentence�����s����D
	// ��i�����̌v�Z��ɂ͂��̕ύX�����ɖ߂��D
	// ��F commnad = "rox 9 3.4; rox 9 2.55; rox 9 1.7; rox 9 0"
	//      �Ƃ���΁C�X�L���i���w�n�̉������}������D

	// AbeTheory ���^�̂Ƃ��́C�����W�����狁�߂��ߎ����ʂɂ�����O���t�����D
	// AbeTheory==5 �̂Ƃ���5�������W���������D����ȊO�̂Ƃ���3�������̂݁D

	int i,j,ii,iii,iiii;
	const int M= BothSide ? 7 : 4;  // ���̐�
	const int N=8;                  // �����v�Z�_�� =Nx2+1
	const double H=12, W=26;
	double ****X,****Y;
	double Xmax=0,Ymax=0; int Xdiv,Ydiv;
	double *yObj,*xObj;
	int *found;
	int DYIs1DXIs2;
	double y1pr,x1pr;
	point pr, max,min, *ymax,*ymin,*xmax,*xmin, p, height;
	int write_hpr, write_axis_title;
	std::string text,com;

	cLens1 buf=*this;    // ����command�ɂ��f�[�^�ύX�̉\�������邽��
	cShapes* ps=CurrentShapes(DELTAH_GRAPH);
	
	yObj=new double[M+1]; xObj=new double[M+1];
	ymax=new point[M+1];
	ymin=new point[M+1];
	xmax=new point[M+1];
	xmin=new point[M+1];
	found=new int [M+1];
	X=new double***[M+1]; Y=new double***[M+1];

	for(iii=0; iii<=M; ++iii){
		X[iii]=new double**[3];
		Y[iii]=new double**[3];
		for(ii=0; ii<=2; ++ii){
			X[iii][ii]=new double*[N*2+2]; Y[iii][ii]=new double*[N*2+2];
			for(i=0; i<=N*2+1; ++i) {
				X[iii][ii][i]=new double[cn+1]; Y[iii][ii][i]=new double[cn+1];
			}
		}
	}
	
	if(yObjMax==0 && xObjMax==0){
		yObjMax=yObjectMax;
		xObjMax=xObjectMax;
	}
	if(PupilMax==0)  PupilMax= EPDx==0 ? EPD/2:Max(EPD,EPDx)/2;

	for(iii=1; iii<=M; ++iii){
		*this=buf; // iii���Ɍ��ɖ߂��D

		// "command" ��sentence�̈�K�̓���q��F�߂�D
		// ��F
		//    command="rox 18 5; roy 21 5;"; "rox 18 2.5; roy 18 2.5"; ....
		//    -> M=1 �� roy 18 5; rox 21 5; �����s���C
		//       M=2 �� roy 18 2.5; rox 18 2.5; �����s����D
		com=sentence(command,iii);
		for(iiii=1; iiii<=sentences(com); ++iiii){
			scmd(sentence(com,iiii),0);
		}

		switch(iii){
		case(1):
			yObj[iii]=yObjMax;       xObj[iii]=xObjMax;
			break;
		case(2):
			yObj[iii]=yObjMax*0.75;  xObj[iii]=xObjMax*0.75;
			break;
		case(3):
			yObj[iii]=yObjMax*0.5;   xObj[iii]=xObjMax*0.5;
			break;
		case(4):
			yObj[iii]=0;             xObj[iii]=0;
			break;
		case(5):
			yObj[iii]=-yObjMax*0.5;  xObj[iii]=-xObjMax*0.5;
			break;
		case(6):
			yObj[iii]=-yObjMax*0.75; xObj[iii]=-xObjMax*0.75;
			break;
		case(7):
			yObj[iii]=-yObjMax;      xObj[iii]=-xObjMax;
			break;
		}

		if( FindPupil==0 ){
			ymax[iii]=point(0,PupilMax); ymin[iii]=point(0,-PupilMax);
			xmax[iii]=point(PupilMax,0); xmin[iii]=point(-PupilMax,0);
			found[iii]=1;
		}
		else{
			found[iii]=FindPupilEdge(ymax[iii],ymin[iii],xmax[iii],xmin[iii],yObj[iii],xObj[iii],1,1);
		}

		for(ii=1; ii<=2; ++ii){
			if( ii==2 && yObj[iii]==0 && xObj[iii]==0 && IsRotationallySymmetric() ){
				break;
			}
			else{				
				switch(ii){
				case(1):
					DYIs1DXIs2=1;
					max=ymax[iii];
					min=ymin[iii];
					break;
				case(2):
					DYIs1DXIs2=2;
					max=xmax[iii];
					min=xmin[iii];
					break;
				}
			}
			for(i=0; i<=N*2; ++i) for(j=1; j<=cn; ++j) {
				p=min+(max-min)*i/(N*2);
				pr=(max+min)/2;
				if( found[iii] ){
					Y[iii][ii][i][j]= 
					DeltaH(DYIs1DXIs2,yObj[iii],xObj[iii], p.y,p.x, pr.y,pr.x, defocus,j,1,i1,i2,AbeTheory);
				}
				else {
					Y[iii][ii][i][j]=0;
				}
				switch(DYIs1DXIs2){
				case 1: X[iii][ii][i][j]=(p-pr).y; break;
				case 2:	X[iii][ii][i][j]=(p-pr).x; break;
				}
				if( fabs(Y[iii][ii][i][j])>Ymax ) Ymax=fabs(Y[iii][ii][i][j]);
				if( fabs(X[iii][ii][i][j])>Xmax ) Xmax=fabs(X[iii][ii][i][j]);
			}
		}
	}
	ps->Reset();
    ps->SetWindow(-35,30,120,-(M-1)*26-30);

	if( Xmax!=0 ){
		for(iii=1; iii<=M; ++iii){
			switch(iii){
			case(1):
				ps->OffsetY=0;
				write_hpr=1;
				write_axis_title=1;
				break;
			case(2):
				ps->OffsetY=-26;
				write_hpr=0;
				write_axis_title=0;
				break;
			case(3):
				ps->OffsetY=-52;
				write_hpr=0;
				write_axis_title=0;
				break;
			case(4):
				ps->OffsetY=-78;
				write_hpr=0;
				write_axis_title=0;
				break;
			case(5):
				ps->OffsetY=-104;
				write_hpr=0;
				write_axis_title=0;
				break;
			case(6):
				ps->OffsetY=-130;
				write_hpr=0;
				write_axis_title=0;
				break;
			case(7):
				ps->OffsetY=-156;
				write_hpr=1;
				write_axis_title=0;
				break;
			}
			for(ii=1; ii<=2; ++ii){
				if( ii==2 && yObj[iii]==0 && xObj[iii]==0 && IsRotationallySymmetric() ){
					break;
				}
				else {
					switch(ii){
					case(1):
						if( yObj[iii]==0 && xObj[iii]==0 && IsRotationallySymmetric() ){
							ps->OffsetX=37.5;
						}
						else{
							ps->OffsetX=0;
						}
						break;
					case(2):
						ps->OffsetX=75;
						break;
					}
				}
				Xdiv=4;

				if(FullScale==0){
					if(Ymax==0) Ymax=0.000001;  // ���z�����Y�̌n�ȂǂŔ�������/0�G���[��h��
					Ymax=fceil(Ymax);
					Ydiv=(int)fceil_mantissa(Ymax);
				}
				else{
					Ymax=FullScale;  // �ۂ߂Ȃ�
					if(FullScale==fceil(FullScale)){
						Ydiv=(int)fceil_mantissa(FullScale);
					}
					else{
						Ydiv=4;
					}
				}
				
				ps->AddLine(-W,0, W,0);
				for(i=-Xdiv; i<=Xdiv; ++i) {
					ps->AddLine((W/Xdiv)*i,0.7, (W/Xdiv)*i,-0.7);
				}
				ps->AddLine(0,-H, 0,H);
				for(i=-Ydiv; i<=Ydiv; ++i) {
					ps->AddLine(-0.7,(H/Ydiv)*i, 0.7,(H/Ydiv)*i);
				}
				if(found[iii]){
					for(j=cn; j>=1; --j) {
						// for(j=1; j<=cn; ++j) �Ƃ���ƁC�O���t���d�Ȃ����Ƃ��C
						// ��1�g���̐F����cn�g���̐F�ɉB����Ă��܂��D
						long col=cSpectrum::ApproximateColor(wl[j]);
						for(i=0; i<=N*2-1; ++i) {
							ps->AddLine
							(X[iii][ii][i][j]*W/Xmax,Y[iii][ii][i][j]*H/Ymax, 
							 X[iii][ii][i+1][j]*W/Xmax,Y[iii][ii][i+1][j]*H/Ymax,col);
						}
					}
				}
				if(ii==1) {
					double h;
					if(write_axis_title){
						ps->AddText("DY", -2,H+5, 3);
						ps->AddText("Y", W+1,3.6, 3);
						ps->AddText(str(Ymax)+AfocalImageUnitStr(), 1,H+2, 3);
						ps->AddText(str(Round(Xmax,-3)), W,-0.4, 3);
					}
					if(write_hpr){
						h= xObjMax==0 ? H : H+4;
					}
					else{
						h=H;
					}
					text="Y="+str(yObj[iii]);
					if(yObj[iii]!=0) text+="("+str(Round(AngleOfView(yObj[iii]),-1))+"�)";
					ps->AddText(text,W+2,h, 3);
					if(xObjMax!=0){
						text="X="+str(xObj[iii]);
						if(xObj[iii]!=0) text+="("+str(Round(AngleOfView(xObj[iii]),-1))+"�)";
						ps->AddText(text,W+2,h-=3, 3);
					}
					if(write_hpr){
						pr=(ymax[iii]+ymin[iii])/2;
						if(AbeTheory==0){
							if( ImageHeight(height,point(xObj[iii],yObj[iii]),point(pr.x,pr.y),defocus,1,0,0,1) ){
								y1pr=height.y;
								x1pr=height.x;
							}
							else{
								y1pr=x1pr=0;
							}
						}
						else{
							if( ImageHeightAbe(height.y,height.x, yObj[iii],xObj[iii],"",pr.y,pr.x,defocus,1,AbeTheory) ){
								y1pr=height.y;
								x1pr=height.x;
							}
							else{
								y1pr=x1pr=0;
							}
						}
						ps->AddText("Y1pr="+str(Round(y1pr,-4)), W+2,h-=3, 3);
						if(xObjMax!=0){
							ps->AddText("X1pr="+str(Round(x1pr,-4)), W+2,h-=3, 3);
						}
					}
				}
				if(ii==2) {
					if(write_axis_title){
						ps->AddText("DX", -2,H+5, 3);
						ps->AddText("X", W+1,3.6, 3);
						ps->AddText(str(Ymax)+AfocalImageUnitStr(), 1,H+2, 3);
						ps->AddText(str(Round(Xmax,-3)), W,-0.4, 3);
					}
				}
			}
		}
	}

	delete[] yObj; delete[] xObj;
	delete[] ymax; delete[] ymin; delete[] xmax; delete[] xmin;
	delete[] found;
	for(iii=0; iii<=M; ++iii) {
		for(ii=0; ii<=2; ++ii){
			for(i=0; i<=N*2+1; ++i)  {
				delete[] X[iii][ii][i]; delete[] Y[iii][ii][i];
			}
			delete[] X[iii][ii]; delete[] Y[iii][ii];
		}
		delete[] X[iii]; delete[] Y[iii];
	}
	delete[] X; delete[] Y;

	*this=buf;  // shapes[DELTAH_GRAPH], current_shapes �͑�����Z�q�Ɋ܂܂�Ȃ��̂ł��̂܂܂ƂȂ�D
}

int cLens1::LensView() { return LENSVIEW; }

void cLens1::MakeLensView(std::string command,int FindPupil){
	double defocus=0;
	int i, ii, iii, j;
	const vector<double> Ez=vector<double>(0,0,1);
	list<ray_data> ray_list;
	list< vector<double> > sp,ep;
	list<long> col;
	double zo,zi;
	int color_start=1, color_end=1;
	const double COSTH=0.3420;  // TH=70deg
	matrix<double> T;
	int rays=5;   // �������̏����l
	int ray_from=0, ray_to=(int)LN; // ������`���͈͗��[�̖ʔԍ�. LN�͐�����݂��Ȃ����Ƃ��Ӗ�����D
	                                // �Ⴆ�΁Cray_to=k+1 �Ƃ���� for(iii=1; ...) �̃��[�v���� Turn, DoubleTurn
	                                // ���ɂ��ʐ����������Ƃ��Ɍ����`�悪�r���ʂŏI����Ă��܂��D
	int marginal_only=0, principal_only=0, rays_only=0, draw_virtual=0;
	int image_surface=1, object_surface=1;
	int with_ea=0;
	int draw_normal=0; 
	double normal_length=0.1;
	int window_fix=1;
	int distinguish_ref=0;
	int clip_enable=0;  // �f�t�H���g�ł� findpupil=0 �̂Ƃ��L���a���ɂ�邯���͍l�����Ȃ��Ō�����`���D
	                    // findpupil=0 �ł��������������Ƃ��� clip_enable=true �Ƃ���D
	cShapes* ps=CurrentShapes(LENSVIEW);
	int size_ini=ps->Size();
	std::string s1,s2,s3;
	bool b1,b2,b3;
	list<shape> centerlines;

	ps->Clear();
	if( !ea_is_defined() ) return;

	push(); // ����command�ɂ�背���Y�f�[�^�̕ύX���\�Ȃ��ߌ��f�[�^���X�^�b�N�ɕۑ�����D
	        // ******* �֐��̍Ō��pop����̂ŁC�ȉ��r����return�������Ȃ����� *******

	ps->xView(); ps->BackView();  // �f�t�H���g�̎���������cLens1�����W��+x����

	for(iii=1; iii<=sentences(command); ++iii){
		std::string com,s;
		com=sentence(command,iii);
		
		// �܂��Cscmd()�Ɋ܂܂��R�}���h����������D
		// ����ɂ��C�Ⴄ�������̌����𓯎��ɕ`���Ȃǂł���D
		scmd(com,0);
		s=arg(com,0);

		if(s=="col"){
			s1=arg(com,1); b1=is_numeric(s1);
			s2=arg(com,2); b2=is_numeric(s2);
			if(b1 && b2){
				int j1,j2;
				j1=atoi(s1.c_str());
				j2=atoi(s2.c_str());
				if( 1<=j1 && j1<=j2 && j2<=cn ){
					color_start=j1;
					color_end=j2;
				}
			}
			else if(b1){
				int j;
				j=atoi(s1.c_str());
				if( 1<=j && j<=cn ){
					color_start=color_end=j;
				}
			}
		}

		if(s=="clip_enable"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				clip_enable=atoi(s1.c_str());
			}
			else{
				clip_enable=1;
			}
		}

		if(s=="distinguish_ref"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				distinguish_ref=atoi(s1.c_str());
			}
			else{
				distinguish_ref=1;
			}
		}

		if(s=="marginal_only"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				marginal_only=atoi(s1.c_str());
			}
			else{
				marginal_only=1;
			}
			if(marginal_only) principal_only=0;
		}

		if(s=="principal_only"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				principal_only=atoi(s1.c_str());
			}
			else{
				principal_only=1;
			}
			if(principal_only) marginal_only=0;
		}

		if(s=="rays_only"){  // CAD�ŕ`�������w�z�u�}�ɁC���H���d�˂���ł���D
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				rays_only=atoi(s1.c_str());
			}
			else{
				rays_only=1;
			}
		}

		if(s=="draw_virtual"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				draw_virtual=atoi(s1.c_str());
			}
			else{
				draw_virtual=1;
			}
		}

		if(s=="rays"){  // ��������ݒ�( �������C"yfan(s)" "xfan(s)" �ł�17�{�܂łɐ�������� �j
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				rays=atoi(s1.c_str());
			}
		}

		if(s=="ray_from"){ // ������`���͈͂�ݒ�i�܂�Ԃ��ŉ��H���H�̌������d�Ȃ��Č��ɂ����Ƃ����Ɏg�p�j
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				ray_from=atoi(s1.c_str());
			}
		}
		if(s=="ray_to"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				ray_to=atoi(s1.c_str());
			}
		}

		if(s=="yfan" || s=="xfan" || s=="yfans" || s=="xfans" || s=="xyfan" || s=="xyfans"){
			double yObj0,xObj0,yObj,xObj,a;
			point p[18];  // �����W
			int i1,i2, n0, fans,l,ll;
			point ymax,ymin,xmax,xmin;
			ray_data ray;
			int errcode,n;
			std::string s0;

			ray_list.RemoveAll();
			make_coordinate(defocus); // ��̃R�}���h�ɂ��f�[�^�ύX�̉\��������K�v

			s1=arg(com,1); b1=is_numeric(s1);
			s2=arg(com,2); b2=is_numeric(s2);

			if(b1 && b2){
					yObj0=atof(s1.c_str());
					xObj0=atof(s2.c_str());
			}
			else if(b1){
				if(s=="yfan" || s=="yfans" || s=="xyfan" || s=="xyfans"){
					yObj0=atof(s1.c_str());
					xObj0=0;
				}
				if(s=="xfan" || s=="xfans"){
					yObj0=0;
					xObj0=atof(s1.c_str());
				}
			}
			else{
				yObj0=yObjectMax; xObj0=xObjectMax;
			}

			if( (s=="yfans" || s=="xfans" || s=="xyfans") && (yObj0!=0 || xObj0!=0) ){
				fans=4; 
			}
			else{
				fans=1;
				// yObj0=xObj0=0 �̂Ƃ��͌������d�Ȃ�̂�fans=1�Ƃ���D
			}

			for(l=1; l<=fans; l++){

				switch(l){
				case 1: a=1;    break;   // 100%����
				case 2: a=0.75; break;   // 75%
				case 3: a=0.5;  break;   // 50%
				case 4: a=0;    break;   // 0%
				}

				yObj=yObj0*a; xObj=xObj0*a;

				for(j=color_end; j>=color_start; j--){ // ������`�悪���g���̉��ɂȂ��Ă��܂�Ȃ��悤�C
				                                       // color_end���珈������D
					s0=s;
					for(ll=1; ll<=2; ll++){
						
						if(s0=="xyfan"){
							if(ll==1) s="yfan";
							if(ll==2) s="xfan";
						}
						else if(s0=="xyfans"){
							if(ll==1) s="yfans";
							if(ll==2) s="xfans";
						}
						else{
							if(ll==2) break;
						}

						if(FindPupilEdge(ymax,ymin,xmax,xmin,yObj,xObj,j,FindPupil,1)){
							if(s=="yfan" || s=="yfans"){
								p[1]=(ymax+ymin)/2;
								p[2]=ymax;
								p[3]=ymin;
							}
							if(s=="xfan" || s=="xfans"){
								p[1]=(xmax+xmin)/2;
								p[2]=xmax;
								p[3]=xmin;
							}
							p[4]=(p[1]+p[2])/2;
							p[5]=(p[1]+p[3])/2;
							
							p[6]=(p[2]+p[4])/2;
							p[7]=(p[4]+p[1])/2;
							p[8]=(p[1]+p[5])/2;
							p[9]=(p[5]+p[3])/2;

							p[10]=(p[2]+p[6])/2;
							p[11]=(p[6]+p[4])/2;
							p[12]=(p[4]+p[7])/2;
							p[13]=(p[7]+p[1])/2;
							p[14]=(p[1]+p[8])/2;
							p[15]=(p[8]+p[5])/2;
							p[16]=(p[5]+p[9])/2;
							p[17]=(p[9]+p[3])/2;
							
							i1=1; i2=17;
							if(marginal_only ){i1=2; i2=3;}  // ���������̂�
							if(principal_only){i1=1; i2=1;}  // ������̂�
							
							// ���˓��ɈÕ��̓�������Ƃ��ɔ�����p[17]�܂ŗp�ӂ��Ă��邪�C
							// �`��������rays�{�܂ŁD
							n0=ray_list.GetSize();
							for(ii=i1; ii<=i2 && ray_list.GetSize()-n0<rays; ++ii){
								errcode=RayTrace(yObj,xObj, p[ii].y,p[ii].x, 0,j,
								                 FindPupil ? 1 : clip_enable, FindPupil ? 1 : clip_enable,
									             0,0,0,vector<complex>(),0);
								if(errcode==0 || FindPupil==0){
									// FindPupil���U�̏ꍇ�́C�r���Ō���������ꂽ�肵���ꍇ�ł��C
									// ���̖ʂ܂ł͌�����`���D
									// �������C�������ʂƌ����Ȃ��Ƃ��͌��������s��Ȃ̂ŁC
									// �O�̖ʂ܂łƂ���D
									if(errcode==0){
										n=k+1;
									}
									else{
										n=NotThruSurf(errcode);
										if(ErrorMessage(errcode)=="not intersect"){
											n--;
										}
									}
									ray=ray_data(n);
									for(int i=0; i<=n; ++i){
										ray.color=this->color[j];
										ray.x[i]=this->x[i];
										ray.y[i]=this->y[i];
										ray.z[i]=this->z[i];
										ray.X[i]=this->X[i];
										ray.Y[i]=this->Y[i];
										ray.Z[i]=this->Z[i];
										ray.X1[i]=this->X1[i]; 
										ray.Y1[i]=this->Y1[i];
										ray.Z1[i]=this->Z1[i];
										ray.Ref[i]=this->N(i,1)>0 ? 0 : 1;
										ray.GrinRay_i=this->GrinRay_i;
										ray.GrinRay=this->GrinRay;
									}
									ray_list.AddTail(ray);
								}
							}
						}
					} // (ll=1; ll<=2; ll++)
					s=s0;
				} // for(j=color_end; j>=color_start; j--)
			} // for(l=1; l<=fans; l++)

			// zo,zi=region of drawing rays ////////
			if( fabs(this->s)<LN && N(0,1)*this->s<0 ){
				zo=this->s*1.1;
			}
			else if( fabs(t)<LN && N(0,1)*t<0 ){
				zo=t*1.1;
			}
			else{
				zo=-fabs( TotalThickness()+ea_max(1) )*sgn(N(0,1))*0.1;
			}

			if( dk(defocus)*N(k,1)>0 ){
					zi=dk(defocus)*1.1;
			}
			else{
				zi= fabs( TotalThickness()+ea_max(1) )*sgn(N(k,1))*0.1; // ea_max(1)���Ȃ��Ɩʐ�1�̂Ƃ�0
			}

			// zi�̏C��  ////////
			if( dk(defocus)*N(k,1)>0 ){
				if(ray_list.GetSize()>0){
					// Scheimpflug�n�ȂǑ��ʂ��X�΂��Ă���ꍇ�ŏI�ʂ��瑜�܂ł̋�����
					// �e���_�ő傫���قȂ�̂Ō������l������D
					ray_data ray;
					vector<double> a;
					double l,max;
					max=0;
					for(i=1; i<=ray_list.GetSize(); i++){
						ray=ray_list[i];
						if(ray.n==this->k+1){  // ���ʂ܂Ńf�[�^������ꍇ
							// a = k�ʒ��_���瑜�_�Ɍ������x�N�g��
							a=(o[k+1]-o[k])+ray.x[k+1]*ex[k+1]+ray.y[k+1]*ey[k+1]+ray.z[k+1]*ez[k+1];
							l=fabs(sProduct(a,ez[k]));
							if(l>max) max=l;
						}
					}
					if( max<fabs(TotalThickness(0,k+1))*0.3 ){
						zi=sgn(dk(defocus))*max*1.3;	
					}
					else{
						zi=sgn(dk(defocus))*max*1.1;
					}
				}
			}

			// �O���[�o�����W�����f�[�^���X�g sp,ep,col�𐶐�����  ///////////////////////
			{
				int ii,i,ig;
				vector<double> *v=new vector<double>[k+2];
				vector<double> Q, v1,v2, vn,vn1,vn2;
				long color0,color;					

				for(ii=1; ii<=ray_list.GetSize(); ++ii){
					ray_data ray;
					ray=ray_list[ii];
					color0=cSpectrum::ApproximateColor(Wavelength(ray.color));
					for(i=0; i<=ray.n; ++i){
						v[i]=o[i]+ray.x[i]*ex[i]+ray.y[i]*ey[i]+ray.z[i]*ez[i];
					}
					if( !( fabs(this->s)<LN && N(0,1)*this->s<0 || ExcludeVirtualObject!=0 ) ) {
						Q=ray.X[1]*ex[1]+ray.Y[1]*ey[1]+ray.Z[1]*ez[1];
						v[0]=v[1]-Q*(sProduct(v[1]-o0[1],ez0[1])-zo)/sProduct(Q,ez0[1]);
					}
					if(ray.n==k+1){   // ���ʂ܂Ńf�[�^������ꍇ
						Q=ray.X1[k]*ex[k]+ray.Y1[k]*ey[k]+ray.Z1[k]*ez[k];
						v[k+1]=v[k]+Q*(zi-sProduct(v[k]-o[k],ez[k]))/sProduct(Q,ez[k]);
					}
					for(i=1; i<=ray.n; ++i){
						if(this->N(i-1,1)<0 && distinguish_ref){
							color=RGBComplement(color0);
						}
						else{
							color=color0;
						}
						if(i==1) v1=v[0];
						if( is_acting_surf(i) || i==ray.n || draw_virtual ){
							// ���܁C���˂��Ȃ��ʂ͔�΂����Ƃɂ��C
							// �܂�Ȃ��~���[�ɗאڂ̉ˋ�ʂ�����ꍇ�Ƀ~���[�̌��܂Ō�����
							// �`���Ă��܂��C�Ȃǂ������D
							// �`�������̕��������̐i�s����Q�Ɠ������ǂ����ł����ʂł�������
							// �v���邪�C�܂�Ȃ��~���[�ɗאڂ̉ˋ�ʂ�����ꍇ�ȂǁC
							// Q�Ɠ������ł����Ă����̌����̎n�_����I�_�̑S�Ă����̌���
							// (�G�l���M�[���𔺂�����)�Ƃ͌���Ȃ��D
							// �������Ci==ray.n(�Ⴆ��ray.n�Ԗڂ̖ʂŌ����������Ă���)
							// �̏ꍇ�͔�΂��Ȃ��D

							if(IsGrin(i-1)){
								for(ig=1; ig<=ray.GrinRay_i.GetSize(); ig++){
									if(ray.GrinRay_i[ig]==i-1){
										v2=o0[i]+ray.GrinRay[ig].x*ex0[i]
												+ray.GrinRay[ig].y*ey0[i]
												+ray.GrinRay[ig].z*ez0[i];
										if( ray_from<=i-1 && i<=ray_to ){
											sp.AddTail(v1); ep.AddTail(v2); col.AddTail(color);
										}
										v1=v2;
									}
								}
							}
							else{
								v2=v[i];
								if( ray_from<=i-1 && i<=ray_to){
									sp.AddTail(v1); ep.AddTail(v2); col.AddTail(color);
								}
								v1=v2;
							}
						}
					}
					if(draw_normal){
						// �����Ɩʂ̌�_�ɖʖ@����`��
						// �ʂ̗L���a��normal_length�{�̒����ŕ`��
						for(i=1; i<=k; ++i){
							if(gname(i-1)!=gname(i) || asph_type(i)==IDEAL ){
								vn=surface_normal(i,ray.y[i],ray.x[i],0);
								vn=vn.x*ex[i]+vn.y*ey[i]+vn.z*ez[i];
								vn1=v[i]+vn*ea_max(i)*normal_length*0.5;
								vn2=v[i]-vn*ea_max(i)*normal_length*0.5;
								sp.AddTail(vn1); ep.AddTail(vn2); col.AddTail(0);
							}
						}
					}
				}

				delete []v;
			}
		}

		if(s=="view"){
			s1=arg(com,1); b1=is_numeric(s1);
			s2=arg(com,2); b2=is_numeric(s2);
			s3=arg(com,3); b3=is_numeric(s3);
			if(b1 && b2 && b3){
				// �摜���W��rox,roy,roz��]����D
				// ������-z�������C�}�̏オy�������C�E��x�����ł���D
				double rox,roy,roz;
				rox=atof(s1.c_str());
				roy=atof(s2.c_str());
				roz=atof(s3.c_str());
				ps->ViewAngle(rox,roy,roz);
				// ��F
				//  - "view ��x ��y 0" �� "view ��x 0 0; view 0 ��y 0" �Ɠ���
				//  - "view -90 0 0" ("top-view") ��ǉ�����ƌ��݂̕\�����ォ�猩��D
				//  - "view 0 90 0" ("right-view") ����ǉ�����ƌ��݂̕\���𗠂��猩��D
				// �Q�l�F
				//  �܂�Ȃ��̂���n�ŏ����l�ł̌����������Â炢���C
				//  ���� "view" ���w�肷��̂��ʓ|�ȂƂ��́C
				//  ���̖ʂ̉�]�ΐS�̐ݒ�ł��������𒲐��ł���D
			}
		}
		
		// ���݂̉摜����,��,��, ... ���猩��D
		if(s=="topview"   ) ps->TopView();
		if(s=="bottomview") ps->BottomView();
		if(s=="leftview"  ) ps->LeftView();
		if(s=="rightview" ) ps->RightView();
		if(s=="backview"  ) ps->BackView();

		// ���݂̎��������ɂ�����炸�C�����W��x,y,z�������猩��D
		if(s=="xview"     ) ps->xView();
		if(s=="yview"     ) { 
			ps->yView();
			ps->ViewAngle(0,0,-90);   // �����W+z�������E�ɂ���D
		}
		if(s=="zview"     ) ps->zView();
		
		if(s=="image_surface"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				image_surface=atoi(s1.c_str());
			}
			else{
				image_surface=1;
			}
		}

		if(s=="object_surface"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				object_surface=atoi(s1.c_str());
			}
			else{
				object_surface=1;
			}
		}

		if(s=="with_ea"){
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				with_ea=atoi(s1.c_str());
			}
			else{
				with_ea=1;
			}
		}

		if(s=="draw_normal"){
			s1=arg(com,1); b1=is_numeric(s1);
			draw_normal=1;
			if(b1){
				normal_length=atof(s1.c_str());
			}
		}

		if(s=="window_fix"){
			// �E�C���h�E�̌Œ�^��Œ�D�Ǐ��g����ێ��������Ƃ��ȂǂɎg���D
			s1=arg(com,1); b1=is_numeric(s1);
			if(b1){
				window_fix=atoi(s1.c_str());
			}
			else{
				window_fix=1;
			}
		}
	}  // next iii

	// T�s�񂪒�܂����̂Ŏ擾���� /////////
	T=ps->Tmatrix();

	make_coordinate(defocus);

	if(rays_only==0){
		// draw surface and stop ////////////////////
		for(i=1; i<=k; ++i){
			if( !is_dummy(i) ){  //  �_�~�[�ʂ͕`���Ȃ�
				if( fabs(sProduct(T*ez[i],Ez))<COSTH ){
					if(is_stop(i) && r(i)==0){  // �i��ʂŁC����������
												// (���ʂłȂ��i��̗� : OPD�̃`���b�p�[�͉~����)
						vector<double> v1,v2,v3,v4, r_in,r_out, oo;
						point p;
						const double R_OUT_MIN=2.5;   //  �Œᒼ�a5(mm)

						p=edge( ea_x(i),ea_y(i),EAtype(i),
								arg(complex(sProduct(T*ex[i],Ez),sProduct(T*ey[i],Ez)))+PI/2 );
						r_in=p.x*ex[i]+p.y*ey[i];
						r_out=1.4*r_in;
						if(abs(r_out)<R_OUT_MIN){
							r_out=r_out*(R_OUT_MIN/abs(r_out));  // �a���������Ƃ��C�O�a��傫�����Č��₷������D
						}
						oo=o[i]+EAdx(i)*ex[i]+EAdy(i)*ey[i];
						v1=oo+r_in; v2=oo+r_out;
						v3=oo-r_in; v4=oo-r_out;
						ps->Add3DLine(v1,v2);
						ps->Add3DLine(v3,v4);
					}
					else{
						const int DIV_MAX=40; vector<double> v[DIV_MAX+1]; int err[DIV_MAX+1];
						int div;
						point p;
						double x,y, phix,phiy;
						
						if(with_ea || is_stop(i) || is_mask(i)){
							phix=ea_x(i); phiy=ea_y(i);
						}
						else{
							phix=phi_x(i);phiy=phi_y(i);
						}

						p=edge( phix,phiy,EAtype(i),
							  arg(complex(sProduct(T*ex[i],Ez),sProduct(T*ey[i],Ez)))+PI/2 );
						div= (asph_type(i)>0 && asph_type(i)!=IDEAL) ? 40 : 20;   // �񋅖ʂׂ͍�����������
						for(ii=0; ii<=div; ++ii){
							x=-p.x+2*p.x*ii/div+EAdx(i);
							y=-p.y+2*p.y*ii/div+EAdy(i);
							v[ii]=o[i]+x*ex[i]+y*ey[i]+surface_sag(i,y,x,0,err[ii])*ez[i];
						}
						for(ii=0; ii<=div-1; ++ii) {
							if(err[ii]==0 && err[ii+1]==0) ps->Add3DLine(v[ii],v[ii+1]);
						}
					}
				}
				else{
					if(EAtype(i)==1){
						vector<double> v1,v2,v3,v4,v5,v6,v7,v8;
						double x1,y1, x2,y2;

						// (x1,y1),(x2,y2) �͎l���ʎ��̂����`�̑�1�ی���2�_�D
						if(with_ea){
							x1=fabs(ea_x(i)/2)-CHMx(i); y1=fabs(ea_y(i)/2);
							x2=fabs(ea_x(i)/2);         y2=fabs(ea_y(i)/2)-CHMy(i);	
						}
						else{
							double chmx,chmy;
							
							if(CHMx(i)==0 || CHMy(i)==0){
								chmx=chmy=0;
							}
							else{
								chmx=CHMx(i)+fabs(ea_y(i)-phi_y(i))/2*CHMx(i)/CHMy(i);
								chmy=CHMy(i)+fabs(ea_x(i)-phi_x(i))/2*CHMy(i)/CHMx(i);
							}
							x1=fabs(phi_x(i)/2)-chmx+EAdx(i);  y1=fabs(phi_y(i)/2)+EAdy(i);
							x2=fabs(phi_x(i)/2)+EAdx(i);       y2=fabs(phi_y(i)/2)-chmy+EAdy(i);
						}
						
						v1=o[i] +x1*ex[i] +y1*ey[i] +surface_sag(i, y1, x1,0)*ez[i];
						v2=o[i] +x2*ex[i] +y2*ey[i] +surface_sag(i, y2, x2,0)*ez[i];
						v3=o[i] +x2*ex[i] -y2*ey[i] +surface_sag(i,-y2, x2,0)*ez[i];
						v4=o[i] +x1*ex[i] -y1*ey[i] +surface_sag(i,-y1, x1,0)*ez[i];
						v5=o[i] -x1*ex[i] -y1*ey[i] +surface_sag(i,-y1,-x1,0)*ez[i];
						v6=o[i] -x2*ex[i] -y2*ey[i] +surface_sag(i,-y2,-x2,0)*ez[i];
						v7=o[i] -x2*ex[i] +y2*ey[i] +surface_sag(i, y2,-x2,0)*ez[i];
						v8=o[i] -x1*ex[i] +y1*ey[i] +surface_sag(i, y1,-x1,0)*ez[i];
						
						ps->Add3DLine(v2,v3);
						ps->Add3DLine(v4,v5);
						ps->Add3DLine(v6,v7);
						ps->Add3DLine(v8,v1);

						if(CHMx(i)>0 && CHMy(i)>0){
							ps->Add3DLine(v1,v2);
							ps->Add3DLine(v3,v4);
							ps->Add3DLine(v5,v6);
							ps->Add3DLine(v7,v8);
						}
					}
					else if(EAtype(i)==0){
						const int DIV=20;
						vector<double> v[DIV+1];
						double a,b, x,y;
						if(with_ea){
							a=ea_x(i)/2;
							b=ea_y(i)/2;
						}
						else{
							a=phi_x(i)/2;
							b=phi_y(i)/2;
						}
						for(ii=1; ii<=DIV; ii++){
							x=a*cos(PI*2/DIV*ii)+EAdx(i);
							y=b*sin(PI*2/DIV*ii)+EAdy(i);
							v[ii]=o[i]+x*ex[i]+y*ey[i]+surface_sag(i,y,x,0)*ez[i];
						}
						for(ii=2; ii<=DIV; ii++){
							ps->Add3DLine(v[ii-1],v[ii]);
						}
						ps->Add3DLine(v[DIV],v[1]);
					}
				}
			}
		}  // next i

		// draw image surface ///////////////////////
		if( image_surface && dk(defocus)*N(k,1)>0 ){  // ������N����(��)�̂Ƃ�z�̐�(��)�̕����֐i�s����Ƃ���
			if( fabs(sProduct(T*ez[k+1],Ez))<COSTH ){
				const int DIV=20; vector<double> v[DIV+1]; 
				double x,y;
				double xmax=0,ymax=0;
				point p;

				for(ii=1; ii<=ray_list.GetSize(); ++ii){
					ray_data ray;
					ray=ray_list[ii];
					if(ray.n==this->k+1){
						x=fabs(ray.x[k+1]);  if( x>xmax ) xmax=x;
						y=fabs(ray.y[k+1]);  if( y>ymax ) ymax=y;
					}
				}
				xmax*=1.2; ymax*=1.2;
				if( xmax<fabs(TotalThickness(0,k+1))*0.1 ) xmax=fabs(TotalThickness(0,k+1))*0.1;
				if( ymax<fabs(TotalThickness(0,k+1))*0.1 ) ymax=fabs(TotalThickness(0,k+1))*0.1;
				if(this->r(k+1)!=0){
					if( xmax>fabs(this->r(k+1)) ) xmax=fabs(this->r(k+1))*0.95;
					if( ymax>fabs(this->r(k+1)) ) ymax=fabs(this->r(k+1))*0.95;
				}
				p=edge( xmax*2,ymax*2,0,
					  arg(complex(sProduct(T*ex[k+1],Ez),sProduct(T*ey[k+1],Ez)))+PI/2 );
				for(ii=0; ii<=DIV; ++ii) {
					x=-p.x+2*p.x*ii/DIV;
					y=-p.y+2*p.y*ii/DIV;
					v[ii]=o[k+1]+x*ex[k+1]+y*ey[k+1]+surface_sag(k+1,y,x,0)*ez[k+1];
				}
				for(ii=0; ii<=DIV-1; ++ii) {
					ps->Add3DLine(v[ii],v[ii+1]);
				}
			}
			else{

			}
		}

		// draw object surface ///////////////////////
		if( object_surface && fabs(s)<LN && s*N(0,1)<0 ){
			if( fabs(sProduct(T*ez[0],Ez))<COSTH ){
				const int DIV=20; vector<double> v[DIV+1]; 
				double x,y;
				double xmax=0,ymax=0;
				point p;

				for(ii=1; ii<=ray_list.GetSize(); ++ii){
					ray_data ray;
					ray=ray_list[ii];
					x=fabs(ray.x[0]);  if( x>xmax ) xmax=x;
					y=fabs(ray.y[0]);  if( y>ymax ) ymax=y;
				}
				if(xmax==0) xmax=ymax*0.1;   // �Ⴆ�� ymax=1,xmax=0�̂Ƃ��C����p�̎���arg(..)��
				if(ymax==0) ymax=xmax*0.1;   // ���m��PI/4�ɂȂ�Ȃ���p��0�ɂȂ�C�`����Ȃ��Ȃ�D
				p=edge( xmax*2,ymax*2,0,
					  arg(complex(sProduct(T*ex[0],Ez),sProduct(T*ey[0],Ez)))+PI/2 );
				if(p!=point(0,0)){
					for(ii=0; ii<=DIV; ++ii) {
						x=-p.x+2*p.x*ii/DIV;
						y=-p.y+2*p.y*ii/DIV;
						v[ii]=o[0]+x*ex[0]+y*ey[0]+surface_sag(0,y,x,0)*ez[0];
					}
					for(ii=0; ii<=DIV-1; ++ii) {
						ps->Add3DLine(v[ii],v[ii+1]);
					}	
				}
			}
			else{

			}
		}

		// draw center line /////////////////////////
		{
			vector<double> v1,v2;

			// zo,zi= region of drawing center line  ////////
			if( fabs(this->s)<LN && N(0,1)*this->s<0 ){
				zo=this->s*1.1;
			}
			else if( fabs(t)<LN && N(0,1)*t<0 ){
				zo=t*1.1;
			}
			else{
				zo=-fabs( TotalThickness()+ea_max(1) )*sgn(N(0,1))*0.1;
			}

			if( dk(defocus)*N(k,1)>0 ){
				zi=dk(defocus)*1.1;
			}
			else{
				zi= fabs( TotalThickness()+ea_max(1) )*sgn(N(k,1))*0.1; // ea_max(1)���Ȃ��Ɩʐ�1�̂Ƃ�0
			}

			v1= o0[1]+ez0[1]*zo; 
			v2= o0[1];
			ps->Add3DLine(v1,v2);
			for(i=1; i<=k; ++i){
				switch(decenter_type(i)){
				case 1:
					v1=o0[i+1]-ez0[i+1]*d(i);
					break;
				case 2:
					v1=o0[i];
					break;
				default:
					v1=o[i];
					break;
				}
				v2= o0[i+1];
				if(i==k) v2=v1+((v2-v1)/abs(v2-v1))*zi*sgn(sProduct(v2-v1,ez[k]));
				centerlines.AddTail(ps->Add3DLine(v1,v2));  // centerlines�� "draw edge" �ŗ��p����
			}
		}

		// draw edge ////////////////////////////////
		for(i=1; i<=k-1; ++i){
			if(  fabs(sProduct(T*ez[i],  Ez))<COSTH && fabs(sProduct(T*ez[i+1],Ez))<COSTH ){
				if(is_solid(i) && !is_mask(i) && !is_mask(i+1)){
					vector<double> v11,v12,v21,v22;
					point p1,p2,p;
					double phix,phiy;
					if(with_ea){phix=ea_x(i); phiy=ea_y(i); }
					else       {phix=phi_x(i);phiy=phi_y(i);}
					p1=edge( phix,phiy,EAtype(i),
						   arg(complex(sProduct(T*ex[i],Ez),sProduct(T*ey[i],Ez)))+PI/2 );
					if(with_ea){phix=ea_x(i+1); phiy=ea_y(i+1); }
					else       {phix=phi_x(i+1);phiy=phi_y(i+1);}
					p2=edge( phix,phiy,EAtype(i+1),
						   arg(complex(sProduct(T*ex[i+1],Ez),sProduct(T*ey[i+1],Ez)))+PI/2 );
					if( p1.x*p2.x+p1.y*p2.y<0 ) p2=-p2;    // p1��p2���������ɂ��邱�Ƃ͕ۏ؂���Ă��Ȃ�����
					p.x=p1.x+EAdx(i); p.y=p1.y+EAdy(i);
					v11=o[i]+p.x*ex[i]+p.y*ey[i]+surface_sag(i,p.y,p.x,0)*ez[i];
					p.x=-p1.x+EAdx(i); p.y=-p1.y+EAdy(i);
					v12=o[i]+p.x*ex[i]+p.y*ey[i]+surface_sag(i,p.y,p.x,0)*ez[i];
					p.x=p2.x+EAdx(i+1); p.y=p2.y+EAdy(i+1);
					v21=o[i+1]+p.x*ex[i+1]+p.y*ey[i+1]+surface_sag(i+1,p.y,p.x,0)*ez[i+1];
					p.x=-p2.x+EAdx(i+1); p.y=-p2.y+EAdy(i+1);
					v22=o[i+1]+p.x*ex[i+1]+p.y*ey[i+1]+surface_sag(i+1,p.y,p.x,0)*ez[i+1];

					{
						double x1,y1,x2,y2,x3,y3,x4,y4;
						shape s;

						s=ps->Add3DLine(v11,v21);
						x1=s.x1; y1=s.y1;
						x2=s.x2; y2=s.y2;
						for(ii=1; ii<centerlines.GetSize(); ++ii){							
							x3=centerlines[ii].x1; y3=centerlines[ii].y1;
							x4=centerlines[ii].x2; y4=centerlines[ii].y2;
							if(IsCrossLines(x1,y1,x2,y2,x3,y3,x4,y4)) ps->RemoveTail();
							// ���S���ƌ�������G�b�W���͍폜����D
							// ����ɂ��C�Ⴆ��FT�̎Q�ƃv���Y���̕`�悪���Ղ��Ȃ�D
						}

						s=ps->Add3DLine(v12,v22);
						x1=s.x1; y1=s.y1;
						x2=s.x2; y2=s.y2;
						for(ii=1; ii<centerlines.GetSize(); ++ii){							
							x3=centerlines[ii].x1; y3=centerlines[ii].y1;
							x4=centerlines[ii].x2; y4=centerlines[ii].y2;
							if(IsCrossLines(x1,y1,x2,y2,x3,y3,x4,y4)) ps->RemoveTail();
						}
					}
				}
			}
			else{

			}
		}
	} // if(rays_only==0)

	// draw ray /////////////////////////////////
	{
		int i;
		for(i=1; i<=sp.GetSize(); i++){
			ps->Add3DLine(sp[i],ep[i],col[i]);
		}
	}

	if(window_fix==0){
		// �E�B���h�E���Œ肵�Ȃ��i= SetWindow()�����s����j
		//   ���̂Ƃ��ł��C���ύX��̍Ď��s�ł̓E�B���h�E(�g��k����)���ۑ������.
		//   �������C���XWindow���ɉ����Ȃ��Ƃ��͖�������SetWindow()�����s�D
		//     ��)
		//       �ŏ��̕`��ł��� -> size_ini==0
		//       ���̋��������ŕ��̂������Y���E�ɏo�Ă��� -> �f�t�H���g��Window�̒��ɉ����Ȃ��Ȃ� ps->Size()==0
		//       ���������� "size_ini==0" ���Ȃ��Ɖ����`�悳��Ȃ��D
		if(size_ini==0 || (ps->Size()!=size_ini)) ps->SetWindow(4.0/3.0);
	}
	else{
		// �E�B���h�E���Œ肷��i= SetWindow()�����s���Ȃ��j
		//   �������C�ŏ��̕`��(size_ini==0)�ł͖������ɃE�B���h�E�̐ݒ���s�Ȃ��D
		//   (�Ⴆ�΃t�@�C�����J�����Ƃ��C �R���X�g���N�^->initialize() �ɂ��shapes[]���N���A�����.)
		if(size_ini==0) ps->SetWindow(4.0/3.0);
	}

	pop();
}

int cLens1::SelectedView() { return SELECTEDVIEW; }

void cLens1::MakeSelectedView(int WhichDraw,int FindPupil,double FullScale,
								std::string DHCommand,std::string LVCommand){
	cShapes *ps;

	ps=CurrentShapes(SELECTEDVIEW);
	ps->Clear();

	switch(WhichDraw){
	case 1: // SA�O���t
		MakeSAGraph(FindPupil,0,FullScale,0);
		ps=CurrentShapes(SELECTEDVIEW);
		*ps=shapes[SA_GRAPH];
		break;
	case 2: // AS�O���t
		MakeASGraph(0,FindPupil,FullScale,0);
		ps=CurrentShapes(SELECTEDVIEW);
		*ps=shapes[AS_GRAPH];
		break;
	case 3: // DH�O���t
		MakeDeltaHGraph(0,0,FindPupil,0,0,0,FullScale,DHCommand,0,0);
		ps=CurrentShapes(SELECTEDVIEW);
		*ps=shapes[DELTAH_GRAPH];
		break;
	case 10: // LensView
		MakeLensView(LVCommand,FindPupil);
		ps=CurrentShapes(SELECTEDVIEW);
		*ps=shapes[LENSVIEW];
		break;
	case 11: // SA�O���t��LensView
		MakeSAGraph(FindPupil,0,FullScale,0);
		MakeLensView(LVCommand,FindPupil);
		ps=CurrentShapes(SELECTEDVIEW);
		*ps=shapes[SA_GRAPH];
		*ps=ps->AddInLine(shapes[LENSVIEW]);
		break;
	case 12: // AS�O���t��LensView
		MakeASGraph(0,FindPupil,FullScale,0);
		MakeLensView(LVCommand,FindPupil);
		ps=CurrentShapes(SELECTEDVIEW);
		*ps=shapes[AS_GRAPH];
		*ps=ps->AddInLine(shapes[LENSVIEW]);
		break;
	case 13: // DH�O���t��LensView
		MakeDeltaHGraph(0,0,FindPupil,0,0,0,FullScale,DHCommand,0,0);
		MakeLensView(LVCommand,FindPupil);
		ps=CurrentShapes(SELECTEDVIEW);
		*ps=shapes[DELTAH_GRAPH];
		*ps=ps->AddInLine(shapes[LENSVIEW]);
		break;
	}
}

int cLens1::SpotDiagram() { return SPOT; }

void cLens1::MakeSpotDiagram(cShapes& shapes,cSpot& sp,double FullScale,int WithFrame) {
	// WithFrame = 0(�g�Ȃ�), 1(�����`�g), 2(�~�`�g)
	int i;
	cPoint p;
	long col,r,g,b,back;
	double x,y;
	double div;

	shapes.Reset();
	if(sp.GetSize()==0) return;
	div= FullScale==0 ? fceil(sp.XYAbsMax()*2) : FullScale;
	
	if(WithFrame==1){
		shapes.AddLine(-div/2, div/2,  div/2, div/2,  0, SOLID);
		shapes.AddLine( div/2, div/2,  div/2,-div/2,  0, SOLID);
		shapes.AddLine( div/2,-div/2, -div/2,-div/2,  0, SOLID);
		shapes.AddLine(-div/2,-div/2, -div/2, div/2,  0, SOLID);
		shapes.AddLine(-div/2,0, -div/2+div/16,0,  0, SOLID);
		shapes.AddLine( div/2,0,  div/2-div/16,0,  0, SOLID);
		shapes.AddLine( 0,-div/2, 0,-div/2+div/16, 0, SOLID);
		shapes.AddLine( 0, div/2, 0, div/2-div/16, 0, SOLID);
		shapes.AddText( str(div) + sp.Unit + "/a side", div*0.05,div/2+div*0.1, div*0.07);
	}
	else if(WithFrame==2){
		const int N=24;
		double th1,th2;

		for(i=1; i<=N; ++i){
			th1=2.0*PI*(i-1)/N;
			th2=2.0*PI*i/N;
			shapes.AddLine((div/2)*cos(th1),(div/2)*sin(th1),(div/2)*cos(th2),(div/2)*sin(th2),0,SOLID);
			shapes.AddText( "��" + str(div) + sp.Unit, div*0.05,div/2+div*0.1, div*0.07);
		}
	}

	for(i=1; i<=sp.GetSize(); ++i){
		sp.GetData(p,i);
		x=p.p.x;
		y=p.p.y;
		col=cSpectrum::ApproximateColor(p.wl);
		if(p.weight<1){
			back=255;  // ���w�i�D���w�i�̂Ƃ���0�Ƃ���D
			r=back+(long)((double)(Rrgb(col)-back)*p.weight);
			g=back+(long)((double)(Grgb(col)-back)*p.weight);
			b=back+(long)((double)(Brgb(col)-back)*p.weight);
			col=rgb(r,g,b);
		}
		shapes.AddLine(x-div*0.005,y, x+div*0.005,y, col, SOLID);
	}
}

void cLens1::MakeSpotDiagram(double yObj,double xObj,int FindPupil,double defocus,int footprint,
						   int ColorStart,int ColorEnd,
						   int IsAreaSource,int IsLambert,double OriginIsGravityCenter,int Add,double FullScale,
						   int defocusN,double defocusStep,std::string commands,
						   int FieldDiagram,int DIVy,int DIVx,int UseObjectList,int Scan/*=0*/){
	// Scan==true �̂Ƃ��́CyScan(yObj),xScan(xObu),yObj=xObj=0 �Ƃ��ăX�|�b�g�_�C�A�O�������쐬����
	// <�Q�l>
	// �V���ɃX�|�b�g����������this->spot����cShapes���쐬����ꍇ�C
	// ColorStart>ColorEnd, Add!=0 �Ƃ���D

	// OriginIsGravityCenter<0 �̂Ƃ��́C���a�� |OriginIsGravityCenter| �̉~��`���i��܌��E�̕\�����ɗp����j

	int m,n,mc,nc, i,j,jj;
	double **yh,**xh, **defo;
	double rx,ry,r, max,temp, div, beta;
	int Randomize,CircularFrame;
	cShapes s;
	cShapes& shapes=*CurrentShapes(SPOT);
	list<point> object_list;
	cSpot *ptr;
	cSpot spot;
	cSpot** sp;
	std::string str1,str2;

	shapes.Reset();

	// FullScale < 0 �̂Ƃ��́C���a�� |FullScale| �̉~�`�g��t����
	if(FullScale<0){
		FullScale=-FullScale;
		CircularFrame=1;
	}
	else{
		CircularFrame=0;
	}

	// m=�s��, n=��
	if(FieldDiagram){  // FieldDiagram��UseObjectList�ɗD�悷��
		if(yObj==0 && xObj==0) { yObj=this->yObjectMax; xObj=this->xObjectMax; }
		if(yObj==0 && xObj==0) return;
		// yObj,xObj�̂����ꂩ��0�̂Ƃ�
		if(yObj==0) { yObj=xObj; }
		if(xObj==0) { xObj=yObj; }

		if(footprint==0){
			// ���ʂ̐}�̏㉺���E�𑜖ʂ̏㉺���E�ɍ��킹��D�t�b�g�v�����g�̂Ƃ��͉������Ȃ��D
			beta = M()==0 ? f() : M();
			yObj=sgn(beta)*fabs(yObj);
			xObj=sgn(beta)*fabs(xObj);
		}

		m= DIVy%2==0 ? DIVy+1 : DIVy;  // m,n����ɂ���
		n= DIVx%2==0 ? DIVx+1 : DIVx;
		mc=m/2+1;
		nc=n/2+1;
	}
	else{
		if(UseObjectList){
			object_list=this->object_list; object_list.EraseDuplicates();
			m=object_list.GetSize();
		}
		else{
			object_list.AddTail(point(xObj,yObj));  // ����yObj,xObj���g��
			m=1;
		}
		n=defocusN*2+1;  // �f�t�H�[�J�X��n
		nc=defocusN+1;   // �f�t�H�[�J�X�̒���
	}

	if( m>1 || n>1 ){
		// ���[�J���� cSpot spot ���g���̂ŁCAdd�͖���
		Add=0;
	}

	// ���������m��
	sp=new cSpot* [m+1];    for(i=1; i<=m; i++) sp[i]=new cSpot [n+1];
	xh=new double* [m+1];   for(i=1; i<=m; i++) xh[i]=new double [n+1];
	yh=new double* [m+1];   for(i=1; i<=m; i++) yh[i]=new double [n+1];
	defo=new double* [m+1]; for(i=1; i<=m; i++) defo[i]=new double [n+1];

	// �X�|�b�g���������Ƃ���Randomize=true�Ƃ���
	// �摜�̎����I�ȃm�C�Y����������D
	// ���Ȃ��Ƃ��͎����`��������悤��Radomize=false�Ƃ��Ă����D
	// (������Randomize�������悤�Ƃ������CClassWizard�ł͈�����16�ȏ���͂ł��Ȃ������D)
	Randomize= this->nSpot>=100 ? 1 : 0;
	
	for(i=1; i<=m; i++) for(j=1; j<=n; j++){
		if(FieldDiagram){
			xh[i][j]=-xObj+2*xObj/(n-1)*(j-1);
			yh[i][j]=-yObj+2*yObj/(m-1)*(i-1);
			defo[i][j]=defocus;
		}
		else{
			xh[i][j]=object_list[i].x;
			yh[i][j]=object_list[i].y;
			defo[i][j]=defocus-defocusStep*defocusN+defocusStep*(j-1);
		}
	}

	// �e���� i�C�f�t�H�[�J�X j ��cSpot���쐬
	for(i=1; i<=m; i++) for(j=1; j<=n; j++){

		if(FieldDiagram){
			if(i==mc && j==nc){
				ptr=&(this->spot);  // �������_�ł�*this->spot�ɍ쐬
			}
			else{
				spot.RemoveAll();
				ptr=&spot;
			}
		}
		else{
			if(i==1 && j==nc){
				ptr=&(this->spot);  // ��1���_�Ńf�t�H�[�J�X�̒����ł�*this->spot�ɍ쐬
			}
			else{
				spot.RemoveAll();
				ptr=&spot;
			}
		}
		
		if(words(commands)==0){   // command����̏ꍇ
			if(Scan==0){
				MakeSpot(*ptr, yh[i][j],xh[i][j],FindPupil,defo[i][j],footprint,
						 ColorStart,ColorEnd,IsAreaSource,IsLambert,OriginIsGravityCenter!=0,
						 Add,Randomize);
			}
			else{
				push();           // y(x)Scan �Ōn���ύX�����̂ŕۑ����Ă���
				yScan(yh[i][j]);
				xScan(xh[i][j]);
				MakeSpot(*ptr, 0,0,FindPupil,defo[i][j],footprint,
						 ColorStart,ColorEnd,IsAreaSource,IsLambert,OriginIsGravityCenter!=0,
						 Add,Randomize);
				pop();            // ��: pop���Ă�this->spot�͏㏑������Ȃ��iassignment���Q�Ɓj
			}
		}
		else{
			// ��j�f�t�H�[�J�X�ʒu�ɑ΂��� j�Ԗڂ�word�����s�Dword����1�̂Ƃ��͂��ׂẴf�t�H�[�J�X�ʒu�œ���word�����s�D
			push();
			// ��Fcommands="(rox 19 .1; add; rox 25 .1) (rox 19 .15; add; rox 25 .15)"
			if(words(commands)==1){
				str1=word(commands,1,1);  // str1="rox 19 .1; add; rox 25 .1"
			}
			else{
				str1=word(commands,j,1);  // j=1�̂Ƃ� str1="rox 19 .1; add; rox 25 .1"
			}
			if(!Add) ptr->RemoveAll();
			for(jj=1; jj<=sentences(str1); ++jj){
				str2=sentence(str1,jj);            // j=1,jj=1�̂Ƃ��Cstr2="rox 19 .1"
				cmd(str2,1);
				if(str2=="Add" || str2=="add" || jj==sentences(str1)){    // "Add (�����Ȃ�)"�͓��ʂȖ��߁D���̎��_�ŃX�|�b�g�_�C�A�O������ǉ�����D
					if(Scan==0){
						MakeSpot(*ptr, yh[i][j],xh[i][j],FindPupil,defo[i][j],footprint,
								 ColorStart,ColorEnd,IsAreaSource,IsLambert,OriginIsGravityCenter!=0,
								 1,Randomize);
					}
					else{
						yScan(yh[i][j]);
						xScan(xh[i][j]);
						MakeSpot(*ptr, 0,0,FindPupil,defo[i][j],footprint,
								 ColorStart,ColorEnd,IsAreaSource,IsLambert,OriginIsGravityCenter!=0,
								 1,Randomize);
					}
				}
			}
			pop();   // ��: pop���Ă�this->spot�͏㏑������Ȃ��iassignment���Q�Ɓj�D
		}

		sp[i][j]=*ptr;
	}

	// �X�|�b�g�̍L����̍ő�max������, �Ή�����X�P�[��div���v�Z����D
	max=0;
	for(i=1; i<=m; i++) for(j=1; j<=n; j++){
		if( (temp=sp[i][j].XYAbsMax()) > max ) max=temp;
	}
	div= FullScale==0 ? fceil(max*2) : FullScale;

	// sp[i][j]���}���쐬
	footprint= footprint>=0 ? footprint : -footprint;
	for(i=1; i<=m; i++) for(j=1; j<=n; j++){
		if(1<=footprint && footprint<=k && OriginIsGravityCenter==0){
		    // �t�b�g�v�����g�ŏd�S���_�ł͂Ȃ��Ƃ�
			//   �E�����Y�̗L���͈͂܂��͎Օ��͈͂�`��
			//   �E���ɁC�Օ��łȂ����(ry>0)EA�ӂ�\���e�L�X�g�������C
			//     �Օ��̂Ƃ��͘g��t����D

			rx=ea_x(footprint);   ry=ea_y(footprint);
			
			if(ry>0){  // �Օ��łȂ��Ƃ�
				r=Max(rx,ry);
				div= FullScale==0 ? r : div;  
				// ����FullScale�̎w�肪�Ȃ����r���g���C�����div(==FullScale)���g��
			}
			else{
				r=Max(-rx,-ry);
				div= FullScale==0 ? Max(r,div) : div;
				// ����FullScale�̎w�肪�Ȃ���΁C
				//   div����������΁Cdiv=r�Ƃ��ĎՕ����̗֊s��Window�̊O�ɏo�Ȃ��悤�ɂ���D
				//   div���傫����΁Cdiv=div�Ƃ��ĎՕ����̊O�̌����iFootPrint<0)
				//   ��Window�̊O�ɏo�Ȃ��悤�ɂ���D
			}

			MakeSpotDiagram(s,sp[i][j],div, ry>0 ? 0:1);
			switch(EAtype(footprint)){
			case 0:
				s.AddEllipse(0,0,rx/2,ry/2,24,0,SOLID);	
				if(ry>0){
					if(rx==ry){
						s.AddText("EA��"+str(r),0,ry/2+r*0.1,r*0.07);
					}
					else{
						s.AddText("EA��"+str(rx)+"x"+str(ry),0,ry/2+r*0.1,r*0.07);
					}
				}
				break;
			case 1:
				s.AddBox(rx/2,ry/2,-rx/2,-ry/2,0,SOLID,CHMx(footprint),CHMy(footprint));
				if(ry>0){
					s.AddText("EA"+str(rx)+"x"+str(ry),rx*0.1,ry/2+r*0.1,r*0.07);
				}
				break;
			}
		}
		else if(m==1 && n==1){
			MakeSpotDiagram(s,sp[i][j],div,CircularFrame ? 2 : 1);
		}
		else{
			MakeSpotDiagram(s,sp[i][j],div,0);
			
			// �X�P�[��������
			if(i==1 && j==n){  // �E��
				s.AddLine(-div/2,-div/2,         div/2,-div/2, 0,SOLID);
				s.AddLine(-div/2,-div/2+div/40, -div/2,-div/2-div/40);
				s.AddLine( div/2,-div/2+div/40,  div/2,-div/2-div/40);
				s.AddText( str(div) + this->spot.Unit + "/a side", -div*0.125,-div*0.525, div*0.1);
			}

			// ���̍���rms�a����������������
			if(FieldDiagram){
				if( (i==1 || i==m) && (j==1 || j==n) || (i==mc && j==nc) ){   // �l���ƒ���
					s.AddText("Y="+str(yh[i][j])+" X="+str(xh[i][j]), -div/2,div/2,div*0.1);
				}
			}
			else{
				if(j==1){
					s.AddText("Y="+str(yh[i][j])+" X="+str(xh[i][j]), -div/2,div/2,div*0.1);
				}
				s.AddText("rms="+str(Round(sp[i][j].RmsPhi(),-5)), -div/2,div*0.4,div*0.1);
			}
		}

		if(OriginIsGravityCenter<0){
			// ���a�� |OriginIsGravityCenter| �̉~��`��
			s.AddCircle(0,0,-OriginIsGravityCenter/2,24);
		}

		shapes.OffsetY=div*(i-1);
		shapes.OffsetX=div*(j-1);
		shapes.Add(s);	
	} // next i,j

	shapes.SetWindow(-div,div+div*(m-1), div+div*(n-1),-div);
	// �Ăяo�����v���O�����Ńr���[�̍������K�肷��Ƃ��C�E�C���h�E�������̏ꍇ��
	// �������̒[����ʂ̊O�֏o�Ă��܂��̂ŁC�E�C���h�E�̃A�X�y�N�g��𒲐�����D
	if(shapes.WindowAspectRatio()>1.5) shapes.SetWindow(1.5);
	
	// �������̊J��
	for(i=1; i<=m; i++){ delete [] sp[i];   } delete [] sp;
	for(i=1; i<=m; i++){ delete [] xh[i];   } delete [] xh;
	for(i=1; i<=m; i++){ delete [] yh[i];   } delete [] yh;
	for(i=1; i<=m; i++){ delete [] defo[i]; } delete [] defo;
}

void cLens1::MakeSpotDiagram(double yObj,double xObj,int FindPupil,double defocus,int footprint,
                 int ColorStart,int ColorEnd,
                 int IsAreaSource,int IsLambert,double OriginIsGravityCenter,int Add,double FullScale,
                 int defocusN,double defocusStep,int UseObjectList){
	// VC6 ClassWizard->�I�[�g���[�V����->���\�b�h�̒ǉ� �ł͈����̐��̏����15�ł���l�q
	// (����ȏ������ǉ��ł��Ȃ�)
	// MakeSpotDiagram2()�� commands �̑���� Add �������ɓ����Ă���D
	MakeSpotDiagram(yObj,xObj,FindPupil,defocus,footprint,
		            ColorStart,ColorEnd,
					IsAreaSource,IsLambert,OriginIsGravityCenter,Add,FullScale,
					defocusN,defocusStep,"",0,0,0,UseObjectList);
}

void cLens1::MakeSpotDiagram2(double yObj,double xObj,int FindPupil,double defocus,int footprint,
                 int ColorStart,int ColorEnd,
                 int IsAreaSource,int IsLambert,double OriginIsGravityCenter,double FullScale,
                 int defocusN,double defocusStep,std::string commands,int UseObjectList){
	// VC6 ClassWizard->�I�[�g���[�V����->���\�b�h�̒ǉ� �ł͈����̐��̏����15�ł���l�q
	// (����ȏ������ǉ��ł��Ȃ�)
	// MakeSpotDiagram()�� Add �̑���� commands �������ɓ����Ă��� (Add��commands�� "add" ���܂܂��邱�Ƃő�ւł���)
	MakeSpotDiagram(yObj,xObj,FindPupil,defocus,footprint,
		            ColorStart,ColorEnd,
				    IsAreaSource,IsLambert,OriginIsGravityCenter,0,FullScale,
				    defocusN,defocusStep,commands,0,0,0,UseObjectList);
}

void cLens1::MakeSpotDiagram3(int objects,double *yObj,double *xObj,int FindPupil,double defocus,int footprint,
	             int ColorStart,int ColorEnd,
	             int IsAreaSource,int IsLambert,double OriginIsGravityCenter,int Add,double FullScale,
				 int defocusN,double defocusStep){
	// VC6 ClassWizard->�I�[�g���[�V����->���\�b�h�̒ǉ� �ł͈����̐��̏����15�ł���l�q
	// (����ȏ������ǉ��ł��Ȃ�)
	// MakeSpotDiagram(..)�� int UseObjectList �̑���� int objects �������ɓ���CyObj,xObj��
	// �|�C���^�ɂȂ��Ă���D

	// (yObj[0],xObj[0]), ... (yObj[objects-1],xObj[objects-1]) �̕��_�ɂ�������X�|�b�g�_�C�A�O����
	// ���d�˂ĕ\������D
	int i, add;	

	for(i=0; i<=objects-1; i++){
		add = i==0 ? Add : 1;
		MakeSpotDiagram(yObj[i],xObj[i],FindPupil,defocus,footprint,
		                ColorStart,ColorEnd,
		                IsAreaSource,IsLambert,OriginIsGravityCenter,add,FullScale,
		                defocusN,defocusStep,0);
	}
}

void cLens1::MakeSpotFieldDiagram(double yObj,double xObj,int DIVy,int DIVx,
                 int FindPupil,double defocus,int footprint,
                 int ColorStart,int ColorEnd,
                 int IsAreaSource,int IsLambert,double OriginIsGravityCenter,double FullScale,std::string commands,int Scan){
	MakeSpotDiagram(yObj,xObj,
		            FindPupil,defocus,footprint,ColorStart,ColorEnd,
		            IsAreaSource,IsLambert,OriginIsGravityCenter,0,FullScale,0,0,commands,1,DIVy,DIVx,0,Scan);
}

int cLens1::AddSpotToImage(cImage& image,double weight){
	// this->spot �� image �ɉ�����
	int i,n,size, coherent;
	cPoint p;
	double ave,rms, r,th;

	size=spot.GetSize();
	
	ave=0;
	for(i=1; i<=size; ++i){
		spot.GetData(p,i);
		ave+=p.opl;
	}
	ave=ave/size;

	rms=0;
	for(i=1; i<=size; ++i){
		spot.GetData(p,i);
		rms+=(p.opl-ave)*(p.opl-ave);
	}
	rms=rms/(size-1);
	rms=sqrt(rms);              // ���H���̕W���΍�
	coherent = rms<0.1 ? 1:0;   // �W���΍���100��m(�`200��)�̂Ƃ��̂݃R�q�[�����g�Ƃ��Ĉ����D
	                            // �W���΍����傫���Ƃ��ɃR�q�[�����g�Ƃ���ƁCimage�̊e��f���̕����̌�����
	                            // ���H�����傫���Ȃ邽�߉摜���X�y�b�N����ƂȂ��Ă��܂��C�Ӑ}����摜�������Ȃ��D
	n=0;
	for(i=1; i<=size; i++){
		spot.GetData(p,i);
		r=sqrt(weight);
		th=2*PI*p.opl/(p.wl/1000000);
		if(coherent){
			if( image.AddComplexAmplitude(p.p.x,p.p.y,complex(r*cos(th),r*sin(th)))) ++n;  // ���f�U�������Z
		}
		else{
			if( image.Add(p.p.x,p.p.y,weight)) ++n;  // ���x�����Z
		}
	}

	return n;
}

int cLens1::AddSpotToImage(double weight){
	// this->spot �� this->image �ɉ�����
	return AddSpotToImage(this->image,weight);
}

int cLens1::CopySpotToImage(int xPixels,int yPixels,double xWidth) {
	// xPixels,yPixels : image�̈�ӂ̉�f��(xPixels�̂ݎw��̏ꍇ�͐����`�摜�Ƃ���j
	// xWidth          : image�̉����̎�����(0���w��̏ꍇ��|x|,|y|�̍ő��4�{�Ƃ���)
	int n;
	const int DEFAULT_SIZE=256;
	
	ImageClear();
	if(xPixels<=0) xPixels=DEFAULT_SIZE;
	if(yPixels<=0) yPixels=xPixels;
	SetImageXPixels(xPixels);
	SetImageYPixels(yPixels);

	if(xWidth<=0) xWidth=spot.XYAbsMax()*4;
	SetImagePitch( xWidth/xPixels );
	n=AddSpotToImage(1);
	ImageNormalize(1);  // �ő�l��1�ɋK�i������ 
	return n;
}

void cLens1::Image0ToImage(int FindPupil,double defocus,int ColorStart,int ColorEnd,
						   double ImageWidth/*=0*/,int ImageXPixels/*=0*/,int ImageYPixels/*=0*/,int RemoveMoire/*=0*/){
	// this->image0�𕨑̂Ƃ��āC�X�|�b�g�_�C�A�O�����̂����ݍ��݂ɂ�葜�摜���v�Z���C
	// this->image�ɐݒ肷��D����S�̂�P������Ƃ͂��Ȃ��D
	
	// �����̏ꍇ�C�����̓f�t�H���g�ł悢�D(�摜���͎����v�Z�D��f�����͌��摜�Ɠ����D���A�������͖����j

	int i,j, m0,n0,m,n;
	double p0,p, x0,y0,weight, l0,l;
	cImage temp;

	m0=image0.GetYpixels();
	n0=image0.GetXpixels();
	p0=image0.GetPitch();

	if(m0==0 || n0==0 || p0==0) return;

	if(ImageXPixels==0 || ImageYPixels==0){
		SetImageYPixels(m0);
		SetImageXPixels(n0);
	}
	else{
		SetImageXPixels(ImageXPixels);
		SetImageYPixels(ImageYPixels);
	}

	if(ImageWidth==0){
		l0=p0*n0;
		l=l0*( fabs(s)>LN ? fabs(f()) : fabs(M()) );
		n=GetImageXPixels();
		p=l/n;
		SetImagePitch(p);
	}
	else{
		SetImagePitch(ImageWidth/GetImageXPixels());
	}

	ImageClear();
	temp=this->image;

	for(i=1; i<=m0; i++) for(j=1; j<=n0; j++){
		weight=image0.GetIntensity(i,j);
		if(weight>0){
			// image0�̒��������_�Ƃ��Ĉʒu(x0,y0)���v�Z����
			y0=-((double)i-((double)m0+1)/2)*p0;
			x0= ((double)j-((double)n0+1)/2)*p0;
			MakeSpot(y0,x0,FindPupil,defocus,0,ColorStart,ColorEnd,0,0,0,0);
			AddSpotToImage(weight);
			if(RemoveMoire!=0) AddSpotToImage(temp,1);  // ���x����l�ȕ��̂̑��i���A�������p�j
		}
	}

	// ���A������
	//   ��������Ȃ��Ɗi�q��̃m�C�Y���o�邱�Ƃ�����D���A���ƍl������D
	//   ImageX(Y)Pixels,�܂���ImagePitch��ω�������Ɗi�q�̊Ԋu���ς��D
	if(RemoveMoire!=0){
		m=GetImageYPixels();
		n=GetImageXPixels();
		for(i=1; i<=m; i++) for(j=1; j<=n; j++){
			if(temp.GetIntensity(i,j)!=0){
				SetImageIntensity(i,j,GetImageIntensity(i,j)/temp.GetIntensity(i,j));
			}
		}
	}

	ImageNormalize(1);  // �ő�l��1�ɋK�i������
}

void cLens1::Image0ToImageScan(int FindPupil,double defocus,int ColorStart,int ColorEnd,
	                           double xScanStart,double xScanEnd,double yScanStart,double yScanEnd,
							   double xImageStart,double xImageEnd,double yImageStart,double yImageEnd,
	                           int xScanPoints,int yScanPoints,int RemoveMoire/*=0*/){
	// �{�֐���Image0ToImage�֐��𑖍����w�n�ɑΉ����������́D
	// this->image0�𕨑̂Ƃ��āC�X�|�b�g�_�C�A�O�����̂����ݍ��݂ɂ�葜�摜���v�Z���C
	// this->image�ɐݒ肷��D

	// xScanStart�`xScanEnd ����� yScanStart�`yScanEnd �͈̔͂𑖍�����D���̓��C
	// xImageStart�`xImageEnd ����� yImageStart�`yImageEnd �͈̔͂ŉ摜���쐬����D
	// �摜�쐬�͈͂𑖍��͈͂�苷�����邱�Ƃɂ��C�v�Z���Ԃ�Z�k�ł���D
	// ����͓��ɁC�摜�����ɂ�艡�F������␳����Ƃ������͈͎��̂��k���ł��Ȃ�
	// (�k�����Ă��܂��ƕ␳�W���̕ϊ����K�v�ɂȂ��Ă��܂�)�̂ŗL���ł���D
	// x(y)ImageStart=x(y)ImageEnd=0 �Ƃ���Ή摜�쐬�͈͂͑����͈͂Ɠ����ł���D

	int i0,j0, i,j, m,n, ii;
	double intensity,intensity0, rimg,x,y, th,th_max,image0_r, r,r0, thx,thy;
	cPoint p;
	cImage temp;

	if(image0.GetYpixels()==0 || image0.GetXpixels()==0 || image0.GetPitch()==0) return;
	push();  // xScan,yScan�ŕύX�����̂ŏ����l��ۑ�����D

	SetImageYPixels(yScanPoints);  // �o�͉摜�̉�f�� (= �����_���j
	SetImageXPixels(xScanPoints);
	ImageClear();
	temp=this->image;
	rimg=this->rImage();

	image0_r=image0.GetHeight()/2;                               // image0�̍����̔����D���ʂ����ʂłȂ��Ƃ��p
	th_max= fabs(image0_r/rimg)>1 ? PI/2 : asin(image0_r/rimg);  // ���ʋȗ����S���猩��image0�̒[�̋p�D���ʂ����ʂłȂ��Ƃ��p

	if(xImageStart==0 && xImageEnd==0){ xImageStart=xScanStart; xImageEnd=xScanEnd; }
	if(yImageStart==0 && yImageEnd==0){ yImageStart=yScanStart; yImageEnd=yScanEnd; }

	for(i=1; i<=yScanPoints; ++i) for(j=1; j<=xScanPoints; ++j){
		thy=yScanStart+(yScanEnd-yScanStart)*(i-1)/(yScanPoints-1);
		thx=xScanStart+(xScanEnd-xScanStart)*(j-1)/(xScanPoints-1);
		
		if( (thy-yImageStart)*(thy-yImageEnd)<=0 && (thx-xImageStart)*(thx-xImageEnd)<=0 ){ // �����͈͂̓��̉摜�쐬�͈͂̂ݏ������s��
			yScan(thy);
			xScan(thx);

			n=MakeSpot(0,0,FindPupil,defocus,0,ColorStart,ColorEnd,0,0,0,0);  // n=�X�|�b�g���D 
																			  // �ȉ������̓s����C�ڕ��ʓ��e�͂��Ȃ�
			intensity=intensity0=0;
			for(ii=1; ii<=n; ++ii){
				spot.GetData(p,ii);
				if(is_plane(this->k+1)){ // ���ʂ����ʂ̂Ƃ�
					x=p.p.x;
					y=p.p.y;
				}
				else{  // ���ʂ����ʂłȂ��Ƃ��DImage0�͑��ʂ̋ȗ����S���猩���}�ł���Ƃ���D
					   // ���Ȃ킿�C�ȗ����S���猩���Image0��������z���ɐ����ȕ��ʂɓ\���Ă���悤�Ɍ�����D
					x=p.p.x;
					y=p.p.y;
					r0=sqrt(x*x+y*y);
					th=asin(fabs(r0/rimg));  // �ȗ����S���猩��(x,y)�������ƂȂ��p�x
					r=image0_r*(th/th_max);
					x*=(r/r0);
					y*=(r/r0);
				}
				image0.PixelLocation(i0,j0,x,y);
				intensity +=p.weight*image0.GetIntensity(i0,j0);
				if(RemoveMoire!=0) intensity0+=p.weight;
			}
			
			image.SetIntensity(i,j,intensity);
			if(RemoveMoire!=0) temp.SetIntensity(i,j,intensity0);
		}
	}

	if(RemoveMoire!=0){
		m=GetImageYPixels();
		n=GetImageXPixels();
		for(i=1; i<=m; i++) for(j=1; j<=n; j++){
			if(temp.GetIntensity(i,j)!=0){
				SetImageIntensity(i,j,GetImageIntensity(i,j)/temp.GetIntensity(i,j));
			}
		}
	}

	ImageNormalize(1);
	pop();
}

double cLens1::GetImage0Width() { return image0.GetWidth(); }
void   cLens1::SetImage0Width(double value) {
	image0.SetWidth(value);
	image0.SetYpitch(image0.GetXpitch());  // xy�����ŉ�f�s�b�`�����������Ƃ�O��Ƃ���
}
double cLens1::GetImageWidth() { return image.GetWidth(); }
void   cLens1::SetImageWidth(double value) { image.SetWidth(value); }
int  cLens1::GetImageXPixels() { return image.GetXpixels(); }
void cLens1::SetImageXPixels(int value) { image.SetXpixels(value); }
int  cLens1::GetImageYPixels() { return image.GetYpixels(); }
void cLens1::SetImageYPixels(int value) { image.SetYpixels(value); }
double cLens1::GetImagePitch() { return image.GetPitch(); }
void   cLens1::SetImagePitch(double value) { image.SetPitch(value); }
void cLens1::ImageClear() { image.Clear(); }
int cLens1::ImageAdd(double x,double y,double weight) { return image.Add(x,y,weight); }
void cLens1::ImageNormalize(double new_max) { image.Normalize(new_max); }
void cLens1::ImageReverseXY(){ image.ReverseXY(); }
double cLens1::GetImageIntensity(int i,int j) { return image.GetIntensity(i,j); }
void   cLens1::SetImageIntensity(int i,int j,double value) { image.SetIntensity(i,j,value); }
int cLens1::LoadImage0(std::string filename){
	return image0.Open(filename);
}
int cLens1::LoadImage0FromBmp(std::string filename){
	return image0.OpenFromBmp(filename);
}
int cLens1::LoadImage(std::string filename){
	return image.Open(filename);
}
int cLens1::LoadImageFromBmp(std::string filename){
	return image.OpenFromBmp(filename);
}
int cLens1::SaveImage(std::string filename){
	return image.Save(filename);
}
int cLens1::SaveImageAsBmp(std::string filename,double gamma){
	image.Normalize(255);
	return image.SaveAsBmp(filename,gamma);
}


double cLens1::EyeRefractivePower(double y,double VD,double SAm,int UseCorneaHeight,double PupilLocation){
	// ��̋��ܗ͂��v�Z����D��1�g���ɂ��Čv�Z����D�f�B�I�v�^�[�P�ʂŕԂ��D
	//   ��̃f�[�^��*this�Ɏw�肷��D
	//     ��1�ʂ��p���\�ʁD
	//     ���ʂ����D���̈ʒu��s1fix�Ŏw�肷��D
	//   Y = ���ܗ͂��v�Z�����������
	//   VD = �p�����_�Ɗዾ�ʒu�̋���(���̂Ƃ��ዾ���p���̃}�C�i�X���D�ʏ�͐�)�D
	//        �ዾ�ʒu�͋��ܗ͂̌��_�ƂȂ�D
	//   SAm = ���ʎ����ʂɑ΂���ŗǑ��ʈʒu
	//   UseCorneaHeight = y���p���ƌ����̌�_�Ƃ��邩�ǂ����D�U�̂Ƃ��͓����ʏ�̌��������g���D
	//                     (���ʂ͋U�ɂ���)
	//   PupilLocation = �p������̓��ʒu

	cLens1 X;
	double s,s10,s1,sa,da;

	X=*this;

	// s���ߎ����_�ɉ��ݒ肷��
	s=-X.Reversed(1).s1();

	// ���������������Ƃ��Đݒ肷�邽�߂̗L���a�̐ݒ�
	if(UseCorneaHeight){
		X.EAy(1)=y*2;
	}
	else{
		X.Add(0,2,0);
		X.EAy(1)=100;
		X.EAy(2)=y*2;
		X.EAy(3)=100;
		X.d(1)=PupilLocation;
		X.d(2)=-PupilLocation;
		X.gname(1)=X.gname(0);
		X.gname(2)=X.gname(0);
	}

	do{
		// �����ߎ��ōŗǑ��ʂ����ƈ�v���镨�_�����߂�D
		//   SAm==1�ł���Ί�ꂩ��̋t�ǐՂɂ���x�ŋ��܂邪�C
		//   SAm!=1�̂Ƃ��ɐ��m�Ɍv�Z����ɂ͒����ߎ������Ȃ��͗l�D
		X.s=s;
		s10=X.s1();    // �ߎ����_
		sa=X.LSA();    // ���ʎ���
		s1=s10+sa*SAm; // ���ʎ������l���������ʒu

		// s1��dk()�ɋ߂Â��邽�߁C�j���[�g���� zz'=-NN'f^2 �̔����C
		//   da = -dz'/NN'/f^2, a=1/z
		// �ɂ�蒀���ߎ�����D
		da=-(X.dk(0)-s1)/X.Get_N(X.k,1)/X.f()/X.f();
		s=1/(1/s+da);
	} while(fabs(s1-X.dk(0))>=0.00001);

	return 1000/(s+VD);
}

void cLens1::EyeToAxial(double diopter,double VD){
	// ���ʒus1fix��ς��āC���ܗ͂�diopter�ɂȂ�悤�ɂ���D
	cLens1 X;

	X=*this;
	X.s= diopter==0 ? LN : 1000/diopter-VD;
	this->s1fix=X.s1();
}

void cLens1::EyeToRefractive(double diopter,double VD){
	// �p�����ܗ�(��1�ʂ̋��ܗ�)��ς��āC���ܗ͂�diopter�ɂȂ�悤�ɂ���D
	cLens1 X;
	double k,k_target, dp,p0,p;

	X=*this;
	k_target= diopter==0 ? 0 : 1/(1000/diopter-VD); // VD=0���ܗ͖ڕW
	X.Reverse();
	k=-1/X.s1();                  // ����VD=0���ܗ�
	dp=-k_target-(-k);            // �p���ʃp���[��ω��������
	p0=1/X.f(X.k,X.k);            // ���̊p���ʃp���[
	p=p0+dp;                      // �C�������p���ʃp���[
	X.AdjustFocalLength(1/p,X.k,X.k,0,0,0);
	X.Reverse();
	this->surf[1]=X.surf[1];
}

double cLens1::EyeZernikeCoefficient(double PupilZ,double RefSurfZ,double PupilPhi,int j){
	// ���Zernike�W�����v�Z����D��1�g���ɂ��Čv�Z����D�߂�l��0.001(=1��m)�P�ʂƂ���D
	//   ��̃f�[�^��*this�Ɏw�肷��D
	//     ��1�ʂ��p���\�ʁD
	//     ���ʂ����D���̈ʒu��s1fix�Ŏw�肷��D
	//   PupilZ   = �p�����_������˓��܂ł̋����iOPD�ł�3.7)
	//   RefSurfZ = �p�����_����Q�ƕ��ʂ܂ł̋��� (OPD�ł�3.7)
	//   PupilPhi = ���˓����a�i�����a�����肷��j
	//   j = �v�Z����Zernike�W���̔ԍ��i�X�^���_�[�h�I�[�_�[�Ƃ��C�s�X�g����0�Ƃ���)

	point p;
	cLens1 X=*this;
	double dummy;

	X.Add(0,2,0);            // �n�̍ŏ��ɎQ�ƕ���(��1��)�Ɠ��˓���(��2��)��ǉ�����
	X.EAy(1)=100;            // �Q�ƕ���
	X.EAy(2)=PupilPhi;       // ���˓�
	X.d(1)=PupilZ-RefSurfZ;  // �Q�ƕ��ʂ�����˓�
	X.d(2)=-PupilZ;          // ���˓�����p��
	X.s=LN;

	X.ImageHeight(p.y,p.x,dummy, 0,0,"",0,0,0,1,0,0,0,0);  // �����ɉ��������ˌ����Ɗ��̌�_p�����߂�
	                                                 // �ʏ�C�����v�̌�����ɌŎ��������邽�߁C
	                                                 // p�������ƍl������D
	X.Reverse();
	X.SetOPDMapZernikeR0Fix(PupilPhi/2);       // ���˓����a��1�Ƃ���
	X.SetOPDMapZernikeJBase(0);                // �s�X�g�����0���Ƃ���
	X.SetOPDMapZernikeIsFringeOrder(0);        // �X�^���_�[�h�I�[�_�[�Ƃ���
	X.MakeOPDMap(p.y,p.x,0,1,1,0,0,cZernike::nNumber(j,0,0),-1,X.k-1,0,0,""); // p�𕨓_�Ƃ��C
	                                                                          // ���˓���̍��W�œW�J����
	return X.OPDMapZernikeCoefficient(j);
}

double cLens1::ModelEyeS1(double diopter,double VD){
	// �O���X�g�����h�����͌^��i���U�t�j�̐����̌�ʂ�����܂ł̋�����Ԃ�
	cLens1 X(7,1);
	X.Set_color(1,"e");   
	                               X.gname(0)="1";
	X.r(1)= 7.7;   X.d(1)=0.5;     X.gname(1)="376450";  // cornea
	X.r(2)= 6.8;   X.d(2)=3.1;     X.gname(2)="336500";  // cornea
	X.r(3)= 0;     X.d(3)=0;       X.gname(3)="336500";  // iris
	X.r(4)=10;     X.d(4)=0.546;   X.gname(4)="386490";  // cortex
	X.r(5)= 7.911; X.d(5)=2.419;   X.gname(5)="406465";  // nuclear
	X.r(6)=-5.76;  X.d(6)=0.635;   X.gname(6)="386490";  // nuclear
	X.r(7)=-6;                     X.gname(7)="336500";  // cortex

	X.EyeToAxial(diopter,VD);
	return X.s1fix;
}

double cLens1::ModelEyeFundusR(double diopter,double VD){
	// �O���X�g�����h�����͌^��ɐݒ肷����ȗ����a��Ԃ��D
	// 0diopter��ł̊��ȗ����a��12�o�Ƃ��C�������܈ُ��ł�
	// ���̈ʒu��(e���ߎ��l)�����ዅ���a���ω������Ƒz�肵�C
	// ���̈ʒu���̔����������ȗ����a���ω�����Ƃ���D

	return 12+(ModelEyeS1(diopter,VD)-ModelEyeS1(0,VD))/2;
}

cLens1 cLens1::LMTestLens(double sph,double cyl){
	//   sph = ���ʓx��(��ʒ��_���ܗ�) -25<=sph<=+25D
	//   cyl = �~���x��  cyl<=0
	//
	// ISO 9342-1 �� Table 1 �Ɋ�Â��ă����Y���[�^�p�e�X�g�����Y�𐶐�����D
	//    �ގ��� S-NSL5
	//    ��ʃJ�[�u�� Table 1 �̎w��l�̒�����Ԃɂ��
	//    ���S���� Table 1 �̉����l�Ƃ���
	//    �~������x���Ƃ���	

	cXYList BSPList, TcList;
	cLens1 X(2,3);
	double BSP,FSP, Tc,N;

	if( sph<-25 || 25<sph || cyl>0 ) return X;

	const std::string GNAME="S-NSL5";

	BSPList.AddData(-25,-25);  TcList.AddData(-25,2); // Table 1 ���
	BSPList.AddData(-20,-20);  TcList.AddData(-20,2);
	BSPList.AddData(-15,-15);  TcList.AddData(-15,2);
	BSPList.AddData(-10,-12);  TcList.AddData(-10,2);
	BSPList.AddData( -5, -9);  TcList.AddData( -5,2);
	BSPList.AddData(  5, -5);  TcList.AddData(  5,3);
	BSPList.AddData( 10, -3);  TcList.AddData( 10,3);
	BSPList.AddData( 15, -1);  TcList.AddData( 15,5);
	BSPList.AddData( 20,  0);  TcList.AddData( 20,7);
	BSPList.AddData( 25,  0);  TcList.AddData( 25,9);

	BSP=BSPList.y(sph);         // ��ʋ��ܗ�(dptr)
	Tc=TcList.y(sph);           // ���S��
	N=Re(Index(GNAME,546.07));  // �����Y���ܗ�(e��)
	FSP=(sph/1000-BSP/1000)/((Tc/N)*(sph/1000-BSP/1000)+1) *1000;   // �O�ʋ��ܗ�(dptr, ��-h�ǐՎ����瓱�o)

	X.r(1)= FSP==0 ? 0 : 1000*(N-1)/FSP;
	X.r(2)= BSP==0 ? 0 : 1000*(1-N)/BSP;
	if(cyl<0){
		X.asph_type(2)=TOROID;
		X.ry(2)=X.r(2);
		X.rx(2)= BSP+cyl==0 ? 0 : 1000*(1-N)/(BSP+cyl);
	}
	X.d(1)=Tc;
	X.gname(1)=GNAME;
	X.EAy(1)=X.EAy(2)=30;  // �L���a�͂Ƃ肠������30�Ƃ���

	return X;
}

void cLens1::ToLMTestLens(double sph,double cyl){
	*this=LMTestLens(sph,cyl);
}

double cLens1::LensCenterPowerWorn(std::string SCA,double So,double Co,double Ao,
							       double WrapAngle,double TiltAngle){
	cLens1 X;
	double S,C,A, delta,xp,yp;

	X=LMTestLens(So,Co);
	X.Add(2,1,0);    // �Q�Ɩʂ�ǉ�
	X.Afocal=1;

	X.decenter_type(1)=1; X.rox(1)=TiltAngle; X.roy(1)=-WrapAngle; X.roz(1)=-Ao;
	X.decenter_type(2)=1; X.ret(2)=1;

	delta=X.delta(1,2,1);  // �����Y�̕\�ʂ���O����_�܂ł̋���
	// (xp,yp) = �����Y���X�����Ƃ��̑O����_�̈ʒu�i�����ɓ�������̕Ίp�̓[���j
	xp= delta*sin(X.roy(1)*PI/180);
	yp=-delta*cos(X.roy(1)*PI/180)*sin(X.rox(1)*PI/180);

	X.SCA(S,C,A,0,0,yp,xp,0,0,1);

	if     (SCA=="S" || SCA=="s") return S;
	else if(SCA=="C" || SCA=="c") return C;
	else if(SCA=="A" || SCA=="a") return A;
	else                          return 0;
}

