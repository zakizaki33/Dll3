#include "stdafx.h"
#include "cFilm.h"

///// protected members //////////////////////////////////////////////////////////

void cFilm::create(){
	mname=new std::string [k+2];
	t=new double [k+2];
	tVariable= new int [k+2];
	ndTol=new double[k+2];
	dTol=new double[k+2];
}

void cFilm::initialize(){
	int i;
	th_med_deg=0;
	wl0=550;
	tMode=OPTICAL;
	mname[0]="1";
	mname[k+1]="1";
	for(i=0; i<=k+1; i++) {
		t[i]=0;
		tVariable[i]=0;
		ndTol[i]=0.03;     // 光学膜厚誤差初期値 3%
		dTol[i]=0;         // 物理膜厚誤差初期値 0nm
	}
	SubThickness=0;
	M.redim(2,2);
	Ms.redim(2,1);
}

void cFilm::erase(){
	delete [] mname;
	delete [] t;
	delete [] tVariable;
	delete [] ndTol;
	delete [] dTol;
}

void cFilm::assignment(const cFilm& x){
	int i;

	filename=x.filename;
	wl0=x.wl0;
	tMode=x.tMode;
	th_med_deg=x.th_med_deg;
	for(i=0; i<=k && i<=x.k; i++) {
		mname[i]=x.mname[i];
		t[i]=x.t[i];
		tVariable[i]=x.tVariable[i];
		ndTol[i]=x.ndTol[i];
		dTol[i]=x.dTol[i];
	}
	mname[k+1]=x.mname[x.k+1];
	t[k+1]=x.t[x.k+1];
	SubThickness=x.SubThickness;
}

complex cFilm::costh(int i,double th_med_deg,double wl_nm){
	complex sinth_med,costh, N,Nmed;
	N   =Index(mname[i],  wl_nm);
	Nmed=Index(mname[k+1],wl_nm);
	sinth_med=sin(th_med_deg*PI/180);
	costh=sqrt( 1-(Nmed*sinth_med/N)*(Nmed*sinth_med/N) );
	
	// cosθ=±√(1-(sinθ)^2) のどちらの解をとるかは次による．
	//   “光学薄膜と成膜技術 (李正中 アルバック)” 2.2.5より，
	//    cosθは第4象限にあることが要求される．
	//
	// 以前は引数zが負の実数のとき，
	//     sqrt(z) = ( √|z|) * exp(i*arg(z)/2))
	// のarg(z)が-PIでなく，PIを返していた．したがって，sqrt(z)の偏角はPI/2であった．
	// 全反射では 1-(sinθ)^2は負の実数なので
	//     cosθ= √(1-(sinθ^2))
	// の偏角はPI/2，すなわちcosθは第2象限にあった．
	// 全反射を含むときに反射波の位相計算がうまくいかなかったのはこのためだと考えられる．
	// 2011.03.15

	return Im(costh)>0 ? -costh : costh;  // x=±sprt(...) の2つの解のうち第4象限側の解をとる
}

complex cFilm::Y(int i,double th_med_deg,double wl_nm,char s_or_p){
	// アドミタンスを計算する
	complex y,N,cs;

	N=Index(mname[i],wl_nm);
	cs=costh(i,th_med_deg,wl_nm);
	switch(s_or_p){
		case 's': y=N*cs; break;
		case 'p': y=N/cs; break;
	}
	return y;
}

void cFilm::calc_M(double th_med_deg,double wl_nm,char s_or_p){
	int i;
	matrix<complex> m(2,2);
	double d;
	complex N, cs, g, y, cosg,sing;
	const complex j(0,1);

	M=unit(M);
	for(i=1; i<=k; i++){
		d=this->d(i);
		N=Index(mname[i],wl_nm);
		cs=costh(i,th_med_deg,wl_nm);
		if(s_or_p=='s') y=N*cs;
		if(s_or_p=='p') y=N/cs;
		g=2*PI/wl_nm*(N*d*cs);
		cosg=cos(g);
		sing=sin(g);
		m[1][1]=cosg;      m[1][2]=j*sing/y;
		m[2][1]=j*sing*y;  m[2][2]=cosg;
		M=m*M;
	}
}

void cFilm::calc_Ms(double th_med_deg,double wl_nm,char s_or_p){
	// 基板も含めた特性行列
	// C/B = Ms[2][1]/Ms[1][1] が媒質から見たアドミタンスとなる
	complex y;

	calc_M(th_med_deg,wl_nm,s_or_p);
	y=Y(0,th_med_deg,wl_nm,s_or_p);
	Ms[1][1]=M[1][1]+M[1][2]*y;
	Ms[2][1]=M[2][1]+M[2][2]*y;
}

double cFilm::phase_thickness(double th_med_deg,double wl_nm){
	// 主波長での位相膜厚(rad)
	// ArgTs(),ArgTp()で使用する．
	int i;
	complex cs,N;
	// complex test_comp;
	double d,p;

	p=0;
	for(i=1; i<=k; i++){
		cs=costh(i,th_med_deg,this->wl0);
		d=this->d(i);
		N=Index(mname[i],this->wl0);
		// test_comp = (2.0 * PI / wl_nm) * N * d * cs;
		p+=Re((2.0*PI/wl_nm)*N*d*cs);
		// p += Re(test_comp);
	}
	return p;
}

void cFilm::make_xylist(std::string AboutWhat,double wl1,double wl2,double dwl){
	double wl;
	std::string s;
	char buf[100];

	xylist.RemoveAll();
	dwl= dwl*(wl2-wl1)>=0 ? dwl : -dwl;
	for(wl=wl1; (wl-wl1)*(wl-wl2)<=1e-10; wl+=dwl) {
		s=AboutWhat;
		sprintf(buf," %g",wl); s+=buf;
		xylist.AddData( wl,atof(cmd(s,1).c_str()) );
	}
}

void cFilm::make_spectrum(std::string AboutWhat,double wl1,double wl2,double dwl){
	double wl;
	std::string s;
	char buf[100];

	spectrum.RemoveAllData();
	for(wl=wl1; (wl-wl1)*(wl-wl2)<=1e-10; wl+=dwl) {
		s=AboutWhat;
		sprintf(buf," %g",wl); s+=buf;
		spectrum.AddData( wl,atof(cmd(s,1).c_str()) );
	}
}


///// public members /////////////////////////////////////////////////////
double cFilm::NdOfSingleLayerARCoat(std::string material,double wl_nm,double th_deg) {
	// wl_nm  : wavelength in air
	// th_deg : incident angle in air 
	double N=Re(Index(material,wl_nm));
	double sinth=sin(th_deg*PI/180);
	return N==0 ? 0 : wl_nm/4/sqrt(1-sinth*sinth/N/N);
}

double cFilm::EquiFilmN(std::string SideMaterial1,std::string CenterMaterial2,
                        double TotalPhi,double Phi2Phi1Ratio)
{
	double Phi1,Phi2,Nd1,Nd2;
	complex N;
	cFilm eq(3); eq.wl0=550;
	eq.set_mname(1,SideMaterial1); eq.set_mname(2,CenterMaterial2); eq.set_mname(3,SideMaterial1);
	Phi1=TotalPhi/(2+Phi2Phi1Ratio); Phi2=Phi1*Phi2Phi1Ratio;
	Nd1=(Phi1*PI/180)*eq.wl0/2/PI; Nd2=(Phi2*PI/180)*eq.wl0/2/PI;
	eq.set_t(1,Nd1); eq.set_t(2,Nd2); eq.set_t(3,Nd1);
	eq.calc_M(0,eq.wl0,'s');
	N=sqrt( eq.M[2][1]/eq.M[1][2] );
	return Im(N)==0 ? Re(N) : 999;
}
double cFilm::EquiFilmN(std::string SideMaterial1,std::string CenterMaterial2,
                        double TotalNd,double Nd2Nd1Ratio,double wl_nm)
{
	double Nd1,Nd2;
	complex N;
	cFilm eq(3);
	eq.set_mname(1,SideMaterial1); eq.set_mname(2,CenterMaterial2); eq.set_mname(3,SideMaterial1);
	Nd1=TotalNd/(2+Nd2Nd1Ratio); Nd2=Nd1*Nd2Nd1Ratio;
	eq.set_t(1,Nd1); eq.set_t(2,Nd2); eq.set_t(3,Nd1);
	eq.calc_M(0,wl_nm,'s');
	N=sqrt( eq.M[2][1]/eq.M[1][2] );
	return Im(N)==0 ? Re(N) : 999;
}

double cFilm::EquiFilmPhi(std::string SideMaterial1,std::string CenterMaterial2,
                         double TotalPhi,double Phi2Phi1Ratio)
{
	double Phi1,Phi2,Nd1,Nd2;
	cFilm eq(3); eq.wl0=550;
	eq.set_mname(1,SideMaterial1); eq.set_mname(2,CenterMaterial2); eq.set_mname(3,SideMaterial1);
	Phi1=TotalPhi/(2+Phi2Phi1Ratio); Phi2=Phi1*Phi2Phi1Ratio;
	Nd1=(Phi1*PI/180)*eq.wl0/2/PI; Nd2=(Phi2*PI/180)*eq.wl0/2/PI;
	eq.set_t(1,Nd1); eq.set_t(2,Nd2); eq.set_t(3,Nd1);
	eq.calc_M(0,eq.wl0,'s');
	return acos( Re(eq.M[1][1]) )*180/PI;   // acosの引数が複素数は未対応
}
double cFilm::EquiFilmPhi(std::string SideMaterial1,std::string CenterMaterial2,
                         double TotalNd,double Nd2Nd1Ratio,double wl_nm)
{
	double Nd1,Nd2;
	cFilm eq(3);
	eq.set_mname(1,SideMaterial1); eq.set_mname(2,CenterMaterial2); eq.set_mname(3,SideMaterial1);
	Nd1=TotalNd/(2+Nd2Nd1Ratio); Nd2=Nd1*Nd2Nd1Ratio;
	eq.set_t(1,Nd1); eq.set_t(2,Nd2); eq.set_t(3,Nd1);
	eq.calc_M(0,wl_nm,'s');
	return acos( Re(eq.M[1][1]) )*180/PI;   // acosの引数が複素数は未対応
}

cFilm::cFilm(int k/*=0*/){
	this->k=k;
	create();
	initialize();
}

cFilm::cFilm(const cFilm &x){
	k=x.k;
	create();
	assignment(x);
}

cFilm::~cFilm(){
	erase();
}

void cFilm::Reset(){
	*this=cFilm(0);
}

cFilm& cFilm::operator=(const cFilm& x){ 
	erase();
	k=x.k;
	create();
	assignment(x);
	return *this;
}

int cFilm::push() { Stack.push(*this); return 1; }
int cFilm::pop()  { return Stack.pop(*this);     }
int cFilm::undo() { return Stack.undo(*this);    }
int cFilm::redo() { return Stack.redo(*this);    }

void cFilm::AddToList() { List.AddTail(*this); }
void cFilm::ClearList() { List.RemoveAll(); }
int  cFilm::GetListData(int num) { return List.GetData(*this,num); }

int cFilm::get_k(){
	return k;
}
void cFilm::set_k(int value){
	cFilm buf;
	buf=*this;
	erase();
	k=value;
	create();
	initialize();
	assignment(buf);
}

double cFilm::get_wl0(){
	return wl0;
}
void cFilm::set_wl0(double value_nm){
	wl0=value_nm;
}

std::string cFilm::get_tMode(){
	switch(tMode){
		case OPTICAL:   return "optical";   break;
		case PHYSICAL:  return "physical";  break;
		case OPTICALQW: return "opticalqw"; break;
		default:        return "";
	}
}
void cFilm::set_tMode(std::string value){
	if( value=="optical" || value=="Optical" || value=="OPTICAL" ){
		tMode=OPTICAL;
	}
	else if( value=="physical" || value=="Physical" || value=="PHYSICAL" ){
		tMode=PHYSICAL;
	}
	else if( value=="opticalqw" || value=="OpticalQw" || value=="OpticalQW" || value=="OPTICALQW" ){
		tMode=OPTICALQW;
	}
}

void cFilm::tToOptical(){
	int i;
	for(i=1; i<=k; i++) t[i]=Round(nd(i),-2);
	tMode=OPTICAL;
}

void cFilm::tToPhysical(){
	int i;
	for(i=1; i<=k; i++) t[i]=Round(d(i),-2);
	tMode=PHYSICAL;
}

void cFilm::tToOpticalQW(){
	int i;
	for(i=1; i<=k; i++) t[i]=Round(nd(i)/(wl0/4),-3);
	tMode=OPTICALQW;
}

double cFilm::d(int i){
	// 物理膜厚
	if(0<=i && i<=k){
		switch(tMode){
			case OPTICAL:
				return get_t(i)/Re(Index(mname[i],wl0));	
				break;
			case PHYSICAL:
				return get_t(i);
				break;
			case OPTICALQW:
				return get_t(i)*(wl0/4)/Re(Index(mname[i],wl0));
				break;
			default:
				return 0;
		}
	}
	else{
		return 0;
	}
}

void cFilm::set_d(int i,double value){
	if(0<=i && i<=k){
		switch(tMode){
			case OPTICAL:
				set_t(i,value*Re(Index(mname[i],wl0)));
				break;
			case PHYSICAL:
				set_t(i,value);
				break;
			case OPTICALQW:
				set_t(i,value*Re(Index(mname[i],wl0))/(wl0/4));
				break;
		}
	}
}

double cFilm::nd(int i,double wl){
	if(0<=i && i<=k){
		return Re(Index(mname[i],wl))*d(i);
	}
	else{
		return 0;
	}
}

double cFilm::nd(int i){
	return nd(i,wl0);
}

double cFilm::get_th_med_deg(){
	return th_med_deg;
}
void cFilm::set_th_med_deg(double th_deg){
	th_med_deg=th_deg;
}

std::string cFilm::get_mname(int i){
	if(0<=i && i<=k+1) return mname[i];
	else return "";
}
void cFilm::set_mname(int i,std::string name){
	if(0<=i && i<=k+1)  mname[i]=name;
}

void cFilm::ReplaceMaterial(std::string from,std::string to){
	int i;
	
	if(from!="" && to!=""){
		for(i=1; i<=k; ++i){
			if(get_mname(i)==from) set_mname(i,to);
		}
	}
}

double cFilm::get_t(int i){
	if(0<=i && i<=k) return t[i];
	else return 0;
}
void cFilm::set_t(int i,double value){
	if(0<=i && i<=k) t[i]=value;
}

int cFilm::get_tVariable(int i){
	if(0<=i && i<=k) return tVariable[i];
	else return 0;
}
void cFilm::set_tVariable(int i,int value){
	if(0<=i && i<=k) tVariable[i]=value;
}

double cFilm::get_ndTol(int i){
	if(0<=i && i<=k) return ndTol[i];
	else return 0;
}
void cFilm::set_ndTol(int i,double value){
	if(0<=i && i<=k) ndTol[i]=value;
}

double cFilm::get_dTol(int i){
	if(0<=i && i<=k) return dTol[i];
	else return 0;
}
void cFilm::set_dTol(int i,double value){
	if(0<=i && i<=k) dTol[i]=value;
}

double cFilm::TotalNd(){
	int i;
	double sum=0;
	for(i=1; i<=k; i++) sum+=nd(i,wl0);
	return sum;
}

double cFilm::TotalD(){
	int i;
	double sum=0;
	for(i=1; i<=k; i++) sum+=d(i);
	return sum;
}

int cFilm::SubstrateAbsorptionEnable=0;

double cFilm::getSubThickness(){
	return SubThickness;
}
void cFilm::setSubThickness(double value){
	SubThickness=value;
}

complex cFilm::rs(double th_med_deg,double wl_nm){
	complex B,C, y_med,r;
	
	calc_Ms(th_med_deg,wl_nm,'s');
	B=Ms[1][1];
	C=Ms[2][1];
	y_med=Y(k+1,th_med_deg,wl_nm,'s');
	r=(y_med-C/B)/(y_med+C/B) *sqrt(Tfilters(wl_nm));
	return r;
}
complex cFilm::rs(double wl_nm){
	return rs(th_med_deg,wl_nm);
}

complex cFilm::rp(double th_med_deg,double wl_nm){
	complex B,C, y_med,r;
	
	calc_Ms(th_med_deg,wl_nm,'p');
	B=Ms[1][1];
	C=Ms[2][1];
	y_med=Y(k+1,th_med_deg,wl_nm,'p');
	r=(y_med-C/B)/(y_med+C/B) *sqrt(Tfilters(wl_nm));
	return r;

	// 何冊かの教科書では，位相変化のない状態が arg(rs)=0, arg(rp)=180° となるが，
	// ここでは arg(rs)=arg(rp)=0 が位相変化のない状態である．この方が直感的に把握し易い．
	// これは，s,p方向の正方向の取り方の違いによる．
	// Cを面の法線方向，Qを入射光線方向，Q1を透過あるいは反射光線方向とすると，
	// s,pの正方向は，
	//   S=(CxQ)/|CxQ|
	//   P=SxQ
	//   P1= SxQ1 (透過)
	//      -SxQ1 (反射)
	// とする．
	// これは，例えば，
	//    (1)“光学薄膜と成膜技術 李正中 アルバック”
	//    (2)“光学薄膜フィルターデザイン 小檜山光信 オプトロニクス”
	// と同じ取り方である．
	//    (3)“光応用技術1988 波動光学 田中俊一 JOEM”
	//    (4)“光学薄膜 藤原史郎編 共立出版”
	// とは異なる．
	// (3)はP1の取り方が透過反射共 P1=SxQ1 となり，この点簡単だが，
	// 反射においてs,p位相差が与えられない状態が，arg(rs)=0，arg(rp)=180°
	// となり，直感的に把握しにくい．（2010.4.8まではこの方法であった．）
	//
	// 例: 空気からガラスへ垂直入射するとき，rs,rp<0（位相反転）．
	//     (位相差を混同しないこと．s,p位相差は0である．)
}
complex cFilm::rp(double wl_nm){
	return rp(th_med_deg,wl_nm);
}

complex cFilm::ts(double th_med_deg,double wl_nm){
	complex B,C, y_med, N_sub, t;
	
	calc_Ms(th_med_deg,wl_nm,'s');
	B=Ms[1][1];
	C=Ms[2][1];
	N_sub=Index(mname[0],wl_nm);
	y_med=Y(k+1,th_med_deg,wl_nm,'s');
	t=2*y_med/(y_med*B+C);
	if(Im(N_sub)==0) return t*sqrt( Tsub(th_med_deg,wl_nm)*Tfilters(wl_nm) ); else return 0;
}
complex cFilm::ts(double wl_nm){
	return ts(th_med_deg,wl_nm);
}

complex cFilm::tp(double th_med_deg,double wl_nm){
	complex B,C, y_med, N_sub, t;
	
	calc_Ms(th_med_deg,wl_nm,'p');
	B=Ms[1][1];
	C=Ms[2][1];
	N_sub=Index(mname[0],wl_nm);
	y_med=Y(k+1,th_med_deg,wl_nm,'p');
	t=2*y_med/(y_med*B+C);
	if(Im(N_sub)==0) return t*sqrt( Tsub(th_med_deg,wl_nm)*Tfilters(wl_nm) ); else return 0;
}
complex cFilm::tp(double wl_nm){
	return tp(th_med_deg,wl_nm);
}

double cFilm::Rs(double th_med_deg,double wl_nm){
	complex r;

	r=rs(th_med_deg,wl_nm);
	return sqabs(r);
}
double cFilm::Rs(double wl_nm){
	return Rs(th_med_deg,wl_nm);
}

double cFilm::Rp(double th_med_deg,double wl_nm){
	complex r;

	r=rp(th_med_deg,wl_nm);
	return sqabs(r);
}
double cFilm::Rp(double wl_nm){
	return Rp(th_med_deg,wl_nm);
}

double cFilm::Rave(double th_med_deg,double wl_nm){
	return ( Rs(th_med_deg,wl_nm)+Rp(th_med_deg,wl_nm) )/2;
}
double cFilm::Rave(double wl_nm){
	return Rave(th_med_deg,wl_nm);
}	

double cFilm::Ts(double th_med_deg,double wl_nm){
	complex y_med,y_sub, t;
	
	t=ts(th_med_deg,wl_nm);
	y_med=Y(k+1,th_med_deg,wl_nm,'s');
	y_sub=Y(0,  th_med_deg,wl_nm,'s');
	return Re(y_sub/y_med)*sqabs(t) *Tsub(th_med_deg,wl_nm);
}
double cFilm::Ts(double wl_nm){
	return Ts(th_med_deg,wl_nm);
}

double cFilm::Tp(double th_med_deg,double wl_nm){
	complex y_med,y_sub, t;

	t=tp(th_med_deg,wl_nm);
	y_med=Y(k+1,th_med_deg,wl_nm,'p');
	y_sub=Y(0,  th_med_deg,wl_nm,'p');
	return Re(y_sub/y_med)*sqabs(t) *Tsub(th_med_deg,wl_nm);
}
double cFilm::Tp(double wl_nm){
	return Tp(th_med_deg,wl_nm);
}

double cFilm::Tave(double th_med_deg,double wl_nm){
	return ( Ts(th_med_deg,wl_nm)+Tp(th_med_deg,wl_nm) )/2;
}
double cFilm::Tave(double wl_nm){
	return Tave(th_med_deg,wl_nm);
}

double cFilm::As(double th_med_deg,double wl_nm){
	return 1-Ts(th_med_deg,wl_nm)-Rs(th_med_deg,wl_nm);
}
double cFilm::As(double wl_nm){
	return As(th_med_deg,wl_nm);
}

double cFilm::Ap(double th_med_deg,double wl_nm){
	return 1-Tp(th_med_deg,wl_nm)-Rp(th_med_deg,wl_nm);
}
double cFilm::Ap(double wl_nm){
	return Ap(th_med_deg,wl_nm);
}

double cFilm::Aave(double th_med_deg,double wl_nm){
	return ( As(th_med_deg,wl_nm)+Ap(th_med_deg,wl_nm) )/2;
}
double cFilm::Aave(double wl_nm){
	return Aave(th_med_deg,wl_nm);
}

double cFilm::RsBothSide(double th_med_deg,double wl_nm){
	cFilm x; double R,Rinv, T,Tinv, A;
	R=Rs(th_med_deg,wl_nm);
	T=Ts(th_med_deg,wl_nm);
	A=Tsub(th_med_deg,wl_nm);
	x=*this;
	x.reverse(wl_nm);
	Rinv=x.Rs(wl_nm);
	Tinv=x.Ts(wl_nm);
	return R + (T*Tinv*R)/(1-Rinv*Rinv*A*A);
}
double cFilm::RsBothSide(double wl_nm){
	return RsBothSide(th_med_deg,wl_nm);
}

double cFilm::RpBothSide(double th_med_deg,double wl_nm){
	cFilm x; double R,Rinv, T,Tinv, A;
	R=Rp(th_med_deg,wl_nm);
	T=Tp(th_med_deg,wl_nm);
	A=Tsub(th_med_deg,wl_nm);
	x=*this;
	x.reverse(wl_nm);
	Rinv=x.Rp(wl_nm);
	Tinv=x.Tp(wl_nm);
	return R + (T*Tinv*R)/(1-Rinv*Rinv*A*A);
}
double cFilm::RpBothSide(double wl_nm){
	return RpBothSide(th_med_deg,wl_nm);
}

double cFilm::RaveBothSide(double th_med_deg,double wl_nm){
	return ( RsBothSide(th_med_deg,wl_nm)+RpBothSide(th_med_deg,wl_nm) )/2;
}
double cFilm::RaveBothSide(double wl_nm){
	return RaveBothSide(th_med_deg,wl_nm);
}

double cFilm::TsBothSide(double th_med_deg,double wl_nm){
	cFilm x; double Rinv, T,Tinv, A;
	T=Ts(th_med_deg,wl_nm);
	A=Tsub(th_med_deg,wl_nm);
	x=*this;
	x.reverse(wl_nm);
	Rinv=x.Rs(wl_nm);
	Tinv=x.Ts(wl_nm);
	return (T*Tinv*A)/(1-Rinv*Rinv*A*A);
}
double cFilm::TsBothSide(double wl_nm){
	return TsBothSide(th_med_deg,wl_nm);
}

double cFilm::TpBothSide(double th_med_deg,double wl_nm){
	cFilm x; double Rinv, T,Tinv, A;
	T=Tp(th_med_deg,wl_nm);
	A=Tsub(th_med_deg,wl_nm);
	x=*this;
	x.reverse(wl_nm);
	Rinv=x.Rp(wl_nm);
	Tinv=x.Tp(wl_nm);
	return (T*Tinv*A)/(1-Rinv*Rinv*A*A);
}
double cFilm::TpBothSide(double wl_nm){
	return TpBothSide(th_med_deg,wl_nm);
}

double cFilm::TaveBothSide(double th_med_deg,double wl_nm){
	return ( TsBothSide(th_med_deg,wl_nm)+TpBothSide(th_med_deg,wl_nm ))/2;
}
double cFilm::TaveBothSide(double wl_nm){
	return TaveBothSide(th_med_deg,wl_nm);
}

// ArgRs(),ArgRp()は0次反射点の位相を返す．
//
//  【例1】
//   媒質物質からなる単層膜でも入射角により位相が変わってしまうが，
//   凹面鏡で言えば面位置がずれて焦点ずれが起こっている状態だから当然といえる．
//  【例2】
//   ある点からみて光線進行方向の点の位相はその点の過去の位相であるから，
//   “遅れ”ており，位相はより小さい値になる．
//  【例3】
//   媒質物質からなる単層膜での反射主成分(基板表面反射：光路差 2ndcosθ)の位相は，
//      波長が長くなると遅れが少なくなり，進む.
//      入射角が大きくなると光路差が短くなって遅れが少なくなり，進む．
double cFilm::ArgRs(double th_med_deg,double wl_nm){
	return arg(rs(th_med_deg,wl_nm))*180/PI;
}
double cFilm::ArgRs(double wl_nm){
	return ArgRs(th_med_deg,wl_nm);
}

double cFilm::ArgRp(double th_med_deg,double wl_nm){
	return arg(rp(th_med_deg,wl_nm))*180/PI;
}
double cFilm::ArgRp(double wl_nm){
	return ArgRp(th_med_deg,wl_nm);
}

double cFilm::DArgR(double th_med_deh,double wl_nm){
	double x;
	x=ArgRp(th_med_deg,wl_nm)-ArgRs(th_med_deg,wl_nm);
	while(x<=-180 || x>180){
		if(x<=-180) x+=360;
		if(x>180)   x-=360;
	}
	return x;
}
double cFilm::DArgR(double wl_nm){
	return DArgR(th_med_deg,wl_nm);
}

// ArgTs(),ArgTp()は0次透過点Bではなく，入射点Aから出射面に下した垂線の足A0のAを基準とした位相を返す．
// 光線追跡は薄膜を考慮しておらず，知りたいのはあくまで光線追跡計算による光線の経路，
// すなわち，AとA0を通る光線に沿った位相であるからA0の位相を見るのは妥当である．
// また，このままでは光学系の面間隔が薄膜の厚さ分だけ長くなることになるから，これを補正するという考えもある．
// （例えば "渋谷眞人，大木裕史 回折と結像の光学（朝倉書店)”の付録）
// しかし，薄膜の厚さが媒質を食うか，基板を食うかによって補正式が異なり，さらに，
// 通常の光学系では薄膜の厚さ分程度の面間隔のずれは製造上制御できていないこともあり，補正はしない．
//
//  【例1】
//   A0の位相はAの過去の位相であるから“遅れ”ている．
//  【例2】
//   単層膜での0次透過光(光路差 ndcosθ)の位相は，
//      波長が長くなると遅れが少なくなり，進む.
//      入射角が大きくなると光路差が短くなって遅れが少なくなり，進む．
//  【例3】
//   収束光に対する平行平面板は補正過剰の球面収差を持つから周縁光線の波面は参照球から遅れている．
//   したがって，周縁光線の位相は主光線と参照球交点の過去の位相であり，遅れている．
//   これは，例2の“入射角が大きくなると光路差が短くなって遅れが少なくなり，進む．”に反するように見えるが，
//   光線経路A-Bに沿って計算するか，A-A0に沿って計算するかの違いによっている．
//
// <追記: 2013.04.23>
//   AとA0が離れている，すなわち膜厚が厚いとArgTs(),ArgTp()の波長に対する変化が激しすぎるため，
//   主波長wl0での位相膜厚，すなわち膜物質に分散がないとした位相膜厚，による位相を減じる．
//   この位相膜厚は波数に比例する量であるのでGDDなど本質には影響しないと考えられる．
//   これにより例えば，OCTで分散補償板，眼媒質など厚い“膜”の分散を表すのに使いやすくなる．
double cFilm::ArgTs(double th_med_deg,double wl_nm){
	double x;
	
	x=( arg(ts(th_med_deg,wl_nm))+phase_thickness(th_med_deg,wl_nm) )*180/PI;
	while(x<=-180 || x>180){
		if(x<=-180) x+=360;
		if(x>180)   x-=360;
	}
	return x;
}
double cFilm::ArgTs(double wl_nm){
	return ArgTs(th_med_deg,wl_nm);
}

double cFilm::ArgTp(double th_med_deg,double wl_nm){
	double x;

	x=( arg(tp(th_med_deg,wl_nm))+phase_thickness(th_med_deg,wl_nm) )*180/PI;
	while(x<=-180 || x>180){
		if(x<=-180) x+=360;
		if(x>180)   x-=360;
	}
	return x;
}
double cFilm::ArgTp(double wl_nm){
	return ArgTp(th_med_deg,wl_nm);
}

double cFilm::DArgT(double th_med_deg,double wl_nm){
	double x;
	x=ArgTp(th_med_deg,wl_nm)-ArgTs(th_med_deg,wl_nm);
	while(x<=-180 || x>180){
		if(x<=-180) x+=360;
		if(x>180)   x-=360;
	}
	return x;
}
double cFilm::DArgT(double wl_nm){
	return DArgT(th_med_deg,wl_nm);
}

double cFilm::Tsub(double th_deg_med,double wl_nm){
	double l;
	if(SubstrateAbsorptionEnable==0) return 1;
	if(SubThickness==0) return 1;
	l=SubThickness/Re(costh(0,th_med_deg,wl_nm));
	return cMaterial::ColorGlassInnerTransmittance(mname[0],wl_nm,l);
}
double cFilm::Tsub(double wl_nm){
	return Tsub(th_med_deg,wl_nm);
}

double cFilm::GD(double th_med_deg,std::string AboutWhat,double wl_nm,double dwl_nm){
	// 群遅延 dφ/dω を fs(femto second)単位で返す．ω=2π/T=2πc/λ(c,λは真空中の光速度，波長)
	const double c=299.7925;   // 真空光速度[nm/fs]
	double p,p1,p2,domega;

	if(dwl_nm==0) dwl_nm=0.1;
	
	if(AboutWhat=="Rs" || AboutWhat=="rs"){
		p=ArgRs(th_med_deg,wl_nm);
		p1=ArgRs(th_med_deg,wl_nm-dwl_nm);
		p2=ArgRs(th_med_deg,wl_nm+dwl_nm);

	}
	else if(AboutWhat=="Rp" || AboutWhat=="rp"){
		p=ArgRp(th_med_deg,wl_nm);
		p1=ArgRp(th_med_deg,wl_nm-dwl_nm);
		p2=ArgRp(th_med_deg,wl_nm+dwl_nm);
	}
	else if(AboutWhat=="Ts" || AboutWhat=="ts"){
		p=ArgTs(th_med_deg,wl_nm);
		p1=ArgTs(th_med_deg,wl_nm-dwl_nm);
		p2=ArgTs(th_med_deg,wl_nm+dwl_nm);
	}
	else if(AboutWhat=="Tp" || AboutWhat=="tp"){
		p=ArgTp(th_med_deg,wl_nm);
		p1=ArgTp(th_med_deg,wl_nm-dwl_nm);
		p2=ArgTp(th_med_deg,wl_nm+dwl_nm);
	}
	else{
		p=p1=p2=0;
	}

	Unwrap(p1,p,p2);
	domega=-2*PI*c*(-dwl_nm/wl_nm/wl_nm);  // dω[1/fs] =2πc(-dλ/λ^2)
	return (p2-p1)/2 *(PI/180)/domega;     // GD =dφ/dω ={φ(λ+dλ)-φ(λ)}/dω
}
double cFilm::GD(std::string AboutWhat,double wl_nm,double dwl_nm){
	return GD(this->th_med_deg,AboutWhat,wl_nm,dwl_nm);
}

double cFilm::GDD(double th_med_deg,std::string AboutWhat,double wl_nm,double dwl_nm){
	// 群遅延分散 (d/dω)(dφ/dω) を fs^2 単位で返す．
	const double c=299.7925;   // 真空光速度[nm/fs]
	double domega;

	if(dwl_nm==0) dwl_nm=0.1;

	domega=-2*PI*c*(-dwl_nm/wl_nm/wl_nm);  // dω[1/fs]
	return (GD(th_med_deg,AboutWhat,wl_nm+dwl_nm,dwl_nm)-GD(th_med_deg,AboutWhat,wl_nm-dwl_nm,dwl_nm))
		   /2/domega;
}
double cFilm::GDD(std::string AboutWhat,double wl_nm,double dwl_nm){
	return GDD(this->th_med_deg,AboutWhat,wl_nm,dwl_nm);
}

double cFilm::GDRs(double th_med_deg,double wl_nm,double dwl_nm){
	return GD(th_med_deg,"Rs",wl_nm,dwl_nm);
}
double cFilm::GDRs(double wl_nm,double dwl_nm){
	return GD("Rs",wl_nm,dwl_nm);
}
double cFilm::GDRp(double th_med_deg,double wl_nm,double dwl_nm){
	return GD(th_med_deg,"Rp",wl_nm,dwl_nm);
}
double cFilm::GDRp(double wl_nm,double dwl_nm){
	return GD("Rp",wl_nm,dwl_nm);
}
double cFilm::GDTs(double th_med_deg,double wl_nm,double dwl_nm){
	return GD(th_med_deg,"Ts",wl_nm,dwl_nm);
}
double cFilm::GDTs(double wl_nm,double dwl_nm){
	return GD("Ts",wl_nm,dwl_nm);
}
double cFilm::GDTp(double th_med_deg,double wl_nm,double dwl_nm){
	return GD(th_med_deg,"Tp",wl_nm,dwl_nm);
}
double cFilm::GDTp(double wl_nm,double dwl_nm){
	return GD("Tp",wl_nm,dwl_nm);
}

double cFilm::GDDRs(double th_med_deg,double wl_nm,double dwl_nm){
	return GDD(th_med_deg,"Rs",wl_nm,dwl_nm);
}
double cFilm::GDDRs(double wl_nm,double dwl_nm){
	return GDD("Rs",wl_nm,dwl_nm);
}
double cFilm::GDDRp(double th_med_deg,double wl_nm,double dwl_nm){
	return GDD(th_med_deg,"Rp",wl_nm,dwl_nm);
}
double cFilm::GDDRp(double wl_nm,double dwl_nm){
	return GDD("Rp",wl_nm,dwl_nm);
}
double cFilm::GDDTs(double th_med_deg,double wl_nm,double dwl_nm){
	return GDD(th_med_deg,"Ts",wl_nm,dwl_nm);
}
double cFilm::GDDTs(double wl_nm,double dwl_nm){
	return GDD("Ts",wl_nm,dwl_nm);
}
double cFilm::GDDTp(double th_med_deg,double wl_nm,double dwl_nm){
	return GDD(th_med_deg,"Tp",wl_nm,dwl_nm);
}
double cFilm::GDDTp(double wl_nm,double dwl_nm){
	return GDD("Tp",wl_nm,dwl_nm);
}

std::string cFilm::DispersionTable(double wl_start,double wl_end,int wl_points){
	// 位相の分散の表を作成する．
	// 表は波数（1/wl_end～1/wl_start)で等間隔に作成する．位相はdegで表す．
	// また，定数成分，波数に対する1次式成分は除去したものを出力する．
	// これにより，分散のOCTへの影響を評価できる．
	cFitting X[5];
	std::string s;
	int i,j;
	double k0,dk;
	double *wl, *ts,*ts0,*tp,*tp0,*rs,*rs0,*rp,*rp0;
	char buf[1000];

	wl =new double [wl_points+1];
	ts =new double [wl_points+1];  // 分散
	ts0=new double [wl_points+1];  // 一次式成分を除去した分散
	tp =new double [wl_points+1];
	tp0=new double [wl_points+1];
	rs =new double [wl_points+1];
	rs0=new double [wl_points+1];
	rp =new double [wl_points+1];
	rp0=new double [wl_points+1];

	k0=1/wl_end;
	dk=((1/wl_start)-(1/wl_end))/(wl_points-1);
	
	for(i=1; i<=wl_points; i++){
		wl[i]=1/(k0+dk*(i-1));
		ts[i]=ArgTs(wl[i]);
		tp[i]=ArgTp(wl[i]);
		rs[i]=ArgRs(wl[i]);
		rp[i]=ArgRp(wl[i]);
	}

	Unwrap(ts+1,wl_points);  // ArgTs()は-180～180°なので位相接続が必要
	Unwrap(tp+1,wl_points);
	Unwrap(rs+1,wl_points);
	Unwrap(rp+1,wl_points);

	for(i=1; i<=wl_points; i++){
		X[1].AddData(wl[i],0,ts[i]);
		X[2].AddData(wl[i],0,tp[i]);
		X[3].AddData(wl[i],0,rs[i]);
		X[4].AddData(wl[i],0,rp[i]);
	}

	for(j=1; j<=4; j++){
		X[j].SetNumberOfTerms(2);
		X[j].dimensionSet(1,0,0);
		X[j].dimensionSet(2,1,0);
		X[j].RemoveTerms();
	}

	for(i=1; i<=wl_points; i++){
		ts[i]=X[1].ydataGet(i);
		tp[i]=X[2].ydataGet(i);
		rs[i]=X[3].ydataGet(i);
		rp[i]=X[4].ydataGet(i);
	}

	sprintf(buf,"  i\t   wl   \t    k   \t  ArgTs \t  ArgTp \t  ArgRs \t  ArgRp \n");
	s+=buf;
	for(i=1; i<=wl_points; i++){
		sprintf(buf,"%3d\t%8.3f\t%8.2f\t%8.5f\t%8.5f\t%8.5f\t%8.5f\n",
		        i,wl[i],1000000/wl[i]*2*PI,ts[i],tp[i],rs[i],rp[i]);
		s+=buf;
	}

	delete [] wl;
	delete [] ts; delete [] ts0;
	delete [] tp; delete [] tp0;
	delete [] rs; delete [] rs0;
	delete [] rp; delete [] rp0;

	return s;
}

std::string cFilm::FileName(){
	return remove_path(filename);
}

std::string cFilm::FilmData(double monitor_wl_nm){
	int i, count;
	std::string s;
	char buf[1000];
	std::string s1,s2,s3,s4;
	double t1,t2,t3,t4;
	double wl;

	wl= monitor_wl_nm==0 ? wl0 : monitor_wl_nm;
	
	sprintf(buf, "                n       k         d      nd    nd/(λ/4)\n");	
	//           "xxx ????????? ##.### ###.### #####.#######.#####.## 
	s+=buf;
	sprintf(buf, "sub %-9s %6.3f %7.3f", 
		    mname[0].c_str(),Re(Index(mname[0],wl)),Im(Index(mname[0],wl)) );
	s+=buf;
	if( SubstrateAbsorptionEnable && SubThickness!=0 ){
		sprintf(buf, "   (t=%.2f)\n", SubThickness ); 
	}
	else{
		sprintf(buf, "\n");
	}
	s+=buf;
	for(i=1; i<=k; i++){
		count=0;
		while(true){
			s1=get_mname(i);      t1=get_t(i);
			s2=get_mname(i+1);    t2=get_t(i+1);
			s3=get_mname(i+2);    t3=get_t(i+2);
			s4=get_mname(i+3);    t4=get_t(i+3);
			if( s1!=s3 || s2!=s4 || t1!=t3 || t2!=t4 ) break;
			i+=2;
			count++;
		}
		if(count==0){
			sprintf(buf, "    %-9s %6.3f %7.3f %8.2f%8.2f%7.3f\n", 
			        mname[i].c_str(),Re(Index(mname[i],wl)),Im(Index(mname[i],wl)),
			        d(i),nd(i,wl),nd(i,wl)/wl*4 );
			s+=buf;
		}
		else{
			sprintf(buf, "    %-9s %6.3f %7.3f %8.2f%8.2f%7.3f  x%d\n", 
			        mname[i].c_str(),Re(Index(mname[i],wl)),Im(Index(mname[i],wl)),
			        d(i),nd(i,wl),nd(i,wl)/wl*4,count+1 );
			s+=buf;
			i++;
			sprintf(buf, "    %-9s %6.3f %7.3f %8.2f%8.2f%7.3f  x%d\n",
			        mname[i].c_str(),Re(Index(mname[i],wl)),Im(Index(mname[i],wl)),
			        d(i),nd(i,wl),nd(i,wl)/wl*4,count+1 );
			s+=buf;
		}
	}
	sprintf(buf, "med %-9s %6.3f %7.3f\n", 
	        mname[k+1].c_str(),Re(Index(mname[k+1],wl)),Im(Index(mname[k+1],wl)) );
	s+=buf;
	s+='\n';
	sprintf(buf, "λo=%gnm\n", wl );
	s+=buf;
	sprintf(buf, "total layers=%d\n", k);
	s+=buf;
	sprintf(buf, "Σnd=%gnm Σd=%gnm\n", TotalNd(),TotalD());
	s+=buf;

	return s;
}

std::string cFilm::FilmDataZEMAX(){
	// ZEMAXの形式での膜構成データ
	int i;
	std::string s;
	char buf[100];

	for(i=k; i>=1; i--){
		sprintf(buf,"%-9s %9.6f 1\n",mname[i].c_str(),d(i)/1000);
		s+=buf;
	}
	return s;
}

std::ostream& operator<<(std::ostream& to,cFilm& x){
	int file_ver=102;
	int i;
	to<<file_ver<<std::endl;
	to<<x.k<<std::endl;
	to<<x.wl0<<std::endl;
	to<<x.tMode<<std::endl;
	to<<x.th_med_deg<<std::endl;
	for(i=0; i<=x.k+1; i++){
		to<<blank_filled(x.mname[i])<<' '<<x.t[i]<<' '<<x.tVariable[i]<<std::endl;
	}
	to<<x.SubThickness<<std::endl;
	return to;
}

std::istream& operator>>(std::istream& from,cFilm& x){
	// <fstream>
	int i,file_ver;
	from>>file_ver;
	switch(file_ver){
	case 102:
		// tVariableを追加(20120822)
		from>>i;
		x=cFilm(i);
		x.k=i;
		from>>x.wl0;
		from>>x.tMode;
		from>>x.th_med_deg;
		for(i=0; i<=x.k+1; i++){
			from>>x.mname[i]; x.mname[i]=inv_blank_filled(x.mname[i]);
			from>>x.t[i]>>x.tVariable[i];
		}
		from>>x.SubThickness;
		return from;
		break;
	case 101:
		// tModeを追加(20100625)
		from>>i;
		x=cFilm(i);
		x.k=i;
		from>>x.wl0;
		from>>x.tMode;
		from>>x.th_med_deg;
		for(i=0; i<=x.k+1; i++){
			from>>x.mname[i]; x.mname[i]=inv_blank_filled(x.mname[i]);
			from>>x.t[i];
		}
		from>>x.SubThickness;
		return from;
		break;
	case 100:
		from>>i;
		x=cFilm(i);
		x.k=i;
		from>>x.wl0;
		from>>x.th_med_deg;
		for(i=0; i<=x.k+1; i++){
			from>>x.mname[i]; x.mname[i]=inv_blank_filled(x.mname[i]);
			from>>x.t[i];
		}
		from>>x.SubThickness;
		return from;
		break;
	default:
		// versionのない旧形式に対応する．
		// 厳密には層数=100～現在のバージョン数 のとき誤動作となるが多分そのようなデータはないとする．
		i=file_ver;
		x=cFilm(i);
		x.k=i;
		from>>x.wl0;
		from>>x.th_med_deg;
		for(i=0; i<=x.k+1; i++){
			from>>x.mname[i];
			from>>x.t[i];
		}
		from>>x.SubThickness;
		return from;
		break;
	}
}

int cFilm::open(std::string filename){
	std::ifstream from(filename.c_str());
	if(from){
		from >> *this;
		this->filename=filename;
		return 1;
	}
	else {
		return 0;
	}
}

int cFilm::save(std::string filename){
	std::ofstream to(filename.c_str());
	if(to){
		to << *this;
		this->filename=filename;
		return 1;
	}
	else{
		return 0;
	}
}

void cFilm::reverse(double wl_nm){
	int i;

	for(i=0; i<=k+1; i++){
		if(k+1-i>i){
			std::swap(mname[i],mname[k+1-i]);
			std::swap(t[i],t[k+1-i]);
		}
	}
	SubThickness=0;

	// calc th_med
	complex Nsub=Index(mname[0],wl_nm);
	complex Nmed=Index(mname[k+1],wl_nm);
	th_med_deg=asin( sin(th_med_deg*PI/180)*Re(Nsub)/Re(Nmed) )*180/PI;
}

void cFilm::reverse(){
	reverse(this->wl0);
}

void cFilm::AdjustThickness(double multiplier){
	int i;
	if(multiplier==0) return;
	for(i=1; i<=k; i++) t[i]*=multiplier;
}

void cFilm::AdjustForIncidentAngle(double new_th_med_deg,double wl_nm){
	int i;
	for(i=1; i<=k; i++){
		t[i] *= Re( costh(i,th_med_deg,wl_nm) )/Re( costh(i,new_th_med_deg,wl_nm) );
	}
	th_med_deg=new_th_med_deg;
}

void cFilm::ndPerturbe(){
	int i;
	double ratio;

	for(i=1; i<=k; i++){
		ratio=(wl0/4)*ndTol[i] /nd(i);    // 変化量の膜厚に対する割合
		t[i]+=Random( t[i]*ratio ,0 );
		// nd を (wl0/4) * ndTol[i] だけ変化させる．
		// (薄い膜ほど，膜厚に対する変動量の割合を大きくする）
		// nd * ndTol[i] よりも膜厚制御誤差を反映すると考えられる．
	}
}

cFilm cFilm::ndPerturbed(){
	cFilm x=*this;
	x.ndPerturbe();
	return x;
}

void cFilm::dPerturbe(){
	int i;
	double ratio;

	for(i=1; i<=k; i++){
		ratio=dTol[i]/d(i);              // 変化量の膜厚に対する割合
		t[i]+=Random( t[i]*ratio ,0 );
		// d を dTol[i](nm) だけ変化させる．
	}
}

cFilm cFilm::dPerturbed(){
	cFilm x=*this;
	x.dPerturbe();
	return x;
}

double cFilm::xOfXYZ(std::string AboutWhat,std::string IlluminantName){
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.x(IlluminantName,0);
}

double cFilm::yOfXYZ(std::string AboutWhat,std::string IlluminantName){
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.y(IlluminantName,0);
}

double cFilm::uOfUCS(std::string AboutWhat,std::string IlluminantName){
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.u(IlluminantName,0);
}

double cFilm::vOfUCS(std::string AboutWhat,std::string IlluminantName){
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.v(IlluminantName,0);
}

double cFilm::Tcp(std::string AboutWhat,std::string IlluminantName){
	// 相関色温度 (JIS Z 8725:2015 による)
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.Tcp(IlluminantName,0);
}

double cFilm::duv(std::string AboutWhat,std::string IlluminantName){
	// 黒体からの偏差 (JIS Z 8725:2015 による)
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.duv(IlluminantName,0);
}

double cFilm::DominantWavelength(std::string AboutWhat,std::string IlluminantName){
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.XYToDominantWavelength(spectrum.x(IlluminantName,0),spectrum.y(IlluminantName,0));
}

double cFilm::ExcitationPurity(std::string AboutWhat,std::string IlluminantName){
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.ExcitationPurity(spectrum.x(IlluminantName,0),spectrum.y(IlluminantName,0));
}

double cFilm::LuminousFluxRatio(std::string AboutWhat,std::string IlluminantName){
	make_spectrum(AboutWhat,380,780,10);
	return spectrum.LuminousFluxRatio(IlluminantName,0);
}

double cFilm::Average(std::string AboutWhat,double wl1,double wl2,double dwl){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.Average(wl1,wl2);
}

double cFilm::Max(std::string AboutWhat,double wl1,double wl2,double dwl){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.Max(wl1,wl2);
}

double cFilm::Min(std::string AboutWhat,double wl1,double wl2,double dwl){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.Min(wl1,wl2);
}

double cFilm::FWHM(std::string AboutWhat,double wl1,double wl2,double dwl){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.FWHM();
}

double cFilm::FWHMCenter(std::string AboutWhat,double wl1,double wl2,double dwl){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.FWHMCenter();
}

double cFilm::FW(std::string AboutWhat,double wl1,double wl2,double dwl,double val){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.FW(val);
}

double cFilm::FWCenter(std::string AboutWhat,double wl1,double wl2,double dwl,double val){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.FWCenter(val);
}

double cFilm::TransitionPoint(std::string AboutWhat,double wl1,double wl2,double dwl,double val){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.TransitionPoint(val);
}

double cFilm::TransitionInterval(std::string AboutWhat,double wl1,double wl2,double dwl,
                                 double val1,double val2){
	double x1,x2;
	x1=TransitionPoint(AboutWhat,wl1,wl2,dwl,val1);
	x2=TransitionPoint(AboutWhat,wl1,wl2,dwl,val2);
	return x1==0 || x2==0 ? 0 : fabs(x1-x2);
}

double cFilm::Integral(std::string AboutWhat,double wl1,double wl2,double dwl){
	make_xylist(AboutWhat,wl1,wl2,dwl);
	return xylist.Integral(wl1,wl2);
}

double cFilm::FiniteSlit(std::string AboutWhat,double wl,double dwl,int n){
	// 波長wlを中心とするdwlの幅の中のn個の波長での平均値を返す．
	// 例えば厚さがmm程度のガラス板の透過率は激しく振動するが，
	// 本関数により分光光度計のスリット幅をdwl(nm)としたときの測定値をシミュレーションする．
	// （ただしnを相当大きくしないとなめらかにならない）
	int i;
	double sum;

	make_xylist(AboutWhat,wl-dwl/2+dwl/n/2,wl+dwl/2,dwl/n);
	sum=0;
	for(i=1; i<=n; i++){
		sum+=xylist[i].y;
	}
	return sum/n;
}

double cFilm::tDerivative(std::string command,int i,double dt){
	// をdt変化させて微分係数df/dtを計算する
	if(1<=i && i<=k){
		return Derivative(command,&t[i],dt);
	}
	else return 0;
}

double cFilm::optimize(std::string command,double dt_ratio,double rho) {
	const int M=500;        // 評価関数の最大数
	int i,j,m,n, ii;  
	cLeneq E;  
	double MeritFunction, pf;
	std::string sign, s;
	std::string *com;
	double *target,*weight;
	double val,tol;
	double dt;
	cFilm buf_ini,buf_min;
	double ini,min;
	
	if(dt_ratio==0) return 0;

	n=0;
	for( j=1; j<=k; j++ ) if( tVariable[j] ) n++;
	E.SetNumberOfVar(n);
	
	E.SetNumberOfEq(M);
	com=new std::string [M+1];
	target=new double [M+1];
	weight=new double [M+1];

	push();  // 変更を伴う場合もあるため保存する
	m=sentences(command);

	for(i=1; i<=m; ++i){
		interpret_target(sentence(command,i),com[i],sign,tol,target[i],weight[i]);

		if(weight[i]==0){
			cmd(com[i],1);
			E.SetConstrKind(i,cLeneq::DLS);
			E.SetWeight(i,0);
		}
		else if(weight[i]<0){
			// Lagrange未定乗数法(等号)，有効制約法(不等号) による拘束
			if(sign=="<"){
				E.SetConstrKind(i,cLeneq::LT);
			}
			else if(sign==">"){
				E.SetConstrKind(i,cLeneq::GT);
			}
			else{
				E.SetConstrKind(i,cLeneq::EQ);
			}
		}
		else{
			// weight[i]>0 のときはDLS法による最適化
			// (不等号にも対応．すなわち不等号はDLS法でも有効制約法でも処理できる．)
			E.SetConstrKind(i,cLeneq::DLS);
			E.SetWeight(i,weight[i]);
			
			if(sign==">"){
				val=atof(cmd(com[i],1).c_str());
				if(val>target[i]) E.SetWeight(i,0);   // 条件を満たしているときはウエイトを開放する
			}
			if(sign=="<"){
				val=atof(cmd(com[i],1).c_str());
				if(val<target[i]) E.SetWeight(i,0);
			}
			if(sign=="+-" || sign=="±"){
				val=atof(cmd(com[i],1).c_str());
				if(fabs(val-target[i])<tol) E.SetWeight(i,0);
			}

			// 注：最小二乗法計算でのウエイトはこのように解放するが，
			//     評価関数計算ではウエイトを0としない
			//    （すなわち rhoの最適化中に不等式が侵害されると評価関数は増加する）．
		}

		if( ! (E.GetConstrKind(i)==cLeneq::DLS && E.GetWeight(i)<=0) ){

			n=0;
			for( j=1; j<=k; j++ ) if( tVariable[j] ){
				n++;
				dt=t[j]*dt_ratio;
				E.SetA( i,n,tDerivative(com[i],j,dt)*dt ); 
			}

			E.SetB( i,target[i]-atof(cmd(com[i],1).c_str()) );
		}
	}

	pop();

	buf_ini=*this;
	for(ii=1; true; ii++){

		if(ii!=1){
			n=0;
			for( j=1; j<=k; j++ ) if( tVariable[j] ){
				n++;
				t[j]*=1+dt_ratio*E.GetDampedX(n,rho);
			}
		}
		
		// 評価関数の計算
		push();   // 変更を伴う場合もあるため保存する
		MeritFunction=0;
		for(i=1; i<=m; ++i) {
			if(weight[i]==0){
				cmd(com[i],1);   // 命令のみ実行
			}
			else{
				pf=atof( cmd(com[i],1).c_str() );
				MeritFunction+=E.GetWeight(i)*E.GetWeight(i)*(pf-target[i])*(pf-target[i]);
			}
		}
		pop();

		if(ii==1){
			ini=min=MeritFunction;
			buf_min=*this;
		}
		else{
			if(MeritFunction<min){
				// 最良の状態を保存する．
				min=MeritFunction;
				buf_min=*this;
				// 減衰係数を緩めて評価関数が減少するか試みる．
				*this=buf_ini;
				rho/=10;
			}
			else{
				if(min==ini && ii<10){
					// 初期状態から一度も改良されていない場合は減衰係数を強めて再試行する(回数制限有り)．
					*this=buf_ini;
					rho*=10;	
				}
				else{
					// 評価関数が増加または不変のときは最良の状態に戻し終了する
					*this=buf_min;
					break;
				}
			}
		}

	}

	delete [] com;
	delete [] target; delete [] weight;
	return min;
}

void cFilm::RemoveMinusLayers(double wl_nm){
	// 特性行列Mにおける換算光路長gを膜厚が正になるまで2πステップで増やす．
	// (すなわち，wl_nmにおけるMを一定に保ちつつ膜厚を正にする．）
	// 自動設計で膜厚が負になってしまうことへの対策として本関数を作成した（2017.03.02)．
	int i;
	double dd;
	
	for(i=1; i<=k; ++i){
		if(d(i)<0){
			dd=wl_nm/  Re(Index(mname[i],wl_nm)) / Re(costh(i,this->th_med_deg,wl_nm));
			do{
				set_d(i,d(i)+dd);
			} while(d(i)<0);
		}
	}
}

void cFilm::RemoveMinusLayers(){
	RemoveMinusLayers(wl0);
}

list<cSpectrum> cFilm::filters;
list<std::string> cFilm::filternames;

std::string cFilm::FilterDataFileLocation;

int cFilm::AddFilter(std::string filename){
	cSpectrum x;
	std::string path = FilterDataFileLocation=="" ? "" : FilterDataFileLocation+'/';
	if(x.Open(path+filename)){
		filters.AddTail(x);
		filternames.AddTail(filename);
		return 1;
	}
	else{
		return 0;
	}
}

void cFilm::RemoveFilters(){
	filters.RemoveAll();
	filternames.RemoveAll();
}

std::string cFilm::FilterNames(){
	int i;
	std::string x,s;
	for(i=1; i<=filternames.GetSize(); i++){
		filternames.GetData(x,i);
		s+=x+'\n';
	}
	return s=="" ? "no filter is installed\n" : s;
}

double cFilm::Tfilters(double wl_nm){
	int i;
	double T=1;
	for(i=1; i<=filters.GetSize(); i++){
		cSpectrum x;
		filters.GetData(x,i);
		T*=x.Data(wl_nm);
	}
	return fabs(T); 
	// 例えば分光測定実測値では0%に近い場合負になることがある．
	// fabs()がないとts,tp,rs,rpの計算でsqrt()の引数が負となってしまう．
}


std::string cFilm::scmd(std::string com,int val){
	// val=tureのとき，一部のコマンドはBasic側にてval関数で処理するために
	// 数値を表す文字列のみ返す．
	std::string s;
	char buf[1000];
	std::string s0,s1,s2,s3,s4,s5,s6;
	bool b1,b2,b3,b4,b5,b6;

	s0=arg(com,0);

	s+=scmd_general_func(com,val); if(s0!="??" && s!="") return s;
	s+=cOptics::scmd(com,val);     if(s0!="??" && s!="") return s;
	
	if(s0=="??"){
		s+="Aave AddFilter AdjustForIncidentAngle Ap ArgRp ArgRs ArgTp ArgTs As ";
		s+="Average DArgR DArgT DispersionTable ";
		s+="FilmData FilmDataZEMAX FilterDataFileLocation FilterNames FiniteSlit FW FWCenter FWHM FWHMCenter ";
		s+="GDDRp GDDRs GDDTp GDDTs GDRp GDRs GDTp GDTs ";
		s+="Integral Max Min nd Optimize ";
		s+="pop push ";
		s+="Rave RaveBothSide RemoveFilters RemoveMinusLayers ReplaceMaterial Reset Reverse ";
		s+="Rp RpBothSide Rs RsBothSide ";
		s+="Tave TaveBothSide th_med_deg Tp TpBothSide TransitionInterval TransitionPoint Ts TsBothSide ";
		s+="Tsub tToOptical tToPhysical tToOpticalQW ";
		s+="Values\n";
		return s;
	}

	{   // 0変数プロパティ /////////////////////////////
		double *p=0;
		int *ip=0;
		std::string *sp=0;
		double val;
		int ival;
		std::string sval;

		if(s0=="th_med_deg" || s0=="afocal") p=&th_med_deg;

		if(p!=0 || ip!=0 || sp!=0){
			s1=arg(com,1); b1=is_numeric(s1);

			if(s1=="?"){
				s+=s0+" (no arguments)\n";
				s+=s0+" new_value\n";
			}
			else{
				if(sp!=0){
					if(args(com)==1){
						sval=s1;
						*sp=sval;
						sprintf(buf,"%s\n",(*sp).c_str()); s=buf;
					}
					else{
						sprintf(buf,"%s\n",(*sp).c_str()); s=buf;
					}
				}
				else if( b1 ) {
					if(p!=0){
						val=atof(s1.c_str());
						*p=val;
						sprintf(buf,"%.15g\n",*p); s=buf;
					}
					else if(ip!=0){
						ival=atoi(s1.c_str());
						*ip=ival;
						sprintf(buf,"%d\n",*ip); s=buf;
					}
				}
				else{
					if(p!=0){
						sprintf(buf,"%.15g\n",*p); s=buf;
					}
					else if(ip!=0){
						sprintf(buf,"%d\n",*ip); s=buf;
					}
				}
			}

			return s;
		}
	}


	if(s0=="Aave" || s0=="aave"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Aave [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Aave(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Aave(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="AddFilter" || s0=="addfilter"){
		s1=arg(com,1);
		if(s1=="?"){
			s="AddFilter filename\n";
		}
		else{
			AddFilter(s1);
			s=FilterNames();
		}
		return s;
	}
	if(s0=="AdjustForIncidentAngle" || s0=="adjustforincidentangle"){
		double new_th_med_deg,wl_nm;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="AdjustForIncidentAngle new_th_med_deg wl\n";
		}
		else if( b1 && b2 ) {
			new_th_med_deg=atof(s1.c_str());
			wl_nm=atof(s2.c_str());
			AdjustForIncidentAngle(new_th_med_deg,wl_nm);
		}
		return s;
	}
	if(s0=="Ap" || s0=="ap"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Ap [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Ap(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Ap(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="ArgRp" || s0=="argrp"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ArgRp [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",ArgRp(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",ArgRp(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="ArgRs" || s0=="argrs"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ArgRs [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",ArgRs(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",ArgRs(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="As" || s0=="as"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="As [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",As(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",As(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="ArgTp" || s0=="argtp"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ArgTp [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",ArgTp(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",ArgTp(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="ArgTs" || s0=="argts"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ArgTs [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",ArgTs(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",ArgTs(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="Average" || s0=="average"){
		double wl1,wl2,dwl;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Average AboutWhat wl1 wl2 dwl\n";
		}
		else if( b2 && b3 && b4 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			if(val){
				sprintf(buf,"%.15g\n",Average(s1,wl1,wl2,dwl));
				s=buf;
			}
			else{
				sprintf(buf,("Average "+s1).c_str());
				s+=buf;
				sprintf(buf,"(%g-%g)=%g\n", wl1,wl2,Average(s1,wl1,wl2,dwl));
				s+=buf;
			}
		}
		return s;
	}
	if(s0=="DArgR" || s0=="dargr"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DArgR [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",DArgR(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",DArgR(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="DArgT" || s0=="dargt"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DArgT [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",DArgT(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",DArgT(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="DispersionTable" || s0=="dispersiontable"){
		double wl_start,wl_end;
		int wl_points;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="DispersionTable wl_start wl_end wl_points\n";
		}
		else if(b1 && b2 && b3){
			wl_start=atof(s1.c_str());
			wl_end=atof(s2.c_str());
			wl_points=atoi(s3.c_str());
			s=DispersionTable(wl_start,wl_end,wl_points);
		}
		return s;
	}
	if(s0=="FilmData" || s0=="filmdata"){
		double monitor_wl_nm;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="FilmData [monitor_wl=wl0]\n";
		}
		else if( b1 ) {
			monitor_wl_nm=atof(s1.c_str());
			s=FilmData(monitor_wl_nm);
		}
		else{
			s=FilmData(0);
		}
		return s;
	}
	if(s0=="FilmDataZEMAX" || s0=="filmdatazemax"){
		s1=arg(com,1);
		if(s1=="?"){
			s="FilmDataZEMAX (no arguments)\n";
		}
		else{
			s=FilmDataZEMAX();
		}
		return s;
	}
	if(s0=="FilterDataFileLocation" || s0=="filterdatafilelocation"){
		s1=arg(com,1);
		if(s1=="?"){
			s="FilterDataFileLocation foldername\n";
		}
		else{
			FilterDataFileLocation=s1;
		}
		return s;
	}
	if(s0=="FilterNames" || s0=="filternames"){
		s1=arg(com,1);
		if(s1=="?"){
			s="FilterNames (no arguments)\n";
		}
		else{
			s=FilterNames();
		}
		return s;
	}
	if(s0=="FiniteSlit" || s0=="finiteslit"){
		double wl,dwl;
		int n;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="FiniteSlit AboutWhat wl dwl n\n";
		}
		else if( b2 && b3 && b4 ) {
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			n=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",FiniteSlit(s1,wl,dwl,n)); s=buf;
		}
		return s;
	}
	if(s0=="FW" || s0=="fw"){
		double wl1,wl2,dwl,val;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s="FW AboutWhat wl1 wl2 dwl val\n";
		}
		else if( b2 && b3 && b4 && b5 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			val=atof(s5.c_str());
			sprintf(buf,"%.15g\n",FW(s1,wl1,wl2,dwl,val));
			s=buf;
		}
		return s;
	}
	if(s0=="FWCenter" || s0=="fwcenter"){
		double wl1,wl2,dwl,val;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s="FWCenter AboutWhat wl1 wl2 dwl val\n";
		}
		else if( b2 && b3 && b4 && b5 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			val=atof(s5.c_str());
			sprintf(buf,"%.15g\n",FWCenter(s1,wl1,wl2,dwl,val));
			s=buf;
		}
		return s;
	}
	if(s0=="FWHM" || s0=="fwhm"){
		double wl1,wl2,dwl;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="FWHM AboutWhat wl1 wl2 dwl\n";
		}
		else if( b2 && b3 && b4 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			sprintf(buf,"%.15g\n",FWHM(s1,wl1,wl2,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="FWHMCenter" || s0=="fwhmcenter"){
		double wl1,wl2,dwl;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="FWHMCenter AboutWhat wl1 wl2 dwl\n";
		}
		else if( b2 && b3 && b4 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			sprintf(buf,"%.15g\n",FWHMCenter(s1,wl1,wl2,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="GDDRp" || s0=="gddrp"){
		double th,wl,dwl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDDRp [th_med_deg=this->th_med_deg] wl_nm dwl_nm\n";
		}
		else if( b1 && b2 && b3 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			sprintf(buf,"%.15g\n",GDDRp(th,wl,dwl));
			s=buf;
		}
		else if( b1 && b2 ) {
			wl=atof(s1.c_str());
			dwl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",GDDRp(wl,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="GDDRs" || s0=="gddrs"){
		double th,wl,dwl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDDRs [th_med_deg=this->th_med_deg] wl_nm dwl_nm\n";
		}
		else if( b1 && b2 && b3  ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			sprintf(buf,"%.15g\n",GDDRs(th,wl,dwl));
			s=buf;
		}
		else if( b1 && b2 ) {
			wl=atof(s1.c_str());
			dwl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",GDDRs(wl,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="GDDTp" || s0=="gddtp"){
		double th,wl,dwl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDDTp [th_med_deg=this->th_med_deg] wl_nm dwl_nm\n";
		}
		else if( b1 && b2 && b3  ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			sprintf(buf,"%.15g\n",GDDTp(th,wl,dwl));
			s=buf;
		}
		else if( b1 && b2 ) {
			wl=atof(s1.c_str());
			dwl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",GDDTp(wl,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="GDDTs" || s0=="gddts"){
		double th,wl,dwl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDDTs [th_med_deg=this->th_med_deg] wl_nm dwl_nm\n";
		}
		else if( b1 && b2 && b3  ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			sprintf(buf,"%.15g\n",GDDTs(th,wl,dwl));
			s=buf;
		}
		else if( b1 && b2 ) {
			wl=atof(s1.c_str());
			dwl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",GDDTs(wl,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="GDRp" || s0=="gdrp"){
		double th,wl,dwl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDRp [th_med_deg=this->th_med_deg] wl_nm dwl_nm\n";
		}
		else if( b1 && b2 && b3  ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			sprintf(buf,"%.15g\n",GDRp(th,wl,dwl));
			s=buf;
		}
		else if( b1 && b2 ) {
			wl=atof(s1.c_str());
			dwl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",GDRp(wl,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="GDRs" || s0=="gdrs"){
		double th,wl,dwl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDRs [th_med_deg=this->th_med_deg] wl_nm dwl_nm\n";
		}
		else if( b1 && b2 && b3  ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			sprintf(buf,"%.15g\n",GDRs(th,wl,dwl));
			s=buf;
		}
		else if( b1 && b2 ) {
			wl=atof(s1.c_str());
			dwl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",GDRs(wl,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="GDTp" || s0=="gdtp"){
		double th,wl,dwl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDTp [th_med_deg=this->th_med_deg] wl_nm dwl_nm\n";
		}
		else if( b1 && b2 && b3  ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			sprintf(buf,"%.15g\n",GDTp(th,wl,dwl));
			s=buf;
		}
		else if( b1 && b2 ) {
			wl=atof(s1.c_str());
			dwl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",GDTp(wl,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="GDTs" || s0=="gdts"){
		double th,wl,dwl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDTs [th_med_deg=this->th_med_deg] wl_nm dwl_nm\n";
		}
		else if( b1 && b2 && b3  ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			dwl=atof(s3.c_str());
			sprintf(buf,"%.15g\n",GDTs(th,wl,dwl));
			s=buf;
		}
		else if( b1 && b2 ) {
			wl=atof(s1.c_str());
			dwl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",GDTs(wl,dwl));
			s=buf;
		}
		return s;
	}
	if(s0=="Integral" || s0=="integral"){
		double wl1,wl2,dwl;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Integral AboutWhat wl1 wl2 dwl\n";
		}
		else if( b2 && b3 && b4 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			if(val){
				sprintf(buf,"%.15g\n",Integral(s1,wl1,wl2,dwl)); s=buf;
			}
			else{
				sprintf(buf,("Integral "+s1).c_str());
				s=buf;
				sprintf(buf,"(%g-%g) : %g\n", wl1,wl2,Integral(s1,wl1,wl2,dwl));
				s=buf;	
			}
		}
		return s;
	}
	if(s0=="Max" || s0=="max"){
		double wl1,wl2,dwl;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Max AboutWhat wl1 wl2 dwl\n";
		}
		else if( b2 && b3 && b4 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			if(val){
				sprintf(buf,"%.15g\n",Max(s1,wl1,wl2,dwl));
				s=buf;
			}
			else{
				sprintf(buf,("Max "+s1).c_str());
				s=buf;
				sprintf(buf,"(%g-%g)=%g\n", wl1,wl2,Max(s1,wl1,wl2,dwl));
				s=buf;
			}
		}
		return s;
	}
	if(s0=="Min" || s0=="min"){
		double wl1,wl2,dwl;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Min AboutWhat wl1 wl2 dwl\n";
		}
		else if( b2 && b3 && b4 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			if(val){
				sprintf(buf,"%.15g\n",Min(s1,wl1,wl2,dwl));
				s=buf;
			}
			else{
				sprintf(buf,("Min "+s1).c_str());
				s+=buf;
				sprintf(buf,"(%g-%g)=%g\n", wl1,wl2,Min(s1,wl1,wl2,dwl));
				s+=buf;
			}
		}
		return s;
	}
	if(s0=="nd"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="nd i\n";
		}
		else if(b1) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",nd(i));
			s=buf;
		}
		return s;
	}
	if(s0=="Optimize" || s0=="optimize"){
		double dd_ratio,rho,meritfunc;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s+="optimize command dd_ratio rho\n";
			s+="(ex) optimize \"\'Average Tave 450 550 2\' [=|<|>] .7 10;";
			s+=" \'Average Tave 550 650 2\' [=|<|>] .7 10\" .01 .001\n";
			// "'Average Tave 550 650 2' [=|<|>] .7 10; 'Average Tave 550 650' [=|<|>] .7 10" .01 .001
		}
		else if( b2 && b3 ) {
			dd_ratio=atof(s2.c_str());
			rho=atof(s3.c_str());
			meritfunc=optimize(s1,dd_ratio,rho);
			if(val){
				sprintf(buf,"%f\n", meritfunc);
				s=buf;
			}
			else{
				sprintf(buf,"merit function=%g\n", meritfunc);
				s=buf;
			}
		}
		return s;
	}
	if(s0=="pop"){
		s1=arg(com,1);
		if(s1=="?"){
			s="pop (no argument)\n";
		}
		else{
			pop();
		}
	}
	if(s0=="push"){
		s1=arg(com,1);
		if(s1=="?"){
			s="push (no argument)\n";
		}
		else{
			push();
		}
	}
	if(s0=="Rave" || s0=="rave"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Rave [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Rave(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Rave(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="RaveBothSide" || s0=="ravebothside"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="RaveBothSide [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",RaveBothSide(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",RaveBothSide(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="RemoveFilters" || s0=="removefilters"){
		s1=arg(com,1);
		if(s1=="?"){
			s="RemoveFilters (no arguments)\n";
		}
		else{
			RemoveFilters();
			s=FilterNames();
		}
		return s;
	}
	if(s0=="RemoveMinusLayers" || s0=="removeminuslayers"){
		double wl_nm;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="RemoveMinusLayers [wl_nm=wl0]\n";
		}
		else if(b1){
			wl_nm=atof(s1.c_str());
			RemoveMinusLayers(wl_nm);
		}
		else{
			RemoveMinusLayers();
		}
		return s;
	}
	if(s0=="ReplaceMaterial" || s0=="replacematerial"){
		std::string from,to;
		s1=arg(com,1);
		s2=arg(com,2);
		if(s1=="?"){
			s="ReplaceMaterial from to\n";
		}
		else{
			from=s1;
			to=s2;
			ReplaceMaterial(from,to);
		}
		return s;
	}
	if(s0=="Reset" || s0=="reset"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Reset (no arguments)\n";
		}
		else{
			Reset();
		}
		return s;
	}
	if(s0=="Reverse" || s0=="reverse"){
		double wl_nm;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Reverse [wl_nm=wl0]\n";
		}
		else if( b1 ) {
			wl_nm=atof(s1.c_str());
			reverse(wl_nm);
		}
		else{
			reverse();
		}
		return s;
	}
	if(s0=="Rp" || s0=="rp"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Rp [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Rp(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Rp(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="RpBothSide" || s0=="rpbothside"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="RpBothSide [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",RpBothSide(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",RpBothSide(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="Rs" || s0=="rs"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Rs [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Rs(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Rs(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="RsBothSide" || s0=="rsbothside"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="RsBothSide [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",RsBothSide(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",RsBothSide(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="Tave" || s0=="tave"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Tave [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Tave(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Tave(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="TaveBothSide" || s0=="tavebothside"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TaveBothSide [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",TaveBothSide(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",TaveBothSide(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="Tp" || s0=="tp"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Tp [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Tp(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Tp(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="TpBothSide" || s0=="tpbothside"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TpBothSide [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",TpBothSide(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",TpBothSide(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="TransitionPoint" || s0=="transitionpoint"){
		double wl1,wl2,dwl,val;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s="TransitionPoint AboutWhat wl1 wl2 dwl val\n";
		}
		else if( b2 && b3 && b4 && b5 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			val=atof(s5.c_str());
			if(val){
				sprintf(buf,"%.15g\n",TransitionPoint(s1,wl1,wl2,dwl,val));
				s=buf;
			}
			else{
				sprintf(buf,("TransitionPoint "+s1).c_str());
				s+=buf;
				sprintf(buf,"(%g-%g @%g) : %gnm\n", wl1,wl2,val,TransitionPoint(s1,wl1,wl2,dwl,val));
				s+=buf;
			}
		}
		return s;
	}
	if(s0=="TransitionInterval" || s0=="transitioninterval"){
		double wl1,wl2,dwl,val1,val2;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="TransitionInterval AboutWhat wl1 wl2 dwl val1 val2\n";
		}
		else if( b2 && b3 && b4 && b5 && b6 ) {
			wl1=atof(s2.c_str());
			wl2=atof(s3.c_str());
			dwl=atof(s4.c_str());
			val1=atof(s5.c_str());
			val2=atof(s6.c_str());
			if(val){
				sprintf(buf,"%.15g\n",TransitionInterval(s1,wl1,wl2,dwl,val1,val2));
				s=buf;
			}
			else{
				sprintf(buf,("TransitionInterval "+s1).c_str());
				s+=buf;
				sprintf(buf,"(%g-%g %g-%g) : %gnm\n", 
					   wl1,wl2,val1,val2,TransitionInterval(s1,wl1,wl2,dwl,val1,val2));
				s+=buf;
			}
		}
		return s;
	}
	if(s0=="Ts" || s0=="ts"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Ts [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Ts(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Ts(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="TsBothSide" || s0=="tsbothside"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TsBothSide [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",TsBothSide(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",TsBothSide(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="Tsub" || s0=="tsub"){
		double th,wl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Tsub [th_med_deg=this->th_med_deg] wl_nm\n";
		}
		else if( b1 && b2 ) {
			th=atof(s1.c_str());
			wl=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Tsub(th,wl));
			s=buf;
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"%.15g\n",Tsub(wl));
			s=buf;
		}
		return s;
	}
	if(s0=="tToOptical" || s0=="ttooptical"){
		s1=arg(com,1);
		if(s1=="?"){
			s="tToOptical (no arguments)\n";
		}
		else{
			tToOptical();
		}
		return s;
	}
	if(s0=="tToOpticalQW" || s0=="ttoopticalqw"){
		s1=arg(com,1);
		if(s1=="?"){
			s="tToOpticalQW (no arguments)\n";
		}
		else{
			tToOpticalQW();
		}
		return s;
	}
	if(s0=="tToPhysical" || s0=="ttophysical"){
		s1=arg(com,1);
		if(s1=="?"){
			s="tToPhysical (no arguments)\n";
		}
		else{
			tToPhysical();
		}
		return s;
	}
	if(s0=="Values" || s0=="values"){
		double wl;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Values wl\n";
		}
		else if( b1 ) {
			wl=atof(s1.c_str());
			sprintf(buf,"               both sides\n"                   ); s+=buf;
			sprintf(buf,"Ts  =%8.6f  %8.6f\n", Ts(wl),  TsBothSide(wl)  ); s+=buf;
			sprintf(buf,"Tp  =%8.6f  %8.6f\n", Tp(wl),  TpBothSide(wl)  ); s+=buf;
			sprintf(buf,"Tave=%8.6f  %8.6f\n", Tave(wl),TaveBothSide(wl)); s+=buf;
			sprintf(buf,"Rs  =%8.6f  %8.6f\n", Rs(wl),  RsBothSide(wl)  ); s+=buf;
			sprintf(buf,"Rp  =%8.6f  %8.6f\n", Rp(wl),  RpBothSide(wl)  ); s+=buf;
			sprintf(buf,"Rave=%8.6f  %8.6f\n", Rave(wl),RaveBothSide(wl)); s+=buf;
		}
		return s;
	}

	return s="";
}



