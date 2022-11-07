#include "stdafx.h"
#include "cLens1.h"

// group members //////////////////////////////////////////////////////////////

group::group(int i1/*=0*/,int i2/*=0*/) { 
	this->i1=i1; this->i2=i2; 
}

bool group::operator==(const group& x) const {
	return i1==x.i1 && i2==x.i2;
}

bool group::operator>(const group& x) const {
	if     (i2==0 && x.i2==0) return i1>x.i1;
	else if(i2==0 && x.i2!=0) return i1>x.i2;
	else if(i2!=0 && x.i2==0) return i2>=x.i1;
	else {
		if(i2==x.i2) return i1>x.i1;
		else         return i1>x.i2;
	}
}

std::string group::str() const {
	char buf[100];
	
	if(i2==0) sprintf(buf,"%2d",i1);
	else      sprintf(buf,"%2d-%2d",i1,i2);
	return buf;
}

std::string group::str(int k) const {
	char buf[100];
	
	if(i2==0){
		if(i1==0)        sprintf(buf,"obj");
		else if(i1==k+1) sprintf(buf,"img");
		else             sprintf(buf,"%2d",i1);
	}
	else{
		sprintf(buf,"%2d-%2d",i1,i2);
	}
	return buf;
}


// ray_data members ///////////////////////////////////////////////////////////

ray_data::ray_data(){
	n=1; 
	color="";
	x=new double[2]; y=new double[2]; z=new double[2];
	X=new double[2]; Y=new double[2]; Z=new double[2];
	X1=new double[2]; Y1=new double[2]; Z1=new double[2];
	Ref=new int[2];
}

ray_data::ray_data(int n){
	this->n=n; 
	color="";
	x=new double[n+2]; y=new double[n+2]; z=new double[n+2];
	X=new double[n+2]; Y=new double[n+2]; Z=new double[n+2];
	X1=new double[n+2]; Y1=new double[n+2]; Z1=new double[n+2];
	Ref=new int[n+2];
}

ray_data::ray_data(const ray_data& a){
	n=a.n;
	color=a.color;
	x=new double[n+2]; y=new double[n+2]; z=new double[n+2];
	X=new double[n+2]; Y=new double[n+2]; Z=new double[n+2];
	X1=new double[n+2]; Y1=new double[n+2]; Z1=new double[n+2];
	Ref=new int[n+2];
	for(int i=0; i<=n+1; ++i){
		x[i]=a.x[i]; y[i]=a.y[i]; z[i]=a.z[i];
		X[i]=a.X[i]; Y[i]=a.Y[i]; Z[i]=a.Z[i];
		X1[i]=a.X1[i]; Y1[i]=a.Y1[i]; Z1[i]=a.Z1[i];
		Ref[i]=a.Ref[i];
	}
	GrinRay_i=a.GrinRay_i;
	GrinRay=a.GrinRay;
}

ray_data::~ray_data(){
	delete []x; delete []y; delete []z;
	delete []X; delete []Y; delete []Z;
	delete []X1; delete []Y1; delete []Z1;
	delete []Ref;
}

ray_data& ray_data::operator=(const ray_data& a){
	delete []x; delete []y; delete []z;
	delete []X; delete []Y; delete []Z;	
	delete []X1; delete []Y1; delete []Z1;
	delete []Ref;
	n=a.n;
	color=a.color;
	x=new double[n+2]; y=new double[n+2]; z=new double[n+2];
	X=new double[n+2]; Y=new double[n+2]; Z=new double[n+2];
	X1=new double[n+2]; Y1=new double[n+2]; Z1=new double[n+2];
	Ref=new int[n+2];
	for(int i=0; i<=n+1; ++i){
		x[i]=a.x[i]; y[i]=a.y[i]; z[i]=a.z[i];
		X[i]=a.X[i]; Y[i]=a.Y[i]; Z[i]=a.Z[i];
		X1[i]=a.X1[i]; Y1[i]=a.Y1[i]; Z1[i]=a.Z1[i];
		Ref[i]=a.Ref[i];
	}
	GrinRay_i=a.GrinRay_i;
	GrinRay=a.GrinRay;
	return *this;
}


// cPoint members /////////////////////////////////////////////////////////////

cPoint::cPoint(){
	p.x=p.y=0;
	wl=0;
	weight=1;
	opl=0;
}

cPoint::cPoint(double x,double y,double wl,double weight,double opl){
	p.x=x; p.y=y;
	this->wl=wl;
	this->weight=weight;
	this->opl=opl;
}

cPoint::cPoint(double x,double y,double wl,double weight){
	p.x=x; p.y=y;
	this->wl=wl;
	this->weight=weight;
	opl=0;
}

cPoint::cPoint(double x,double y,double wl){
	p.x=x; p.y=y;
	this->wl=wl;
	weight=1;
	opl=0;
}


// cSpot members //////////////////////////////////////////////////////////////

cSpot::cSpot(){
	origin=point(0,0);
}

int cSpot::GetSize() const {
	return spot.GetSize();
}

void cSpot::AddTail(const cPoint& p) {
	spot.AddTail(p);
}  

int cSpot::GetData(cPoint& p,int i) {
	int a;
	a=spot.GetData(p,i);
	p.p-=origin;
	return a;
}

void cSpot::RemoveAll() { spot.RemoveAll(); }

point cSpot::GravityCenter() {
	int i,n;
	cPoint p;
	point g;
	double nw=0;
	n=spot.GetSize();
	if(n>0){
		for(i=1; i<=n; ++i){
			spot.GetData(p,i);
			g+=p.p*p.weight;
			nw+=p.weight;
		}
	}
	return nw==0 ? point(0,0) : g/nw;
}

double cSpot::InertiaPrincipalAxis(){
	// 慣性主軸とX軸のなす角度(deg)を返す．
	point g;
	double vx,vy,vxy;
	cPoint p;
	int i,n;

	g=GravityCenter();
	n=spot.GetSize();
	vx=vy=vxy=0;

	for(i=1; i<=n; ++i){
		p=spot[i];
		vx +=(p.p.x-g.x)*(p.p.x-g.x)*p.weight;
		vy +=(p.p.y-g.y)*(p.p.y-g.y)*p.weight;
		vxy+=(p.p.y-g.y)*(p.p.x-g.x)*p.weight;
	}

	return vxy==0 ? 90 : atan((vx-vy)/vxy)*180/PI;
}

void cSpot::PrincipalWidth(double &Wx,double &Wy){
	int i,n;
	double th,sn,cs, x,y, xmax,xmin,ymax,ymin;
	vector<double> ex,ey,r;
	cPoint p;
	point g;
	
	n=spot.GetSize();
	th=InertiaPrincipalAxis(); sn=sin(th*PI/180); cs=cos(th*PI/180);
	
	ex=vector<double>(cs,sn,0);    // 主軸方向の単位ベクトル
	ey=vector<double>(-sn,cs,0);   // 主軸に垂直な単位ベクトル

	xmax=-1e30; xmin=1e30;
	ymax=-1e30; ymin=1e30;

	for(i=1; i<=n; ++i){
		p=spot[i];
		r=vector<double>(p.p.x-g.x, p.p.y-g.y, 0);   // gからpへのベクトル
		x=sProduct(r,ex);                            // rのex方向の成分
		y=sProduct(r,ey);                            // rのey方向の成分
		
		if(x>xmax) xmax=x;
		if(x<xmin) xmin=x;
		if(y>ymax) ymax=y;
		if(y<ymin) ymin=y;
	}

	Wx=xmax-xmin;
	Wy=ymax-ymin;
}

double cSpot::PrincipalWx(){
	double Wx,Wy;
	PrincipalWidth(Wx,Wy);
	return Wx;
}
double cSpot::PrincipalWy(){
	double Wx,Wy;
	PrincipalWidth(Wx,Wy);
	return Wy;
}

void cSpot::OriginToGravityCenter() {
	origin=GravityCenter();
}

void cSpot::OriginToNewCenter(double x,double y) {
	origin=point(x,y);
}

complex cSpot::OTF(double nu_y,double nu_x) {
	int i;
	cPoint p;
	double Rc,Rs,w;
	Rc=Rs=w=0;
	for(i=1; i<=GetSize(); ++i) {
		GetData(p,i);
		Rc+=cos( 2*PI*(nu_y*p.p.y+nu_x*p.p.x) )*p.weight;
		Rs+=sin( 2*PI*(nu_y*p.p.y+nu_x*p.p.x) )*p.weight;
		w+=p.weight;
	}
	return w==0 ? 0 : complex(Rc/w,Rs/w);
}

double cSpot::MTF(double nu_y,double nu_x){
	return abs( OTF(nu_y,nu_x) );
}

double cSpot::TotalIntensity() {
	cPoint p; double result=0;
	for(int i=1; i<=spot.GetSize(); ++i){
		GetData(p,i);
		result+=p.weight;
	}
	return result;
}

double cSpot::Intensity(double y,double x,double SensorPhi) {
	cPoint p; double result=0;
	for(int i=1; i<=spot.GetSize(); ++i){
		GetData(p,i);
		if( (p.p.y-y)*(p.p.y-y)+(p.p.x-x)*(p.p.x-x) <= SensorPhi*SensorPhi/4 ) result+=p.weight;
	}
	return result;
}

double cSpot::EncircledEnergy(double SensorPhi){
	cSpot buf=*this;
	
	buf.OriginToGravityCenter();
	return buf.Intensity(0,0,SensorPhi)/buf.TotalIntensity();
}

double cSpot::XIntensity(double x,double SensorFW) {
	cPoint p; double result=0;
	for(int i=1; i<=spot.GetSize(); ++i){
		GetData(p,i);
		if( fabs(p.p.x-x) <= SensorFW/2 ) result+=p.weight;
	}
	return result;
}

double cSpot::YIntensity(double y,double SensorFW) {
	cPoint p; double result=0;
	for(int i=1; i<=spot.GetSize(); ++i){
		GetData(p,i);
		if( fabs(p.p.y-y) <= SensorFW/2 ) result+=p.weight;
	}
	return result;
}

double cSpot::RmsPhi() {
	cSpot temp;
	int i; cPoint p;
	double nw=0, a=0;
	temp=*this;
	if( temp.GetSize()==0 ) return 0;
	temp.OriginToGravityCenter();
	for(i=1; i<=temp.GetSize(); ++i){
		temp.GetData(p,i);
		a+=( p.p.x*p.p.x+p.p.y*p.p.y )*p.weight;
		nw+=p.weight;
	}
	return nw==0 ? 0 : 2*sqrt(a/nw);
}

double cSpot::XYAbsMax() {
	int i; cPoint p;
	double max=0;
	for(i=1; i<=GetSize(); ++i){
		GetData(p,i);
		if( fabs(p.p.x)>max ) max=fabs(p.p.x);
		if( fabs(p.p.y)>max ) max=fabs(p.p.y);
	}
	return max;
}

double cSpot::XMax() {
	int i; cPoint p;
	double max=-1e100;
	for(i=1; i<=GetSize(); ++i){
		GetData(p,i);
		if( p.p.x > max ) max=p.p.x;
	}
	return max;
}

double cSpot::YMax() {
	int i; cPoint p;
	double max=-1e100;
	for(i=1; i<=GetSize(); ++i){
		GetData(p,i);
		if( p.p.y > max ) max=p.p.y;
	}
	return max;
}

double cSpot::XMin() {
	int i; cPoint p;
	double min=1e100;
	for(i=1; i<=GetSize(); ++i){
		GetData(p,i);
		if( p.p.x < min ) min=p.p.x;
	}
	return min;
}

double cSpot::YMin() {
	int i; cPoint p;
	double min=1e100;
	for(i=1; i<=GetSize(); ++i){
		GetData(p,i);
		if( p.p.y < min ) min=p.p.y;
	}
	return min;
}

double cSpot::XWidth() {
	int i; cPoint p;
	double min=1e100, max=-1e100;
	for(i=1; i<=GetSize(); ++i){
		GetData(p,i);
		if( p.p.x < min ) min=p.p.x;
		if( p.p.x > max ) max=p.p.x;
	}
	return max-min;
}

double cSpot::YWidth() {
	int i; cPoint p;
	double min=1e100, max=-1e100;
	for(i=1; i<=GetSize(); ++i){
		GetData(p,i);
		if( p.p.y < min ) min=p.p.y;
		if( p.p.y > max ) max=p.p.y;
	}
	return max-min;
}

double cSpot::XGravityCenter() {
	return GravityCenter().x;
}

double cSpot::YGravityCenter() {
	return GravityCenter().y;
}

void cSpot::RemoveTC(){
	// 横色収差を除く．例えば，SLOのソフトウエアでの横色収差の補正を想定している．
	// 具体的には，各波長の重心が一致するようにする．
	int i,j, spots,cols;
	point g,g0;
	list<double> wl_list;
	cSpot *buf;

	spots=GetSize();  // スポットのサイズ

	// 使われている波長のリスト wl_list を作成する．
	for(i=1; i<=spots; ++i){
		wl_list.AddTail(spot[i].wl);
	}
	wl_list.EraseDuplicates();

	cols=wl_list.GetSize();   // 波長数
	buf=new cSpot[cols+1];

	// 波長毎に分ける
	for(i=1; i<=spots; ++i){
		for(j=1; j<=cols; ++j){
			if(spot[i].wl==wl_list[j]) buf[j].AddTail(spot[i]);
		}
	}

	// 各波長の重心を全体の重心に一致させる
	g0=GravityCenter();
	for(j=1; j<=cols; ++j){
		spots=buf[j].GetSize();
		g=buf[j].GravityCenter();
		for(i=1; i<=spots; ++i){
			buf[j].spot[i].p-=g-g0;
		}
	}

	// 波長毎のcSpotを統合する
	RemoveAll();
	origin=point(0,0);
	for(j=1; j<=cols; ++j){
		spots=buf[j].GetSize();
		for(i=1; i<=spots; ++i) spot.AddTail(buf[j].spot[i]);
	}

	delete[] buf;
}
