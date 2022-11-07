#include "stdafx.h"
#include "cXYList.h"


/////////  cXYData members  ////////////////////////////////////////////////////////////////
cXYData::cXYData(){ }

cXYData::cXYData(double x,double y/*=1*/){
	this->x=x;
	this->y=y;
}

bool cXYData::operator>(const cXYData& a) const {
	return x>a.x;
}

bool cXYData::operator<(const cXYData& a) const {
	return x<a.x;
}

std::istream& operator>>(std::istream& from, cXYData& a){
	from >> a.x >> a.y;
	return from;
}

std::ostream& operator<<(std::ostream& to, const cXYData& a){
	to << a.x << ' ' << a.y;
	return to;
}


////////  cXYList members  ////////////////////////////////////////////////////////////////////
void cXYList::AddData(double x,double y) {
	AddTail( cXYData(x,y) );
}

double cXYList::x(int i) { 
	cXYData a;
	return GetData(a,i) ? a.x : 0;
}

double cXYList::y(double x) {
	// 複製liを操作する．
	// そうでないとSort()を含むため次の例で"xxx"が降順のとき
	// 意図する結果が得られない．
	//
	//  cXYlist li;
	//  cXYData a;
	//  li.open("xxx");
	//	for(int i=1; i<=li.GetSize(); i++)
	//  {
	//      std::cout << x(i) << '\t' << y(x(i)) << std::endl;
	//  }

	int i,n;
	cXYList li;
	cXYData a1,a2;
	
	li=*this;
	li.Sort();
	n=li.GetSize();
	if(n==1){
		if(li[1].x==x) return li[1].y;
	}
	else{
		for(i=1; i<=n-1; ++i){
			a1=li[i];
			a2=li[i+1];
			if(a1.x<=x && x<=a2.x){
				return a1.y+(a2.y-a1.y)*(x-a1.x)/(a2.x-a1.x);
			}
		}
	}
	// xがリストのxデータの最小値と最大値の間になければ0を返す．
	return 0;
}

double cXYList::dydx(double x) {
	int i;
	cXYList li;
	cXYData a1,a2;
	
	li=*this;
	li.Sort();
	for(i=1; i<=li.GetSize()-1; ++i){
		a1=li[i];
		a2=li[i+1];
		if(a1.x<=x && x<=a2.x){
			return (a2.y-a1.y)/(a2.x-a1.x);
		}
	}
	return 0;
}

cXYList cXYList::Subset(double x_start,double x_end) {
	int i;
	cXYData a,b;
	cXYList li;

	if(x_start>=x_end) return *this;
	Sort();
	for(i=1; i<=GetSize()-1; ++i){
		GetData(a,i);
		GetData(b,i+1);
		if( a.x<x_start && x_start<b.x){
			li.AddData(x_start,y(x_start));
			break;
		}
	}
	for(i=1; i<=GetSize(); ++i){
		GetData(a,i);
		if( x_start<=a.x && a.x<=x_end){
			li.AddTail(a);
		}
	}
	for(i=1; i<=GetSize()-1; ++i){
		GetData(a,i);
		GetData(b,i+1);
		if( a.x<x_end && x_end<b.x){
			li.AddData(x_end,y(x_end));
			break;
		}
	}
	return li;
}

void cXYList::Apply(cXYList &x) {
	int i;
	cXYData a;
	for(i=1; i<=GetSize(); ++i){
		GetData(a,i);
		a.y*=x.y(a.x);
		SetData(a,i);
	}
}

double cXYList::Integral(double x_start,double x_end) {
	int i;
	cXYList li;
	cXYData a,b;
	double s=0;

	li=this->Subset(x_start,x_end);
	for(i=1; i<=li.GetSize()-1; ++i){
		li.GetData(a,i);
		li.GetData(b,i+1);
		s+=((a.y+b.y)/2)*(b.x-a.x);
	}
	return s;
}

double cXYList::Average(double x_start,double x_end) {
	return Integral(x_start,x_end)/(x_end-x_start);
}

double cXYList::Max(double x_start,double x_end){
	int i;
	cXYList li;
	cXYData a;
	double max=-1e30;

	li=this->Subset(x_start,x_end);
	for(i=1; i<=li.GetSize(); ++i){
		a=li[i];
		if(a.y>max) max=a.y;
	}
	return max;
}

double cXYList::Min(double x_start,double x_end){
	int i;
	cXYList li;
	cXYData a;
	double min=1e30;

	li=this->Subset(x_start,x_end);
	for(i=1; i<=li.GetSize(); ++i){
		a=li[i];
		if(a.y<min) min=a.y;
	}
	return min;
}

double cXYList::Peak() {
	int i;
	cXYData a;
	double max=-1e30, peak=0;
	for(i=1; i<=GetSize(); ++i){
		GetData(a,i);
		if(a.y>max){
			max=a.y;
			peak=a.x;
		}
	}
	return peak;
}

double cXYList::PeakLocal(double x_start,double x_end) {
	cXYList li;
	li=Subset(x_start,x_end);
	return li.Peak();	
}

double cXYList::Valley() {
	int i;
	cXYData a;
	double min=1e30, valley=0;
	for(i=1; i<=GetSize(); ++i){
		GetData(a,i);
		if(a.y<min){
			min=a.y;
			valley=a.x;
		}
	}
	return valley;
}

double cXYList::ValleyLocal(double x_start,double x_end) {
	cXYList li;
	li=Subset(x_start,x_end);
	return li.Valley();	
}

void cXYList::Normalize(double new_max) {
	int i;
	cXYData a;
	double max=0;
	for(i=1; i<=GetSize(); ++i){
		GetData(a,i);
		if(a.y>max) max=a.y;
	}
	for(i=1; i<=GetSize(); ++i){
		GetData(a,i);
		if(max!=0) a.y*=new_max/max;
		SetData(a,i);
	}
}

void cXYList::NormalizeByLocal(double x_start,double x_end,double new_max) {
	int i;
	cXYData a;
	cXYList li;
	double max=0;
	li=Subset(x_start,x_end);
	for(i=1; i<=li.GetSize(); ++i){
		li.GetData(a,i);
		if(a.y>max) max=a.y;
	}
	for(i=1; i<=GetSize(); ++i){
		GetData(a,i);
		if(max!=0) a.y=new_max/max;
		SetData(a,i);
	}
}

void cXYList::xShift(double dx){
	// dxシフトさせる．
	int i;
	cXYList li;
	double x,y;

	if(dx!=0){
		li=*this;
		this->RemoveAll();
		li.Sort();
		for(i=1; i<=li.GetSize(); i++){
			x=li.x(i);
			y=li.y(x-dx);
			this->AddData(x,y);
		}
	}
}

void cXYList::Unwrap(){
	// 位相接続を行なう
	int i,n;
	double *th_deg;
	cXYData a;

	n=GetSize(); th_deg=new double[n];
	Sort();
	for(i=1; i<=n; ++i){
		GetData(a,i);
		th_deg[i-1]=a.y;
	}
	::Unwrap(th_deg,n);
	for(i=1; i<=n; ++i){
		GetData(a,i);
		a.y=th_deg[i-1];
		SetData(a,i);
	}
	delete [] th_deg;
}

void cXYList::RemoveLinearTerm(double x_start,double x_end){
	// y=Ax+B にフィッティングした後，yからAx+Bを減ずる．
	// フィッティングは x=x_start～x_end の範囲で行う．
	// x_start=x_end=0 のときは全てのデータでフィッティングする．
	cFitting F;
	cXYData a;
	int i;
	double A,B;

	F.x10VariableSet(0); F.x20VariableSet(0);
	F.SetNumberOfTerms(2); F.dimensionSet(1,0,0); F.dimensionSet(2,1,0);

	for(i=1; i<=GetSize(); ++i){
		if( (x_start==0 && x_end==0) || (x(i)>=x_start && x(i)<=x_end) ){
			GetData(a,i);
			F.AddData(a.x,0,a.y);
		}
	}

	A=F.coefficientGet(1,0,1);
	B=F.coefficientGet(0,0,0);

	for(i=1; i<=GetSize(); ++i){
		GetData(a,i);
		a.y=a.y-A*a.x-B;
		SetData(a,i);
	}
}


int cXYList::FW(double& x1,double & x2,double threshold){
	int i;
	cXYList li;
	cXYData e, e1,e2;
	double max; int max_i;
	
	if( GetSize()<=2 ) return 0;
	li=*this;
	li.Sort();
	li.GetData(e,1); max=e.y; max_i=1;
	for(i=2; i<li.GetSize(); i++) {
		li.GetData(e,i);
		if( e.y>max ) { max=e.y; max_i=i; }
	}
/*
	// 外側から探す
	for(i=1; i<=max_i-1; ++i) {
		li.GetData(e1,i);
		li.GetData(e2,i+1);
		if( e1.y<=max*threshold && max*threshold<=e2.y ) {
			x1=e1.x + ( max*threshold-e1.y )/( e2.y-e1.y )*( e2.x-e1.x );
			for(i=li.GetSize()-1; i>=max_i; --i) {
				li.GetData(e1,i);
				li.GetData(e2,i+1);
				if( e1.y>=max*threshold && max*threshold>=e2.y ) {
					x2=e1.x + ( max*threshold-e1.y )/( e2.y-e1.y )*( e2.x-e1.x );
					return 1;
				}
			}
			return 0;
		}
	}
	return 0;
*/
	// 内側から探す
	for(i=max_i-1; i>=1; --i) {
		li.GetData(e1,i);
		li.GetData(e2,i+1);
		if( e1.y<=max*threshold && max*threshold<=e2.y ) {
			x1=e1.x + ( max*threshold-e1.y )/( e2.y-e1.y )*( e2.x-e1.x );
			for(i=max_i; i<=li.GetSize()-1; ++i) {
				li.GetData(e1,i);
				li.GetData(e2,i+1);
				if( e1.y>=max*threshold && max*threshold>=e2.y ) {
					x2=e1.x + ( max*threshold-e1.y )/( e2.y-e1.y )*( e2.x-e1.x );
					return 1;
				}
			}
			return 0;
		}
	}
	return 0;
}

double cXYList::FW(double threshold){
	double x1,x2;
	if(FW(x1,x2,threshold)) return x2-x1; else return 0;
}

double cXYList::FWStart(double threshold){
	double x1,x2;
	if(FW(x1,x2,threshold)) return x1; else return 0;	
}

double cXYList::FWEnd(double threshold){
	double x1,x2;
	if(FW(x1,x2,threshold)) return x2; else return 0;
}

double cXYList::FWCenter(double threshold){
	double x1,x2;
	if(FW(x1,x2,threshold)) return (x1+x2)/2; else return 0;
}

double cXYList::FWLocal(double x_start,double x_end,double threshold){
	cXYList li;
	li=Subset(x_start,x_end);
	return li.FW(threshold);
}

double cXYList::FWStartLocal(double x_start,double x_end,double threshold){
	cXYList li;
	li=Subset(x_start,x_end);
	return li.FWStart(threshold);
}

double cXYList::FWEndLocal(double x_start,double x_end,double threshold){
	cXYList li;
	li=Subset(x_start,x_end);
	return li.FWEnd(threshold);
}

double cXYList::FWCenterLocal(double x_start,double x_end,double threshold){
	cXYList li;
	li=Subset(x_start,x_end);
	return li.FWCenter(threshold);
}

double cXYList::FWHM(){
	return FW(0.5);
}

double cXYList::FWHMCenter(){
	return FWCenter(0.5);
}

double cXYList::FWHMLocal(double x_start,double x_end){
	return FWLocal(x_start,x_end,0.5);
}

double cXYList::TransitionPoint(double val){
	// y が val になる (最小の)x を直線補間で求める．
	int i;
	cXYList li;
	cXYData e1,e2;
	if( GetSize()<=1 ) return 0;
	li=*this;
	li.Sort();
	for(i=1; i<=li.GetSize()-1; i++){
		li.GetData(e1,i);
		li.GetData(e2,i+1);
		if( (val-e1.y)*(val-e2.y)<=0 ){
			return e2.y==e1.y ? 0 : e1.x+(val-e1.y)/(e2.y-e1.y)*(e2.x-e1.x);
		}
	}
	return 0;
}

double cXYList::TransitionPointLocal(double x_start,double x_end,double val){
	cXYList li;
	li=Subset(x_start,x_end);
	return li.TransitionPoint(val);
}

double cXYList::TransitionInterval(double val1,double val2){
	double x1,x2;
	x1=TransitionPoint(val1);
	x2=TransitionPoint(val2);
	return fabs(x1-x2);
}

double cXYList::TransitionIntervalLocal(double x_start,double x_end,double val1,double val2){
	double x1,x2;
	x1=TransitionPointLocal(x_start,x_end,val1);
	x2=TransitionPointLocal(x_start,x_end,val2);
	return fabs(x1-x2);
}

double cXYList::Correl(){
	// 相関係数 Σ{(Xi-Xave)(Yi-Yave)} / √Σ{(Xi-Xave)^2} / √Σ{(Xi-Xave)^2} を計算する．
	int i;
	list<double> x,y;
	double xave,yave, dx,dy, A,B,C;

	for(i=1; i<=GetSize(); i++){
		x.AddTail((*this)[i].x);
		y.AddTail((*this)[i].y);
	}

	xave=x.Ave();
	yave=y.Ave();

	A=B=C=0;

	for(i=1; i<=GetSize(); i++){
		dx=x[i]-xave;
		dy=y[i]-yave;

		A+=dx*dy;
		B+=dx*dx;
		C+=dy*dy;
	}

	return A/sqrt(B)/sqrt(C);
}

double cXYList::R2(){
	double r;

	r=Correl();
	return r*r;
}