#include "stdafx.h"
#include "cSelfoc.h"

double cSelfoc::N(vector<double> r){
	return N(no,rA,r.x,r.y);					
}

vector<double> cSelfoc::gradN(vector<double> r){
	return vector<double>(-no*rA*rA*r.x, -no*rA*rA*r.y, 0);
}

vector<double> cSelfoc::NgradN(vector<double> r){
	return N(r)*gradN(r);
}

double cSelfoc::N(double no,double rA,double x,double y){
	return no*(1-rA*rA/2*(x*x+y*y));		
}

cSelfoc::cSelfoc(){}

cSelfoc::cSelfoc(double no,double rA,double phi){
	Set_no(no);
	Set_rA(rA);
	Set_phi(phi);
}

cSelfoc::~cSelfoc(){}

double cSelfoc::Get_no(){ return no; }
void   cSelfoc::Set_no(double value){ no=value; }

double cSelfoc::Get_rA(){ return rA; }
void   cSelfoc::Set_rA(double value){ rA=value; }

double cSelfoc::Get_phi(){ return phi; }
void   cSelfoc::Set_phi(double value){ phi=value; }

double cSelfoc::Get_L(){ return L; }
void   cSelfoc::Set_L(double value){
	L=value;
	Pn=L/(2*PI/rA);
}

double cSelfoc::Get_Pn(){ return Pn; }
void   cSelfoc::Set_Pn(double value){
	Pn=value;
	L=(2*PI/rA)*Pn;
}

int cSelfoc::trace(vector<double> &r,vector<double> &Q,double ds){
	//   r(x,y,z) : 位置座標
	//   Q(X,Y,Z) : 光線方向余弦
	// を与え，光線に沿ってrからdsだけ離れた点における,r,Qを計算する．

	// 光線方程式  (d/ds){n(r)(d/ds)r(s)} = grad(n(r)) による．
	// これに置換 dt=ds/n を使うと，(d/dt)(d/dt)r = n*grad(n)
	// となる．これに初期値r,dr/dtを与えルンゲクッタ法で解く．

	vector<double> v;
	double dt;

	v=Q/abs(Q)*N(r);   // v=dr/dt=N(dr/ds)=NQ の初期値
	dt=ds/N(r);

	{
		// 2階常微分方程式 y''(x) =f(x,y(x),y'(x)) をルンゲクッタ法で解く．
		// ここでは， y''(x) = (d/dt)(d/dt)r
		//                x  = t
		//              y(x) = r
		//             y'(x) = dr/dt = v
		//           f(x,y(x),y'(x))=NgradN(r) (すなわち第2引数のみに依存)
		// 参考書：“FORTRANによる演習数値計算 洲之内治男他(サイエンス社)” 9.4章
		double dx;
		vector<double> y,y1;
		vector<double> k11,k21,k12,k22,k13,k23,k14,k24;

		y=r;
		y1=v;
		dx=dt;
		
		k11=y1;
		k21=NgradN(y);
		
		k12=y1+(dx/2)*k21;
		k22=NgradN(y+(dx/2)*k11);

		k13=y1+(dx/2)*k22;
		k23=NgradN(y+(dx/2)*k12);

		k14=y1+dx*k23;
		k24=NgradN(y+dx*k13);

		y +=(dx/6)*(k11+2.0*k12+2.0*k13+k14);
		y1+=(dx/6)*(k21+2.0*k22+2.0*k23+k24);

		r=y;
		v=y1;
	}
	
	Q=v/abs(v);

	return r.x*r.x+r.y*r.y<=phi*phi/4 ? 1 : 0;  // 戻り値でレンズ径を超えたかどうかを返す．
}

int cSelfoc::Trace(double &x,double &y,double &Qx,double &Qy,double &Qz,double N,double N1,double ds){
	// 入射面における x,y,Q(入射側媒質内)により光線追跡をし，出射面のx,y,Qで上書きする．
	//   N  : 入射側媒質屈折率
	//   N1 : 出射側媒質屈折率
	//   ds : ルンゲクッタ法のステップ(光線に沿った実距離)
	double xi,xi1,n,G;
	vector<double> r,r0,Q;
	
	Q=vector<double>(Qx,Qy,Qz);
	xi=Q.z/abs(Q);                                        // (3.20) "松居 レンズ設計法"
	r=vector<double>(x,y,0);
	n=this->N(r);
	xi1=N*n*xi/fabs(N*n*xi)*sqrt(1-(N*N/n/n)*(1-xi*xi));  // (3.21)
	G=xi1-fabs(N/n)*xi;                                   // (3.22)
	Q=fabs(N/n)*Q+vector<double>(0,0,G);                  // (3.23)

	do{
		r0=r;
		trace(r,Q,ds);
		if (r.x*r.x+r.y*r.y>phi*phi/4) return 0;
	} while( fabs(r.z)<L );

	r=r0+(r-r0)*(L-fabs(r0.z))/fabs(r.z-r0.z);
	x=r.x;
	y=r.y;

	xi=Q.z/abs(Q);                                            // (3.20)
	r=vector<double>(r.x,r.y,0);
	n=this->N(r);
	xi1=n*N1*xi/fabs(n*N1*xi)*sqrt(1-(n*n/N1/N1)*(1-xi*xi));  // (3.21)
	G=xi1-fabs(n/N1)*xi;                                      // (3.22)
	Q=fabs(n/N1)*Q+vector<double>(0,0,G);                     // (3.23)
	Qx=Q.x; Qy=Q.y; Qz=Q.z;

	return 1;
}

