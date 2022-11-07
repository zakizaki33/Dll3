#include "stdafx.h"
#include "cPolarization.h"

void cPolarization::rotate(complex &Ex,complex &Ey,double xAngleDeg){
	// Ex,Eyをx座標の方向がxAngleDegの新座標に変換する．
	complex Ex0,Ey0;
	
	xAngleDeg*=PI/180;
	Ex0=Ex; Ey0=Ey;
	Ex= cos(xAngleDeg)*Ex0+sin(xAngleDeg)*Ey0;
	Ey=-sin(xAngleDeg)*Ex0+cos(xAngleDeg)*Ey0;
}


///////////////  public members  ///////////////////////////////////////////////////////////////////////

void cPolarization::JonesToStokes(double& S0,double& S1,double& S2,double& S3,complex Ex,complex Ey){
	// ジョーンズベクトル(Ex,Ey)から完全偏光S0^2=S1~2+S2~2+S3~2として
	// ストークスパラメータ(S0,S1,S2,S3)を求める．
	double Ix,Iy,I45, S3S3;
	
	Ix=sqabs(Ex);
	Iy=sqabs(Ey);
	I45=sqabs(Ex/sqrt(2)+Ey/sqrt(2));  // 45°方向に偏光子透過軸を置いたときの強度
	S0=Ix+Iy;
	S1=Ix-Iy;
	S2=2*I45-S0;
	S3S3=S0*S0-S1*S1-S2*S2;
	S3= S3S3<=0 ? 0 : sqrt(S3S3);  // 計算誤差で根号内が負になりエラーが発生するのを防ぐ．
}

void cPolarization::StokesToJones(complex &Ex,complex &Ey,double S1,double S2,double S3){
	// ストークスパラメータの完全偏光部分(S1,S2,S3)よりジョーンズベクトル(Ex,Ey)を求める．
	// Ex=ax*exp(argEx), Ey=ay*exp(argEy)
	// ax^2 = Ix = (S0+S1)/2
	// ay^2 = Iy = (S0-S1)/2
	// S2 = 2*I45-Ix-Iy
	//    = 2*| Ex/√2 + Ex/√2 |^2 -Ix-Iy = |Ex+Ey|^2 -Ix-Iy
	//    = ( Ix+Iy+2*ax*ay*cos(argEx-argEy) ) -Ix-Iy
	//    = 2*ax*ay*cos(argEx-argEy)

	double S0p, ax,ay,px;

	S0p=sqrt(S1*S1+S2*S2+S3*S3);   // 完全偏光の強度 != this->S0
	if(S0p==0){
		Ex=Ey=complex(0,0);
		return;
	}
	ax=sqrt( (S0p+S1)/2 );
	ay=sqrt( (S0p-S1)/2 );
	px=acos(S2/2/ax/ay);
	Ex=complex(ax*cos(px),ax*sin(px));
	Ey=ay;
}

matrix<double> cPolarization::MullerMatrix
    (double xAngleDeg,double Tx,double ArgExDeg,double Ty,double ArgEyEdeg){

	matrix<double> A(4,4),B(4,4);
	double t,p;
	
	t=xAngleDeg*PI/180;
	p=(ArgExDeg-ArgEyEdeg)*PI/180;

	A[1][1]=1; A[1][2]=0;         A[1][3]=0;        A[1][4]=0;
	A[2][1]=0; A[2][2]= cos(t*2); A[2][3]=sin(t*2); A[2][4]=0;
	A[3][1]=0; A[3][2]=-sin(t*2); A[3][3]=cos(t*2); A[3][4]=0;
	A[4][1]=0; A[4][2]=0;         A[4][3]=0;        A[4][4]=1;

	B[1][1]=(Tx+Ty)/2; B[1][2]=(Tx-Ty)/2; B[1][3]=0;                   B[1][4]=0;
	B[2][1]=(Tx-Ty)/2; B[2][2]=(Tx+Ty)/2; B[2][3]=0;                   B[2][4]=0;
	B[3][1]=0;         B[3][2]=0;         B[3][3]= sqrt(Tx*Ty)*cos(p); B[3][4]=sqrt(Tx*Ty)*sin(p);
	B[4][1]=0;         B[4][2]=0;         B[4][3]=-sqrt(Tx*Ty)*sin(p); B[4][4]=sqrt(Tx*Ty)*cos(p);

	return inv(A)*B*A;
}

double cPolarization::DepolarizerEffect
	(double Angle1Deg, double dPhiRangeY1Deg,double dPhiRangeX1Deg,
	 double Angle2Deg, double dPhiRangeY2Deg,double dPhiRangeX2Deg,
	 double DetectorEfficiencyXYRatio,
	 int IsGaussNotRect)
	 //	                                                                       typical value
	 // Angle1Deg       crystal plate1 crystal axis direction                  45
	 // dPhiRangeY1Deg  crystal plate1 y-direction phase difference range      360
	 // dPhiRangeX1Deg  crystal plate1 x-direction phase difference range      0
	 // Angle2Deg       crystal plate2 crystal axis direction                  -45
	 // dPhiRangeY2Deg  crystal plate2 y-direction phase difference range      -360
	 // dPhiRangeX2Deg  crystal plate2 x-direction phase difference range      0
	 // DetectorEfficiencyXYRatio  (ex.) monitor glass Rs/Rp
{
	cPolarization pol;
	list<double> li1,li2,li3;
	double Ix_ratio, phi_y;
	double x,y, step=0.1;
	double phi1,phi2;
	double dPhi2Base;
	for( dPhi2Base=0; dPhi2Base<359.999; dPhi2Base+=20 )
	{
		for( Ix_ratio=0; Ix_ratio<1.001; Ix_ratio+=0.1 ) for( phi_y=0; phi_y<359.999; phi_y+=20 )
		{
			if( dPhiRangeX1Deg==0 && dPhiRangeX2Deg==0 )
			{
				for( y=-1; y<1.001; y+=step )
				{
					pol.Initialize(0,1,Ix_ratio,0,phi_y,1);
					phi1=dPhiRangeY1Deg*y/2;
					pol.Apply(Angle1Deg,1,phi1,1,0);
					//phi2=dPhiRangeY2Deg*y/2;
					phi2=dPhiRangeY2Deg*y/2+dPhi2Base;
					pol.Apply(Angle2Deg,1,phi2,1,0);
					pol.Apply(0,1,0,DetectorEfficiencyXYRatio,0);
					if(IsGaussNotRect){
						// gaussian distribution : 1/e^2 at y=±1
						li1.AddTail(pol.Intensity() * exp(-2*y*y));
					}
					else{
						li1.AddTail(pol.Intensity() * sqrt(1.001-y*y));
					}
				}
				li2.AddTail(li1.Ave());
				li1.RemoveAll();
			}
			else
			{
				for( x=-1; x<1.001; x+=step ) for( y=-1; y<1.001; y+=step ) if( x*x+y*y<=1.001*1.001 ) 
				{
					pol.Initialize(0,1,Ix_ratio,0,phi_y,1);
					phi1=dPhiRangeX1Deg*x/2 + dPhiRangeY1Deg*y/2;
					pol.Apply(Angle1Deg,1,phi1,1,0);
					phi2=dPhiRangeX2Deg*x/2 + dPhiRangeY2Deg*y/2;
					pol.Apply(Angle2Deg,1,phi2,1,0);
					pol.Apply(0,1,0,DetectorEfficiencyXYRatio,0);
					if(IsGaussNotRect){
						li1.AddTail(pol.Intensity() * exp(-2*(x*x+y*y)));
					}
					else{
						li1.AddTail(pol.Intensity());
					}
				}
				li2.AddTail(li1.Ave());
				li1.RemoveAll();
			}
		}
		li3.AddTail( li2.Min()/li2.Max() );
		li2.RemoveAll();
	}
	return li3.Min();
}

/*
void cPolarization::EllipseShape
	(double& a,double& b,double& phi_deg,double a1,double a2,double delta_deg) {
	// Ex=a1cosθ
	// Ey=a2cos(θ+δ)
	// より，楕円長軸半径a，短軸半径b，長軸とx軸のなす角φを求める．a,bは負でもよい．
	// b>0は右回りの偏光，b<0は左回りの偏光を表す．
	// 参考書：光学の原理Ⅰ 1.4
	double delta=delta_deg*PI/180;
	double phi,phi2, chi;
	if(a1==a2 || a1==-a2){
		phi= a1*a2*cos(delta)>=0 ? PI/4 : PI*3/4;
	}
	else {
		phi2=atan( 2*a1*a2*cos(delta)/(a1*a1-a2*a2) );
		// この式だけでは，phi2が対角関係の象限のどちらにあるか区別できないため，
		// 次の処理が必要となる(図1.7参考)．
		if(phi2>=0){
			phi2= a1*a2*cos(delta)>=0 ? phi2 : phi2+PI;
		}
		else{
			phi2= a1*a2*cos(delta)>=0 ? phi2+PI : phi2+PI*2;
		}
		phi=phi2/2;
	}
	phi_deg=phi*180/PI;
	chi=asin( 2*a1*a2/(a1*a1+a2*a2)*sin(delta) )/2;         // 4.1 (27)(29)
	a=sqrt( (a1*a1+a2*a2)/(1+tan(chi)*tan(chi)) );          // 4.1 (23)(28)
	b=sqrt( fabs(a1*a1+a2*a2-a*a) );                        // 4.1 (23)
	b= chi>=0 ? b : -b;
}
*/
void cPolarization::EllipseShape
	(double& a,double& b,double& phi_deg,double a1,double a2,double delta_deg) {
	// Ex=a1cosθ
	// Ey=a2cos(θ+δ)
	// より，楕円長軸半径a，短軸半径b，長軸とx軸のなす角φを求める．a1,a2は負でもよい．
	// b>0は右回りの偏光，b<0は左回りの偏光を表す．
	// 参考書：光学の原理Ⅰ 1.4
	double delta;
	const complex i=complex(0,1);
	complex Ex,Ey, Eplus,Eminus;
	double chi;

	delta=delta_deg*PI/180;
	Ex=a1;
	Ey=a2*exp(i*delta);
	Eplus =Ex+i*Ey;
	Eminus=Ex-i*Ey;
	if( Eplus==0 || Eminus==0 ){
		phi_deg=0; // 円偏光
	}
	else{
		phi_deg=arg(Eplus/Eminus)/2 *180/PI;
	}
	chi=asin( 2*a1*a2/(a1*a1+a2*a2)*sin(delta) )/2;         // 4.1 (27)(29)
	a=sqrt( (a1*a1+a2*a2)/(1+tan(chi)*tan(chi)) );          // 4.1 (23)(28)
	b=sqrt( fabs(a1*a1+a2*a2-a*a) );                        // 4.1 (23)
	b= chi>=0 ? b : -b;
}

void cPolarization::Analyze
	(double& a1,double& a2,double& delta_deg,double a,double b,double phi_deg) {
	// 楕円軸半径a,b，半径aの軸とx軸のなす角φより，
	// Ex=a1cosθ
	// Ey=a2cos(θ+δ)
	// における，a1>0,a2>0,δを求める．|a|>|b|は要請しない．
	// ab>0は右回りの偏光，ab<0は左回りの偏光を表す．
	// 参考書：光学の原理Ⅰ 1.4
	//
	// 楕円の軸方向に，ξ軸，η軸をとると，
    //     Eξ=acosθ		(3)
    //     Eη=bsinθ		(4)
	// また，1.4(18)より，
    //     Ex=Eξcosφ-Eηsinφ		(5)
    //     Ey=Eξsinφ+Eηcosφ		(6)
	// (5)に(3)(4)を代入すると，
    //     Ex=acosθcosφ-bsinθsinφ		(7)
	// Ex=a1のときdEx/dθ=0だから，
    //     dEx/dθ=-asinθcosφ-bcosθsinφ=0		(8)
	// より，
    //     tanθ=-(bsinφ)/(acosφ)		(9)
	// だから，1つの解は，
    //     sinθ=-bsinφ/√{(acosφ)^2+(bsinφ)^2}	(10)
    //     cosθ= acosφ/√{(acosφ)^2+(bsinφ)^2}	(11)
	// (10),(11)を(7)に代入して，
    //     a1=√{(acosφ)^2+ (bsinφ)^2}		(12)
	// となる．同様にして，
    //     b1=√{(asinφ)^2+ (bcosφ)^2}		(13)

	double phi=phi_deg*PI/180;
	double x,y;
	a1=sqrt( a*a*cos(phi)*cos(phi)+b*b*sin(phi)*sin(phi) );
	a2=sqrt( a*a*sin(phi)*sin(phi)+b*b*cos(phi)*cos(phi) );
	if(a1==0 || a2==0) {
		delta_deg=0;
		return;
	}
	else{
		y=a*b/a1/a2;         // 1.4 (24)
		if( fabs(sin(2*phi))<0.001 ){
			x=tan(2*phi)*(a1*a1-a2*a2)/2/a1/a2;    
		}
		else{
			x=( a*a-b*b+cos(2*phi)*(a2+a2-a1*a1) )/( 2*a1*a2*sin(2*phi) ); // 1.4 (22)
		}
		delta_deg=arg( complex(x,y) )*180/PI;
	return;
	}
}

cPolarization::cPolarization(){
	S0=1; S1=S2=S3=0;  // 自然光
}

cPolarization::cPolarization(const complex& Ex,const complex& Ey){
	JonesToStokes(S0,S1,S2,S3,Ex,Ey);  // 完全偏光
}

cPolarization::cPolarization(const vector<complex>& E){
	JonesToStokes(S0,S1,S2,S3,E.x,E.y);  // 完全偏光
}

void cPolarization::Initialize
     (double xAngleDeg,double I,double IxRatio,double ArgExDeg,double ArgEyDeg,double PolDeg){
	// xAngleDeg ローカル座標のx軸の方向(deg)
	// I         全パワー
	// IxRatio   =Ix/I
	// ArgExDeg  x方向位相(deg)
	// ArgEyDeg  y方向位相(deg)
	// PolDeg    偏光度 = 完全偏光成分強度／強度
	const complex i=complex(0,1);
	complex Ex,Ey;

	if(PolDeg==0){
		S0=I; S1=S2=S3=0;
	}
	else{
		if(IxRatio<0) IxRatio=0;
		if(IxRatio>1) IxRatio=1;
		Ex=sqrt( I*IxRatio     )*exp(i*ArgExDeg*PI/180);
		Ey=sqrt( I*(1-IxRatio) )*exp(i*ArgEyDeg*PI/180);
		rotate(Ex,Ey,-xAngleDeg);
		JonesToStokes(S0,S1,S2,S3,Ex,Ey);
		if(0<PolDeg && PolDeg<1){
			S1*=sqrt(PolDeg);
			S2*=sqrt(PolDeg);
			S3*=sqrt(PolDeg);
		}
	}
}

void cPolarization::Apply(double xAngleDeg, double Tx,double ArgExDeg, double Ty,double ArgEyDeg){
	// 偏光子，波長板などを作用させる
	//     xAngleDeg Retarderのx軸の位置
	//     Tx        x軸方向エネルギー透過率
	//     ArgExDeg  x軸方向位相変化(deg)
	//     Ty        y軸方向エネルギー透過率
	//     ArgEyDeg  y軸方向位相変化(deg)
	matrix<double> S(4,1);

	S[1][1]=S0; S[2][1]=S1; S[3][1]=S2; S[4][1]=S3;
	S=MullerMatrix(xAngleDeg,Tx,ArgExDeg,Ty,ArgEyDeg)*S;
	S0=S[1][1]; S1=S[2][1]; S2=S[3][1]; S3=S[4][1];
}

double cPolarization::Intensity() const {
	return S0;
}

double cPolarization::PolarizedIntensity(double AngleDeg) const {
	// 透過軸方向がAngleDegの偏光子をかけた強度．自然光成分も含む．
	matrix<double> S(4,1);

	S[1][1]=S0; S[2][1]=S1; S[3][1]=S2; S[4][1]=S3;
	S=MullerMatrix(AngleDeg,1,0,0,0)*S;
	return S[1][1];
}

double cPolarization::Amplitude(double AngleDeg) const {
	complex Ex,Ey,E;

	StokesToJones(Ex,Ey,S1,S2,S3);
	E=Ex*cos(AngleDeg*PI/180)+Ey*sin(AngleDeg*PI/180);
	return abs(E);	
}

double cPolarization::EllipseLongAxisDirection() const {
	complex Ex,Ey;
	double a,b,phi;

	StokesToJones(Ex,Ey,S1,S2,S3);
	EllipseShape(a,b,phi, abs(Ex),abs(Ey),(arg(Ey)-arg(Ex))*180/PI);
	return phi;
}

double cPolarization::EllipseLongShortRatio() const {
	complex Ex,Ey;
	double a,b,phi;

	StokesToJones(Ex,Ey,S1,S2,S3);
	EllipseShape(a,b,phi, abs(Ex),abs(Ey),(arg(Ey)-arg(Ex))*180/PI);
	return a==0 ? 0 : b/a;
}

double cPolarization::EllipseLongShortIntensityRatio() const {
	return EllipseLongShortRatio()*EllipseLongShortRatio();
}

double cPolarization::Stokes(int i) const {
	// 規格化(S0=1)されたストークスパラメータ
	if(S0==0) return 0;
	switch(i){
		case 1:  return S1/S0; break;
		case 2:  return S2/S0; break;
		case 3:  return S3/S0; break;
		default: return 0;
	}
}

double cPolarization::PolDeg() const {
	// 偏光度
	return S0==0 ? 0 : sqrt(S1*S1+S2*S2+S3*S3)/S0;
}

double cPolarization::PoincareLat() const {
	// ポアンカレ球緯度(deg)
	double r;
	
	r=sqrt(S1*S1+S2*S2+S3*S3);
	if(r==0){
		return 0;
	}
	else if(S3/r>1){
		return PI;    // 計算誤差によりasin(S3/r)でエラーが発生するのを防ぐ．
	}
	else if(S3/r<-1){
		return -PI;
	}
	else{
		return asin(S3/r)*180/PI;
	}
}

double cPolarization::PoincareLng() const {
	// ポアンカレ球経度(deg)
	return S1*S1+S2*S2+S3*S3==0 ? 0 : arg( complex(S1,S2) )*180/PI;
}

