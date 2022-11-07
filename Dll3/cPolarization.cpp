#include "stdafx.h"
#include "cPolarization.h"

void cPolarization::rotate(complex &Ex,complex &Ey,double xAngleDeg){
	// Ex,Ey��x���W�̕�����xAngleDeg�̐V���W�ɕϊ�����D
	complex Ex0,Ey0;
	
	xAngleDeg*=PI/180;
	Ex0=Ex; Ey0=Ey;
	Ex= cos(xAngleDeg)*Ex0+sin(xAngleDeg)*Ey0;
	Ey=-sin(xAngleDeg)*Ex0+cos(xAngleDeg)*Ey0;
}


///////////////  public members  ///////////////////////////////////////////////////////////////////////

void cPolarization::JonesToStokes(double& S0,double& S1,double& S2,double& S3,complex Ex,complex Ey){
	// �W���[���Y�x�N�g��(Ex,Ey)���犮�S�Ό�S0^2=S1~2+S2~2+S3~2�Ƃ���
	// �X�g�[�N�X�p�����[�^(S0,S1,S2,S3)�����߂�D
	double Ix,Iy,I45, S3S3;
	
	Ix=sqabs(Ex);
	Iy=sqabs(Ey);
	I45=sqabs(Ex/sqrt(2)+Ey/sqrt(2));  // 45�������ɕΌ��q���ߎ���u�����Ƃ��̋��x
	S0=Ix+Iy;
	S1=Ix-Iy;
	S2=2*I45-S0;
	S3S3=S0*S0-S1*S1-S2*S2;
	S3= S3S3<=0 ? 0 : sqrt(S3S3);  // �v�Z�덷�ō����������ɂȂ�G���[����������̂�h���D
}

void cPolarization::StokesToJones(complex &Ex,complex &Ey,double S1,double S2,double S3){
	// �X�g�[�N�X�p�����[�^�̊��S�Ό�����(S1,S2,S3)���W���[���Y�x�N�g��(Ex,Ey)�����߂�D
	// Ex=ax*exp(argEx), Ey=ay*exp(argEy)
	// ax^2 = Ix = (S0+S1)/2
	// ay^2 = Iy = (S0-S1)/2
	// S2 = 2*I45-Ix-Iy
	//    = 2*| Ex/��2 + Ex/��2 |^2 -Ix-Iy = |Ex+Ey|^2 -Ix-Iy
	//    = ( Ix+Iy+2*ax*ay*cos(argEx-argEy) ) -Ix-Iy
	//    = 2*ax*ay*cos(argEx-argEy)

	double S0p, ax,ay,px;

	S0p=sqrt(S1*S1+S2*S2+S3*S3);   // ���S�Ό��̋��x != this->S0
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
						// gaussian distribution : 1/e^2 at y=�}1
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
	// Ex=a1cos��
	// Ey=a2cos(��+��)
	// ���C�ȉ~�������aa�C�Z�����ab�C������x���̂Ȃ��p�ӂ����߂�Da,b�͕��ł��悢�D
	// b>0�͉E���̕Ό��Cb<0�͍����̕Ό���\���D
	// �Q�l���F���w�̌����T 1.4
	double delta=delta_deg*PI/180;
	double phi,phi2, chi;
	if(a1==a2 || a1==-a2){
		phi= a1*a2*cos(delta)>=0 ? PI/4 : PI*3/4;
	}
	else {
		phi2=atan( 2*a1*a2*cos(delta)/(a1*a1-a2*a2) );
		// ���̎������ł́Cphi2���Ίp�֌W�̏ی��̂ǂ���ɂ��邩��ʂł��Ȃ����߁C
		// ���̏������K�v�ƂȂ�(�}1.7�Q�l)�D
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
	// Ex=a1cos��
	// Ey=a2cos(��+��)
	// ���C�ȉ~�������aa�C�Z�����ab�C������x���̂Ȃ��p�ӂ����߂�Da1,a2�͕��ł��悢�D
	// b>0�͉E���̕Ό��Cb<0�͍����̕Ό���\���D
	// �Q�l���F���w�̌����T 1.4
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
		phi_deg=0; // �~�Ό�
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
	// �ȉ~�����aa,b�C���aa�̎���x���̂Ȃ��p�ӂ��C
	// Ex=a1cos��
	// Ey=a2cos(��+��)
	// �ɂ�����Ca1>0,a2>0,�����߂�D|a|>|b|�͗v�����Ȃ��D
	// ab>0�͉E���̕Ό��Cab<0�͍����̕Ό���\���D
	// �Q�l���F���w�̌����T 1.4
	//
	// �ȉ~�̎������ɁC�̎��C�Ŏ����Ƃ�ƁC
    //     E��=acos��		(3)
    //     E��=bsin��		(4)
	// �܂��C1.4(18)���C
    //     Ex=E��cos��-E��sin��		(5)
    //     Ey=E��sin��+E��cos��		(6)
	// (5)��(3)(4)��������ƁC
    //     Ex=acos��cos��-bsin��sin��		(7)
	// Ex=a1�̂Ƃ�dEx/d��=0������C
    //     dEx/d��=-asin��cos��-bcos��sin��=0		(8)
	// ���C
    //     tan��=-(bsin��)/(acos��)		(9)
	// ������C1�̉��́C
    //     sin��=-bsin��/��{(acos��)^2+(bsin��)^2}	(10)
    //     cos��= acos��/��{(acos��)^2+(bsin��)^2}	(11)
	// (10),(11)��(7)�ɑ�����āC
    //     a1=��{(acos��)^2+ (bsin��)^2}		(12)
	// �ƂȂ�D���l�ɂ��āC
    //     b1=��{(asin��)^2+ (bcos��)^2}		(13)

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
	S0=1; S1=S2=S3=0;  // ���R��
}

cPolarization::cPolarization(const complex& Ex,const complex& Ey){
	JonesToStokes(S0,S1,S2,S3,Ex,Ey);  // ���S�Ό�
}

cPolarization::cPolarization(const vector<complex>& E){
	JonesToStokes(S0,S1,S2,S3,E.x,E.y);  // ���S�Ό�
}

void cPolarization::Initialize
     (double xAngleDeg,double I,double IxRatio,double ArgExDeg,double ArgEyDeg,double PolDeg){
	// xAngleDeg ���[�J�����W��x���̕���(deg)
	// I         �S�p���[
	// IxRatio   =Ix/I
	// ArgExDeg  x�����ʑ�(deg)
	// ArgEyDeg  y�����ʑ�(deg)
	// PolDeg    �Ό��x = ���S�Ό��������x�^���x
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
	// �Ό��q�C�g���Ȃǂ���p������
	//     xAngleDeg Retarder��x���̈ʒu
	//     Tx        x�������G�l���M�[���ߗ�
	//     ArgExDeg  x�������ʑ��ω�(deg)
	//     Ty        y�������G�l���M�[���ߗ�
	//     ArgEyDeg  y�������ʑ��ω�(deg)
	matrix<double> S(4,1);

	S[1][1]=S0; S[2][1]=S1; S[3][1]=S2; S[4][1]=S3;
	S=MullerMatrix(xAngleDeg,Tx,ArgExDeg,Ty,ArgEyDeg)*S;
	S0=S[1][1]; S1=S[2][1]; S2=S[3][1]; S3=S[4][1];
}

double cPolarization::Intensity() const {
	return S0;
}

double cPolarization::PolarizedIntensity(double AngleDeg) const {
	// ���ߎ�������AngleDeg�̕Ό��q�����������x�D���R���������܂ށD
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
	// �K�i��(S0=1)���ꂽ�X�g�[�N�X�p�����[�^
	if(S0==0) return 0;
	switch(i){
		case 1:  return S1/S0; break;
		case 2:  return S2/S0; break;
		case 3:  return S3/S0; break;
		default: return 0;
	}
}

double cPolarization::PolDeg() const {
	// �Ό��x
	return S0==0 ? 0 : sqrt(S1*S1+S2*S2+S3*S3)/S0;
}

double cPolarization::PoincareLat() const {
	// �|�A���J�����ܓx(deg)
	double r;
	
	r=sqrt(S1*S1+S2*S2+S3*S3);
	if(r==0){
		return 0;
	}
	else if(S3/r>1){
		return PI;    // �v�Z�덷�ɂ��asin(S3/r)�ŃG���[����������̂�h���D
	}
	else if(S3/r<-1){
		return -PI;
	}
	else{
		return asin(S3/r)*180/PI;
	}
}

double cPolarization::PoincareLng() const {
	// �|�A���J�����o�x(deg)
	return S1*S1+S2*S2+S3*S3==0 ? 0 : arg( complex(S1,S2) )*180/PI;
}

