#include "stdafx.h"
#include "cOptics.h"

cOptics::cOptics(){};

cOptics::cOptics(const cOptics &x){
	*this=x;
}

cOptics& cOptics::operator =(const cOptics& x){
	this->symbol=x.symbol;
	this->value=x.value;
	return *this;
}

std::string cOptics::MaterialName(const std::string& gname){
	return cMaterial::MaterialName(gname);
}

double cOptics::HerzEq(double Nd,double nud,const std::string& color){
	return cMaterial::HerzEq(Nd,nud,Wavelength(color));
}

complex cOptics::Index(const std::string& s,double wl_nm){
	return cMaterial::Index(s,wl_nm);
}

double cOptics::GroupIndex(const std::string& s,double wl_nm){
	return cMaterial::GroupIndex(s,wl_nm);
}

complex cOptics::IndexDerivative2(const std::string& s,double wl_nm){
	return cMaterial::IndexDerivative2(s,wl_nm);
}
double cOptics::ReIndex(const std::string& s,double wl_nm) { return Re(Index(s,wl_nm)); }
double cOptics::ImIndex(const std::string& s,double wl_nm) { return Im(Index(s,wl_nm)); }

double cOptics::AbbeNumber(const std::string& glass,
						   const std::string& color0,const std::string& color1,const std::string& color2){
	// AbbeNumber = (n0-1)/(n1-n2) 
	double n0,n1,n2;
	n0=fabs(ReIndex(glass,Wavelength(color0)));
	n1=fabs(ReIndex(glass,Wavelength(color1)));
	n2=fabs(ReIndex(glass,Wavelength(color2)));
	return n1==n2 ? 0 : (n0-1)/(n1-n2);
}

double cOptics::AbbeNumber(const std::string& glass){
	// AbbeNumber = (nd-1)/(nF-nC) 
	return AbbeNumber(glass,"d","F","C");
}

std::string cOptics::GlassCode(double nd,double nud){
	// nd�̏����_�ȉ��ƃA�b�x���ɂ��R�[�h
	// �Ⴆ��TiO2(nd=2.5)��nd=1.5�̃K���X�̃R�[�h�̍ŏ���3���͓����ƂȂ�D
	// �����������̌��w�K���X�Ɍ���΃R�[�h��nd�̑Ή��͈�ӂł���D
	char buf[100]; std::string s;

	while(nd>=1) nd-=1;
	sprintf(buf, "%03.0f", (Round(nd,-3))*1000);
	s=buf;
	sprintf(buf, "%03.0f", Round(nud,-1)*10);
	s+=buf;
	return s;
}

std::string cOptics::GlassCode(const std::string& s){
	return GlassCode(ReIndex(s,Wavelength("d")),AbbeNumber(s));
}

double cOptics::Wavelength(const std::string& color){
	// �g����nm�P�ʂŕԂ�
	double wl;

	wl=atof(color.c_str());
	
	if     (color=="i" ) wl=365.015;
	else if(color=="h" ) wl=404.656;
	else if(color=="g" ) wl=435.835;
	else if(color=="F'") wl=479.99;
	else if(color=="F" ) wl=486.13;
	else if(color=="e" ) wl=546.07;
	else if(color=="d" ) wl=587.56;
	else if(color=="D" ) wl=589.29;
	else if(color=="C'") wl=643.85;
	else if(color=="C" ) wl=656.27;
	else if(color=="r" ) wl=706.52;
	else if(color=="A'") wl=768.19;
	else if(color=="s" ) wl=852.11;
	else if(color=="t" ) wl=1013.98;
	
	return wl;
}

double cOptics::PolygonDutyFactor(double Facets,double ScanFullAngle){
	double C=Facets*ScanFullAngle/720;
	return C>1 ? 0 : C;
}

double cOptics::PolygonMinDiaNotIntensityVary
                (double Facets,double BeamDia,double ScanFullAngle,double FeedAngle,double RollOffLength)
{
	// Facets        : the number of facets
	// BeamDia       : the beam dia in a plane perpendicular to the optical axis
	// ScanFullAngle : the full angle used
	// FeedAngle     : the angle between the input beam and the center of scan
	// RollOffLength : effective length = length - RollOffLength
	double C,d1,l;
	C=Facets*ScanFullAngle/720;
	if(C>=1) return 0;
	d1=BeamDia/cos(FeedAngle/2*PI/180);
	l=(d1+RollOffLength)/(1-C);
	return l/tan( PI/Facets );
}

double cOptics::BeamMaxDiaNotIntensityVary
	            (double Facets,double PolygonDia,double ScanFullAngle,double FeedAngle,double RollOffLength)
{
	double C,l,d1;
	C=Facets*ScanFullAngle/720;
	if(C>=1) return 0;
	l=PolygonDia*tan( PI/Facets );
	d1=l*(1-C)-RollOffLength;
	return d1*cos(FeedAngle/2*PI/180);
}

double cOptics::BeamMaxDiaCenter(double Facets,double PolygonDia,double FeedAngle,double RollOffLength)
{
	double l;
	l=PolygonDia*tan( PI/Facets );
	return (l-RollOffLength)*cos(FeedAngle/2*PI/180);
}

double cOptics::PolygonEfficiency
	            (double Facets,double PolygonDia,double BeamDia,double FeedAngle,double RollOffLength,double ScanAngle)
{
	double l, th,th1,th2, r,x1,x2, R, TH1,TH2;
	th=FeedAngle/2 +ScanAngle/2;
	if( th<=0 || 90<th ) return 0;
	l=PolygonDia*tan(PI/Facets);
	l=l-RollOffLength;
	th1=th-atan(l/PolygonDia)*180/PI; 
	th2=th+atan(l/PolygonDia)*180/PI;
	r=sqrt( PolygonDia*PolygonDia/4 + l*l/4 );
	x1=r*sin(th1*PI/180) -PolygonDia/2*sin(FeedAngle/2*PI/180);
	x2=r*sin(th2*PI/180) -PolygonDia/2*sin(FeedAngle/2*PI/180);
	R=BeamDia/2;
	if(x1<=-R && R<=x2){
		return 1;
	}
	else if(x1<=-R && x2<=-R){
		return 0;
	}
	else if(x1<=-R && -R<x2 && x2<R){
		TH2=acos(x2/R);
		return ( PI -( TH2-sin(2*TH2)/2 ) )/PI;
	}
	else if(-R<x1 && x1<R && -R<x2 && x2<R){
		TH1=acos(x1/R); TH2=acos(x2/R);
		return ( ( TH1-sin(2*TH1)/2 )-( TH2-sin(2*TH2)/2 ) )/PI;
	}
	else if(-R<x1 && x1<R && R<=x2){
		TH1=acos(x1/R);
		return ( TH1-sin(2*TH1)/2 )/PI;
	}
	else {	
		return 0;
	}
}

point cOptics::edge(double phi_x,double phi_y,int is_rect,double th_rad){
	// ��x=phi_x,��y=phi_y�̋�`�܂��͑ȉ~��
	// ����	y=tan(th_rad)x �Ƃ̌�_�̈����Ԃ��D
	//
	// th_rad���}PI/2�̋߂��ł����Ȃ��悤�ł���D
	double x,y;
	double tn=tan(th_rad);
	if(phi_x==0 || phi_y==0) return point(0,0);
	if(is_rect){
		if(tn==0) return point(phi_x/2,0);
		x=fabs(phi_y/2)/tn;
		if(x>=0){
			x=  x<fabs(phi_x)/2 ? x : fabs(phi_x)/2;
		}
		else{
			x= -fabs(phi_x)/2<x ? x : -fabs(phi_x)/2;
		}
		y=tn*x;
	}
	else{
		x=sqrt( 1/( 1/(phi_x*phi_x/4) + tn*tn/(phi_y*phi_y/4) ) );
		y=tn*x;
	}
	return point(x,y);
}

std::string cOptics::Spherometer(double D,double H,double Rref/*=0*/,double Href/*=0*/){
	// ���ʌv�ŋ��ʋȗ����a�𑪒肷��D
	//     D    = ����w�b�h�̒��a
	//     H    = ����ʂɓ��Ă��Ƃ��̃_�C�����Q�[�W�̓ǂ�(�ʖʁC���ʋ��ɐ��Ƃ���j
	//     Rref = �W���ʂ̋ȗ����a�i0�̂Ƃ��́C�W���ʂɂ��␳�͍s��Ȃ��D �ʖʁC���ʋ��ɐ��Ƃ���j
	//     Href = �W���ʂɂ��Ă��Ƃ��̃_�C�����Q�[�W�̓ǂ�
	double r,r1,r2, h_err;
	double dh=0.01;  // ��F�_�C�����Q�[�W�̍ŏ��ڐ�
	char buf[1000];

	r=(H*H+D*D/4)/(H*2);
	if(Rref>0){
		h_err=Href - (D*D/4)/Rref/(1+sqrt(1-(D*D/4)/Rref/Rref));
		H-=h_err;
		r=(H*H+D*D/4)/(H*2);
	}
	
	r1= ((H-dh)*(H-dh)+D*D/4)/((H-dh)*2);
	r2= ((H+dh)*(H+dh)+D*D/4)/((H+dh)*2);
	sprintf(buf,"r=%g (%g�`%g for dh=�}%g)\n",r,r1,r2,dh);
	return buf;
}

std::string cOptics::Plate(std::string glass,std::string color,double thickness,double incident_angle,
                           double dndtemp,double alpha,std::string what){
	std::string s;
	double wl,n,t,th,th1, a,b,x, ds,dp,dsp,dave,shift,dz,shift01,shift1,dnl,dtanth1,l,nl;
	double dwldf,dtdf,dthdf,dtempdf;
	double kp,ks;
	char buf[1000];
	const double DNAIRDTEMP=-0.94e-6; // 20���ɂ������C���ܗ��̉��x�W���i�I�n���J�^���O���j
	
	wl=Wavelength(color);
	n=Re(Index(glass,wl));
	t=thickness;
	th=incident_angle*PI/180;
	x=sin(th)/n; 
	th1=atan( x/sqrt(-x*x+1) );
	a=t/cos(th1)*( cos(th-th1)-1/n );
	b=t/n/cos(th1)*( 1-(cos(th)/cos(th1))*(cos(th)/cos(th1)) );

	sprintf(buf,"�ɍ�= %s\n"                     ,glass.c_str());                         s+=buf;
	sprintf(buf,"�g��= %s\n"                     ,color.c_str());                         s+=buf;
	sprintf(buf,"�� t= %g\n"                   ,t);                                     s+=buf;
	sprintf(buf,"���ˊp th= %g��\n"              ,incident_angle);                        s+=buf;
	sprintf(buf,"  ��' th1= %.6f��\n"                     ,th1*180/PI);                   s+=buf;
	sprintf(buf,"  S�ʓ��̂� ds= %.6f\n"                  ,ds=a);                         s+=buf;
	sprintf(buf,"  P�ʓ��̂� dp= %.6f\n"                  ,dp=a+b);                       s+=buf;
	sprintf(buf,"  ��_�u��P-S dsp= %.6f\n"               ,dsp=b);                        s+=buf;
	sprintf(buf,"  SP�̂ѕ��� dave= %.6f\n"               ,dave=a+b/2);                   s+=buf;
	sprintf(buf,"  �������� ���ˁ`0������  shift= %.6f\n" ,shift=t/cos(th1)*sin(th-th1)); s+=buf;
	sprintf(buf,"  ��z ���˓_�`0�����ߓ_= %.6f\n"         ,dz=t/cos(th1)*cos(th-th1));    s+=buf;
	sprintf(buf,"  �������� 0���`1��     shift01= %.6f\n" ,shift01=2*t*tan(th1)*cos(th)); s+=buf;
	sprintf(buf,"  �������� ���ˁ`1������ shift1= %.6f\n" ,shift1=shift01-shift );        s+=buf;
	sprintf(buf,"  ���H���� 0���`1�����ߔ��� dnl=%.6f\n"  ,dnl=2*n*t*cos(th1));           s+=buf;
	sprintf(buf,"  �� x tan��' dtanth1= %.6f\n"         ,dtanth1=t*tan(th1));           s+=buf;
	sprintf(buf,"  �����H ����    l= %.6f\n"            ,l=t/cos(th1));                 s+=buf;
	sprintf(buf,"  �����H ���H�� nl= %.6f\n"            ,nl=n*t/cos(th1));              s+=buf;
	
	// �g�􉽌��w (�O��a�v�C�����o��)�h��2�͂̍Ō�̐}(�}2.34)�́C
	//  ���ˌ��C�o�ˌ��̕����P�ʃx�N�g����s,s'�Ƃ���Ƃ��C���˓_���n�_�Ƃ����ns,n's'�̏I�_�́C
	//  ���˓_�𒆐S�Ƃ��锼�an,n'�̋���ɂ���Cn's'-ns�͖ʖ@��p�ƕ��s�ł���C���Ƃ������Ă���D
	//  ������n=1�Cn'=n�Ƃ��Đ}�I�ɍl����ƁC
	//  �ψʂ�P�ʓ��̂Ƃ��C
	//      ncos(th1)��(th1) = cos��(th1)��(th1)  ...(1)
	//  �ψʂ�S�ʓ��̂Ƃ��C
	//      n��(th1) = ��(th)  ...(2)
	//  �ƂȂ�D
	// ���s�x�����Ƃ���D
	// �܂�P�ʓ��ł́C
	// 0��1�����˂̂Ȃ��p�̏ꍇ�͐}��`���čl����ƕ\�ʂɂ����āC
	//     1�� ��(th1)=2��
	// �����炱���(1)�֑�����C��(th)���ƂȂ�D
	// �܂��C0��1�����߂̂Ȃ��p�̏ꍇ�͗��ʂɂ����āC
	//     0�� ��(th1)=��
	//     1�� ��(th1)=3��
	// ������C������(1)�ɑ�����C���ꂼ��̃�(th)�̍����ƂȂ�D
	// ����S�ʓ��ł́C
	// 0��1�����˂̂Ȃ��p�̏ꍇ��
	//     1�� ��(th1)=2��cos(th1)   ( <- cos(th1)�������邱�Ƃɒ��ӁD��L�}�ōl����ƕ�����Ղ��D)
	// �����炱���(2)�֑�����C��(th)���ƂȂ�D
	// �܂��C0��1�����߂̂Ȃ��p�̏ꍇ�͗��ʂɂ����āC
	//     0�� ��(th1)=��cos(th1)
	//     1�� ��(th1)=3��cos(th1)
	// ������C������(2)�ɑ�����C���ꂼ��̃�(th)�̍����ƂȂ�D
	// ���[�U���g�������s���ʔ̕��s�x����(���ˌ��������ɓ��e�������𑪒�)�ɗL�p�D
	//�i���ۂ̑���ł͌덷�v���ƂȂ�"�������� 0���`1��"�ɒ��ӁD)
	sprintf(buf,"  �ꎟ�ߎ�\n");                                                        s+=buf;
	sprintf(buf,"    0��1�����˂܂��͓��߂̂Ȃ��p= k x ���s�x(�\���ʖ@���̂Ȃ��p)\n");  s+=buf;
	sprintf(buf,"      ���s�x�̕ψʂ�P�ʓ� kp= %g\n", kp=2*n*cos(th1)/cos(th));         s+=buf;
	sprintf(buf,"      ���s�x�̕ψʂ�S�ʓ� ks= %g\n", ks=2*n*cos(th1));                 s+=buf;
	
	// fringes = 2*n*thickness*cos(th1)/��o ���D
	dwldf=wl*wl/2/n/(t*1000000)/cos(th1);
	dtdf=(wl/1000000)/2/n/cos(th1);
	dthdf= th1==0 ? 0 : (wl/1000000)*cos(th1)/2/t/sin(th1)/cos(th)*180/PI;
	if(dndtemp==0 || alpha==0){
		dtempdf=0;
	}
	else{
		// dndtemp = ���΋��ܗ��̉��x�W��
		// alpha   = ���c���W��
		dtempdf= 1/( 2*t*cos(th1)/(wl/1000000) )/( dndtemp+n*DNAIRDTEMP+n*alpha );
	}
	
	sprintf(buf,"  ����[nm]/fringe  dwldf= %g\n" ,dwldf  ); s+=buf;
	sprintf(buf,"  ��t [mm]/fringe   dtdf= %g\n" ,dtdf   ); s+=buf;
	sprintf(buf,"  ��th[deg]/fringe dthdf= %g\n" ,dthdf  ); s+=buf;
	sprintf(buf,"  ��T [K]/firnge dtempdf= %g\n" ,dtempdf); s+=buf;
	
	if(what=="th1"     ) { sprintf(buf,"%g\n",th1*180/PI); s=buf; }
	if(what=="ds"      ) { sprintf(buf,"%g\n",ds);         s=buf; }
	if(what=="dp"      ) { sprintf(buf,"%g\n",dp);         s=buf; }
	if(what=="dsp"     ) { sprintf(buf,"%g\n",dsp);        s=buf; }
	if(what=="dave"    ) { sprintf(buf,"%g\n",dave);       s=buf; }
	if(what=="shift"   ) { sprintf(buf,"%g\n",shift);      s=buf; }
	if(what=="dz"      ) { sprintf(buf,"%g\n",dz);         s=buf; }
	if(what=="shift01" ) { sprintf(buf,"%g\n",shift01);    s=buf; }
	if(what=="shift1"  ) { sprintf(buf,"%g\n",shift1);     s=buf; }
	if(what=="dnl"     ) { sprintf(buf,"%g\n",dnl);        s=buf; }
	if(what=="dtanth1" ) { sprintf(buf,"%g\n",dtanth1);    s=buf; }
	if(what=="l"       ) { sprintf(buf,"%g\n",l);          s=buf; }
	if(what=="nl"      ) { sprintf(buf,"%g\n",nl);         s=buf; }
	if(what=="kp"      ) { sprintf(buf,"%g\n",kp);         s=buf; }
	if(what=="ks"      ) { sprintf(buf,"%g\n",ks);         s=buf; }
	if(what=="dwldf"   ) { sprintf(buf,"%g\n",dwldf);      s=buf; }
	if(what=="dtdf"    ) { sprintf(buf,"%g\n",dtdf);       s=buf; }
	if(what=="dthdf"   ) { sprintf(buf,"%g\n",dthdf);      s=buf; }
	if(what=="dtempdf" ) { sprintf(buf,"%g\n",dtempdf);    s=buf; }
	
	return s;
}

double cOptics::Deviation(std::string glass,std::string glass1,std::string color,double incident_angle){
	// ���ܕΊp�̌v�Z
	// glass,glass1 : ���ˑ��C�o�ˑ��̔}��
	double th,n,n1;
	th=incident_angle*PI/180;
	n =Re(Index(glass, Wavelength(color)));
	n1=Re(Index(glass1,Wavelength(color)));
	return ( th-asin(n*sin(th)/n1) )*180/PI;
}

std::string cOptics::PrismDeviation(std::string glass,std::string color,double apex_angle,double incident_angle){
	std::string s;
	double n, al, th,th1,th2,th3, gamma1,gamma2;
	char buf[1000];

	n=Re(Index(glass,Wavelength(color)));
	al=apex_angle*PI/180;
	th=incident_angle*PI/180;
	th1=asin( sin(th)/n );
	th2=al-th1;
	th3=asin( sin(th2)*n );
	gamma1=cos(th)/(n*cos(th1));
	gamma2=(n*cos(th2))/cos(th3);

	sprintf(buf,"�ɍ�= %s\n"        ,glass.c_str());       s+=buf;
	sprintf(buf,"�g��= %s\n"        ,color.c_str());       s+=buf;
	sprintf(buf,"���p= %g\n"        ,apex_angle);          s+=buf;
	sprintf(buf,"���ˊp= %g��\n"    ,incident_angle);      s+=buf;
	sprintf(buf,"  �o�ˊp= %g��\n"  ,th3*180/PI);          s+=buf;  
	sprintf(buf,"  �ӂ�p= %g��\n"  ,(th+th3-al)*180/PI ); s+=buf;       
	sprintf(buf,"  ��p= %g\n"       ,100*tan(th+th3-al) ); s+=buf;
	sprintf(buf,"  �p�{��= %g\n"    ,gamma1*gamma2);       s+=buf;
	sprintf(buf,"  ���{��= %g\n"    ,1/gamma1/gamma2);     s+=buf;
	return s;
}

double cOptics::MeasuredDistortion(double x1,double y1,double x2,double y2,double x3,double y3){
	// ���� (x1,y1)-(x3,y3) �� �_(x2,y2) ������Ƃ��C
	// �����Ɠ_�̋����𒼐��̒����Ŋ������l��Ԃ��D
	// ���̒l�́C���ɁC
	//   (x1,y1),(x3,y3),(x2,y2)�������𒆐S�Ƃ��鐳���`(2a x 2a)�̕ӂ̑��̗��[�y�ђ��_�ŁC
	//   ���w�n�̘c�Ȃ�3�������̂�(�����t��)
	// �̂Ƃ��́Ch=a�ɂ�����c�� (y'-y'0)/y'0 �̂قڔ����ƂȂ�D
	point p1(x1,y1), p2(x2,y2), p3(x3,y3);
	return distance(p2,line(p1,p3))/distance(p1,p3);
}

double cOptics::BerekDOF(std::string &s,double wl_nm,double NAo,double TotalM,double VAmin){
	// Berek�̎��ɂ�镨�̑�DOF����m�P�ʂŕԂ�
	// wl_nm  = nm�P�ʔg��
	// NAo    = ���̑�NA
	// TotalM = �����{��
	// VAmin  = ���P�ʂŕ\������̕���\
	double dof1,dof2;
	char buf[1000];
	dof1=(wl_nm/1000)/2/NAo/NAo;                                        // �g�����w��
	dof2= TotalM==0 ? 0 : tan(VAmin/60*PI/180)*(250*1000)/TotalM/NAo;   // �􉽌��w��
	// �􉽌��w���́g�􉽌��w�i�O��j�h(5.2)���邢��(5.3)���ŁC
	//     ��= (250mm/TotalM)*VAmin) = �œ_����*VAmin
	// �Ƃ������̂Ɉ�v�D
	sprintf(buf,"DOF= ��/(2*NAo^2) + (1/NAo)*(250mm/TotalM)*VA\n" ); s=buf;
	sprintf(buf,"��(nm)=%g\n"      , wl_nm);     s+=buf;
	sprintf(buf,"NAo=%g\n"         , NAo);       s+=buf;
	sprintf(buf,"TotalM=%g\n"      , TotalM);    s+=buf;
	sprintf(buf,"VA(min)=%g\n"     , VAmin);     s+=buf;
	sprintf(buf," DOF(��m)=%g\n" , dof1+dof2);   s+=buf;
	sprintf(buf,"  �g�����w��(��m)=%g\n",dof1);  s+=buf;
	sprintf(buf,"  �􉽌��w��(��m)=%g\n",dof2);  s+=buf;
	return dof1+dof2;
}

double cOptics::ScheimpflugTH1(double ms,double th_deg,double th_l_deg){
	// ms       = s�����{��
	// th_deg   = ���̖ʖ@���ƌ����̂Ȃ��p(deg) ���FScheimpflug�J�����̋p��90-th_deg
	// th_deg_l = �����ƃ����Y�����̂Ȃ��p(deg)
	//
	// �߂�l   = ���̖ʁC�����Y�C���ʂ�Scheimpflug�����𖞂����Ƃ���
	//          = ���ʖ@���ƌ����̂Ȃ��p(deg)
	double th,th_l, a;

	th=th_deg*PI/180;
	th_l=th_l_deg*PI/180;
	a=cos(th)/sin(th-th_l);
	return acos( a*cos(th_l)/sqrt(ms*ms+a*a+2*ms*a*sin(th_l)) )*180/PI;
}

std::string cOptics::Scheimpflug(double ms,double th_deg,double th_l_deg){
	std::string s;
	char buf[1000];
	double th1_deg, k;

	th1_deg=ScheimpflugTH1(ms,th_deg,th_l_deg);
	k=cos(th_deg*PI/180)/cos(th1_deg*PI/180);
	sprintf(buf,"s�����{��ms     = %g\n" ,ms);          s+=buf;
	sprintf(buf,"���̖ʓ��ˊpth  = %gdeg\n" ,th_deg);   s+=buf;
	sprintf(buf,"�����Y���ˊpth_l= %gdeg\n" ,th_l_deg); s+=buf;
	sprintf(buf,"  ���ʓ��ˊpth1 = %gdeg\n" ,th1_deg);  s+=buf;
	sprintf(buf,"  p�����{��mp   = %g\n"    ,ms*k);     s+=buf;
	sprintf(buf,"  mp/ms         = %g\n"    ,k);        s+=buf;
	return s;
}

double cOptics::BrewsterAngle(std::string glass_in,std::string glass_out,std::string color){
	complex n_in,n_out;

	n_in =Index(glass_in,Wavelength(color));
	n_out=Index(glass_out,Wavelength(color));
	return atan2(Re(n_out),Re(n_in))*180/PI;
}

double cOptics::ChirpedPulseWidth(double dt0_fs,double GDD_fs2){
	//  �K�E�X�^�p���X���Q�x�����U���󂯂���̃p���X����tout(fs)��Ԃ��D
	//      dt0_fs  : �ŏ��̃p���X����t(fs)
	//      GDD_fs2 : �Q�x�����UGDD(fs^2)
	// 
	//  ��tout = ��(��t^4+16*(ln2)^2*GDD^2) / ��t (http://www.newport-japan.jp/pdf/0321.pdf)

	return sqrt( dt0_fs*dt0_fs*dt0_fs*dt0_fs + 7.68725*GDD_fs2*GDD_fs2 )/dt0_fs;
}

double cOptics::ODToT(double od){
	// ���w�Z�x�ɑ΂��铧�ߗ������߂�
	return pow(10,-od);
}

double cOptics::CoherenceLength(double wl_nm,double dwl_fwhm_nm){
	// �g�����z�����l�����ɂ̃K�E�X���z�̂Ƃ��C
	//    �R�q�[�����X�� = (4ln2/��)(��o^2/����)
	// �y���ӁzOCT�̕���\�͌��̉����ɂ�肳���2�Ŋ���
	return 0.8825424*wl_nm*wl_nm/dwl_fwhm_nm /1000000;  // �R�q�[�����X��(mm�P��)
}

double cOptics::Derivative(std::string command,double *x,double dx){
	std::ofstream to;
	int i, nc=0;
	double x0, y0,y1,y2,y11,y22, a,a_pre, result;
	const int MakeLog=0;
	const int MAXTIMES=10;
	const int METHOD_OLD=0; // METHOD_OLD=true �ł͋��A���S���Y���ɂ��D���̕����Ȃ���������������ꍇ������D
	                        // METHOD_OLD��true�ɂ����ꍇ�̕ύX�_�͈ȉ��D
							//  �E5�_�����͎g��Ȃ�
							//  �E�I������ 3,4 �͎g��Ȃ�
	clock_t start;
	// �Ⴆ�� cyl�l(�}�C�i�X�\��)��0�Ő܂�Ԃ�������D
	// *x���܂�Ԃ��_�̋ߖT�ɂ���Cdx���傫���� *x+dx��*x-dx�ł̌��z�͕������t�ł���C
	// dx���������Ȃ��Č��z���������ɂȂ�܂Ŕ����W�������܂�Ȃ��D
	// �����ł�MAXTIME=5�ł͕s���ŁC�����݌v�̎������x���Ȃ邱�Ƃ��������D

	if(MakeLog){
		to.open("DerivativeLog.txt",std::ios::app);
		to.precision(6);
		to << "command=" << command << std::endl;
		to << "x=" << *x;
		to.precision(16);
		to << " y=" << atof(cmd(command,1).c_str()) << std::endl;
		to.precision(6);
		start=clock();
	}

	if(dx==0){
		result=0;
	}
	else{
		x0=*x;
		y0=atof(cmd(command,1).c_str());

		for(i=1; i<=MAXTIMES; i++){
			*x=x0+dx; y1=atof(cmd(command,1).c_str());
			*x=x0-dx; y2=atof(cmd(command,1).c_str());

			if(MakeLog){
				to << " i=" << i << std::endl;
				to << " x0=" << x0 << " y0=" << y0 << std::endl;
				to.precision(16);
				to << " dx=" << dx << " y1=" << y1 << " y2=" << y2 << std::endl;
				to.precision(6);
			}

			if(METHOD_OLD){
				double dy1,dy2,a1,a2;

				dy1=y1-y0; a1=dy1/dx;
				dy2=y0-y2; a2=dy2/dx;
				a=(a1+a2)/2;  // ����� a=(y1-y2)/(dx*2) �Ə����ƂȂ����������x���Ȃ�ꍇ����i�������̉e���H�j
				
				if(MakeLog){
					to.precision(16);
					to << "  dy1=" << dy1 << " dy2=" << dy2 << std::endl;
					to.precision(6);
					to << " a=" << a << std::endl;
				}

				if(a==0){result=0; nc=0; break;}
				if(fabs((a1-a2)/a)<0.01){result=a; nc=1; break;}     // �����E���̌X���̍��ɂ��Da�̐��x�͂����炭1%�ȓ� (�I������1)
				if(i>=2 && (fabs((a-a_pre)/a)<0.01)){result=a; nc=2; break;} // �O��Ƃ̍��ɂ��Da�̐��x�͂����炭1%�ȓ� (�I������2)

				// �I������1 �����ł��悢���C����2 �������邱�ƂŁC���������Ȃ� (2016.06.29)
				//  (����2 �݂̂ł� i=1 �Ń��[�v�𔲂��邱�Ƃ��ł��Ȃ��i2��ȏ�̃��[�v�v�j�̂ŏ���1 ���K�v)

				a_pre=result=a;
				dx*=0.1;
			}
			else{
				if(i==1){
					a=(y1-y2)/(dx*2);  // i=1��ڂ�3�_�����ix��4�ӏ��ł�y�l�v�Z���K�v��5�_�������ŏ�����g��Ȃ��j
				}
				else{
					a=(-y11+8*y1-8*y2+y22)/(dx*12);  // i=2��ڈȍ~��5�_�����iy11,y22�͑O��̒l�𗘗p�j
				}
				
				if(MakeLog){
					to.precision(16);
					to << "  y1-y0=" << y1-y0 << " y0-y2=" << y0-y2 << std::endl; 
					to.precision(6);
					to << " a=" << a << std::endl;
				}

				if(a==0){result=0; nc=0; break;}
				if(fabs((y1+y2-y0*2)/dx/a)<0.01){result=a; nc=1; break;}     // �����E���̌X���̍��ɂ��Da�̐��x�͂����炭1%�ȓ� (�I������1)
				if(i>=2 && (fabs((a-a_pre)/a)<0.01)){result=a; nc=2; break;} // �O��Ƃ̍��ɂ��Da�̐��x�͂����炭1%�ȓ� (�I������2)

				if(y0!=0){
					if( (fabs((y1-y0)/y0)<1e-13) || (fabs((y0-y2)/y0)<1e-13) ) { // y�̕ω��ʂ̗L�������Ȃ��Ȃ��Ă��邩������Ȃ� (�I������3)
						result= i==1 ? a : a_pre; nc=3; break; 
					} 
				}

				if(x0!=0){
					if( fabs(dx/x0)<1e-13 ){   // x�}dx�̗L�������Ȃ��Ȃ��Ă��邩������Ȃ� (�I������4)
						result= i==1 ? a : a_pre; nc=4; break;
					}
				}

				// �I������1 �����ł��悢���C����2 �������邱�ƂŁC���������Ȃ� (2016.06.29)
				//  (����2 �݂̂ł� i=1 �Ń��[�v�𔲂��邱�Ƃ��ł��Ȃ��i2��ȏ�̃��[�v�v�j�̂ŏ���1 ���K�v)
				// �I������3,4��ǉ��i2018.11.01)

				a_pre=result=a;
				y11=y1;
				y22=y2;
				dx*=0.5;  // 5�_�����őO��̒l�𗘗p�iy11=y1, y22=y2�j���Ă���̂ŐV����dx�� *0.5 �ł���K�v������
			}
		}

		*x=x0; // push(), pop() �ł��ۑ��C�������ł��邪�C���s���x�����ɒx���Ȃ�D
		       // �����Cpush(),pop()�ōs�Ȃ��Ƃ��� virtual int cOptics::push(){ /* ��̊֐� */ }
		       // ��ǉ�����D"virtual" ��t���Ȃ��ꍇ�CcLens�^�C���X�^���X��cOptics::Derivative()���Ă�ł�
		       // cOptics::Derivative�֐��̒���push()��cOptics::push() (��̊֐�) �ł���D
		       // virtual������� cLens::push()���Ă΂��D

		if(MakeLog){
			to << "Lap time=" << LapTime(start) << std::endl;
			to << "Broken by condition " << nc << std::endl;
			to << std::endl;
		}
	}

	return result;
}

std::string cOptics::interpret_target(std::string sentence,
							 std::string &com,std::string &sign,double &tol,double &target,double &weight,
							 double set_weight/*=1*/)
{
	// s�̏��� : com [sign | �} tol] target weight �����߂��āCcom,sgn,tol,target,weight��ݒ肷��D
	// �߂�l��weight��set_weight���悶��i�֐� preprocess_targets �Ŏg���j�D
	com=arg(sentence,0);
	sign=trim(arg(sentence,1),1);

	if(sign=="=" || sign=="<" || sign==">"){
		target=atof(arg(sentence,2).c_str());
		weight=atof(arg(sentence,3).c_str());
		sentence=replace_arg(sentence,3,str(weight*set_weight));
	}
	else if(sign=="+-" || sign=="�}"){
		// �G�N�Z���̃Z����"+- ..." �Ƃ���ƃG���[���o�Ĕς킵���̂�
		// "�}"���g����悤�ɂ����D
		tol=atof(arg(sentence,2).c_str());
		target=atof(arg(sentence,3).c_str());
		weight=atof(arg(sentence,4).c_str());
		sentence=replace_arg(sentence,4,str(weight*set_weight));
	}
	else{
		// "=","<",">","+-" ���ȗ����ꂽ�ꍇ
		sign="=";
		target=atof(arg(sentence,1).c_str());
		weight=atof(arg(sentence,2).c_str());
		sentence=replace_arg(sentence,2,str(weight*set_weight));
	}

	return sentence;
}


void cOptics::symboltoval(std::string &s) const{
	int i;

	if(symbol.GetSize()>0){
		for(i=1; i<=symbol.GetSize(); i++){
			if(s==symbol[i]) s=value[i];
		}
	}
}

double cOptics::arg(const complex &z){
	// �{�֐����Ȃ��ƁC
	// cOptics::arg(const std::string& sentence, int n)
	// �����邽�߁Ccomplex.h�Ő錾����� arg(const complex& z)
	// ���ĂԂ̂ɉ��Z�q::���K�v(::arg(z))�ɂȂ��Ă��܂��D

	return ::arg(z);
}

std::string cOptics::arg(const std::string& sentence, int n){
	std::string s;

	s=::arg(sentence,n);
	symboltoval(s);
	return s;
}

bool cOptics::is_numeric(std::string& s){
	//  s�܂���s=cmd(s,1)������\��������̂Ƃ�true��Ԃ��D
	std::string buf;

	if(s==""){
		return false;
	}
	else if(::is_numeric(s)){   // string.h
		return true;
	}
	else{
		push();   //  ���z�֐��D�q�N���X(cLens�Ȃ�)��push���ĂԁD
		buf=cmd(s,1);  // �y���Ӂz ��2�������^�Ȃ̂ŁC�߂�l���񐔎��̕�����ł���R�}���h�� buf="" �ɂȂ��Ă��܂��D
		pop();    //  ���z�֐��D�q�N���X(cLens�Ȃ�)��pop���ĂԁD
		//�y�d�v�zpush(),pop()�ɂ��C*this�͕ω������Ȃ��D 2016.03.14
		//   ��F
		//   s="Optimize (r 8 .0002) (f 20 10) 10  // r(8)��ω������āCf��20mm�ɂ���D
		//   ��cOptics::cmd�̈������Ƃ����
		//   cLnes1::propery��s�ɑ΂��Ď��s�����D
		//   ���̒��� is_numeric("r 8 .0002") �����s����邪�C
		//   ��L�� ::is_numeric(buf=cmd(s,1))  �Ƃ��Ă��܂��ƁC
		//   �{���� "r(8)��.0002�Âω������čœK������" �Ƃ����Ӑ}�ɔ����āC
		//   r(8)��0.0002���������Ă��܂��D
		if(::is_numeric(buf)){
			s=buf;
			return true;
		}
		else{
			// --- �߂�l���񐔎��̕�����ł���R�}���h���W�J���� ------------------
			buf=cmd(s,0);
			if(buf!="") s=buf;
			// ---------------------------------------------------------------------
			return false;
		}
	}
}

std::string cOptics::cmd(const std::string& command,int val){
	// val=ture�̂Ƃ��C�ꕔ�̃R�}���h��Basic���ɂ�val�֐��ŏ������邽�߂�
	// ���l��\��������̂ݕԂ��D
	// command��������sentence����Ȃ�Ƃ��́C
	//     val=false : �������ʂ𕶎���Ƃ��ĉ����ĕԂ��D
	//     val=true  : �Ō��sentence(���ʂ�""�ł���sentence�͏���)
	//                 �̌��ʐ��l��\���������Ԃ��D
	//                 (�ȑO�́C�esentence�̌��ʂ̍��v��Ԃ��Ă������C
	//                  ����� "+" �R�}���h�ŉ\�ł��邵�C
	//                  �Ō�̌��ʂ����m�肽���ꍇ������̂ł��̂悤�ɂ����D 2014.10.08 )
	//
	//                 �n�̕ύX���N����Ȃ��悤�ɁCcommand�̍Ō�� pop ���K�v�Ȃ��Ƃ�����D
	//                 pop �̖߂�l��""�Ȃ̂ŕK�v�Ȑ��l��Ԃ����߂ɁC
	//                 ��L�g���ʂ�""�ł���sentence�͏����h���ƂƂ����D 2016.05.11
	std::string s,ss,sv,com;
	double v1,v2,v3;
	char buf[1000];
	int i,n;

	s=ss="";
	n=sentences(command);

	for(i=1; i<=n; ++i){
		com=sentence(command,i);
		symboltoval(com);
		// scmd()�͕�����com��P��sentence�Ƃ��ď�������
		s=scmd(com,val);  // ���z�֐��D�q�N���X��scmd���Ă΂��D
		ss+=s;
		if(s!="") sv=s;   // val=true�̏ꍇ�̖߂�l�́C�Ō��""�łȂ� s �ɂ��D
		s="";
	}

	if(val){
		switch(words(sv)){
		case 2:  // 2�����x�N�g����z��
			v1=atof(word(sv,1,1).c_str());
			v2=atof(word(sv,2,1).c_str());
			sprintf(buf,"%.15g %.15g\n",v1,v2);
			break;
		case 3:  // 3�����x�N�g����z��
			v1=atof(word(sv,1,1).c_str());
			v2=atof(word(sv,2,1).c_str());
			v3=atof(word(sv,3,1).c_str());
			sprintf(buf,"%.15g %.15g %.15g\n",v1,v2,v3);
			break;
		default:
			v1=atof(sv.c_str());
			sprintf(buf,"%.15g\n", v1);
			break;
		}
		// ���ӁF����%f�͏����_�ȉ��̌������Œ�Ȃ̂ŁC�l���������ƗL�������Ȃ��Ȃ�D
		//       %.15g���悢�D
		//       �P��%g���ƗL�����̓f�t�H���g��6���ƂȂ�C�����݌v�̔����W���v�Z��
		//       �������Ƃ�Ƃ��L�������s���C���ɂȂ�D
		//       scmd�֐��ł����̂悤�ɂ��邱�ƁD
		return ::is_numeric(ss)  ? buf : "";  // ::is_numeric(ss)�łȂ��Ƃ���"0"�łȂȂ�""��Ԃ��D
		                                      // �������Ȃ��ƁCcOptics::is_numeric(command)��
		                                      // �K��true�ƂȂ��Ă��܂��D
	}
	else{
		return ss;
	}
}

std::string cOptics::scmd(std::string com,int val){
	// val=ture�̂Ƃ��C�ꕔ�̃R�}���h��Basic���ɂ�val�֐��ŏ������邽�߂�
	// ���l��\��������̂ݕԂ��D

	std::string s;
	char buf[1000];
	std::string s0,s1,s2,s3,s4,s5,s6,s7;
	bool b1,b2,b3,b4,b5,b6;
	
	s0=arg(com,0);

	if(s0=="??"){
		s+="+ - * / ^ == > < Abs Ave Call If Rem Set SetStr ";
		s+="\n";
		s+="AbbeNumber AcosD AsinD AtanD Atan2D BerekDOF BrewsterAngle ";
		s+="CoherenceLength CosD Deviation DistanceX DistanceY DistanceZ ";
		s+="GlassCode GroupIndex HerzEq Index ";
		s+="LumenToWatt LuminousEfficiency ";
		s+="ODToT Plate PrismDeviation ";
		s+="Scheimpflug ScheimpflugTH1 SinD Spherometer TanD ";
		s+="WattToLumen Wavelength XYToDominantWavelength ";
		s+="\n";

		return s;
	}

	if(::is_numeric(s0)){
		// s0�̖`��������\���Ƃ�
		sprintf(buf,"%.15g\n",atof(s0.c_str()));
		s=buf;
		return s;
	}

	if(s0=="+" || s0=="-" || s0=="*" || s0=="/" || s0=="^" || s0=="==" || s0==">" || s0=="<"){
		double v1,v2,v3,v4,v5, x;
		s1=arg(com,1);
		s2=arg(com,2);
		s3=arg(com,3);
		s4=arg(com,4);
		s5=arg(com,5);
		if(s1=="?"){
			if(s0=="+" || s0=="*"){
				s=s0+" com1 com2 [com3 [com4 [com5]]]\n";
			}
			else if(s0=="-"){
				s=s0+" com1 [com2]\n";
			}
			else if(s0=="/" || s0=="^" || s0=="==" || s0==">" || s0=="<"){
				s=s0+" com1 com2\n";
			}
		}
		else{
			v1=atof(cmd(s1,1).c_str());
			v2=atof(cmd(s2,1).c_str());
			v3=atof(cmd(s3,1).c_str());
			v4=atof(cmd(s4,1).c_str());
			v5=atof(cmd(s5,1).c_str());
			if(s0=="+"){
				switch(args(com)){
					case 2:  x=v1+v2;          break;
					case 3:  x=v1+v2+v3;       break;
					case 4:  x=v1+v2+v3+v4;    break;
					case 5:  x=v1+v2+v3+v4+v5; break;
					default: x=0;
				}
			}
			else if(s0=="-"){
				switch(args(com)){
					case 1:  x=-v1;   break;
					case 2:  x=v1-v2; break;
					default: x=0;
				}
			}
			else if(s0=="*"){
				switch(args(com)){
					case 2:  x=v1*v2;          break;
					case 3:  x=v1*v2*v3;       break;
					case 4:  x=v1*v2*v3*v4;    break;
					case 5:  x=v1*v2*v3*v4*v5; break;
					default: x=0;
				}
			}
			else if(s0=="/"){
				x= v2==0 ? 1e100 : v1/v2;
			}
			else if(s0=="^"){
				x= pow(v1,v2);
			}
			else if(s0=="=="){
				x= v1==v2 ? 1 : 0;
			}
			else if(s0==">"){
				x= v1>v2 ? 1 : 0;
			}
			else if(s0=="<"){
				x= v1<v2 ? 1 : 0;
			}

			sprintf(buf,"%.15g\n",x); s=buf;
		}
		return s;
	}
	if(s0=="Abs" || s0=="abs"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Abs com\n";
		}
		else{
			sprintf(buf,"%.15g\n",fabs(atof(scmd(s1,1).c_str())));
			s=buf;
		}
		return s;
	}
	if(s0=="Ave" || s0=="ave"){
		double v1,v2,v3,v4,v5, x;
		s1=arg(com,1);
		s2=arg(com,2);
		s3=arg(com,3);
		s4=arg(com,4);
		s5=arg(com,5);
		if(s1=="?"){
			s="Ave com1 com2 [com3 [com4 [com5]]]\n";
		}
		else{
			v1=atof(cmd(s1,1).c_str());
			v2=atof(cmd(s2,1).c_str());
			v3=atof(cmd(s3,1).c_str());
			v4=atof(cmd(s4,1).c_str());
			v5=atof(cmd(s5,1).c_str());
			switch(args(com)){
				case 2:  x=(v1+v2)/2;          break;
				case 3:  x=(v1+v2+v3)/3;       break;
				case 4:  x=(v1+v2+v3+v4)/4;    break;
				case 5:  x=(v1+v2+v3+v4+v5)/5; break;
				default: x=0;
			}
			sprintf(buf,"%.15g\n",x); s=buf;
		}
		return s;
	}
	if(s0=="Call" || s0=="call"){
		// �t�@�C���̓��e���R�}���h�Ƃ��Ď��s����D�g���q .txt �͏ȗ��\�D
		// �t�@�C�����e��"arg1","arg2", ... �͈���arg1,arg2, ... �ɒu�������D
		std::ifstream from;
		std::string tmp,com1;
		int i;
		const int ARGS_MAX=9;            // ������9�܂�
		std::string a[ARGS_MAX+1];

		s1=arg(com,1);
		if(s1=="?"){
			s+="Call filename [arg1 arg2 ... ]\n";
			s+="   - extension '.txt' can be dropped from filename\n";
			s+="   - max number of arguments is 9\n";
		}
		else{
			from.open(s1.c_str());
			if(!from){                          // open�����s�����ꍇ
				from.clear();                   // ���ӁFclear()�����Ȃ��ƁC���s��open�͐������Ȃ��D
				from.open((s1+".txt").c_str()); // �g���q ".txt" ��t�����čēxopen�����݂�D
			}
			for(i=1; i<=ARGS_MAX; ++i){
				a[i]=arg(com,1+i);
				if( (tmp=cmd(a[i],1))!="" ) a[i]=tmp;
			}
			if(from){
				while(from>>tmp){
					for(i=1; i<=ARGS_MAX; ++i){
						if(a[i]!="") tmp=replace(tmp,"arg"+trim(str(i),0),a[i]);  // �u��
					}
					com1+=tmp+" ";
				}
				s=cmd(com1,val);
			}			
		}
		return s;
	}
	if(s0=="If" || s0=="if"){
		std::string condition, action;

		s1=arg(com,1);
		s2=arg(com,2);
		if(s1=="?"){
			s+="If condition action\n";
		}
		else{
			if(atof(cmd(s1,1).c_str())!=0) s+=cmd(s2,1);
			
		}
		return s;
	}
	if(s0=="Rem" || s0=="rem" || s0=="//"){
		// ���ߕ�
		s1=arg(com,1);
		if(s1=="?"){
			s="Rem [comment]\n";
		}
		return s;
	}
	if(s0=="Set" || s0=="set"){
		// �ϐ�symbol�ɒlcmd(value)��������i������value�̎��s���� cmd(value)�����j
		// ��F
		//   r 1 10;
		//   set a f;  a�ɂ�cmd(f)����������D
		//   r 1 11;
		//   a;        r(1)=10�ł̏œ_������Ԃ��D("setstr a f" �Ƃ����Ƃ��� r(1)=11�ł̏œ_�������Ԃ����j
		int i;
		std::string tmp;
		s1=::arg(com,1);  // ���ӁF "::"��t���Ȃ���s1��cOptics::arg()�ɂ�菑���������Ă��܂��D
		s2=::arg(com,2);
		if(s1=="?"){
			s+="Set symbol value\n";
		}
		else{
			if(s1=="" ){
				// ���݂̐ݒ�̈ꗗ��\��
				for(i=1; i<=symbol.GetSize(); i++){
					s+=symbol[i]+'\t'+value[i]+'\n';
				}
			}
			else{
				tmp=cmd(s2,1);  // "set x (* x 2)" �̂悤�ɍċA�ďo������Ƃ��́C
				                // Remove�̑O�� "* x 2" �����s����K�v����
				for(i=symbol.GetSize(); i>0; i--){
					if(symbol[i]==s1){
						symbol.Remove(i);
						value.Remove(i);
					}
				}
				if(s2!=""){
					symbol.AddTail(s1);
					value.AddTail(tmp);
					s=tmp;
				}
			}
		}
		return s;
	}
	if(s0=="SetStr" || s0=="setstr"){
		// �ϐ�symbol�ɕ�����str��������
		// ��F
		//   r 1 10;
		//   setstr a f;  a=scmd("f")�łȂ����Ƃɒ��ӁDa��"f"�ŁCr(1)=10�ł̏œ_�����l�ł͂Ȃ��D
		//   r 1 11;
		//   a;        r(1)=11�ł̏œ_������Ԃ��D
		//
		// ���ӁF
		//   setstr x (* x 2) �̂悤�ɍċA�Ăяo��������Ɩ������[�v�ƂȂ�G���[�D

		int i;
		s1=::arg(com,1);  // ���ӁF "::"��t���Ȃ���s1��cOptics::arg()�ɂ�菑���������Ă��܂��D
		s2=::arg(com,2);
		if(s1=="?"){
			s+="SetStr symbol str\n";
		}
		else{
			if(s1=="" ){
				// ���݂̐ݒ�̈ꗗ��\��
				for(i=1; i<=symbol.GetSize(); i++){
					s+=symbol[i]+'\t'+value[i]+'\n';
				}
			}
			else{
				for(i=symbol.GetSize(); i>0; i--){
					if(symbol[i]==s1){
						symbol.Remove(i);
						value.Remove(i);
					}
				}
				if(s2!=""){
					symbol.AddTail(s1);
					value.AddTail(s2);
					s=cmd(s2,val);
				}
			}
		}
		return s;
	}

	if(s0=="AbbeNumber" || s0=="abbenumber"){
		double nu;
		s1=arg(com,1);
		s2=arg(com,2);
		s3=arg(com,3);
		s4=arg(com,4);
		if(s1=="?"){
			s="AbbeNumber glassname [color0=d color1=F color2=C]\n";
		}
		else{
			if(args(com)==4){
				nu=AbbeNumber(s1,s2,s3,s4);
				if(val){
					sprintf(buf,"%.15g\n",nu);
					s=buf;
				}
				else{
					sprintf(buf,"%s (N(%s)-1)/(N(%s)-N(%s))=%g\n",
					        s1.c_str(),s2.c_str(),s3.c_str(),s4.c_str(),nu);
					s=buf;
				}
			}
			else if(args(com)==1){
				nu=AbbeNumber(s1);
				if(val){
					sprintf(buf,"%.15g\n",nu);
					s=buf;
				}
				else{
					sprintf(buf,"%s (N(d)-1)/(N(F)-N(C))=%g\n", s1.c_str(),nu);
					s=buf;
				}
			}
		}
		return s;
	}
	if(s0=="AcosD" || s0=="acosd"){
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="AcosD val (arccosine in degree)\n";
		}
		else if( b1 ) {
			val=atof(s1.c_str());
			sprintf(buf,"%.15g\n",acos(val)*180/PI); s=buf;
		}
		return s;
	}
	if(s0=="AsinD" || s0=="asind"){
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="AsinD val (arcsine in degree)\n";
		}
		else if( b1 ) {
			val=atof(s1.c_str());
			sprintf(buf,"%.15g\n",asin(val)*180/PI); s=buf;
		}
		return s;
	}
	if(s0=="AtanD" || s0=="atand"){
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="AtanD val (arctangent in degree)\n";
		}
		else if( b1 ) {
			val=atof(s1.c_str());
			sprintf(buf,"%.15g\n",atan(val)*180/PI); s=buf;
		}
		return s;
	}
	if(s0=="Atan2D" || s0=="atan2d"){
		double y,x;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Atan2D y x (atan2 in degree)\n";
		}
		else if(b1 && b2) {
			y=atof(s1.c_str());
			x=atof(s2.c_str());
			sprintf(buf,"%.15g\n",atan2(y,x)*180/PI); s=buf;
		}
		return s;
	}
	if(s0=="BerekDOF" || s0=="berekdof"){
		std::string temp;
		double wl_nm,NAo,TotalM,VAmin,dof;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="BerekDOF wl_nm NAo TotalM VAmin\n";
		}
		else if( b1 && b2 && b3 && b4) {
			wl_nm=atof(s1.c_str());
			NAo=atof(s2.c_str());
			TotalM=atof(s3.c_str());
			VAmin=atof(s4.c_str());
			dof=BerekDOF(temp,wl_nm,NAo,TotalM,VAmin);
			if(val){
				sprintf(buf,"%.15g\n",dof); s=buf;
			}
			else{
				s=temp;
			}
		}
		return s;
	}
	if(s0=="BrewsterAngle" || s0=="brewsterangle"){
		std::string glass_in,glass_out,color;
		s1=arg(com,1);
		s2=arg(com,2);
		s3=arg(com,3);
		if(s1=="?"){
			s="BrewsterAngle glass_in glass_out color\n";
		}
		else{
			glass_in=s1;
			glass_out=s2;
			color=s3;
			sprintf(buf,"%.15g\n",BrewsterAngle(glass_in,glass_out,color));
			s=buf;
		}
		return s;
	}
	if(s0=="CoherenceLength" || s0=="coherencelength"){
		double wl_nm,dwl_fwhm_nm;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="CoherenceLength wl_nm dwl_fwhm_nmn\n";
		}
		else if(b1 && b2){
			wl_nm=atof(s1.c_str());
			dwl_fwhm_nm=atof(s2.c_str());
			sprintf(buf,"%.15g\n",CoherenceLength(wl_nm,dwl_fwhm_nm)); s=buf;
		}
		return s;
	}
	if(s0=="CosD" || s0=="cosd"){
		double th_deg;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Cosd th_deg\n";
		}
		else if( b1 ) {
			th_deg=atof(s1.c_str());
			sprintf(buf,"%.15g\n",cos(th_deg*PI/180)); s=buf;
		}
		return s;
	}
	if(s0=="Deviation" || s0=="deviation"){
		std::string glass,glass1,color;
		double incident_angle;
		s1=arg(com,1);
		s2=arg(com,2);
		s3=arg(com,3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Deviation glass glass1 color incident_angle\n";
		}
		else if(b4){
			glass=s1;
			glass1=s2;
			color=s3;
			incident_angle=atof(s4.c_str());
			sprintf(buf,"%.15g\n",Deviation(glass,glass1,color,incident_angle)); s=buf;
		}
		return s;
	}
	if(s0=="DistanceX" || s0=="distancex"){
		double x1,y1,z1,x2,y2,z2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="DistanceX x1 y1 z1 x2 y2 z2\n";
		}
		else if(b1 && b2 && b3 && b5 && b6){
			x1=atof(s1.c_str());
			y1=atof(s2.c_str());
			z1=atof(s3.c_str());
			x2=atof(s4.c_str());
			y2=atof(s5.c_str());
			z2=atof(s6.c_str());
			sprintf(buf,"%.15g\n",DistanceX(x1,y1,z1,x2,y2,z2)); s=buf;
		}
		else if(b1 && b2){
			x1=atof(word(s1,1,1).c_str());
			y1=atof(word(s1,2,1).c_str());
			z1=atof(word(s1,3,1).c_str());
			x2=atof(word(s2,1,1).c_str());
			y2=atof(word(s2,2,1).c_str());
			z2=atof(word(s2,3,1).c_str());
			sprintf(buf,"%.15g\n",DistanceX(x1,y1,z1,x2,y2,z2)); s=buf;
		}
		return s;
	}
	if(s0=="DistanceY" || s0=="distancey"){
		double x1,y1,z1,x2,y2,z2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="DistanceY x1 y1 z1 x2 y2 z2\n";
		}
		else if(b1 && b2 && b3 && b5 && b6){
			x1=atof(s1.c_str());
			y1=atof(s2.c_str());
			z1=atof(s3.c_str());
			x2=atof(s4.c_str());
			y2=atof(s5.c_str());
			z2=atof(s6.c_str());
			sprintf(buf,"%.15g\n",DistanceY(x1,y1,z1,x2,y2,z2)); s=buf;
		}
		else if(b1 && b2){
			x1=atof(word(s1,1,1).c_str());
			y1=atof(word(s1,2,1).c_str());
			z1=atof(word(s1,3,1).c_str());
			x2=atof(word(s2,1,1).c_str());
			y2=atof(word(s2,2,1).c_str());
			z2=atof(word(s2,3,1).c_str());
			sprintf(buf,"%.15g\n",DistanceY(x1,y1,z1,x2,y2,z2)); s=buf;
		}
		return s;
	}
	if(s0=="DistanceZ" || s0=="distancez"){
		double x1,y1,z1,x2,y2,z2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="DistanceZ x1 y1 z1 x2 y2 z2\n";
		}
		else if(b1 && b2 && b3 && b5 && b6){
			x1=atof(s1.c_str());
			y1=atof(s2.c_str());
			z1=atof(s3.c_str());
			x2=atof(s4.c_str());
			y2=atof(s5.c_str());
			z2=atof(s6.c_str());
			sprintf(buf,"%.15g\n",DistanceZ(x1,y1,z1,x2,y2,z2)); s=buf;
		}
		else if(b1 && b2){
			x1=atof(word(s1,1,1).c_str());
			y1=atof(word(s1,2,1).c_str());
			z1=atof(word(s1,3,1).c_str());
			x2=atof(word(s2,1,1).c_str());
			y2=atof(word(s2,2,1).c_str());
			z2=atof(word(s2,3,1).c_str());
			sprintf(buf,"%.15g\n",DistanceZ(x1,y1,z1,x2,y2,z2)); s=buf;
		}
		return s;
	}
	if(s0=="GlassCode" || s0=="glasscode"){
		std::string glassname;
		s1=arg(com,1);
		if(s1=="?"){
			s="GlassCode glassname\n";
		}
		else{
			glassname=s1;
			s=GlassCode(glassname)+'\n';
		}
		return s;
	}
	if(s0=="GroupIndex" || s0=="groupindex"){
		std::string glassname,color;
		double N;
		s1=arg(com,1);
		s2=arg(com,2);
		if(s1=="?"){
			s="GroupIndex glassname color\n";
		}
		else{
			if(args(com)==2){
				glassname=s1;
				color=s2;				
				N=GroupIndex(glassname,Wavelength(color));
				if(val){
					sprintf(buf,"%.15g\n",N);
					s=buf;
				}
				else{
					sprintf(buf,"%s GroupN(%s)=%g\n", glassname.c_str(),color.c_str(),N);
					s=buf;
				}
			}
		}
		return s;
	}
	if(s0=="HerzEq" || s0=="herzeq"){
		double Nd,nud;
		std::string color;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3);
		if(s1=="?"){
			s="HerzEq Nd nud color\n";
		}
		else if( b1 && b2 ) {
			Nd=atof(s1.c_str());
			nud=atof(s2.c_str());
			color=s3;
			sprintf(buf,"%.15g\n",Round(HerzEq(Nd,nud,color),-5)); s=buf;
		}
		return s;
	}
	if(s0=="Index" || s0=="index"){
		std::string glassname,color;
		complex N;
		s1=arg(com,1);
		s2=arg(com,2);
		if(s1=="?"){
			s="Index glassname color\n";
		}
		else{
			if(args(com)==2){
				glassname=s1;
				color=s2;				
				N=Index(glassname,Wavelength(color));
				if(Im(N)==0){
					if(val){
						sprintf(buf,"%.15g\n",Re(N));
						s=buf;
					}
					else{
						sprintf(buf,"%s N(%s)=%g\n", glassname.c_str(),color.c_str(),Re(N));
						s=buf;
					}
				}
				else{
					if(val){
						sprintf(buf,"%.15g %.15g\n",Re(N),Im(N));
						s=buf;
					}
					else{
						sprintf(buf,"%s N(%s)=%g  %g\n", glassname.c_str(),color.c_str(),Re(N),Im(N));
						s=buf;
					}
				}	
			}
		}
		return s;
	}
	if(s0=="LumenToWatt" || s0=="lumentowatt"){
		double wl_nm,flux_lumen;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="LumenToWatt wl_nm flux_lumen\n";
		}
		else if( b1 && b2 ) {
			wl_nm=atof(s1.c_str());
			flux_lumen=atof(s2.c_str());
			sprintf(buf,"%.15g\n",cSpectrum::LumenToWatt(wl_nm,flux_lumen)); s=buf;
		}
		return s;
	}
	if(s0=="LuminousEfficiency" || s0=="luminousefficiency"){
		double wl_nm;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="LuminousEfficiency wl_nm\n";
		}
		else if(b1){
			wl_nm=atof(s1.c_str());
			sprintf(buf,"%.15g\n",cSpectrum::LuminousEfficiency(wl_nm)); s=buf;
		}
		return s;
	}
	if(s0=="Plate" || s0=="plate"){
		std::string glassname,color,what;
		double d,th,dndtemp,alpha;
		s1=arg(com,1);
		s2=arg(com,2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7);
		if(s1=="?"){
			s="Plate glassname color d incident_angle [dndtemp alpha] [what]\n";
		}
		else if( b3 && b4 && b5 && b6){
			glassname=s1;
			color=s2;
			d=atof(s3.c_str());
			th=atof(s4.c_str());
			dndtemp=atof(s5.c_str());
			alpha=atof(s6.c_str());
			what=s7;
			s=Plate(glassname,color,d,th,dndtemp,alpha,what);
		}
		else if( b3 && b4 ){
			glassname=s1;
			color=s2;
			d=atof(s3.c_str());
			th=atof(s4.c_str());
			what=s5;
			s=Plate(glassname,color,d,th,0,0,what);
		}
		return s;
	}
	if(s0=="ODToT" || s0=="odtot"){
		double od;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="ODToT od\n";
		}
		else if(b1){
			od=atof(s1.c_str());
			sprintf(buf,"%.15g\n",ODToT(od)); s=buf;
		}
		return s;
	}
	if(s0=="PrismDeviation" || s0=="prismdeviation"){
		std::string glassname,color;
		double apex_angle,incident_angle;
		s1=arg(com,1);
		s2=arg(com,2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="PrismDeviation glassname color apex_angle incident_angle\n";
		}
		else if( b3 && b4 ){
			glassname=s1;
			color=s2;
			apex_angle=atof(s3.c_str());
			incident_angle=atof(s4.c_str());
			s=PrismDeviation(glassname,color,apex_angle,incident_angle);
		}
		return s;
	}
	if(s0=="Scheimpflug" || s0=="scheimpflug"){
		double ms,th_deg,th_l_deg;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="Scheimpflug ms th_deg th_l_deg\n";
		}
		else if( b1 && b2 && b3 ) {
			ms=atof(s1.c_str());
			th_deg=atof(s2.c_str());
			th_l_deg=atof(s3.c_str());
			s=Scheimpflug(ms,th_deg,th_l_deg);
		}
		return s;
	}
	if(s0=="ScheimpflugTH1" || s0=="scheimpflugth1"){
		double ms,th_deg,th_l_deg;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="ScheimpflugTH1 ms th_deg th_l_deg\n";
		}
		else if( b1 && b2 && b3 ) {
			ms=atof(s1.c_str());
			th_deg=atof(s2.c_str());
			th_l_deg=atof(s3.c_str());
			sprintf(buf,"%.15g\n",ScheimpflugTH1(ms,th_deg,th_l_deg)); s=buf;
		}
		return s;
	}
	if(s0=="SinD" || s0=="sind"){
		double th_deg;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Sind th_deg\n";
		}
		else if( b1 ) {
			th_deg=atof(s1.c_str());
			sprintf(buf,"%.15g\n",sin(th_deg*PI/180)); s=buf;
		}
		return s;
	}
	if(s0=="Spherometer" || s0=="spherometer"){
		double D,H,Rref,Href;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s+="Spherometer D H [Rref=0 Href=0]\n";
			s+="    D    = diameter of measuring head\n";
			s+="    H    = ��h with measuring surface\n";
			s+="    Rref = curvature radius of reference surf(if zero, not compensate the result)\n";
			s+="    Href = ��h with reference surface\n";
			s+="    (H,Rref,Href >0 regardless of whether convex or concave surface)\n";
		}
		else if( b1 && b2 && b3 && b4 ) {
			D   =atof(s1.c_str());
			H   =atof(s2.c_str());
			Rref=atof(s3.c_str());
			Href=atof(s4.c_str());
			s=Spherometer(D,H,Rref,Href);
		}
		else if(b1 && b2){
			D   =atof(s1.c_str());
			H   =atof(s2.c_str());
			s=Spherometer(D,H);
		}
		return s;
	}
	if(s0=="TanD" || s0=="tand"){
		double th_deg;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Tand th_deg\n";
		}
		else if( b1 ) {
			th_deg=atof(s1.c_str());
			sprintf(buf,"%.15g\n",tan(th_deg*PI/180)); s=buf;
		}
		return s;
	}
	if(s0=="WattToLumen" || s0=="watttolumen"){
		double wl_nm,flux_watt;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="WattToLumen wl_nm flux_watt\n";
		}
		else if( b1 && b2 ) {
			wl_nm=atof(s1.c_str());
			flux_watt=atof(s2.c_str());
			sprintf(buf,"%.15g\n",cSpectrum::WattToLumen(wl_nm,flux_watt)); s=buf;
		}
		return s;
	}
	if(s0=="Wavelength" || s0=="wavelength"){
		std::string color;
		s1=arg(com,1); b1=is_numeric(s1);
		// is_numeric����p�����Ȃ��ƁC�Ⴆ�΁C"wavelength (color 1)" �����҂��铮������Ȃ��D
		// �����ŁC "color 1" �� cLens1�̃R�}���h�ő�1�g����\���������Ԃ��D
		// �������Cis_numeric() �̎d�l�ɂ��C"color 1" �̖߂�l�� "632.8" �̂悤�ɐ����̏ꍇ�͖��Ȃ����C
		// "e" �̂悤�ɔ񐔎��̏ꍇ�� "wavelength (color 1)" �� 0 �ɂȂ��Ă��܂��D
		if(s1=="?"){
			s="Wavelength color\n";
		}
		else {
			color=s1;
			sprintf(buf,"%.15g\n",Wavelength(color)); s=buf;
		}
		return s;
	}
	if(s0=="XYToDominantWavelength" || s0=="xytodominantwavelength"){
		double x,y;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="XYToDominantWavelength x y\n";
		}
		else if(b1 && b2){
			x=atof(s1.c_str());
			y=atof(s2.c_str());
			sprintf(buf,"%.15g\n",cSpectrum::XYToDominantWavelength(x,y)); s=buf;
		}
		return s;
	}

	return s;
}



