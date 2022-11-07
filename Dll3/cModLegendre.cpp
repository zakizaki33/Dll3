#include "stdafx.h"
#include "MyDllOptics.h"
#include "cModLegendre.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

double cModLegendre::P(int n,double x,int power/*=-1*/){
	// èCê≥LegendreëΩçÄéÆ Pn(x)
	// äÓíÍ {1,x^3,x^4, ... x^10} Ç Gram-Schmidt Ç≈íºåâªÇµÇΩÇ‡ÇÃÅD
	// ëQâªéÆÇçÏê¨Ç≈Ç´ÇªÇ§Ç…Ç»Ç¢ÇÃÇ≈äeéüêîÇÃéÆÇíºê⁄èëÇ≠ÅD
	//
	// power>=0ÇÃÇ∆Ç´ÇÕÅCÇ◊Ç´ãâêîÇ…ÇµÇΩÇ∆Ç´ÇÃÇªÇÃéüêîê¨ï™ÇÃÇ›Çï‘Ç∑ÅD
	// Ç±ÇÍÇÕÅCcModLegendreÇ…ÇÊÇÈîÒãÖñ ê¨ï™ÇÇ◊Ç´ãâêîï\é¶Ç…ïœä∑Ç∑ÇÈÇ∆Ç´Ç…óòópÇ∑ÇÈÅD
	double result;

	switch(n){
	case 0:
		if(power>=0){
			if(power==0) result=1;
			else         result=0;
		}
		else{
			result=1;
		}
		break;
	case 3:
		if(power>=0){
			if(power==3) result=pow(x,3);
			else         result=0;
		}
		else{
			result=pow(x,3);
		}
		break;
	case 4:
		if(power>=0){
			if     (power==4) result=(1.0/4.0)*5*pow(x,4);
			else if(power==0) result=(1.0/4.0)*(-1);
			else              result=0;
		}
		else{
			result=(1.0/4.0)*( 5*pow(x,4)-1 );
		}
		break;
	case 5:
		if(power>=0){
			if     (power==5) result=(1.0/2.0)*9*pow(x,5);
			else if(power==3) result=(1.0/2.0)*(-7)*pow(x,3);
			else              result=0;
		}
		else{
			result=(1.0/2.0)*( 9*pow(x,5)-7*pow(x,3) );
		}
		break;
	case 6:
		if(power>=0){
			if     (power==6) result=(1.0/24.0)*154*pow(x,6);
			else if(power==4) result=(1.0/24.0)*(-135)*pow(x,4);
			else if(power==0) result=(1.0/24.0)*5;
			else              result=0;
		}
		else{
			result=(1.0/24.0)*( 154*pow(x,6)-135*pow(x,4)+5 );
		}
		break;
	case 7:
		if(power>=0){
			if     (power==7) result=(1.0/8.0)*143*pow(x,7);
			else if(power==5) result=(1.0/8.0)*(-198)*pow(x,5);
			else if(power==3) result=(1.0/8.0)*63*pow(x,3);
			else              result=0;
		}
		else{
			result=(1.0/8.0)*( 143*pow(x,7)-198*pow(x,5)+63*pow(x,3) );
		}
		break;
	case 8:
		if(power>=0){
			if     (power==8) result=(1.0/192.0)*5265*pow(x,8);
			else if(power==6) result=(1.0/192.0)*(-8008)*pow(x,6);
			else if(power==4) result=(1.0/192.0)*2970*pow(x,4);
			else if(power==0) result=(1.0/192.0)*(-35);
			else              result=0;
		}
		else{
			result=(1.0/192.0)*( 5265*pow(x,8)-8008*pow(x,6)+2970*pow(x,4)-35 );
		}
		break;
	case 9:
		if(power>=0){
			if     (power==9) result=(1.0/16.0)*1105*pow(x,9);
			else if(power==7) result=(1.0/16.0)*(-2145)*pow(x,7);
			else if(power==5) result=(1.0/16.0)*1287*pow(x,5);
			else if(power==3) result=(1.0/16.0)*(-231)*pow(x,3);
			else              result=0;
		}
		else{
			result=(1.0/16.0)*( 1105*pow(x,9)-2145*pow(x,7)+1287*pow(x,5)-231*pow(x,3) );
		}
		break;
	case 10:
		if(power>=0){
			if     (power==10) result=(1.0/128.0)*14212*pow(x,10);
			else if(power==8)  result=(1.0/128.0)*(-29835)*pow(x,8);
			else if(power==6)  result=(1.0/128.0)*20020*pow(x,6);
			else if(power==4)  result=(1.0/128.0)*(-4290)*pow(x,4);
			else if(power==0)  result=(1.0/128.0)*21;
			else               result=0;
		}
		else{
			result=(1.0/128.0)*( 14212*pow(x,10)-29835*pow(x,8)+20020*pow(x,6)-4290*pow(x,4)+21 );
		}
		break;
	default:
		result=0;
		break;
	}

	return result*sqrt( (2.0*n+1.0)/2.0 );																																								
}

double cModLegendre::P1(int n,double x){
	// èCê≥LegendreëΩçÄéÆ Pn(x) ÇÃ1äKî˜ï™
	double result;

	switch(n){
	case 0:
		result=0;
		break;
	case 3:
		result=3*x*x;
		break;
	case 4:
		result=(1.0/4.0)*( 20*pow(x,3) );
		break;
	case 5:
		result=(1.0/2.0)*( 45*pow(x,4)-21*x*x );
		break;
	case 6:
		result=(1.0/24.0)*( 924*pow(x,5)-540*pow(x,3) );
		break;
	case 7:
		result=(1.0/8.0)*( 1001*pow(x,6)-990*pow(x,4)+189*x*x );
		break;
	case 8:
		result=(1.0/192.0)*( 42120*pow(x,7)-48048*pow(x,5)+11880*pow(x,3) );
		break;
	case 9:
		result=(1.0/16.0)*( 9945*pow(x,8)-15015*pow(x,6)+6435*pow(x,4)-693*x*x );
		break;
	case 10:
		result=(1.0/128.0)*( 142120*pow(x,9)-238680*pow(x,7)+120120*pow(x,5)-17160*pow(x,3) );
		break;
	default:
		result=0;
		break;
	}

	return result*sqrt( (2.0*n+1.0)/2.0 );
}

double cModLegendre::P2(int n,double x){
	// èCê≥LegendreëΩçÄéÆ Pn(x) ÇÃ2äKî˜ï™
	double result;

	switch(n){
	case 0:
		result=0;
		break;
	case 3:
		result=6*x;
		break;
	case 4:
		result=(1.0/4.0)*( 60*x*x );
		break;
	case 5:
		result=(1.0/2.0)*( 180*pow(x,3)-42*x );
		break;
	case 6:
		result=(1.0/24.0)*( 4620*pow(x,4)-1620*x*x );
		break;
	case 7:
		result=(1.0/8.0)*( 6006*pow(x,5)-3960*pow(x,3)+378*x );
		break;
	case 8:
		result=(1.0/192.0)*( 294840*pow(x,6)-240240*pow(x,4)+35640*x*x );
		break;
	case 9:
		result=(1.0/16.0)*( 79560*pow(x,7)-90090*pow(x,5)+25740*pow(x,3)-1386*x );
		break;
	case 10:
		result=(1.0/128.0)*( 1279080*pow(x,8)-1670760*pow(x,6)+600600*pow(x,4)-51480*x*x );
		break;
	default:
		result=0;
		break;
	}

	return result*sqrt( (2.0*n+1.0)/2.0 );
}

cModLegendre::cModLegendre(){
	int i,j;

	R0=0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j) a[i][j]=0;
}

cModLegendre::cModLegendre(const cModLegendre& x){
	int i,j;

	R0=x.R0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j) a[i][j]=x.a[i][j];
}

cModLegendre& cModLegendre::operator=(const cModLegendre& x){
	int i,j;

	R0=x.R0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j) a[i][j]=x.a[i][j];
	return *this;
}

double& cModLegendre::C(int i,int j){
	static double err;
	
	if(i>=0 && j>=0 && i+j<=N){
		return a[i][j];
	}
	else{
		return err;
	}
}

void cModLegendre::Clear(){
	int i,j;

	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j) a[i][j]=0;
}


int cModLegendre::NotNull(){
	int i,j;
	
	if(R0!=0) return true;
	for(i=0; i<=N; ++i) for(j=0; i+j<=N; ++j){
		if(C(i,j)!=0) return true;
	}
	return false;
}

cModLegendre cModLegendre::reversed(){
	int i,j;
	cModLegendre x=*this;

	for(i=0; i<=N; ++i) for(j=0; i+j<=N; ++j){
		x.C(i,j)=-x.C(i,j);
	}
	return x;
}

void cModLegendre::scale(double m){
	int i,j;
	if(m==0) return;

	R0*=m;
	for(i=0; i<=N; ++i) for(j=0; i+j<=N; ++j){
		C(i,j)*=m;
	}
}

double cModLegendre::Z(double x,double y,int xpower/*=-1*/,int ypower/*=-1*/){
	int i,j;
	double z=0;

	if(R0!=0){ x/=R0; y/=R0; } else return 0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j){
		if(i+j<=N && a[i][j]!=0) z+=a[i][j]*P(i,x,xpower)*P(j,y,ypower);
	}
	return z;
}

double cModLegendre::Zx(double x,double y){
	int i,j;
	double z=0;

	if(R0!=0){ x/=R0; y/=R0; } else return 0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j){
		if(i+j<=N && a[i][j]!=0) z+=a[i][j]*P1(i,x)*P(j,y);
	}
	return z/R0;    // Åyíçà”Åz R0Ç≈äÑÇÈÇÃÇñYÇÍÇ»Ç¢
	                //    xo=x/R0 ->  dZ/dx = (dZ/dxo)(dxo/dx) = (dZ/dxo)/R0
}

double cModLegendre::Zy(double x,double y){
	int i,j;
	double z=0;

	if(R0!=0){ x/=R0; y/=R0; } else return 0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j){
		if(i+j<=N && a[i][j]!=0) z+=a[i][j]*P(i,x)*P1(j,y);
	}
	return z/R0;
}

double cModLegendre::Zxx(double x,double y){
	int i,j;
	double z=0;

	if(R0!=0){ x/=R0; y/=R0; } else return 0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j){
		if(i+j<=N && a[i][j]!=0) z+=a[i][j]*P2(i,x)*P(j,y);
	}
	return z/R0/R0;
}

double cModLegendre::Zxy(double x,double y){
	int i,j;
	double z=0;

	if(R0!=0){ x/=R0; y/=R0; } else return 0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j){
		if(i+j<=N && a[i][j]!=0) z+=a[i][j]*P1(i,x)*P1(j,y);
	}
	return z/R0/R0;
}

double cModLegendre::Zyy(double x,double y){
	int i,j;
	double z=0;

	if(R0!=0){ x/=R0; y/=R0; } else return 0;
	for(i=0; i<=N; ++i) for(j=0; j<=N; ++j){
		if(i+j<=N && a[i][j]!=0) z+=a[i][j]*P(i,x)*P2(j,y);
	}
	return z/R0/R0;
}

std::string cModLegendre::GetCoefficients(){
	int i,j;
	char buf[1000];
	std::string s;

	if(NotNull()){
		sprintf(buf,"R0 %.15g; ", R0); s+=buf;
		for(i=0; i<=N; ++i) for(j=0; i+j<=N; ++j){
			if(C(i,j)!=0){
				sprintf(buf,"C %d %d %.15g; ", i,j,C(i,j));
				s+=buf;
			}
		}
	}

	return s;
}

void cModLegendre::SetCoefficients(std::string com){
	// ó·ÅF com="R0 1.5; C 2 0 0.1; C 1 1 0.1"

	std::string s,s0;
	int i,j;
	std::string s1,s2,s3;
	bool b1,b2,b3;

	for(i=0; i<=N; ++i) for(j=0; i+j<=N; ++j) C(i,j)=0;

	for(i=1; i<=sentences(com); ++i){
		s=sentence(com,i);
		s0=arg(s,0);
		if(s0=="r0" || s0=="R0"){
			double val;
			s1=arg(s,1); b1=is_numeric(s1);
			if(s1=="?"){
				s+="R0 val\n";
			}
			else if(b1){
				val=atof(s1.c_str());
				R0=val;
			}
		}
		if(s0=="C" || s0=="c"){
			int i,j;
			double val;
			s1=arg(s,1); b1=is_numeric(s1);
			s2=arg(s,2); b2=is_numeric(s2);
			s3=arg(s,3); b3=is_numeric(s3);
			if(s1=="?"){
				s+="C i j val\n";
			}
			else if(b1 && b2 && b3){
				i=atoi(s1.c_str());
				j=atoi(s2.c_str());
				val=atof(s3.c_str());
				C(i,j)=val;
			}
		}
	}
}

