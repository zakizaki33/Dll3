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
	// ndの小数点以下とアッベ数によるコード
	// 例えばTiO2(nd=2.5)とnd=1.5のガラスのコードの最初の3桁は同じとなる．
	// しかし現存の光学ガラスに限ればコードとndの対応は一意である．
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
	// 波長をnm単位で返す
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
	// φx=phi_x,φy=phi_yの矩形または楕円と
	// 直線	y=tan(th_rad)x との交点の一方を返す．
	//
	// th_radが±PI/2の近くでも問題ないようである．
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
	// 球面計で球面曲率半径を測定する．
	//     D    = 測定ヘッドの直径
	//     H    = 測定面に当てたときのダイヤルゲージの読み(凸面，凹面共に正とする）
	//     Rref = 標準面の曲率半径（0のときは，標準面による補正は行わない． 凸面，凹面共に正とする）
	//     Href = 標準面にあてたときのダイヤルゲージの読み
	double r,r1,r2, h_err;
	double dh=0.01;  // 例：ダイヤルゲージの最小目盛
	char buf[1000];

	r=(H*H+D*D/4)/(H*2);
	if(Rref>0){
		h_err=Href - (D*D/4)/Rref/(1+sqrt(1-(D*D/4)/Rref/Rref));
		H-=h_err;
		r=(H*H+D*D/4)/(H*2);
	}
	
	r1= ((H-dh)*(H-dh)+D*D/4)/((H-dh)*2);
	r2= ((H+dh)*(H+dh)+D*D/4)/((H+dh)*2);
	sprintf(buf,"r=%g (%g～%g for dh=±%g)\n",r,r1,r2,dh);
	return buf;
}

std::string cOptics::Plate(std::string glass,std::string color,double thickness,double incident_angle,
                           double dndtemp,double alpha,std::string what){
	std::string s;
	double wl,n,t,th,th1, a,b,x, ds,dp,dsp,dave,shift,dz,shift01,shift1,dnl,dtanth1,l,nl;
	double dwldf,dtdf,dthdf,dtempdf;
	double kp,ks;
	char buf[1000];
	const double DNAIRDTEMP=-0.94e-6; // 20℃における空気屈折率の温度係数（オハラカタログより）
	
	wl=Wavelength(color);
	n=Re(Index(glass,wl));
	t=thickness;
	th=incident_angle*PI/180;
	x=sin(th)/n; 
	th1=atan( x/sqrt(-x*x+1) );
	a=t/cos(th1)*( cos(th-th1)-1/n );
	b=t/n/cos(th1)*( 1-(cos(th)/cos(th1))*(cos(th)/cos(th1)) );

	sprintf(buf,"硝材= %s\n"                     ,glass.c_str());                         s+=buf;
	sprintf(buf,"波長= %s\n"                     ,color.c_str());                         s+=buf;
	sprintf(buf,"板厚 t= %g\n"                   ,t);                                     s+=buf;
	sprintf(buf,"入射角 th= %g°\n"              ,incident_angle);                        s+=buf;
	sprintf(buf,"  θ' th1= %.6f°\n"                     ,th1*180/PI);                   s+=buf;
	sprintf(buf,"  S面内のび ds= %.6f\n"                  ,ds=a);                         s+=buf;
	sprintf(buf,"  P面内のび dp= %.6f\n"                  ,dp=a+b);                       s+=buf;
	sprintf(buf,"  非点隔差P-S dsp= %.6f\n"               ,dsp=b);                        s+=buf;
	sprintf(buf,"  SPのび平均 dave= %.6f\n"               ,dave=a+b/2);                   s+=buf;
	sprintf(buf,"  光軸ずれ 入射～0次透過  shift= %.6f\n" ,shift=t/cos(th1)*sin(th-th1)); s+=buf;
	sprintf(buf,"  Δz 入射点～0次透過点= %.6f\n"         ,dz=t/cos(th1)*cos(th-th1));    s+=buf;
	sprintf(buf,"  光軸ずれ 0次～1次     shift01= %.6f\n" ,shift01=2*t*tan(th1)*cos(th)); s+=buf;
	sprintf(buf,"  光軸ずれ 入射～1次透過 shift1= %.6f\n" ,shift1=shift01-shift );        s+=buf;
	sprintf(buf,"  光路長差 0次～1次透過反射 dnl=%.6f\n"  ,dnl=2*n*t*cos(th1));           s+=buf;
	sprintf(buf,"  板厚 x tanθ' dtanth1= %.6f\n"         ,dtanth1=t*tan(th1));           s+=buf;
	sprintf(buf,"  板内光路 実寸    l= %.6f\n"            ,l=t/cos(th1));                 s+=buf;
	sprintf(buf,"  板内光路 光路長 nl= %.6f\n"            ,nl=n*t/cos(th1));              s+=buf;
	
	// “幾何光学 (三宅和夫，共立出版)”の2章の最後の図(図2.34)は，
	//  入射光，出射光の方向単位ベクトルをs,s'とするとき，入射点を始点とするとns,n's'の終点は，
	//  入射点を中心とする半径n,n'の球上にあり，n's'-nsは面法線pと平行である，ことを示している．
	//  これよりn=1，n'=nとして図的に考えると，
	//  変位がP面内のとき，
	//      ncos(th1)Δ(th1) = cosθ(th1)Δ(th1)  ...(1)
	//  変位がS面内のとき，
	//      nΔ(th1) = Δ(th)  ...(2)
	//  となる．
	// 平行度をαとする．
	// まずP面内では，
	// 0次1次反射のなす角δの場合は図を描いて考えると表面において，
	//     1次 Δ(th1)=2α
	// だからこれを(1)へ代入し，Δ(th)がδとなる．
	// また，0次1次透過のなす角δの場合は裏面において，
	//     0次 Δ(th1)=α
	//     1次 Δ(th1)=3α
	// だから，これらを(1)に代入し，それぞれのΔ(th)の差がδとなる．
	// 次にS面内では，
	// 0次1次反射のなす角δの場合は
	//     1次 Δ(th1)=2αcos(th1)   ( <- cos(th1)がかかることに注意．上記図で考えると分かり易い．)
	// だからこれを(2)へ代入し，Δ(th)がδとなる．
	// また，0次1次透過のなす角δの場合は裏面において，
	//     0次 Δ(th1)=αcos(th1)
	//     1次 Δ(th1)=3αcos(th1)
	// だから，これらを(2)に代入し，それぞれのΔ(th)の差がδとなる．
	// レーザを使った平行平面板の平行度測定(反射光を遠方に投影し分離を測定)に有用．
	//（実際の測定では誤差要因となる"光軸ずれ 0次～1次"に注意．)
	sprintf(buf,"  一次近似\n");                                                        s+=buf;
	sprintf(buf,"    0次1次反射または透過のなす角= k x 平行度(表裏面法線のなす角)\n");  s+=buf;
	sprintf(buf,"      平行度の変位がP面内 kp= %g\n", kp=2*n*cos(th1)/cos(th));         s+=buf;
	sprintf(buf,"      平行度の変位がS面内 ks= %g\n", ks=2*n*cos(th1));                 s+=buf;
	
	// fringes = 2*n*thickness*cos(th1)/λo より．
	dwldf=wl*wl/2/n/(t*1000000)/cos(th1);
	dtdf=(wl/1000000)/2/n/cos(th1);
	dthdf= th1==0 ? 0 : (wl/1000000)*cos(th1)/2/t/sin(th1)/cos(th)*180/PI;
	if(dndtemp==0 || alpha==0){
		dtempdf=0;
	}
	else{
		// dndtemp = 相対屈折率の温度係数
		// alpha   = 線膨張係数
		dtempdf= 1/( 2*t*cos(th1)/(wl/1000000) )/( dndtemp+n*DNAIRDTEMP+n*alpha );
	}
	
	sprintf(buf,"  Δλ[nm]/fringe  dwldf= %g\n" ,dwldf  ); s+=buf;
	sprintf(buf,"  Δt [mm]/fringe   dtdf= %g\n" ,dtdf   ); s+=buf;
	sprintf(buf,"  Δth[deg]/fringe dthdf= %g\n" ,dthdf  ); s+=buf;
	sprintf(buf,"  ΔT [K]/firnge dtempdf= %g\n" ,dtempdf); s+=buf;
	
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
	// 屈折偏角の計算
	// glass,glass1 : 入射側，出射側の媒質
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

	sprintf(buf,"硝材= %s\n"        ,glass.c_str());       s+=buf;
	sprintf(buf,"波長= %s\n"        ,color.c_str());       s+=buf;
	sprintf(buf,"頂角= %g\n"        ,apex_angle);          s+=buf;
	sprintf(buf,"入射角= %g°\n"    ,incident_angle);      s+=buf;
	sprintf(buf,"  出射角= %g°\n"  ,th3*180/PI);          s+=buf;  
	sprintf(buf,"  ふれ角= %g°\n"  ,(th+th3-al)*180/PI ); s+=buf;       
	sprintf(buf,"  Δp= %g\n"       ,100*tan(th+th3-al) ); s+=buf;
	sprintf(buf,"  角倍率= %g\n"    ,gamma1*gamma2);       s+=buf;
	sprintf(buf,"  横倍率= %g\n"    ,1/gamma1/gamma2);     s+=buf;
	return s;
}

double cOptics::MeasuredDistortion(double x1,double y1,double x2,double y2,double x3,double y3){
	// 直線 (x1,y1)-(x3,y3) と 点(x2,y2) があるとき，
	// 直線と点の距離を直線の長さで割った値を返す．
	// この値は，特に，
	//   (x1,y1),(x3,y3),(x2,y2)が光軸を中心とする正方形(2a x 2a)の辺の像の両端及び中点で，
	//   光学系の歪曲が3次収差のみ(光軸付近)
	// のときは，h=aにおける歪曲 (y'-y'0)/y'0 のほぼ半分となる．
	point p1(x1,y1), p2(x2,y2), p3(x3,y3);
	return distance(p2,line(p1,p3))/distance(p1,p3);
}

double cOptics::BerekDOF(std::string &s,double wl_nm,double NAo,double TotalM,double VAmin){
	// Berekの式による物体側DOFをμm単位で返す
	// wl_nm  = nm単位波長
	// NAo    = 物体側NA
	// TotalM = 総合倍率
	// VAmin  = 分単位で表した眼の分解能
	double dof1,dof2;
	char buf[1000];
	dof1=(wl_nm/1000)/2/NAo/NAo;                                        // 波動光学項
	dof2= TotalM==0 ? 0 : tan(VAmin/60*PI/180)*(250*1000)/TotalM/NAo;   // 幾何光学項
	// 幾何光学項は“幾何光学（三宅）”(5.2)あるいは(5.3)式で，
	//     ε= (250mm/TotalM)*VAmin) = 焦点距離*VAmin
	// としたものに一致．
	sprintf(buf,"DOF= λ/(2*NAo^2) + (1/NAo)*(250mm/TotalM)*VA\n" ); s=buf;
	sprintf(buf,"λ(nm)=%g\n"      , wl_nm);     s+=buf;
	sprintf(buf,"NAo=%g\n"         , NAo);       s+=buf;
	sprintf(buf,"TotalM=%g\n"      , TotalM);    s+=buf;
	sprintf(buf,"VA(min)=%g\n"     , VAmin);     s+=buf;
	sprintf(buf," DOF(μm)=%g\n" , dof1+dof2);   s+=buf;
	sprintf(buf,"  波動光学項(μm)=%g\n",dof1);  s+=buf;
	sprintf(buf,"  幾何光学項(μm)=%g\n",dof2);  s+=buf;
	return dof1+dof2;
}

double cOptics::ScheimpflugTH1(double ms,double th_deg,double th_l_deg){
	// ms       = s方向倍率
	// th_deg   = 物体面法線と光軸のなす角(deg) 注：Scheimpflugカメラの仰角は90-th_deg
	// th_deg_l = 光軸とレンズ光軸のなす角(deg)
	//
	// 戻り値   = 物体面，レンズ，像面がScheimpflug条件を満たすときの
	//          = 像面法線と光軸のなす角(deg)
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
	sprintf(buf,"s方向倍率ms     = %g\n" ,ms);          s+=buf;
	sprintf(buf,"物体面入射角th  = %gdeg\n" ,th_deg);   s+=buf;
	sprintf(buf,"レンズ入射角th_l= %gdeg\n" ,th_l_deg); s+=buf;
	sprintf(buf,"  像面入射角th1 = %gdeg\n" ,th1_deg);  s+=buf;
	sprintf(buf,"  p方向倍率mp   = %g\n"    ,ms*k);     s+=buf;
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
	//  ガウス型パルスが群遅延分散を受けた後のパルス幅Δtout(fs)を返す．
	//      dt0_fs  : 最初のパルス幅Δt(fs)
	//      GDD_fs2 : 群遅延分散GDD(fs^2)
	// 
	//  Δtout = √(Δt^4+16*(ln2)^2*GDD^2) / Δt (http://www.newport-japan.jp/pdf/0321.pdf)

	return sqrt( dt0_fs*dt0_fs*dt0_fs*dt0_fs + 7.68725*GDD_fs2*GDD_fs2 )/dt0_fs;
}

double cOptics::ODToT(double od){
	// 光学濃度に対する透過率を求める
	return pow(10,-od);
}

double cOptics::CoherenceLength(double wl_nm,double dwl_fwhm_nm){
	// 波長分布が半値幅Δλのガウス分布のとき，
	//    コヒーレンス長 = (4ln2/π)(λo^2/Δλ)
	// 【注意】OCTの分解能は光の往復によりさらに2で割る
	return 0.8825424*wl_nm*wl_nm/dwl_fwhm_nm /1000000;  // コヒーレンス長(mm単位)
}

double cOptics::Derivative(std::string command,double *x,double dx){
	std::ofstream to;
	int i, nc=0;
	double x0, y0,y1,y2,y11,y22, a,a_pre, result;
	const int MakeLog=0;
	const int MAXTIMES=10;
	const int METHOD_OLD=0; // METHOD_OLD=true では旧アルゴリズムによる．この方がなぜか速く収束する場合がある．
	                        // METHOD_OLDをtrueにした場合の変更点は以下．
							//  ・5点公式は使わない
							//  ・終了条件 3,4 は使わない
	clock_t start;
	// 例えば cyl値(マイナス表示)は0で折り返しがある．
	// *xが折り返し点の近傍にあり，dxが大きいと *x+dxと*x-dxでの勾配は符号が逆であり，
	// dxが小さくなって勾配が同符号になるまで微分係数が求まらない．
	// 実験ではMAXTIME=5では不足で，自動設計の収束が遅くなることがあった．

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
				a=(a1+a2)/2;  // これを a=(y1-y2)/(dx*2) と書くとなぜか収束が遅くなる場合あり（桁落ちの影響？）
				
				if(MakeLog){
					to.precision(16);
					to << "  dy1=" << dy1 << " dy2=" << dy2 << std::endl;
					to.precision(6);
					to << " a=" << a << std::endl;
				}

				if(a==0){result=0; nc=0; break;}
				if(fabs((a1-a2)/a)<0.01){result=a; nc=1; break;}     // 左側右側の傾きの差による．aの精度はおそらく1%以内 (終了条件1)
				if(i>=2 && (fabs((a-a_pre)/a)<0.01)){result=a; nc=2; break;} // 前回との差による．aの精度はおそらく1%以内 (終了条件2)

				// 終了条件1 だけでもよいが，条件2 も加えることで，多少速くなる (2016.06.29)
				//  (条件2 のみでは i=1 でループを抜けることができない（2回以上のループ要）ので条件1 も必要)

				a_pre=result=a;
				dx*=0.1;
			}
			else{
				if(i==1){
					a=(y1-y2)/(dx*2);  // i=1回目は3点公式（xの4箇所でのy値計算が必要な5点公式を最初から使わない）
				}
				else{
					a=(-y11+8*y1-8*y2+y22)/(dx*12);  // i=2回目以降は5点公式（y11,y22は前回の値を利用）
				}
				
				if(MakeLog){
					to.precision(16);
					to << "  y1-y0=" << y1-y0 << " y0-y2=" << y0-y2 << std::endl; 
					to.precision(6);
					to << " a=" << a << std::endl;
				}

				if(a==0){result=0; nc=0; break;}
				if(fabs((y1+y2-y0*2)/dx/a)<0.01){result=a; nc=1; break;}     // 左側右側の傾きの差による．aの精度はおそらく1%以内 (終了条件1)
				if(i>=2 && (fabs((a-a_pre)/a)<0.01)){result=a; nc=2; break;} // 前回との差による．aの精度はおそらく1%以内 (終了条件2)

				if(y0!=0){
					if( (fabs((y1-y0)/y0)<1e-13) || (fabs((y0-y2)/y0)<1e-13) ) { // yの変化量の有効桁がなくなっているかもしれない (終了条件3)
						result= i==1 ? a : a_pre; nc=3; break; 
					} 
				}

				if(x0!=0){
					if( fabs(dx/x0)<1e-13 ){   // x±dxの有効桁がなくなっているかもしれない (終了条件4)
						result= i==1 ? a : a_pre; nc=4; break;
					}
				}

				// 終了条件1 だけでもよいが，条件2 も加えることで，多少速くなる (2016.06.29)
				//  (条件2 のみでは i=1 でループを抜けることができない（2回以上のループ要）ので条件1 も必要)
				// 終了条件3,4を追加（2018.11.01)

				a_pre=result=a;
				y11=y1;
				y22=y2;
				dx*=0.5;  // 5点公式で前回の値を利用（y11=y1, y22=y2）しているので新しいdxは *0.5 である必要がある
			}
		}

		*x=x0; // push(), pop() でも保存，復元ができるが，実行速度が非常に遅くなる．
		       // もし，push(),pop()で行なうときは virtual int cOptics::push(){ /* 空の関数 */ }
		       // を追加する．"virtual" を付けない場合，cLens型インスタンスでcOptics::Derivative()を呼んでも
		       // cOptics::Derivative関数の中のpush()はcOptics::push() (空の関数) である．
		       // virtualがあれば cLens::push()が呼ばれる．

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
	// sの書式 : com [sign | ± tol] target weight を解釈して，com,sgn,tol,target,weightを設定する．
	// 戻り値のweightにset_weightを乗じる（関数 preprocess_targets で使う）．
	com=arg(sentence,0);
	sign=trim(arg(sentence,1),1);

	if(sign=="=" || sign=="<" || sign==">"){
		target=atof(arg(sentence,2).c_str());
		weight=atof(arg(sentence,3).c_str());
		sentence=replace_arg(sentence,3,str(weight*set_weight));
	}
	else if(sign=="+-" || sign=="±"){
		// エクセルのセルに"+- ..." とするとエラーが出て煩わしいので
		// "±"も使えるようにした．
		tol=atof(arg(sentence,2).c_str());
		target=atof(arg(sentence,3).c_str());
		weight=atof(arg(sentence,4).c_str());
		sentence=replace_arg(sentence,4,str(weight*set_weight));
	}
	else{
		// "=","<",">","+-" が省略された場合
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
	// 本関数がないと，
	// cOptics::arg(const std::string& sentence, int n)
	// があるため，complex.hで宣言される arg(const complex& z)
	// を呼ぶのに演算子::が必要(::arg(z))になってしまう．

	return ::arg(z);
}

std::string cOptics::arg(const std::string& sentence, int n){
	std::string s;

	s=::arg(sentence,n);
	symboltoval(s);
	return s;
}

bool cOptics::is_numeric(std::string& s){
	//  sまたはs=cmd(s,1)が数を表す文字列のときtrueを返す．
	std::string buf;

	if(s==""){
		return false;
	}
	else if(::is_numeric(s)){   // string.h
		return true;
	}
	else{
		push();   //  仮想関数．子クラス(cLensなど)のpushを呼ぶ．
		buf=cmd(s,1);  // 【注意】 第2引数が真なので，戻り値が非数字の文字列であるコマンドは buf="" になってしまう．
		pop();    //  仮想関数．子クラス(cLensなど)のpopを呼ぶ．
		//【重要】push(),pop()により，*thisは変化させない． 2016.03.14
		//   例：
		//   s="Optimize (r 8 .0002) (f 20 10) 10  // r(8)を変化させて，fを20mmにする．
		//   をcOptics::cmdの引き数とすると
		//   cLnes1::properyがsに対して実行される．
		//   この中で is_numeric("r 8 .0002") が実行されるが，
		//   上記で ::is_numeric(buf=cmd(s,1))  としてしまうと，
		//   本来の "r(8)を.0002づつ変化させて最適化する" という意図に反して，
		//   r(8)に0.0002が代入されてしまう．
		if(::is_numeric(buf)){
			s=buf;
			return true;
		}
		else{
			// --- 戻り値が非数字の文字列であるコマンドも展開する ------------------
			buf=cmd(s,0);
			if(buf!="") s=buf;
			// ---------------------------------------------------------------------
			return false;
		}
	}
}

std::string cOptics::cmd(const std::string& command,int val){
	// val=tureのとき，一部のコマンドはBasic側にてval関数で処理するために
	// 数値を表す文字列のみ返す．
	// commandが複数のsentenceからなるときは，
	//     val=false : 処理結果を文字列として加えて返す．
	//     val=true  : 最後のsentence(結果が""であるsentenceは除く)
	//                 の結果数値を表す文字列を返す．
	//                 (以前は，各sentenceの結果の合計を返していたが，
	//                  これは "+" コマンドで可能であるし，
	//                  最後の結果だけ知りたい場合もあるのでこのようにした． 2014.10.08 )
	//
	//                 系の変更が起こらないように，commandの最後に pop が必要なことがある．
	//                 pop の戻り値は""なので必要な数値を返すために，
	//                 上記“結果が""であるsentenceは除く”こととした． 2016.05.11
	std::string s,ss,sv,com;
	double v1,v2,v3;
	char buf[1000];
	int i,n;

	s=ss="";
	n=sentences(command);

	for(i=1; i<=n; ++i){
		com=sentence(command,i);
		symboltoval(com);
		// scmd()は文字列comを単一sentenceとして処理する
		s=scmd(com,val);  // 仮想関数．子クラスのscmdが呼ばれる．
		ss+=s;
		if(s!="") sv=s;   // val=trueの場合の戻り値は，最後の""でない s による．
		s="";
	}

	if(val){
		switch(words(sv)){
		case 2:  // 2次元ベクトルを想定
			v1=atof(word(sv,1,1).c_str());
			v2=atof(word(sv,2,1).c_str());
			sprintf(buf,"%.15g %.15g\n",v1,v2);
			break;
		case 3:  // 3次元ベクトルを想定
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
		// 注意：書式%fは小数点以下の桁数が固定なので，値が小さいと有効桁がなくなる．
		//       %.15gがよい．
		//       単に%gだと有効桁はデフォルトの6桁となり，自動設計の微分係数計算で
		//       差分をとるとき有効桁が不足気味になる．
		//       scmd関数でもこのようにすること．
		return ::is_numeric(ss)  ? buf : "";  // ::is_numeric(ss)でないときは"0"でななく""を返す．
		                                      // こうしないと，cOptics::is_numeric(command)が
		                                      // 必ずtrueとなってしまう．
	}
	else{
		return ss;
	}
}

std::string cOptics::scmd(std::string com,int val){
	// val=tureのとき，一部のコマンドはBasic側にてval関数で処理するために
	// 数値を表す文字列のみ返す．

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
		// s0の冒頭が数を表すとき
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
		// ファイルの内容をコマンドとして実行する．拡張子 .txt は省略可能．
		// ファイル内容の"arg1","arg2", ... は引数arg1,arg2, ... に置換される．
		std::ifstream from;
		std::string tmp,com1;
		int i;
		const int ARGS_MAX=9;            // 引数は9個まで
		std::string a[ARGS_MAX+1];

		s1=arg(com,1);
		if(s1=="?"){
			s+="Call filename [arg1 arg2 ... ]\n";
			s+="   - extension '.txt' can be dropped from filename\n";
			s+="   - max number of arguments is 9\n";
		}
		else{
			from.open(s1.c_str());
			if(!from){                          // openが失敗した場合
				from.clear();                   // 注意：clear()をしないと，次行のopenは成功しない．
				from.open((s1+".txt").c_str()); // 拡張子 ".txt" を付加して再度openを試みる．
			}
			for(i=1; i<=ARGS_MAX; ++i){
				a[i]=arg(com,1+i);
				if( (tmp=cmd(a[i],1))!="" ) a[i]=tmp;
			}
			if(from){
				while(from>>tmp){
					for(i=1; i<=ARGS_MAX; ++i){
						if(a[i]!="") tmp=replace(tmp,"arg"+trim(str(i),0),a[i]);  // 置換
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
		// 注釈文
		s1=arg(com,1);
		if(s1=="?"){
			s="Rem [comment]\n";
		}
		return s;
	}
	if(s0=="Set" || s0=="set"){
		// 変数symbolに値cmd(value)を代入する（文字列valueの実行結果 cmd(value)を代入）
		// 例：
		//   r 1 10;
		//   set a f;  aにはcmd(f)が代入される．
		//   r 1 11;
		//   a;        r(1)=10での焦点距離を返す．("setstr a f" としたときは r(1)=11での焦点距離が返される）
		int i;
		std::string tmp;
		s1=::arg(com,1);  // 注意： "::"を付けないとs1がcOptics::arg()により書き換えられてしまう．
		s2=::arg(com,2);
		if(s1=="?"){
			s+="Set symbol value\n";
		}
		else{
			if(s1=="" ){
				// 現在の設定の一覧を表示
				for(i=1; i<=symbol.GetSize(); i++){
					s+=symbol[i]+'\t'+value[i]+'\n';
				}
			}
			else{
				tmp=cmd(s2,1);  // "set x (* x 2)" のように再帰呼出があるときは，
				                // Removeの前に "* x 2" を実行する必要あり
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
		// 変数symbolに文字列strを代入する
		// 例：
		//   r 1 10;
		//   setstr a f;  a=scmd("f")でないことに注意．aは"f"で，r(1)=10での焦点距離値ではない．
		//   r 1 11;
		//   a;        r(1)=11での焦点距離を返す．
		//
		// 注意：
		//   setstr x (* x 2) のように再帰呼び出しがあると無限ループとなりエラー．

		int i;
		s1=::arg(com,1);  // 注意： "::"を付けないとs1がcOptics::arg()により書き換えられてしまう．
		s2=::arg(com,2);
		if(s1=="?"){
			s+="SetStr symbol str\n";
		}
		else{
			if(s1=="" ){
				// 現在の設定の一覧を表示
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
			s+="    H    = Δh with measuring surface\n";
			s+="    Rref = curvature radius of reference surf(if zero, not compensate the result)\n";
			s+="    Href = Δh with reference surface\n";
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
		// is_numericを作用させないと，例えば，"wavelength (color 1)" が期待する動作をしない．
		// ここで， "color 1" は cLens1のコマンドで第1波長を表す文字列を返す．
		// ただし，is_numeric() の仕様により，"color 1" の戻り値が "632.8" のように数字の場合は問題ないが，
		// "e" のように非数字の場合は "wavelength (color 1)" は 0 になってしまう．
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



