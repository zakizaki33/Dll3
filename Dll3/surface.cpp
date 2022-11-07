#include "stdafx.h"
#include "surface.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

int surface::file_ver=0;

void surface::match_r_c(){
	// rとcを整合させる．

	// 例えば，
	//     double *a; x=&c(); または，double &a; x=c();
	// のaでccを操作した後，
	//     double *b; b=&r(); または，double &b; b=r();
	// のbにより取得したrrはccと整合されていない．
	// （この場合でも，r()(r()は本関数を呼ぶ)によって取得したrrは整合されている．）

	// rとcの整合については，rまたはcのみをメンバとして持つことも考えられるが次の問題がある．
	//   例えば，メンバとしてはcのみとし，
	//       GetR(){ return 1/c; }
	//       SetR(double r){ c=1/r; }
	//   とすれば，本関数が不要となり単純にはなる．
	//   しかし，rは設計値として小数以下二桁までで表せる値を持ちたい場合が多く，
	//   r=1/c の計算誤差により，小数以下三桁以降が0でなくなる可能性がある．
	//   逆にメンバをrのみとすると，自動設計で有用なcへの参照やポインタが得られなくなる．

	if( rr!=r0 || cc!=c0 ){
		
		if(rr!=r0 && cc==c0){
			cc= rr==0 ? 0 : 1/rr;
		}
		else if(rr==r0 && cc!=c0){
			rr= cc==0 ? 0 : 1/cc;
		}
		else if(rr!=r0 && cc!=c0){
			// 残る可能性，rr!=r0 && cc!=c0 は，本関数が呼ばれない間に
			// r(),c()の両方が呼ばれることはないから,
			// r(),c()を通してrr,ccにアクセスしている限りあり得ないが，
			// 参照やポインタで直接操作すれば起こりうる．
			// rr=cc=0 として警告する．
			rr=cc=0;
		}

		// r0==rr c0==cc であれば整合が保証される．
		r0=rr;
		c0=cc;
	}	
}



void surface::match_fi_pi(){
	// fideal()とpideal()を整合させる．考え方は match_r_c() と同じ．

	if( fi!=fi0 || pi!=pi0 ){
		
		if(fi!=fi0 && pi==pi0){
			pi= fi==0 ? 0 : 1/fi;
		}
		else if(fi==fi0 && pi!=pi0){
			fi= pi==0 ? 0 : 1/pi;
		}
		else if(fi!=fi0 && pi!=pi0){
			fi=pi=0;
		}

		fi0=fi;
		pi0=pi;
	}	
}

void surface::match_acoa_cacoa(){
	// aCOA()とcaCOA()を整合させる．考え方は match_r_c() と同じ．

	if( acoa!=acoa0 || cacoa!=cacoa0 ){
		
		if(acoa!=acoa0 && cacoa==cacoa0){
			cacoa= acoa==0 ? 0 : 1/acoa;
		}
		else if(acoa==acoa0 && cacoa!=cacoa0){
			acoa= cacoa==0 ? 0 : 1/cacoa;
		}
		else if(acoa!=acoa0 && cacoa!=cacoa0){
			acoa=cacoa=0;
		}

		acoa0=acoa;
		cacoa0=cacoa;
	}	
}

void surface::match_bcoa_cbcoa(){
	// bCOA()とcbCOA()を整合させる．考え方は match_r_c() と同じ．

	if( bcoa!=bcoa0 || cbcoa!=cbcoa0 ){
		
		if(bcoa!=bcoa0 && cbcoa==cbcoa0){
			cbcoa= bcoa==0 ? 0 : 1/bcoa;
		}
		else if(bcoa==bcoa0 && cbcoa!=cbcoa0){
			bcoa= cbcoa==0 ? 0 : 1/cbcoa;
		}
		else if(bcoa!=bcoa0 && cbcoa!=cbcoa0){
			bcoa=cbcoa=0;
		}

		bcoa0=bcoa;
		cbcoa0=cbcoa;
	}	
}

double& surface::r(){
	// r(),c()はrr,ccの参照を返すので値が変えられる場合があるが，
	// return rr(cc) の後では r==1/cの整合に関して何もできない．
	// したがって，次にr(),c()が参照を返す前に r==1/c の関係に整合しておかなけらばならない．
	// (本クラスのメンバであってもrr,ccの変更はr(),c()を使うことが望ましい)
	match_r_c();
	return rr;
}

double& surface::c(){
	match_r_c();
	return cc;
}

double& surface::NormH(){	
	if(normh==0) normh=1;  
	// この行がないと0除算エラーが発生することがある．値0は無意味なので失うものはない．
	// 例：
	//   Lens.xlsの "LensData" シートで手動で(Addコマンドを使わずに)レンズ面数を増やす．
	//   "Calc" ボタンを押すと，"Asph" シートの該当面のNormH値は0なのでnormhに0が読み込まれる．
	//   これらの面を非球面とすると0除算エラーが発生する．
	return normh;
}

double& surface::fideal(){
	match_fi_pi();
	return fi;
}

double& surface::pideal(){
	match_fi_pi();
	return pi;
}

double& surface::aCOA(){
	match_acoa_cacoa();
	return acoa;
}

double& surface::caCOA(){
	match_acoa_cacoa();
	return cacoa;
}

double& surface::bCOA(){
	match_bcoa_cbcoa();
	return bcoa;
}

double& surface::cbCOA(){
	match_bcoa_cbcoa();
	return cbcoa;
}

surface::surface()
{
	rr=cc=r0=c0=0;
	Newton=As0=As45=NewtonTol=AsTol=rVariable=0;
	asph_type=0; 
	kp=0;
	normh=1; a1=a2=a3=a4=a5=a6=a7=a8=a9=a10=a11=a12=a13=a14=a15=a16=a18=a20=0;
	zernike.SetR0(0);
	zernike.SetNormalize(1);
	zernike.SetIsFringeOrder(0);
	zernike.SetJBase(0);
	zernike.SetMaxOrder(6);
	legendre.R0=0;
	spline.StartPointDerivativeZero=1;   // 出発点は光軸上を想定するので，でサグの一次微分は0
	spline.EndPointDerivativeZero=0;
	Apv=sfx=sfy=0;
	UserDefSurf=0;
	cylinder=0;
	fresnel=0; rbase=0;
	ry=rx=kpy=kpx=0; IsXToroid=0;
	coneangle=0;
	fi=pi=fi0=pi0=0;
	acoa=cacoa=acoa0=cacoa0=0;
	bcoa=cbcoa=bcoa0=cbcoa0=0;
	tCOA=0;
	SA0=CM0=0;
	grating=difforder=0; gpitch=grx=gry=grz=0;
	Diffusion=0;
	Fish=0;
	EAtype=0; EAy=0; EAx=0; CHMy=CHMx=0; EAdy=EAdx=0;
	decenter_type=0;
	dx=dy=dz=rox=roy=roz=0;
	order=0;
	ret=0;
	dx1=dy1=dz1=rox1=roy1=roz1=0;
	order1=0;
	CoatName=""; CoatReverse=0;
	rem="";
}

surface::surface(const surface &x){
	*this=x;
}

surface& surface::operator=(const surface& x){
	rr=x.rr; cc=x.cc; r0=x.r0; c0=x.c0;
	match_r_c();  // 代入時に match_r_c() を実行する(このために代入演算子を定義した)．
	              // こうすれば代入直後(match_r_c()を実行していない状態)でも，
	              // ポインタ double *x=&c() や，参照 double &x=c() でccを取得できる．
	              // これは，自動設計に有用．
	zernike=x.zernike;
	legendre=x.legendre;

	Newton=x.Newton; As0=x.As0; As45=x.As45; 
	NewtonTol=x.NewtonTol; AsTol=x.AsTol;
	rVariable=x.rVariable;
	asph_type=x.asph_type;
	kp=x.kp; normh=x.normh;
	a1=x.a1; a2=x.a2; a3=x.a3; a4=x.a4; a5=x.a5; a6=x.a6; a7=x.a7; a8=x.a8; a9=x.a9;
	a10=x.a10; a11=x.a11; a12=x.a12; a13=x.a13; a14=x.a14; a15=x.a15; a16=x.a16; a18=x.a18; a20=x.a20;
	b=x.b;
	Dcon=x.Dcon;
	spline=x.spline;
	Apv=x.Apv; sfx=x.sfx; sfy=x.sfy;
	UserDefSurf=x.UserDefSurf;
	cylinder=x.cylinder;
	fresnel=x.fresnel; rbase=x.rbase;
	ry=x.ry; rx=x.rx; kpy=x.kpy; kpx=x.kpx; IsXToroid=x.IsXToroid;
	coneangle=x.coneangle;
	fi=x.fi; pi=x.pi; fi0=x.fi0; pi0=x.pi0; match_fi_pi();
	acoa=x.acoa; cacoa=x.cacoa; acoa0=x.acoa0; cacoa0=x.cacoa0; match_acoa_cacoa();
	bcoa=x.bcoa; cbcoa=x.cbcoa; bcoa0=x.bcoa0; cbcoa0=x.cbcoa0; match_bcoa_cbcoa();
	tCOA=x.tCOA;
	SA0=x.SA0; CM0=x.CM0;
	grating=x.grating; difforder=x.difforder; gpitch=x.gpitch; grx=x.grx; gry=x.gry; grz=x.grz;
	Diffusion=x.Diffusion;
	Fish=x.Fish;
	EAtype=x.EAtype; EAy=x.EAy; EAx=x.EAx; CHMy=x.CHMy; CHMx=x.CHMx; EAdy=x.EAdy; EAdx=x.EAdx;
	decenter_type=x.decenter_type;
	dx=x.dx;dy=x.dy;dz=x.dz; rox=x.rox;roy=x.roy;roz=x.roz; 
	order=x.order;
	ret=x.ret;
	dx1=x.dx1;dy1=x.dy1;dz1=x.dz1; rox1=x.rox1;roy1=x.roy1;roz1=x.roz1;
	order1=x.order1;
	CoatName=x.CoatName; CoatReverse=x.CoatReverse;
	rem=x.rem;

	return *this;
}

std::ostream& operator<<(std::ostream& to,surface& x){
	std::streamsize np;
	
	to<<x.r()<<' '<<x.NewtonTol<<' '<<x.AsTol<<' '<<x.rVariable<<std::endl;
	to<<x.asph_type<<std::endl;
	to<<x.kp<<std::endl;
	to<<x.NormH()<<std::endl;
	
	np=to.precision(16);
	// 非球面べき級数展開の単項 a(n)*h^n が大きくなることがあり，波長程度の精度を確保するには倍精度の精度が必要
	to<<x.a1<<' '<<x.a2<<' '<<x.a3<<std::endl;
	to<<x.a4<<' '<<x.a5<<' '<<x.a6<<' '<<x.a7<<' '<<x.a8<<' '<<x.a9<<std::endl;
	to<<x.a10<<' '<<x.a11<<' '<<x.a12<<' '<<x.a13<<' '<<x.a14<<' '<<x.a15<<' '<<x.a16<<' '<<x.a18<<' '<<x.a20<<std::endl;
	to.precision(np);
	
	to<<blank_filled(x.b.GetTerms())<<std::endl;
	to<<blank_filled(x.zernike.GetCoefficients())<<std::endl;
	to<<blank_filled(x.Dcon.GetTerms())<<std::endl;
	to<<blank_filled(x.legendre.GetCoefficients())<<std::endl;
	to<<blank_filled(x.spline.GetData())<<std::endl;
	to <<x.Apv<<' '<<x.sfx<<' '<<x.sfy<<std::endl;
	to<<x.UserDefSurf<<std::endl;
	to<<x.cylinder<<std::endl;
	to<<x.fresnel<<' '<<x.rbase<<std::endl;
	to<<x.ry<<' '<<x.rx<<' '<<x.kpy<<' '<<x.kpx<<' '<<x.IsXToroid<<std::endl;
	to<<x.coneangle<<std::endl;
	to<<x.fideal()<<' '<<x.aCOA()<<' '<<x.bCOA()<<' '<<x.tCOA<<' '<<x.SA0<<' '<<x.CM0<<std::endl;
	to<<x.grating<<' '<<x.difforder<<std::endl;
	to<<x.gpitch<<' '<<x.grx<<' '<<x.gry<<' '<<x.grz<<std::endl;
	to<<x.Diffusion<<std::endl;
	to<<x.Fish<< std::endl;
	to<<x.EAtype<<' '<<x.EAy<<' '<<x.EAx<<' '<<x.CHMy<<' '<<x.CHMx<<' '<<x.EAdy<<' '<<x.EAdx<<std::endl;
	to<<x.decenter_type<<std::endl;
	to<<x.dx<<' '<<x.dy<<' '<<x.dz<<' '<<x.rox<<' '<<x.roy<<' '<<x.roz<<' '<<std::endl;
	to<<x.order<<std::endl;
	to<<x.ret<<std::endl;
	to<<x.dx1<<' '<<x.dy1<<' '<<x.dz1<<' '<<x.rox1<<' '<<x.roy1<<' '<<x.roz1<<std::endl;
	to<<x.order1<<std::endl;
	to<<blank_filled(x.CoatName)<<' '<<x.CoatReverse<<std::endl;
	to<<blank_filled(x.rem)<<std::endl;
	return to;
}

std::istream& operator>>(std::istream& from,surface& x){
	std::string buf;

	switch(surface::file_ver){
		case 138:
		case 137: // EAdy,EAdxを追加 20.05.08
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6>>x.a7>>x.a8>>x.a9;
			from>>x.a10>>x.a11>>x.a12>>x.a13>>x.a14>>x.a15>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>buf; x.legendre.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.spline.SetData(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.aCOA()>>x.bCOA()>>x.tCOA; x.match_acoa_cacoa(); x.match_bcoa_cbcoa();
			from>>x.SA0>>x.CM0;
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx>>x.CHMy>>x.CHMx>>x.EAdy>>x.EAdx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 136: // a9,a11,a13,a15を追加 20.04.23
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6>>x.a7>>x.a8>>x.a9;
			from>>x.a10>>x.a11>>x.a12>>x.a13>>x.a14>>x.a15>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>buf; x.legendre.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.spline.SetData(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.aCOA()>>x.bCOA()>>x.tCOA; x.match_acoa_cacoa(); x.match_bcoa_cbcoa();
			from>>x.SA0>>x.CM0;
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx>>x.CHMy>>x.CHMx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 135: // splineを追加 20.01.23
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6>>x.a7>>x.a8;
			from>>x.a10>>x.a12>>x.a14>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>buf; x.legendre.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.spline.SetData(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.aCOA()>>x.bCOA()>>x.tCOA; x.match_acoa_cacoa(); x.match_bcoa_cbcoa();
			from>>x.SA0>>x.CM0;
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx>>x.CHMy>>x.CHMx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 134:  // cLens1にvar,var1,..var5を追加 190729
		case 133:  // cLens1でEPphiをEPDに改称，EPDxを追加 18.09.19
		case 132:  // aCOAとbCOAを入れ替え 18.08.26
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6>>x.a7>>x.a8;
			from>>x.a10>>x.a12>>x.a14>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>buf; x.legendre.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.aCOA()>>x.bCOA()>>x.tCOA; x.match_acoa_cacoa(); x.match_bcoa_cbcoa();
			from>>x.SA0>>x.CM0;
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx>>x.CHMy>>x.CHMx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 131:  // aCOA,bCOA,tCOAを追加 18.03.27 (aCOAとbCOAを入れ替え 18.08.26)
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6>>x.a7>>x.a8;
			from>>x.a10>>x.a12>>x.a14>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>buf; x.legendre.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.bCOA()>>x.aCOA()>>x.tCOA; x.match_acoa_cacoa(); x.match_bcoa_cbcoa();
			from>>x.SA0>>x.CM0;
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx>>x.CHMy>>x.CHMx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 130: // a1,a2,a3,a5,a4を追加 180123
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a1>>x.a2>>x.a3>>x.a4>>x.a5>>x.a6>>x.a7>>x.a8;
			from>>x.a10>>x.a12>>x.a14>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>buf; x.legendre.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx>>x.CHMy>>x.CHMx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 129: // CHMy,CHMxを追加 171006
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a4>>x.a6>>x.a8>>x.a10>>x.a12>>x.a14>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>buf; x.legendre.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx>>x.CHMy>>x.CHMx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 128: // legendreを追加 170526
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a4>>x.a6>>x.a8>>x.a10>>x.a12>>x.a14>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>buf; x.legendre.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 127: // NormHを追加 160610
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.NormH()>>x.a4>>x.a6>>x.a8>>x.a10>>x.a12>>x.a14>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 126: // a18,a20を追加 160317
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10>>x.a12>>x.a14>>x.a16>>x.a18>>x.a20;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			x.match_fi_pi();
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 125: // Fishを追加 160223
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10>>x.a12>>x.a14>>x.a16;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.Fish;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 124: // cLens1にAfocalRotateEyeを追加 150925
		case 123: // a16を追加 150421
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10>>x.a12>>x.a14>>x.a16;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 122: // a12,a14を追加 140619
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10>>x.a12>>x.a14;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 121:  // cLens1にStopDominateを追加 130821
		case 120:  // Dconを追加 130317
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>buf; x.Dcon.SetTerms(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 119:  // CoatName,CoatReverseを追加 130213
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.CoatName; x.CoatName=inv_blank_filled(x.CoatName);
			from>>x.CoatReverse;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 118:  // SA0,CM0を追加 120912
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal()>>x.SA0>>x.CM0;  // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();                 // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 117:  // 自由曲面の表現を変更 120906
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>buf; x.b.SetTerms(inv_blank_filled(buf));
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 116:  // IsUserDefSurfを追加 120605
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.b.bb[2][0]>>x.b.bb[1][1]>>x.b.bb[0][2];
			from>>x.b.bb[3][0]>>x.b.bb[2][1]>>x.b.bb[1][2]>>x.b.bb[0][3];
			from>>x.b.bb[4][0]>>x.b.bb[3][1]>>x.b.bb[2][2]>>x.b.bb[1][3]>>x.b.bb[0][4];
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.UserDefSurf;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 115:  // Apv,sfx,sfyを追加 120410
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.b.bb[2][0]>>x.b.bb[1][1]>>x.b.bb[0][2];
			from>>x.b.bb[3][0]>>x.b.bb[2][1]>>x.b.bb[1][2]>>x.b.bb[0][3];
			from>>x.b.bb[4][0]>>x.b.bb[3][1]>>x.b.bb[2][2]>>x.b.bb[1][3]>>x.b.bb[0][4];
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>x.Apv>>x.sfx>>x.sfy;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 114:  // cLens1にNoteを追加 120406
		case 113:
			// IsXToroidを追加 120316
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.b.bb[2][0]>>x.b.bb[1][1]>>x.b.bb[0][2];
			from>>x.b.bb[3][0]>>x.b.bb[2][1]>>x.b.bb[1][2]>>x.b.bb[0][3];
			from>>x.b.bb[4][0]>>x.b.bb[3][1]>>x.b.bb[2][2]>>x.b.bb[1][3]>>x.b.bb[0][4];
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx>>x.IsXToroid;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 112:
			// kpx,kpyを追加 111220
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c();
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.b.bb[2][0]>>x.b.bb[1][1]>>x.b.bb[0][2];
			from>>x.b.bb[3][0]>>x.b.bb[2][1]>>x.b.bb[1][2]>>x.b.bb[0][3];
			from>>x.b.bb[4][0]>>x.b.bb[3][1]>>x.b.bb[2][2]>>x.b.bb[1][3]>>x.b.bb[0][4];
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx>>x.kpy>>x.kpx;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 111:
			// zernikeを追加 110601
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			x.match_r_c(); // r(),c()を通してrr,ccを取得を取得するのであれば，不要だが，
			               // こうしておけば ポインタ double *x=&c() や，参照 double &x=c()
			               // を通して取得することも可能となる．
			               // これは，cLens::optimize() のコードを簡潔にするのに有用．
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.b.bb[2][0]>>x.b.bb[1][1]>>x.b.bb[0][2];
			from>>x.b.bb[3][0]>>x.b.bb[2][1]>>x.b.bb[1][2]>>x.b.bb[0][3];
			from>>x.b.bb[4][0]>>x.b.bb[3][1]>>x.b.bb[2][2]>>x.b.bb[1][3]>>x.b.bb[0][4];
			from>>buf; x.zernike.SetCoefficients(inv_blank_filled(buf));
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 110:
			// cLens1にSourcePhiY等を追加 100805
		case 109:
			// mediumのメンバgnameで空白を許容する 100412
		case 108:
			// order,order1を追加 100202
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.b.bb[2][0]>>x.b.bb[1][1]>>x.b.bb[0][2];
			from>>x.b.bb[3][0]>>x.b.bb[2][1]>>x.b.bb[1][2]>>x.b.bb[0][3];
			from>>x.b.bb[4][0]>>x.b.bb[3][1]>>x.b.bb[2][2]>>x.b.bb[1][3]>>x.b.bb[0][4];
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz;
			from>>x.order;
			from>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.order1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 107:
			// cLens1にExcludeVirtualRay,ExcludeVirtualObjectを追加 091018
		case 106:
			// cLens1にFRWを追加 090108
		case 105:
			// cLens1にAfocalを追加 080912
		case 104:
			// b20,b11,...,b04を追加 080328
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.b.bb[2][0]>>x.b.bb[1][1]>>x.b.bb[0][2];
			from>>x.b.bb[3][0]>>x.b.bb[2][1]>>x.b.bb[1][2]>>x.b.bb[0][3];
			from>>x.b.bb[4][0]>>x.b.bb[3][1]>>x.b.bb[2][2]>>x.b.bb[1][3]>>x.b.bb[0][4];
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 103:
			// rem(コメント)を追加 080312
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			from>>x.rem; x.rem=inv_blank_filled(x.rem);
			break;
		case 102:
			// Diffusionを追加 070717
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.Diffusion;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			break;
		case 101:
			// dx1,dy1,dz1,rox1,roy1,roz1を追加 070416
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz>>x.ret;
			from>>x.dx1>>x.dy1>>x.dz1>>x.rox1>>x.roy1>>x.roz1;
			break;
		case 100:
		default:
			from>>x.r()>>x.NewtonTol>>x.AsTol>>x.rVariable;
			from>>x.asph_type;
			from>>x.kp>>x.a4>>x.a6>>x.a8>>x.a10;
			from>>x.cylinder;
			from>>x.fresnel>>x.rbase;
			from>>x.ry>>x.rx;
			from>>x.coneangle;
			from>>x.fideal();   // pideal() を追加に伴う変更( fideal -> fideal() ) 170427
			x.match_fi_pi();    // pideal() を追加に伴う追記 170427
			from>>x.grating>>x.difforder;
			from>>x.gpitch>>x.grx>>x.gry>>x.grz;
			from>>x.EAtype>>x.EAy>>x.EAx;
			from>>x.decenter_type;
			from>>x.dx>>x.dy>>x.dz>>x.rox>>x.roy>>x.roz>>x.ret;
			break;
	}

	{
		// 090528
		// decenter_type=4 : dx,dy,dz,rox,roy,roz, ret
		// decenter_type=5 : dx,dy,dz,rox,roy,roz, dx1,dy1,dz1,rox1,roy1,roz1
		// であったものを,
		// decenter_type=4 : dx,dy,dz,rox,roy,roz, ret, dx1,dy1,dz1,rox1,roy1,roz1
		// に統合し，decenter_type=5を廃止
		 
		// 101109
		// decenter_type=1 : dx,dy,dz,rox,roy,roz
		// であったものを
		// decenter_type=1 : dx,dy,dz,rox,roy,roz, ret, dx1,dy1,dz1,rox1,roy1,roz1
		// とし，decenter_type=4を廃止

		if(x.decenter_type==5 || x.decenter_type==4 ) x.decenter_type=1;
	}

	return from;
}

bool operator==(surface& a,surface& b){
	bool x;
	x=(a.r()==b.r());
	x=x&&(a.NewtonTol==b.NewtonTol);
	x=x&&(a.AsTol==b.AsTol);
	x=x&&(a.rVariable==b.rVariable);
	x=x&&(a.asph_type==b.asph_type);
	x=x&&(a.kp==b.kp);
	x=x&&(a.NormH()==b.NormH());
	x=x&&(a.a1==b.a1 && a.a2==b.a2 && a.a3==b.a3);
	x=x&&(a.a4==b.a4 && a.a5==b.a5 && a.a6==b.a6 && a.a7==b.a7);
	x=x&&(a.a8==b.a8 && a.a9==b.a9 && a.a10==b.a10 && a.a11==b.a11);
	x=x&&(a.a12==b.a12 && a.a13==b.a13 && a.a14==b.a14 && a.a15==b.a15 && a.a16==b.a16);
	x=x&&(a.a18==b.a18 && a.a20==b.a20);
	x=x&&(a.b==b.b);
	x=x&&(a.Apv==b.Apv && a.sfx==b.sfx && a.sfy==b.sfy);
	x=x&&(a.UserDefSurf==b.UserDefSurf);
	x=x&&(a.Fish==b.Fish);
	x=x&&(a.cylinder==b.cylinder);
	x=x&&(a.fresnel==b.fresnel && a.rbase==b.rbase );
	x=x&&(a.ry==b.ry && a.rx==b.rx && a.kpy==b.kpy && a.kpx==b.kpx && a.IsXToroid==b.IsXToroid);
	x=x&&(a.coneangle==b.coneangle);
	x=x&&(a.fideal()==b.fideal());
	x=x&&(a.aCOA()==b.aCOA() && a.bCOA()==b.bCOA() && a.tCOA==b.tCOA);
	x=x&&(a.grating==b.grating && a.difforder==b.difforder);
	x=x&&(a.gpitch==b.gpitch);
	x=x&&(a.grx==b.grx && a.gry==b.gry && a.grz==b.grz);
	x=x&&(a.EAtype==b.EAtype);
	x=x&&(a.EAy==b.EAy); x=x&&(a.EAx==b.EAx);
	x=x&&(a.CHMy==b.CHMy); x=x&&(a.CHMx==b.CHMx);
	x=x&&(a.EAdy==b.EAdy); x=x&&(a.EAdx==b.EAdx);
	x=x&&(a.decenter_type==b.decenter_type);
	x=x&&(a.dx==b.dx && a.dy==b.dy && a.dz==b.dz );
	x=x&&(a.rox==b.rox && a.roy==b.roy && a.roz==b.roz && a.ret==b.ret);
	x=x&&(a.ret==b.ret);
	x=x&&(a.order==b.order);
	x=x&&(a.dx1==b.dx1 && a.dy1==b.dy1 && a.dz1==b.dz1);
	x=x&&(a.rox1==b.rox1 && a.roy1==b.roy1 && a.roz1==b.roz1);
	x=x&&(a.order1==b.order1);
	x=x&&(a.CoatName==b.CoatName);
	x=x&&(a.CoatReverse==b.CoatReverse);
	return x;
}

surface surface::reversed() {
	surface x=*this;
	x.r()=-x.r();
	x.Newton=-x.Newton; x.As0=-x.As0; x.As45=-x.As45;
	x.a1=-x.a1; x.a2=-x.a2; x.a3=-x.a3;
	x.a4=-x.a4; x.a5=-x.a5; x.a6=-x.a6; x.a7=-x.a7; x.a8=-x.a8; x.a9=-x.a9;
	x.a10=-x.a10; x.a11=-x.a11; x.a12=-x.a12; x.a13=-x.a13; x.a14=-x.a14; x.a15=-x.a15; x.a16=-x.a16;
	x.a18=-x.a18; x.a20=-x.a20;
	x.b=x.b.reversed();
	x.Dcon=x.Dcon.reversed();
	x.legendre=x.legendre.reversed();
	x.spline=x.spline.reversed();
	x.Apv=-x.Apv;
	x.rbase=-x.rbase;
	x.ry=-x.ry; x.rx=-x.rx;
	x.coneangle=-x.coneangle;
	x.aCOA()=-x.aCOA(); x.bCOA()=-x.bCOA();
	x.grx=-x.grx; x.grz=-x.grz;
	x.CoatReverse= x.CoatReverse==0 ? 1 : 0;
	return x;
}

void surface::scale(double m,int with_EA,int with_decenter) {
	double m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15,m17,m19;
	
	if(m==0) return;
	
	m2=m*m;
	m3=m2*m;
	m4=m3*m;
	m5=m4*m;
	m6=m5*m;
	m7=m6*m;
	m8=m7*m;
	m9=m8*m;
	m10=m9*m;
	m11=m10*m;
	m12=m11*m;
	m13=m12*m;
	m14=m13*m;
	m15=m14*m;
	m17=m15*m*m;
	m19=m17*m*m;

	r()*=m;
	if(NormH()==1){
		a2/=m; a3/=m2;
		a4/=m3; a5/=m4; a6/=m5; a7/=m6; a8/=m7; a10/=m9; a12/=m11; a14/=m13; a16/=m15; a18/=m17; a20/=m19;
	}
	else{
		NormH()*=m;
		a1*=m; a2*=m; a3*=m; a4*=m; a5*=m; a6*=m; a7*=m; a8*=m; a9*=m; a10*=m;
		a11*=m; a12*=m; a13*=m; a14*=m; a15*=m; a16*=m; a18*=m; a20*=m;
	}
	b.scale(m);
	Dcon.scale(m);
	legendre.scale(m);
	spline.scale(m);
	Apv*=m; sfx/=m; sfy/=m;
	rbase*=m;
	ry*=m; rx*=m;
	fideal()*=m;
	aCOA()*=m; bCOA()*=m;
	gpitch*=m;
	if(with_EA) { EAy*=m; EAx*=m; CHMy*=m; CHMx*=m; EAdy*=m; EAdx*=m; }
	if(with_decenter) { dx*=m; dy*=m; dz*=m; dx1*=m; dy1*=m; dz1*=m; }
}


////// freeform メンバ ///////////////////////////////////////////////

freeform::freeform(){
	int i,j;
	for(i=0; i<=N; i++) for(j=0; j<=N-i; j++) bb[i][j]=0;
}

freeform::freeform(const freeform &x){
	*this=x;
}

freeform& freeform::operator =(const freeform& x){
	int i,j;
	for(i=0; i<=N; i++) for(j=0; j<=N-i; j++) bb[i][j]=x.bb[i][j];
	return *this;
}

bool operator ==(const freeform &a,const freeform &b){
	int i,j;
	bool x=true;

	for(i=0; i<=a.N; i++) for(j=0; j<=a.N-i; j++){
		x=x&&(a.bb[i][j]==b.bb[i][j]);
	}
	return x;
}

double& freeform::b(int m,int n){
	if(0<=m && 0<=n && m+n>0 && m+n<=N){  // m=n=0 は許さない
		return bb[m][n];
	}
	else{
		return null;
	}
}

int freeform::NotNull(){
	int i,j;

	for(i=0; i<=N; i++) for(j=0; j<=N-i; j++){
		if(bb[i][j]!=0) return true;
	}
	return false;
}

freeform freeform::reversed(){
	int i,j;
	freeform x=*this;

	for(i=0; i<=N; i++) for(j=0; j<=N-i; j++){
		x.bb[i][j]=-x.bb[i][j];
	}
	return x;
}

void freeform::scale(double m){
	int i,j;
	double x;

	if(m==0) return;

	for(i=0; i<=N; i++) for(j=0; j<=N-i; j++){
		if(i+j>=2){
			x=pw(m,i+j);
			bb[i][j]/=x;
		}
	}
}

std::string freeform::GetTerms(){
	int i,j,n;
	char buf[100];
	std::string s;

	for(n=0; n<=N; n++){
		for(j=0; j<=n; j++){
			i=n-j;
			if(bb[i][j]!=0){
				sprintf(buf,"b %d %d %g; ", i,j,bb[i][j]);
				s+=buf;
			}
		}
	}
	return s;
}

void freeform::SetTerms(std::string terms){
	std::string s,s0,s1,s2,s3;
	bool b1,b2,b3;
	int i,j,k;
	double b;

	for(i=0; i<=N; i++) for(j=0; j<=N-i; j++) bb[i][j]=0;
	
	for(k=1; k<=sentences(terms); k++){
		s=sentence(terms,k);
		s0=arg(s,0);
		s1=arg(s,1); b1=is_numeric(s1);
		s2=arg(s,2); b2=is_numeric(s2);
		s3=arg(s,3); b3=is_numeric(s3);

		if(s0=="b" && b1 && b2 && b3){
			i=atoi(s1.c_str());
			j=atoi(s2.c_str());
			b=atof(s3.c_str());
			if(0<=i && 0<=j && i+j<=N){
				this->bb[i][j]=b;
			}
		}
	}
}
