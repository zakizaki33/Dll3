#include "stdafx.h"
#include "cZernike.h"

void cZernike::ConvertR0(){
	double *a;
	int j;

	if(r0fix!=0 && r0!=r0fix){
		// 展開径の変更
		//   展開径r0fixの多項式に改めて最小二乗フィッティングせずに，
		//   変換公式により展開径r0fixの多項式の係数に変換する．
		//   展開式の表す関数は変換前後でexactに一致するが，
		//   r0fixの多項式に改めてフィッティングしたときと一致するとは限らない．
		//   この変換式を評価する等に使える．(10.04.30)
		a=new double[terms+1];
		for(j=1; j<=terms; j++){
			a[j]=X[j][1];
		}
		ConvertR0(a,terms,r0,r0fix,normalize,IsFringeOrder);
		for(j=1; j<=terms; j++){
			X[j][1]=a[j];
		}
		r0=r0fix;
		delete a;
	}
}


///// public members ///////////////////////////

int cZernike::TotalTerms(int order,int IsFringeOrder){
	// order次までのZernike多項式の総項数を返す．
	if(IsFringeOrder){
		// order=n+m;
		return (is_even(order) && order>=0) ? (order/2+1)*(order/2+1) : 0;
	}
	else{
		// order=n
		return order>=0 ? (order+1)*(order+2)/2 : 0;
	}
}

int cZernike::jNumber(int l,int n,int IsFringeOrder,int jBase){
	// (l,n) の項が通し番号で何番目の項かを返す．
	int j;

	if(IsFringeOrder){
		//   Fringe orderの規則
		//     まず n+m=n+|l| の小さい順
		//     次に nの小さい順
		//     次に m>0, m<0 の順
		//
		//      j         n  l
		//    (n+m=0)
		//      1         0  0          
		//    (n+m=2)
		//      2         1  1
		//      3         1 -1
		//      4         2  0
		//    (n+m=4)
		//      5         2  2
		//      6         2 -2
		//      7         3  1
		//      8         3 -1
		//      9         4  0
		//    (n+m=6)
		//     10         3  3
		//     11         3 -3
		//     12         4  2
		//     13         4 -2
		//     14         5  1
		//     ..        .. ..
		//  式はこの表を見ながら考え出した．
		//   ( (n+m)/2+1 )*( (n+m)/2+1 ) はn+mが等しい各グループの最後の項のjである．

		j=( (n+abs(l))/2+1 )*( (n+abs(l))/2+1 ) -2*abs(l) + (l<0 ? 1 : 0);
	}
	else{
		j=(n*(n+2)+l)/2+1;   // ISO 24157 eq.(17)
	}
	return jBase==0 ? j-1 : j;
}

int cZernike::mNumber(int j,int IsFringeOrder,int jBase){
	int l;
	
	l=lNumber(j,IsFringeOrder,jBase);
	return l>=0 ? l : -l;
}

int cZernike::StandardNo(int FringeNo,int jBase){	
	// Finge no. を Standard no. に変換する．
	//
	//   Fringe orderの規則
	//     まず n+m=n+|l| の小さい順
	//     次に nの小さい順
	//     次に m>0, m<0 の順
	//
	//   fring no.  n  l  standard no.
	//  (n+m=0)
	//    1         0  0   1          
	//  (n+m=2)
	//    2         1  1   3
	//    3         1 -1   2
	//    4         2  0   5
	//  (n+m=4)
	//    5         2  2   6
	//    6         2 -2   4
	//    7         3  1   9
	//    8         3 -1   8
	//    9         4  0  13
	//  (n+m=6)
	//   10         3  3  10
	//   11         3 -3   7
	//   12         4  2  14
	//   13         4 -2  12
	//   14         5  1  19
	//   ..        .. ..  ..

	int j,s,m,l;

	if(jBase==0) FringeNo+=1;

	j=0;  // standard no.
	s=0;  // =n+m=n+|l|
	do{
		for(m=s/2; m>=0; m--){
			j++;
			l=m;
			if(j==FringeNo) return jNumber(l,s-m,0,jBase);

			if(m!=0){
				j++;
				l=-m;
				if(j==FringeNo) return jNumber(l,s-m,0,jBase);
			}
		}
		s+=2;
	}while(true);
}

int cZernike::lNumber(int j,int IsFringeOrder,int jBase){
	// 通し番号j番目の項のlを返す．
	if(jBase==0) j+=1;
	if(IsFringeOrder) j=StandardNo(j,1);
	return -nNumber(j,0,1) + ( j-TotalTerms(nNumber(j,0,1)-1,0)-1 )*2;
}

int cZernike::nNumber(int j,int IsFringeOrder,int jBase){
	// j番目の項のn(次数)を返す．
	int n;

	if(jBase==0) j+=1;
	if(IsFringeOrder) j=StandardNo(j,1);
	n=0;
	while( TotalTerms(n,0) < j ){
		n+=1;
	}
	return n;
}

int cZernike::Order(int j,int IsFringeOrder,int jBase){
	if(IsFringeOrder){
		// フリンジオーダーでは order=n+m
		return nNumber(j,IsFringeOrder,jBase)+mNumber(j,IsFringeOrder,jBase);
	}
	else{
		// スタンダードオーダーでは order=n
		return nNumber(j,IsFringeOrder,jBase);
	}
}

void cZernike::ConvertR0(double *a,int terms,double r1,double r2,int Normalized,int IsFringeOrder){
	// 展開半径r1におけるZernike展開係数を展開半径r2におけるそれに変換する．
	// 引数 a[1],a[2], .... a[terms], terms=項数
	// によりピストン項より順番に元の係数を与え，変換後の係数はa[]に上書きされる．
	// Normalized!=0 内積を積の単位円内平均とするとき自己内積(2乗平均)が1になる正規化
	//           ==0 単位円内の最大最小が±1となる正規化
	// IsFringeOrder!=0 フリンジオーダー
	//              ==0 スタンダードオーダー
	//
	// 変換式は文献：Scaling Zernike expansion coefficients to different pupil sizes
	//              (Jim Schwiegerling, Vol.19,No.10/October 2002/J.Opt.Soc.Am.A)
	// の(A8)式による．
	// (A8)式は複雑であるため，ここでは，(A8)式から規格化に関する要素（根号の付いた要素）
	// は取り除いた式を用い，規格化に関する処理は変換式適用前後に行なっている．
	// (こうすることで規格化されていない係数にも簡単に対応できる．）
	// また，“光学の原理Ⅱ(M.Born,E.Wolf)”の表記との整合のため，
	//   ここでのl=文献のm
	//   ここでのm=文献の|m|
	//   係数a,bの添字は文献では 動径次数,角度周波数 の順だが，ここでは 角度周波数,動径次数
	// としている．

	int n,nmax, l,lmin,lmax, m,i2,j,k;
	double *b;
	double sum1,sum2;

	b=new double [terms+1];

	// 規格化に関する処理
	if(Normalized){
		for(j=1; j<=terms; j++){
			a[j]*=sqrt( 2*(nNumber(j,IsFringeOrder,1)+1)/(1+( lNumber(j,IsFringeOrder,1)==0 ? 1:0 )) );
		}
	}

	// lの最小，最大を求める
	lmin=lmax=0;
	for(j=1; j<=terms; j++){
		l=lNumber(j,IsFringeOrder,1);
		if(l<lmin) lmin=l;
		if(l>lmax) lmax=l;
	}

	// 変換式を適用する
	for(l=lmin; l<=lmax; l++){
		m= l>=0 ? l : -l;

		// lの列の最高次数nmaxを求める
		nmax=0;
		for(j=1; j<=terms; j++){
			if( l==lNumber(j,IsFringeOrder,1) ){
				n=nNumber(j,IsFringeOrder,1);
				if(n>nmax) nmax=n;
			}
		}
		
		for(i2=nmax-m-(is_even(nmax-m) ? 0 :1); i2>=0; i2-=2){
			// ピラミッド図のlの列の次数の高い項から変換後の係数b(l,n)=b(l,m+i2)を求めていく．

			// (A8)式{}内の第1項
			sum1=0;
			for(j=0; j<=(nmax-m-i2)/2; j++){
				sum1+=(is_even(j) ? 1 : -1)*a[jNumber(l,m+i2+j*2,IsFringeOrder,1)]
				      *factorial(m+i2+j)/factorial(j)/factorial(m+i2);
			}
			
			// (A8)式{}内の第2項
			sum2=0;
			for(k=m+i2+2; k<=nmax; k+=2){
				sum2+=b[jNumber(l,k,IsFringeOrder,1)]
				      *(is_even((k-m-i2)/2) ? 1 : -1)*factorial((k+m+i2)/2)
					  /factorial((k-m-i2)/2)/factorial(m+i2)
					  *pw(r1,m+i2)/pw(r2,m+i2);
			}

			b[jNumber(l,m+i2,IsFringeOrder,1)]=pw(r2,m+i2)/pw(r1,m+i2)*(sum1-sum2);
		}
	}

	// 規格化に関する処理
	if(Normalized){
		for(j=1; j<=terms; j++){
			b[j]/=sqrt( 2*(nNumber(j,IsFringeOrder,1)+1)/(1+( lNumber(j,IsFringeOrder,1)==0 ? 1:0 )) );
		}
	}

	// 結果を引数に上書きする
	for(j=1; j<=terms; j++){
		a[j]=b[j];
	}
	
	delete [] b;
}

int cZernike::TotalTerms(int order){
	return TotalTerms(order,this->IsFringeOrder);
}

int cZernike::jNumber(int l,int n){
	return jNumber(l,n,this->IsFringeOrder,this->jBase);
}

int cZernike::StandardNo(int FringeNo){
	return StandardNo(FringeNo,this->jBase);
}

int cZernike::lNumber(int j) const{
	return lNumber(j,this->IsFringeOrder,this->jBase);
}

int cZernike::mNumber(int j) const{
	return mNumber(j,this->IsFringeOrder,this->jBase);
}

int cZernike::nNumber(int j) const{
	return nNumber(j,this->IsFringeOrder,this->jBase);
}

int cZernike::Order(int j) const{
	return Order(j,this->IsFringeOrder,this->jBase);
}

cZernike::cZernike() {
	r0=1;
	r0fix=0;
	normalize=1;
	IsFringeOrder=0;
	jBase=0;
	Digits=0;
	fitted=false;
	SetMaxOrder(4);
}

cZernike::~cZernike() {

}

double cZernike::R(int l,int n,double rho){
	// Zernike動径多項式
	// 参考書: “光学の原理Ⅱ(M.Born,E.Wolf)” 9.2.2 Zernikeのcircle polynomial
	// の(5)式
	int s,m;
	double result;
	
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		result=0;
		for(s=0; s<=(n-m)/2; s++){
			result+=(is_even(s) ? 1: -1)
			       *factorial(n-s)/factorial(s)/factorial((n+m)/2-s)/factorial((n-m)/2-s)
			       *pw(rho,n-2*s);
		}
	}
	else{
		result=0;
	}
	return result;
}

double cZernike::Rrho(int l,int n,double rho){
	// Zernike動径多項式のrho偏微分
	int s,m;
	double result;
	
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		result=0;
		for(s=0; s<=(n-m)/2; s++){
			result+=(is_even(s) ? 1: -1)
			       *factorial(n-s)/factorial(s)/factorial((n+m)/2-s)/factorial((n-m)/2-s)
				   *( n-2*s>=1 ? (n-2*s)*pw(rho,n-2*s-1) : 0 );
		}
	}
	else{
		result=0;
	}
	return result;
}

double cZernike::Rrhorho(int l,int n,double rho){
	// Zernike動径多項式のrho2階偏微分
	int s,m;
	double result;
	
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		result=0;
		for(s=0; s<=(n-m)/2; s++){
			result+=(is_even(s) ? 1: -1)
			       *factorial(n-s)/factorial(s)/factorial((n+m)/2-s)/factorial((n-m)/2-s)
				   *( n-2*s>=2 ? (n-2*s)*(n-2*s-1)*pw(rho,n-2*s-2) : 0 );
		}
	}
	else{
		result=0;
	}
	return result;
}

double cZernike::U(int l,int n,double x,double y) const{
	// Zernike実多項式
	// 参考書: “光学の原理Ⅱ(M.Born,E.Wolf)” 9.2.2 Zernikeのcircle polynomial
	
	// n=0  l=0             piston
	// n=1  l=-1,1          tilt
	// n=2  l=-2,0,2        defocus,astigmatism
	// n=3  l=-3,-1,1,3     coma aberration etc.
	// n=4  l=-4,-2,0,2,4   spherical aberration etc.

	// Zernike関数のような，R(r)f(θ)の関数形は，
	// R(r)がr=0で0でなく，f(θ)が一定でなければ
	// xy座標の原点(r=0)でθによって関数値が変化する(つまりxy座標では定義できない)．
	// このことは，xy座標系の原点はｒθ座標のr=0直線に対応し，
	// 両座標系は一対一の対応になっていないことに関係する．
	// しかし，全てのZernike関数はおそらくx=y=0で連続である．
	// したがって，x=y=0でθによって値は変化しない．
	// ここでは，x=y=0でθ=atan(0,0)≡0 を代表値として与えている．
	// しかし，例えば偏微分Urについては明らかにティルト項のそれはxy座標の原点で不連続であるから，
	// 注意しなけらばならない．

	int m;
	double rho,th, R;
	
	if(r0==0) return 0;
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=sqrt(x*x+y*y)/r0;
		th =atan2(y,x); 
		R=this->R(l,n,rho);
		// normalize!=0 内積を積の単位円内平均とするとき自己内積(2乗平均)が1になる正規化
		// normalize==0 単位円内の最大最小が±1となる正規化
		if(normalize) R*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? R*cos(m*th) : R*sin(m*th);
	}
	else{
		return 0;
	}
}

double cZernike::U(int j,double x,double y) const{
	return U(lNumber(j),nNumber(j),x,y);
}

double cZernike::Ur(int l,int n,double x,double y) const{
	// Uのr偏微分
	int m;
	double rho,th, Rr;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // 原点x=y=0では多価関数であるため使えない．Ur_pol()を使うこと．
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=sqrt(x*x+y*y)/r0;
		th =atan2(y,x);
		Rr=Rrho(l,n,rho)/r0;
		if(normalize) Rr*= m==0 ? sqrt(double(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? Rr*cos(m*th) : Rr*sin(m*th);
	}
	else{
		return 0;
	}
}

double cZernike::Ur_pol(int l,int n,double r,double th) const{
	// Uのr偏微分
	int m;
	double rho, Rr;

	if(r0==0) return 0;
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=r/r0;
		Rr=Rrho(l,n,rho)/r0;
		if(normalize) Rr*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? Rr*cos(m*th) : Rr*sin(m*th);
	}
	else{
		return 0;
	}
}

double cZernike::Uth(int l,int n,double x,double y) const{
	// Uのθ偏微分
	int m;
	double rho,th, R;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // x=y=0では0である．Ur()と違ってこれは正しい．
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=sqrt(x*x+y*y)/r0;
		th =atan2(y,x);
		R=this->R(l,n,rho);;
		if(normalize) R*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? R*(-m*sin(m*th)) : R*(m*cos(m*th));
	}
	else{
		return 0;
	}
}

double cZernike::Urr(int l,int n,double x,double y) const{
	// Uのr2階偏微分
	int m;
	double rho,th, Rrr;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // 原点x=y=0では多価関数であるため使えない．Ur_pol()を使うこと．
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=sqrt(x*x+y*y)/r0;
		th =atan2(y,x);
		Rrr=Rrhorho(l,n,rho)/r0/r0;
		if(normalize) Rrr*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? Rrr*cos(m*th) : Rrr*sin(m*th);
	}
	else{
		return 0;
	}
}

double cZernike::Urr_pol(int l,int n,double r,double th) const{
	// Uのr2階偏微分
	int m;
	double rho, Rrr;

	if(r0==0) return 0;
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=r/r0;
		Rrr=Rrhorho(l,n,rho)/r0;
		if(normalize) Rrr*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? Rrr*cos(m*th) : Rrr*sin(m*th);
	}
	else{
		return 0;
	}
}

double cZernike::Uthth(int l,int n,double x,double y) const{
	// Uのθ2階偏微分
	int m;
	double rho,th, R;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // x=y=0では0である．Urr()と違ってこれは正しい．
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=sqrt(x*x+y*y)/r0;
		th =atan2(y,x);
		R=this->R(l,n,rho);;
		if(normalize) R*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? R*(-m*m*cos(m*th)) : R*(-m*m*sin(m*th));
	}
	else{
		return 0;
	}
}

double cZernike::Urth(int l,int n,double x,double y) const{
	// Uのr,th偏微分
	int m;
	double rho,th, Rr;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // 原点x=y=0では多価関数であるため使えない．Urth_pol()を使うこと．
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=sqrt(x*x+y*y)/r0;
		th =atan2(y,x);
		Rr=Rrho(l,n,rho)/r0;
		if(normalize) Rr*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? Rr*(-m*sin(m*th)) : Rr*(m*cos(m*th));
	}
	else{
		return 0;
	}
}

double cZernike::Urth_pol(int l,int n,double r,double th) const{
	// Uのr,th偏微分
	int m;
	double rho, Rr;

	if(r0==0) return 0;
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=r/r0;
		Rr=Rrho(l,n,rho)/r0;
		if(normalize) Rr*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
		return l>=0 ? Rr*(-m*sin(m*th)) : Rr*(m*cos(m*th));
	}
	else{
		return 0;
	}
}

double cZernike::Ux(int l,int n,double x,double y) const{
	// Uのx偏微分
	double r;

	r=sqrt(x*x+y*y);
	if(r==0){
		// dr/dx,dth/dxは原点で定義されない．
		// しかし，Uxはおそらく原点で滑らかな関数である．
		// そこで，dr/dx,dth/dxはx>0での片側微分をとって，
		// dr/dx=1, dth/dx=0とする．
		// Ur(x,y)も原点で定義されない．Ur_pol(r,th)を使う．
		return Ur_pol(l,n,0,0);
	}
	else{		
		// Ux=Ur(dr/dx)+Uth(dth/dx)
		// dr/dx=x/r
		// また，tan(th)=y/x より， d(th)/cos(th)^2 = (-y/x^2)dx+(1/x)dy
		// だから，dth/dx=-y/r^2
		// 【注】 x=rcos(th)よりdx=-rsin(th)dth -> dth/dx=-1/rsin(th)=-1/y は間違い
		//        偏微分ではyが定数であり，rが定数ではない．
		return Ur(l,n,x,y)*(x/r)+Uth(l,n,x,y)*(-y/r/r);
	}
}

double cZernike::Ux(int j,double x,double y) const{
	return Ux(lNumber(j),nNumber(j),x,y);
}

double cZernike::Uy(int l,int n,double x,double y) const{
	// Uのy偏微分
	double r;

	r=sqrt(x*x+y*y);
	if(r==0){
		// dr/dy,dth/dyは原点で定義されない．
		// しかし，Uyはおそらく原点で滑らかな関数である．
		// そこで，dr/dy,dth/dyはy>0での片側微分をとって，
		// dr/dy=1, dth/dy=0とする．
		// Ur(x,y)も原点で定義されない．Ur_pol(r,th)を使う．
		return Ur_pol(l,n,0,PI/2);
	}
	else{
		// Uy=Ur(dr/dy)+Uth(dth/dy)
		// dr/dy=d(xx+yy)/dy=2y
		// また，tan(th)=y/x より， d(th)/cos(th)^2 = (-y/x^2)dx+(1/x)dy
		// だから，dth/dy=x/r^2
		return Ur(l,n,x,y)*(y/r)+Uth(l,n,x,y)*(x/r/r);
	}
}

double cZernike::Uy(int j,double x,double y) const{
	return Uy(lNumber(j),nNumber(j),x,y);
}

double cZernike::Uxx(int l,int n,double x,double y) const{
	double Urx,Uthx,r;

	r=sqrt(x*x+y*y);
	if(r==0){
		return Urr_pol(l,n,0,0);
	}
	else{
		// Uxx=d/dx(Ux)=d/dx(Ur(dr/dx)+Uth(dth/dx))
		//    =Urx(dr/dx)+Ur(d2r/dx2)+Uthx(dth/dx)+Uth(d2th/dx2)
		Urx =Urr(l,n,x,y)*(x/r)+Urth(l,n,x,y)*(-y/r/r);
		Uthx=Urth(l,n,x,y)*(x/r)+Uthth(l,n,x,y)*(-y/r/r);
		return Urx*(x/r)+Ur(l,n,x,y)*((r*r-x*x)/r/r/r)+Uthx*(-y/r/r)+Uth(l,n,x,y)*(2*x*y/r/r/r/r);
	}
}

double cZernike::Uxx(int j,double x,double y) const{
	return Uxx(lNumber(j),nNumber(j),x,y);
}

double cZernike::Uyy(int l,int n,double x,double y) const{
	double Ury,Uthy,r;
	
	r=sqrt(x*x+y*y);
	if(r==0){
		return Urr_pol(l,n,0,PI/2);
	}
	else{
		// Uyy=d/dy(Uy)=d/dy(Ur(dr/dy)+Uth(dth/dy))
		//    =Ury(dr/dy)+Ur(d2r/dy2)+Uthy(dth/dy)+Uth(d2th/dy2)
		Ury =Urr(l,n,x,y)*(y/r)+Urth(l,n,x,y)*(x/r/r);
		Uthy=Urth(l,n,x,y)*(y/r)+Uthth(l,n,x,y)*(x/r/r);
		return Ury*(y/r)+Ur(l,n,x,y)*((r*r-y*y)/r/r/r)+Uthy*(x/r/r)+Uth(l,n,x,y)*(-2*x*y/r/r/r/r);
	}
}

double cZernike::Uyy(int j,double x,double y) const{
	return Uyy(lNumber(j),nNumber(j),x,y);
}

double cZernike::Uxy(int l,int n,double x,double y) const{
	double Ury,Uthy,r;

	r=sqrt(x*x+y*y);
	if(r==0){
		// 原点でのテーラー展開の2次の項は，
		//     z = (1/2)Uxx(dxdx)+Uxy(dxdy)+(1/2)Uyy(dydy)
		// となり，dx=dy=1のとき，
		//     z(1,1) = (1/2)Uxx+Uxy+(1/2)Uyy
		// 一方，Urrを使うと，dr^2=dx^2+dy^2=2 だから，
		//     z(1,1) = (1/2)Urr(θ=PI/4)*2=Urr(θ=PI/4)
		// したがって，
		//     (1/2)Uxx+Uxy+(1/2)Uyy = Urr(θ=PI/4)
		//     Uxy = Urr(θ=PI/4)-(1/2)Uxx-(1/2)Uyy
		//         = Urr(θ=PI/4)-(1/2)Urr(θ=0)-(1/2)Urr(θ=PI/2)
		return Urr_pol(l,n,0,PI/4)-Urr_pol(l,n,0,0)/2-Urr_pol(l,n,0,PI/2)/2;
	}
	else{
		// Uxy=d/dy(Ux)=d/dy(Ur(dr/dx)+Uth(dth/dx))
	    //    =Ury(dr/dx)+Ur(d2r/dxdy)+Uthy(dth/dx)+Uth(d2th/dx/dy)
		Ury =Urr(l,n,x,y)*(y/r)+Urth(l,n,x,y)*(x/r/r);
		Uthy=Urth(l,n,x,y)*(y/r)+Uthth(l,n,x,y)*(x/r/r);
		return Ury*(x/r)+Ur(l,n,x,y)*(-x*y/r/r/r)+Uthy*(-y/r/r)+Uth(l,n,x,y)*((y*y-x*x)/r/r/r/r);
	}
}

double cZernike::Uxy(int j,double x,double y) const{
	return Uxy(lNumber(j),nNumber(j),x,y);
}

double cZernike::GetR0() const{
	return r0;
}
void cZernike::SetR0(double value){
	r0=value;
	fitted=false;
}

void cZernike::CalcR0(){
	int i,m;
	double x,y, rr,max;
	m=NumberOfData();
	max=0;
	for(i=1; i<=m; i++){
		x=this->x[i];
		y=this->y[i];
		rr=x*x+y*y;
		if(rr>max) max=rr;
	}
	r0=sqrt(max);
}

double cZernike::GetR0Fix() const{
	return r0fix;
}
void cZernike::SetR0Fix(double value){
	r0fix=value;
	fitted=false;
}

int  cZernike::GetNormalize() const{
	return normalize;
}
void cZernike::SetNormalize(int value){
	normalize=value;
	fitted=false;
}

int  cZernike::GetIsFringeOrder() const{
	return IsFringeOrder;
}
void cZernike::SetIsFringeOrder(int value){
	IsFringeOrder=value;
	fitted=false;
}

int cZernike::GetJBase() const{
	return jBase==0 ? 0:1;
}
void cZernike::SetJBase(int value){
	jBase=value;
}

void cZernike::DataClear() {
	x.RemoveAll();
	y.RemoveAll();
	z.RemoveAll();
	i.RemoveAll();
	j.RemoveAll();
	fitted=false;
}

int cZernike::NumberOfData() {
	return x.GetSize();
}

int cZernike::GetNumberOfTerms() const{
	return terms;
}
void cZernike::SetNumberOfTerms(int value){
	if(value>=1){
		matrix<double> temp=X;
		terms=value;
		X.redim(terms,1);
		X.datacopy(temp);   // Xのデータは可能な範囲で維持する
		fitted=false;
	}
}

int cZernike::GetMaxOrder() const{
	return Order(jBase==0 ? terms-1 : terms);
}
void cZernike::SetMaxOrder(int value) {
	if(value>=0){
		SetNumberOfTerms(TotalTerms(value,IsFringeOrder));
		fitted=false;
	}
}

void cZernike::SetData(double x,double y,double z,int i,int j) {
	this->x.AddTail(x);
	this->y.AddTail(y);
	this->z.AddTail(z);
	this->i.AddTail(i);
	this->j.AddTail(j);
	fitted=false;
}

double cZernike::GetXData(int i){
	if(1<=i && i<=NumberOfData()){
		return x[i];
	}
	else return 0;
}

double cZernike::GetYData(int i) {
	if(1<=i && i<=NumberOfData()){
		return y[i];
	}
	else return 0;
}

double cZernike::GetZData(int i) {
	if(1<=i && i<=NumberOfData()){
		return z[i];
	}
	else return 0;
}

int cZernike::GetIData(int i) {
	if(1<=i && i<=NumberOfData()){
		return this->i[i];
	}
	else return 0;
}

int cZernike::GetJData(int i) {
	if(1<=i && i<=NumberOfData()){
		return this->j[i];
	}
	else return 0;
}

void cZernike::SetData(double x,double y,double z) {
	SetData(x,y,z,0,0);
}

void cZernike::SetDigits(int value) {
	Digits=value;
	fitted=false;
}

int cZernike::Fit() {
	// this->x,y,zに係数X[j][1]をフィッティングする．jはjBaseにかかわらず1からtermsとする．

	int i,j, m,n, Rank;

	if(fitted==false){
		m=this->x.GetSize();
		n=terms;
		A.redim(m,n);
		F.redim(m,1);
		
		for(i=1; i<=m; i++) for(j=1; j<=n; j++) {
			A[i][j]=U(lNumber(j,IsFringeOrder,1),nNumber(j,IsFringeOrder,1),x[i],y[i]);
		}
		for(i=1; i<=m; i++) {
			F[i][1]=z[i];
		}
		Rank=rank(t(A)*A);
		if(Rank==terms){
			X=inv( t(A)*A )*t(A)*F;
			for(j=1; j<=n; j++){
				if(Digits>0) X[j][1]=Round(X[j][1],Digits);
			}
			ConvertR0();
			fitted=true;
			return 1;
		}
		else{
			X=zero(X);
			return 0;
		}
	}
	else {
		return 1;
	}
}

double cZernike::GetC(int j) {
	if(this->jBase==0) j+=1;
	if(1<=j && j<=terms) {
		return X[j][1];
	}
	else return 0;
}
void cZernike::SetC(int j,double value) {
	if(this->jBase==0) j+=1;
	if(1<=j && j<=terms) {
		X[j][1]=value;
		fitted=false;
	}
}

double& cZernike::C(int j){
	static double err;

	if(this->jBase==0) j+=1;
	if(1<=j && j<=terms){
		fitted=false;
		return X[j][1];
	}
	else{
		return err;
	}
}

int cZernike::NotNull(){
	int j;
	
	if(r0!=0) return 1;
	for(j=1; j<=terms; ++j){
		if(GetC(j)!=0) return 1;
	}
	return 0;
}


std::string cZernike::GetCoefficients(){
	int j;
	char buf[100];
	std::string s;

	if(NotNull()){
		sprintf(buf,"r0 %g; ", GetR0()); s+=buf;
		for(j=1; j<=terms; ++j){
			if(GetC(j)!=0){
				sprintf(buf,"C %d %g; ", j,GetC(j));
				s+=buf;
			}
		}
	}

	return s;
}

void cZernike::SetCoefficients(std::string com){
	// 例： com="Normalize 0; r0 1.5; C 4 0.1; C 12 0.1"

	std::string s,s0;
	int i,j, j_max;
	std::string s1,s2;
	bool b1,b2;

	// 最初にcomに含まれる最大の項番号を調べて，メモリを確保しておく必要がある．
	j_max=0;
	for(i=1; i<=sentences(com); ++i){
		s=sentence(com,i);
		s0=arg(s,0);

		if(s0=="C" || s0=="c"){
			int j;
			s1=arg(s,1); b1=is_numeric(s1);
			s2=arg(s,2); b2=is_numeric(s2);
			if(b1 && b2){
				j=atoi(s1.c_str());
				if(j>j_max) j_max=j;
			}
		}
	}
	if(j_max+1-jBase > terms) SetNumberOfTerms(j_max+1-jBase);

	for(j=1; j<=terms; ++j) SetC(j,0);

	for(i=1; i<=sentences(com); ++i){
		s=sentence(com,i);
		s0=arg(s,0);

		if(s0=="Normalize" || s0=="normalize"){
			int val;
			s1=arg(s,1); b1=is_numeric(s1);
			if(s1=="?"){
				s+="Normalize val\n";
			}
			else if(b1){
				val=atoi(s1.c_str());
				SetNormalize(val);
			}
		}
		if(s0=="r0" || s0=="R0"){
			double val;
			s1=arg(s,1); b1=is_numeric(s1);
			if(s1=="?"){
				s+="r0 val\n";
			}
			else if(b1){
				val=atof(s1.c_str());
				SetR0(val);
			}
		}
		if(s0=="C" || s0=="c"){
			int j;
			double val;
			s1=arg(s,1); b1=is_numeric(s1);
			s2=arg(s,2); b2=is_numeric(s2);
			if(s1=="?"){
				s+="C j val\n";
			}
			else if(b1 && b2){
				j=atoi(s1.c_str());
				val=atof(s2.c_str());
				SetC(j,val);
			}
		}
	}
}

void cZernike::FitZToCoefficient(){
	// zをZernike係数で表される形状とする
	int i;

	for(i=1; i<=NumberOfData(); ++i){
		this->z[i]=Z(GetXData(i),GetYData(i));
	}
	fitted=true;
}

double cZernike::Z(double x,double y,double Rref) const{
	// Rref = 基準球面の曲率半径
	//        例えば，干渉計のZernike係数から偏差ではなく絶対形状を計算するのに使う．
	int j;
	double z,c,rr;

	z=0;
	for(j=1; j<=terms; ++j){
		if(X[j][1]!=0){
			z+=X[j][1]*U(lNumber(j,IsFringeOrder,1),nNumber(j,IsFringeOrder,1),x,y);
		}
	}
	if(Rref!=0){
		c=1/Rref;
		rr=x*x+y*y;
		z+=c*rr/(1+sqrt(1-c*c*rr));
	}
	return z;
}

double cZernike::Z(double x,double y) const{
	return Z(x,y,0);
}

double cZernike::ZApproximate(double x,double y,double Rref) {
	Fit();
	return Z(x,y,Rref);
}

double cZernike::ZApproximate(double x,double y) {
	return ZApproximate(x,y,0);
}

double cZernike::Zr(double x,double y,double Rref,int RemoveTilt) const{
	int j;
	double zra=0, c,r;

	for(j= RemoveTilt==0 ? 2:4 ; j<=terms; ++j){
		if(X[j][1]!=0){
			zra+=X[j][1]*Ur(lNumber(j,IsFringeOrder,1),nNumber(j,IsFringeOrder,1),x,y);
		}
	}
	if(Rref!=0){
		c=1/Rref;
		r=sqrt(x*x+y*y);
		zra+=c*r/(sqrt(1-c*c*r*r));
	}
	return zra;
}

double cZernike::Zx(double x,double y) const{
	int j;
	double zx=0;

	for(j=1 ; j<=terms; ++j){
		if(X[j][1]!=0){
			zx+=X[j][1]*Ux(lNumber(j,IsFringeOrder,1),nNumber(j,IsFringeOrder,1),x,y);
		}
	}
	return zx;
}

double cZernike::Zy(double x,double y) const{
	int j;
	double zy=0;

	for(j=1 ; j<=terms; ++j){
		if(X[j][1]!=0){
			zy+=X[j][1]*Uy(lNumber(j,IsFringeOrder,1),nNumber(j,IsFringeOrder,1),x,y);
		}
	}
	return zy;
}

double cZernike::Zxx(double x,double y) const{
	int j;
	double zxx=0;

	for(j=1 ; j<=terms; ++j){
		if(X[j][1]!=0){
			zxx+=X[j][1]*Uxx(lNumber(j,IsFringeOrder,1),nNumber(j,IsFringeOrder,1),x,y);
		}
	}
	return zxx;
}

double cZernike::Zyy(double x,double y) const{
	int j;
	double zyy=0;

	for(j=1 ; j<=terms; ++j){
		if(X[j][1]!=0){
			zyy+=X[j][1]*Uyy(lNumber(j,IsFringeOrder,1),nNumber(j,IsFringeOrder,1),x,y);
		}
	}
	return zyy;
}


double cZernike::Zxy(double x,double y) const{
	int j;
	double zxy=0;

	for(j=1 ; j<=terms; ++j){
		if(X[j][1]!=0){
			zxy+=X[j][1]*Uxy(lNumber(j,IsFringeOrder,1),nNumber(j,IsFringeOrder,1),x,y);
		}
	}
	return zxy;
}

double cZernike::AxialR(double x,double y,double Rref,int RemoveTilt){
	// 多項式が定義する面のAxial曲率半径Raを計算する．
	// 一般曲面を極座標により，
	//     Z=(r,θ)
	// と表すとき，
	//     Ra=r/[ sin{ atan(dZ/dr) } ]
	// とする．
	double zr;

	zr=Zr(x,y,Rref,RemoveTilt);
	return zr==0 ? 0 : sqrt(x*x+y*y)/sin(atan(zr));
}

void cZernike::ParaxialExpansion(double &cx,double &cxy,double &cy){
	//  zの近軸展開，
	//      z = (cx/2)*x^2 +cxy*x*y +(cy/2)*y^2
	//  の係数を求める．
	//  参考書：“光学の原理Ⅱ(W.Born E.Wolf 草川徹訳)”
	
	// ここではテーラー展開したときの2次の係数を求めている．
	// 最もフィットした（それを除去したときの残余rmsが最小の)
	// 2次曲面とは違うものである．
	// したがって，眼の収差を最も補償する眼鏡のS,C,Aを求める，
	// という問題であれば素直にZernike展開の2次の係数だけを
	// 見た方がよい．

	int s,m,l, j,n;
	double b;

	if(r0==0) return;
	cx=cy=cxy=0;
	
	for(j=1; j<=terms; ++j){
		l=lNumber(j);
		n=nNumber(j);
		m= l>=0 ? l : -l;

		for(s=0; s<=(n-m)/2; ++s){       // 9.2.1 (5)式
			if( n-2*s==2 ){              // n-2*s : rの次数
				b=(is_even(s) ? 1: -1)
				  *factorial(n-s)/factorial(s)/factorial((n+m)/2-s)/factorial((n-m)/2-s);
				if(normalize) b*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
				switch(l){
				case 0:  // r^2=x^2+y^2, C*b*(x^2+y^2) = (cx/2)*x^2+(cy/2)*y^2
					cx+=b*2*GetC(j);
					cy+=b*2*GetC(j);
					break;
				case -2: // r^2*sin(2θ)=2*x*y, C*b*2*x*y = cxy*x*y
					cxy+=b*2*GetC(j);
					break;
				case 2:  // r^2*cos(2θ)=x^2-y^2, C*b*(x^2-y^2)=(cx/2)*x^2+(cy/2)*y^2
					cx+=b*2*GetC(j);
					cy-=b*2*GetC(j);
					break;
				}
			}
		}
	}

	cx *=(1/r0/r0);
	cxy*=(1/r0/r0);
	cy *=(1/r0/r0);
}

void cZernike::ParaxialR(double &r1,double &r2,double &axis_deg,double Rref){
	double cx,cxy,cy,axx,axy,ayy;

	ParaxialExpansion(cx,cxy,cy);
	// (cx/2)*x^2+cxy*x*y+(cy/2)*y^2 = axx*x^2+axy*x*y+ayy*y^2
	axx=cx/2;
	axy=cxy;
	ayy=cy/2;
	axx+= Rref==0 ? 0 : 1.0/2.0/Rref;
	ayy+= Rref==0 ? 0 : 1.0/2.0/Rref;
	::ParaxialR(r1,r2,axis_deg,axx,ayy,axy);
}

double cZernike::ParaxialR1(double Rref){
	double r1,r2,axis;

	ParaxialR(r1,r2,axis,Rref);
	return r1;
}

double cZernike::ParaxialR2(double Rref){
	double r1,r2,axis;
	
	ParaxialR(r1,r2,axis,Rref);
	return r2;
}

double cZernike::ParaxialRaxis(){
	double rx,ry,axis_deg;
	
	ParaxialR(rx,ry,axis_deg,0);
	return axis_deg;
}

double cZernike::Error(int i) {
	if(1<=i && i<=NumberOfData()) {
		double x,y,z;
		x=this->x[i];
		y=this->y[i];
		z=this->z[i];
		return ZApproximate(x,y)-z;
	}
	else return 0;
}

double cZernike::RMSError() {
	int i;
	double x,y,z;
	int m=NumberOfData();
	int n=terms;
	double *za=new double[m+1];
	double error;
	for(i=1; i<=m; i++){
		x=this->x[i];
		y=this->y[i];
		za[i]=ZApproximate(x,y);
	}
	error=0;
	for(i=1; i<=m; ++i){
		z=this->z[i];
		error+=( za[i]-z )*( za[i]-z );
	}
	error=sqrt(error/m);
	delete [] za;
	return error;
}

double cZernike::PVError() {
	int i;
	double x,y,z;
	int m=NumberOfData();
	int n=terms;
	double *za=new double[m+1];
	double peak=-1e30;
	double valley=1e30;
	for(i=1; i<=m; ++i){
		x=this->x[i];
		y=this->y[i];
		za[i]=ZApproximate(x,y);
		z=this->z[i];
		if(za[i]-z>peak) peak=za[i]-z;
		if(za[i]-z<valley) valley=za[i]-z;
	}
	delete [] za;
	return peak-valley;
}

int cZernike::RemoveTerms(int piston,int tilt,int sph) {
	// フィッティング後，this->zよりpiston,tilt,sph項の分を除去する．
	//     piston<>0  piston成分を除去
	//     tilt  <>0  tilt成分を除去
	//     sph   <>0  sph成分を除去
	// 展開次数Ordersが1次以下のときは何もしない．
	if(Fit() && terms>=4) {
		int i;
		double x,y,z;

		for(i=1; i<=NumberOfData(); ++i){
			x=this->x[i];
			y=this->y[i];
			z=this->z[i];
			if(piston){
				z-=X[jNumber( 0,0)-jBase+1][1]*U( 0,0,x,y);
			}
			if(tilt){
				z-=X[jNumber(-1,1)-jBase+1][1]*U(-1,1,x,y);
				z-=X[jNumber( 1,1)-jBase+1][1]*U( 1,1,x,y);
			}
			if(sph){
				z-=X[jNumber( 0,2)-jBase+1][1]*U( 0,2,x,y);
			}
			this->z[i]=z;
		}
		fitted=false;
		return 1;
	}
	else{
		return 0;
	}
}

int cZernike::AdjustTerms(int piston,int tilt,int sph) {
	// フィッティング後，
	// piston,tilt,sph成分を調整して，this->zを分散がx,y分布領域内で最小になるようにする．
	//     tilt=tureならばpistonも，
	//     sph=trueならばpiston,tiltも
	// 調整する．
	// 例えば，x,yの分布領域が円からずれている(楕円など)とき，
	// Zernike多項式はその領域では直交しない．
	// したがって，RemoveTerms()関数で，
	//     zからpiston項の分を除去してもzの分布が正負対称にならない．
	//     zからsph項の分を除去しても却ってrmsが増える．
	// ことがあったため本関数を作成した． (080910)
	//
	// 【注意】
	//    zが波面収差，x,yが瞳座標であるとき，
	//    瞳面と主光線が垂直でないとき(例えばFTのウエハ面を瞳としたとき),
	//    sph成分は光の進行方向から見れば非点収差を含む．
	//    したがって，sphの調整により非点収差も調整され，
	//    収差が実際より良く見えることがある．

	int i,j, m,n;
	matrix<double> A1,X1,F1;
	double x,y,z;

	if(Fit() && terms>=4){
		// データ数m, 項数nの設定．行列用メモリの確保．
		m=this->x.GetSize();	
		n=0;
		if(piston) n=1;
		if(tilt)   n=3;
		if(sph)    n=4;
		if(n==0) return 0;
		A1.redim(m,n);
		F1.redim(m,1);
		// 最小二乗法の偏微分係数行列の設定
		for(i=1; i<=m; ++i){
			x=this->x[i];
			y=this->y[i];
			j=1;
			if(piston || tilt || sph){
				A1[i][j]=U( 0,0,x,y); j++;
			}
			if(tilt || sph){
				A1[i][j]=U(-1,1,x,y); j++;
				A1[i][j]=U( 1,1,x,y); j++;
			}
			if(sph){
				A1[i][j]=U( 0,2,x,y); j++;
			}
		}
		// 最小二乗法の目標値はz=0;
		for(i=1; i<=m; ++i){
			F1[i][1]=-(this->z[i]);
		}
		// 最小二乗法の実行
		if( rank(t(A1)*A1)<n ) return 0; // 非正則ならば何もしないで終了
		X1=inv( t(A1)*A1 )*t(A1)*F1;
		// 係数を修正する．(X[]の添字はピストンが1であることに注意)
		if(piston || tilt || sph){
			X[jNumber( 0,0)-jBase+1][1]+=X1[1][1];
		}
		if(tilt || sph){
			X[jNumber(-1,1)-jBase+1][1]+=X1[2][1];
			X[jNumber( 1,1)-jBase+1][1]+=X1[3][1];
		}
		if(sph){
			X[jNumber( 0,2)-jBase+1][1]+=X1[4][1];
		}
		// zを新しい係数を使った展開式の値にする．
		for(i=1; i<=m; ++i){
			x=this->x[i];
			y=this->y[i];
			// ここで，z=zApproximate(x,y) というのも考えられるが，
			// こうするとRMSError(),PVError()が0になってしまう．（すなわち高次の収差が失われる．）
			// フィッティングエラー量を保存するため次のようにする．
			this->z.GetData(z,i);
			if(piston || tilt || sph){
				z+=X1[1][1]*U( 0,0,x,y);
			}
			if(tilt || sph){
				z+=X1[2][1]*U(-1,1,x,y);
				z+=X1[3][1]*U( 1,1,x,y);
			}
			if(sph){
				z+=X1[4][1]*U( 0,2,x,y);
			}
			this->z[i]=z;
		}

		fitted=false;
		return 1;
	}
	else{
		return 0;
	}
}

double cZernike::RMS(int tilt,int sph,int cyl,int high_order){
	// tilt !=0       のときはtilt項(n==1)を加える
	// sph  !=0       のときはsph項(n==2;l==0)を加える
	// cyl  !=0       のときはcyl項(n==2;l==-2,2)を加える
	// high_order !=0 のときは高次項(n>=3)を加える
	// (直交性を使っていないので，x,yの分布領域が円でなくても有効)
	if(this->x.GetSize()==0) return 0;
	int i,j;
	double x,y,z, zzsum=0;
	for(i=1; i<=this->x.GetSize(); ++i){
		z=0;
		x=this->x[i];
		y=this->y[i];
		if(tilt!=0){
			z+=GetC(2)*U(lNumber(2),nNumber(2),x,y);
			z+=GetC(3)*U(lNumber(3),nNumber(3),x,y);
		}
		if(sph!=0){
			z+=GetC(5)*U(lNumber(5),nNumber(5),x,y);
		}
		if(cyl!=0){
			z+=GetC(4)*U(lNumber(4),nNumber(4),x,y);
			z+=GetC(6)*U(lNumber(6),nNumber(6),x,y);
		}
		if(high_order!=0){
			for(j=7; j<=terms; j++){
				z+=GetC(j)*U(lNumber(j),nNumber(j),x,y);
			}
		}
		zzsum+=z*z;
	}
	return sqrt(zzsum/this->x.GetSize());
}

double cZernike::Shear(double dx,double dy,int order){
	// *thisのｚデータとdx,dyだけずらした場所のzデータとの差をとり，*thisへ代入する．
	// x,yデータは*thisと同じになる．
    // ずらしたzデータにはZernike展開を用いる．この展開のpvエラーを返す．
	// 展開のオーダーを指定するため，コピーを使って計算する．
	// zが波面収差であればShearing干渉計のシミュレーションとなる．
	cZernike a;
	int i;
	double pv;

	a=*this;
	a.SetMaxOrder(order);
	pv=a.PVError();

	for(i=1; i<=a.NumberOfData(); ++i){
		this->z[i]=a.z[i]-a.ZApproximate(a.x[i]+dx,a.y[i]+dy);
	}
	
	AdjustTerms(1,1,0);   // piston,tilt成分は除去する
	this->fitted=false;
	return pv;
}




