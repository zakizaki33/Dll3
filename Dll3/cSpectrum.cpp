#include "stdafx.h"
#include "cSpectrum.h"

double cSpectrum::WLpb1=380, cSpectrum::WLpb2=780;  // 純紫軌跡端の波長

double cSpectrum::Tcp(double& duv){
	// JIS Z 8725(2015) 付録B.2 により相関色温度を計算する
	
	const int N=65;
	double UT,VT, UN[N+1],VN[N+1],RSP[N+1], DT[N+1],DDD,Tm,Tc,PPP,QQQ,UON,VON,DUV;
	int J;

	UT=u(); VT=v();

	UN[ 1]=.180046;  VN[ 1]=.263577;  RSP[ 1]=-4.09562;
	UN[ 2]=.180638;  VN[ 2]=.265948;  RSP[ 2]=-3.91330;
	UN[ 3]=.181309;  VN[ 3]=.268506;  RSP[ 3]=-3.71055;
	UN[ 4]=.182067;  VN[ 4]=.271236;  RSP[ 4]=-3.49518;
	UN[ 5]=.182919;  VN[ 5]=.274118;  RSP[ 5]=-3.72420;

	UN[ 6]=.183872;  VN[ 6]=.277131;  RSP[ 6]=-3.05386;
	UN[ 7]=.184932;  VN[ 7]=.280251;  RSP[ 7]=-2.83890;
	UN[ 8]=.186103;  VN[ 8]=.283452;  RSP[ 8]=-2.63279;
	UN[ 9]=.187389;  VN[ 9]=.286709;  RSP[ 9]=-2.43778;
	UN[10]=.188792;  VN[10]=.289997;  RSP[10]=-2.25517;
	
	UN[11]=.190312;  VN[11]=.293293;  RSP[11]=-2.08544;
	UN[12]=.191949;  VN[12]=.296575;  RSP[12]=-1.92856;
	UN[13]=.193701;  VN[13]=.299825;  RSP[13]=-1.78409;
	UN[14]=.195566;  VN[14]=.303025;  RSP[14]=-1.65136;
	UN[15]=.197540;  VN[15]=.306162;  RSP[15]=-1.52956;
	
	UN[16]=.199619;  VN[16]=.309223;  RSP[16]=-1.41784;
	UN[17]=.201799;  VN[17]=.312199;  RSP[17]=-1.31534;
	UN[18]=.204074;  VN[18]=.315083;  RSP[18]=-1.22121;
	UN[19]=.206440;  VN[19]=.317868;  RSP[19]=-1.13468;
	UN[20]=.208891;  VN[20]=.320550;  RSP[20]=-1.05503;
	
	UN[21]=.211423;  VN[21]=.323126;  RSP[21]=-.98159;
	UN[22]=.214030;  VN[22]=.325595;  RSP[22]=-.91377;
	UN[23]=.216706;  VN[23]=.327956;  RSP[23]=-.85104;
	UN[24]=.219449;  VN[24]=.330208;  RSP[24]=-.79290;
	UN[25]=.222251;  VN[25]=.332354;  RSP[25]=-.73895;
	
	UN[26]=.225110;  VN[26]=.334393;  RSP[26]=-.68880;
	UN[27]=.228020;  VN[27]=.336329;  RSP[27]=-.64211;
	UN[28]=.230978;  VN[28]=.338163;  RSP[28]=-.59859;
	UN[29]=.233979;  VN[29]=.339897;  RSP[29]=-.55795;
	UN[30]=.237020;  VN[30]=.341536;  RSP[30]=-.51998;
	
	UN[31]=.240097;  VN[31]=.343080;  RSP[31]=-.48444;
	UN[32]=.243206;  VN[32]=.344534;  RSP[32]=-.45115;
	UN[33]=.246345;  VN[33]=.345901;  RSP[33]=-.41994;
	UN[34]=.249511;  VN[34]=.347183;  RSP[34]=-.39065;
	UN[35]=.252699;  VN[35]=.348384;  RSP[35]=-.36315;
	
	UN[36]=.255909;  VN[36]=.349508;  RSP[36]=-.33729;
	UN[37]=.259136;  VN[37]=.350557;  RSP[37]=-.31298;
	UN[38]=.262379;  VN[38]=.351534;  RSP[38]=-.29010;
	UN[39]=.265635;  VN[39]=.352443;  RSP[39]=-.26855;
	UN[40]=.268902;  VN[40]=.353287;  RSP[40]=-.24826;
	
	UN[41]=.272179;  VN[41]=.354069;  RSP[41]=-.22914;
	UN[42]=.275462;  VN[42]=.354791;  RSP[42]=-.21111;
	UN[43]=.278750;  VN[43]=.355457;  RSP[43]=-.19410;
	UN[44]=.282042;  VN[44]=.356070;  RSP[44]=-.17806;
	UN[45]=.285335;  VN[45]=.356631;  RSP[45]=-.16293;
	
	UN[46]=.288629;  VN[46]=.357144;  RSP[46]=-.14864;
	UN[47]=.291922;  VN[47]=.357611;  RSP[47]=-.13515;
	UN[48]=.295211;  VN[48]=.358034;  RSP[48]=-.12241;
	UN[49]=.298497;  VN[49]=.358417;  RSP[49]=-.11038;
	UN[50]=.301778;  VN[50]=.358760;  RSP[50]=-.09902;
	
	UN[51]=.305053;  VN[51]=.359066;  RSP[51]=-.08828;
	UN[52]=.308320;  VN[52]=.359338;  RSP[52]=-.07814;
	UN[53]=.311579;  VN[53]=.359577;  RSP[53]=-.06856;
	UN[54]=.314829;  VN[54]=.359785;  RSP[54]=-.05950;
	UN[55]=.318068;  VN[55]=.359964;  RSP[55]=-.05094;
	
	UN[56]=.321297;  VN[56]=.360115;  RSP[56]=-.04285;
	UN[57]=.324514;  VN[57]=.360240;  RSP[57]=-.03520;
	UN[58]=.327718;  VN[58]=.360342;  RSP[58]=-.02797;
	UN[59]=.330909;  VN[59]=.360420;  RSP[59]=-.02114;
	UN[60]=.334087;  VN[60]=.360477;  RSP[60]=-.01467;
	
	UN[61]=.337250;  VN[61]=.360513;  RSP[61]=-.00856;
	UN[62]=.340397;  VN[62]=.360531;  RSP[62]=-.00279;
	UN[63]=.343530;  VN[63]=.360531;  RSP[63]=.00267;
	UN[64]=.346646;  VN[64]=.360515;  RSP[64]=.00784;
	UN[65]=.349746;  VN[65]=.360483;  RSP[65]=.01272;

	J=0;

	do{
		J+=1;
		if(J==N+1) return 0;
		DT[J]=((UN[J]-UT)-RSP[J]*((VN[J]-VT)))/sqrt(1+RSP[J]*RSP[J]);
	} while(DT[J]<0);

	if(J<=2) return 0;

	DDD=-DT[J-1]/(DT[J]-DT[J-1]);
	Tm=J-2+DDD;
	Tc=100000/Tm;

	PPP=(UN[J]+UN[J-2])/2-UN[J-1];
	QQQ=(UN[J]-UN[J-2])/2;
	UON=(PPP*DDD+QQQ)*DDD+UN[J-1];

	PPP=(VN[J]+VN[J-2])/2-VN[J-1];
	QQQ=(VN[J]-VN[J-2])/2;
	VON=(PPP*DDD+QQQ)*DDD+VN[J-1];

	DUV=sqrt((UT-UON)*(UT-UON)+(VT-VON)*(VT-VON));
	if(VT<VON) DUV=-DUV;

	duv=DUV;
	return Tc;
}

// ---- public members ----

double cSpectrum::StandardIlluminant(std::string name,double wl){
		cSpectrum x; cXYList li;
		if(name=="A")     li=x.A;
		if(name=="C")     li=x.C;
		if(name=="D65")   li=x.D65;
		if(name=="WHITE") li=x.WHITE;
		return li.y(wl);
}

double cSpectrum::BlackBodyRadiation(double temp,double wl_nm,double wl0_nm){
	// wl_nm  = wavelength[nm]
	// temp   = temperature[K]
	// return = black body radiant intensity normalized as 1 at wl0_nm
	double wl;                      // wave length[m]
	double s0,s;                    // relative spectral intensity
	const double C2=0.014388;       // [mK]
	wl = wl0_nm/1e9;
	s0 = 1/(wl*wl*wl*wl*wl)/( exp(C2/wl/temp)-1 );
	wl = wl_nm/1e9;
	s  = 1/(wl*wl*wl*wl*wl)/( exp(C2/wl/temp)-1 );
	return s/s0;
}

double cSpectrum::LuminousEfficiency(double wl_nm){
	cSpectrum x;
	return x.y_bar.y(wl_nm);
}

void cSpectrum::SpectrumLocus(double wl_nm,double& xp,double & yp){
	// 波長wl_nmに対応するスペクトル軌跡上の点(xp,yp)を求める．
	// wl_nmが負のときは|wl_nm|に対応するスペクトル軌跡上の点と白色点を結ぶ直線が
	// 純紫軌跡と交わる点(xp,yp)を求める．
	static cSpectrum s;
	double wl, X,Y,Z, xw,yw,x1,y1,x2,y2,x,y, delta;

	wl= wl_nm>=0 ? wl_nm : -wl_nm;
	// 380nm～780nmの外では等色関数のデータがない
	if(wl<WLpb1) wl=WLpb1;
	if(wl>WLpb2) wl=WLpb2;
	X=s.x_bar.y(wl); Y=s.y_bar.y(wl); Z=s.z_bar.y(wl);
	xp=X/(X+Y+Z); yp=Y/(X+Y+Z);
	if(wl_nm<0){
		// 白色点(xw,yw)
		WhitePoint(xw,yw);
		// 純紫軌跡両端点(x1,y1),(x2,y2)
		SpectrumLocus(WLpb1,x1,y1);
		SpectrumLocus(WLpb2,x2,y2);
		// (xp,yp),(xw,yw)を通る直線 (yp-yw)x+(xw-xp)y=xwyp-xpyw
		// (x1,y1),(x2,y2)を通る直線 (y2-y1)x+(x1-x2)y=x1y2-x2y1
		// 両直線の交点(x,y)は，
		// x=(1/delta){(x1-x2)(xwyp-xpyw)+(xp-xw)(x1y2-x2y1)}
		// y=(1/delta){(y1-y2)(xwyp-xpyw)+(yp-yw)(x1y2-x2y1)}
		// delta=(yp-yw)(x1-x2)-(xw-xp)(y2-y1)
		delta=(yp-yw)*(x1-x2)-(xw-xp)*(y2-y1);
		if(delta==0){
			xp=yp=0;
		}
		else{
			x=(1/delta)*( (x1-x2)*(xw*yp-xp*yw)+(xp-xw)*(x1*y2-x2*y1) );
			y=(1/delta)*( (y1-y2)*(xw*yp-xp*yw)+(yp-yw)*(x1*y2-x2*y1) );
			xp=x;
			yp=y;
		}
	}
}

double cSpectrum::SpectrumLocusX(double wl_nm){
	double x,y;
	SpectrumLocus(wl_nm,x,y);
	return x;
}

double cSpectrum::SpectrumLocusY(double wl_nm){
	double x,y;
	SpectrumLocus(wl_nm,x,y);
	return y;
}

void cSpectrum::WhitePoint(double& xw,double& yw){
	// 白色点=分光分布が平坦な光の色度点
	static cSpectrum s;
	static double x=0,y=0;
	if(x==0 && y==0){ // 高速化のため．呼び出し回数の多いプログラムでは効果がある．080828
		s.spectrum=s.WHITE;
		x=s.x();
		y=s.y();
	}
	xw=x;
	yw=y;
}

double cSpectrum::XYToDominantWavelength(double x,double y){
	// 色度座標(x,y)に対応する主波長を返す．
	//   (x,y)が白色点に対し純紫軌跡側のときは補色主波長を負値で返す．
	//   (x,y)がスペクトル軌跡，純紫軌跡の外のときは0を返す．
	// 主波長とは(x,y)と白色点を結ぶ直線と純色軌跡の交点の波長．
	// (x,y)がスペクトル軌跡，純紫軌跡の内側かどうかの判断にも利用できる．
	
	static cSpectrum s;
	int i;
	double wl,wlc, wl1,x1,y1, wl2,x2,y2, l,m,n, xw,yw;
	matrix<double> A(3,3),C(3,1),X(3,1);

	WhitePoint(xw,yw);

	// 明らかにスペクトル軌跡，純紫軌跡より外の場合は0を返して終了
	if(x<0 || y<0 || 1<x+y) return 0;
	
	// (x,y)が純紫色線より下の場合は0を返して終了
	SpectrumLocus(WLpb1,x1,y1);
	SpectrumLocus(WLpb2,x2,y2);
	if( y < (y2-y1)/(x2-x1)*x + y1 -(y2-y1)/(x2-x1)*x1 ) return 0;

	wl=wlc=0;
	for(i=1; i<=s.x_bar.GetSize()-1; i++){
		wl1=s.x_bar.x(i);
		wl2=s.x_bar.x(i+1);
		if(WLpb1<=wl2 && wl1<=WLpb2){
			SpectrumLocus(wl1,x1,y1);
			SpectrumLocus(wl2,x2,y2);
			//    |x1     |    |x2     |    |xw     |   |x    |
			//   l|y1     | + m|y2     | + n|yw     | = |y    |
			//    |1-x1-y1|    |1-x2-y2|    |1-xw-yw|   |1-x-y|
			// として，l,m,nを解く．
			A[1][1]=x1;      A[1][2]=x2;      A[1][3]=xw;       C[1][1]=x;
			A[2][1]=y1;      A[2][2]=y2;      A[2][3]=yw;       C[2][1]=y;
			A[3][1]=1-x1-y1; A[3][2]=1-x2-y2; A[3][3]=1-xw-yw;  C[3][1]=1-x-y;
			X=inv(A)*C;
			l=X[1][1];
			m=X[2][1];
			n=X[3][1];
			if( 0<=l*m ){ // (x,y)と(xw,yw)を結ぶ直線が(x1,y1)-(x2,y2)を内分する
				if( n<0 ){
					return 0; // スペクトル軌跡より外なので0を返して終了
				}
				if( 0<=n && n<1 ){
					wl=(wl1*l+wl2*m)/(l+m);  // 主波長
				}
				if( n==1 ){
					return 0; // 白色点なので0を返して終了
				}
				if( 1<n ){
					wlc=(wl1*l+wl2*m)/(l+m); // 補色主波長
				}
			}
		}
	}
	// 補色主波長のみ求まった場合(白色点に対し純紫色側)，補色主波長を負値で返す
	if( wl==0 && wlc!=0 ) wl=-wlc;

	return wl;
}

double cSpectrum::ExcitationPurity(double x,double y){
	double wl, xw,yw, xp,yp;
	WhitePoint(xw,yw);
	wl=XYToDominantWavelength(x,y);
	SpectrumLocus(wl,xp,yp);
	return sqrt( ((x-xw)*(x-xw)+(y-yw)*(y-yw))/((xp-xw)*(xp-xw)+(yp-yw)*(yp-yw)) );
}

long cSpectrum::XYToApproximateRGB(double x,double y) {
	// 表したい色の色度座標(x,y)に対するRGB値を返す．
	matrix<double> A(3,3),C(3,1),X(3,1);
	static matrix<double> invA(3,3);
	double rx,ry,gx,gy,bx,by,wx,wy, l,m,n, r,g,b, max;
	long R,G,B;
	static bool Calculated=false;
	const double Gamma=2.2;

	if(!Calculated){
		// RGB(255,0,0), RGB(0,255,0), RGB(0,0,255), RGB(255,255,255)
		// のときのモニタ色の色度座標をそれぞれ，
		rx=0.640; gx=0.300; bx=0.150; wx=0.3127;
		ry=0.330; gy=0.600; by=0.060; wy=0.3290;
		// とする(sRGB規格 IEC61966-2-1(1999))．
		// モニタ色が灰色(R=G=B)のとき，XYZ三刺激値の関係，
		//    |rx     |    |gx     |    |bx     |   |wx     |
		//   l|ry     | + m|gy     | + n|by     | = |wy     |
		//    |1-rx-ry|    |1-gx-gy|    |1-bx-by|   |1-wx-wy|
		// より，
		A[1][1]=rx;      A[1][2]=gx;      A[1][3]=bx;       C[1][1]=wx;
		A[2][1]=ry;      A[2][2]=gy;      A[2][3]=by;       C[2][1]=wy;
		A[3][1]=1-rx-ry; A[3][2]=1-gx-gy; A[3][3]=1-bx-by;  C[3][1]=1-wx-wy;
		X=inv(A)*C;
		l=X[1][1];
		m=X[2][1];
		n=X[3][1];
		// のように l,m,n が決まる．
		// したがって，モニタ色色度座標が(x,y)のときのXYZ三刺激値の関係，
		//    |lrx       |    |mgx       |    |nbx       |   |x    |
		//   r|lry       | + g|mgy       | + b|nby       | = |y    |
		//    |l(1-rx-ry)|    |m(1-gx-gy)|    |n(1-bx-by)|   |1-x-y|
		// より，
		A[1][1]=l*rx;        A[1][2]=m*gx;        A[1][3]=n*bx;
		A[2][1]=l*ry;        A[2][2]=m*gy;        A[2][3]=n*by;
		A[3][1]=l*(1-rx-ry); A[3][2]=m*(1-gx-gy); A[3][3]=n*(1-bx-by);
		invA=inv(A);
		Calculated=true;
	}
	C[1][1]=x;
	C[2][1]=y;
	C[3][1]=1-x-y;
	X=invA*C;
	r=X[1][1];
	g=X[2][1];
	b=X[3][1];
	// のように r,g,b が決まる．
	// r,g,bのどれかが負になると色三角形の外である．負の値を0とすることにより圧縮する．
	if(r<0) r=0;
	if(g<0) g=0;
	if(b<0) b=0;
	// r,g,bの最大を1に規格化する．
	max= r>g ? r:g;
	max= max>b ? max:b;
	r=r/max;
	g=g/max;
	b=b/max;
	// ここまでのr,g,bは光量である．ガンマを補正する．
	r=pow(r,1/Gamma);
	g=pow(g,1/Gamma);
	b=pow(b,1/Gamma);
	R=(long)(r*255);
	G=(long)(g*255);
	B=(long)(b*255);
	return B*256*256 +G*256 +R;
}

long cSpectrum::ApproximateColor(double wl_nm) {
	double x,y;
	static double wl_cache=0;
	static long color_cache=0;
	if( wl_nm==wl_cache ) return color_cache;
	SpectrumLocus(wl_nm,x,y);
	wl_cache=wl_nm;
	return color_cache=XYToApproximateRGB(x,y);
}

double cSpectrum::LumenToWatt(double wl_nm,double flux_lumen){
	// 単色光光束をlumen単位からWatt単位に変換する
	cSpectrum s0,s;
	double K;

	// wl_nmにおける，1Wあたりのlumen数は，
	s0.AddData(555,1);
	s.AddData(wl_nm,1);
	K=s.y()/s0.y()*683;   // 1Watt=683lumen, at 555nm

	return flux_lumen/K;
}

double cSpectrum::WattToLumen(double wl_nm,double flux_watt){
	// 単色光光束をwatt単位からlumen単位に変換する

	return flux_watt/cSpectrum::LumenToWatt(wl_nm,1);
}


cSpectrum::cSpectrum() {
	make_constants();
}

cSpectrum::cSpectrum(const cSpectrum& x) {
	make_constants();
	spectrum=x.spectrum;
}

cSpectrum& cSpectrum::operator=(const cSpectrum& x) {
	make_constants();
	spectrum=x.spectrum;
	return *this;
}

int cSpectrum::Size() { 
	return spectrum.GetSize(); 
}

double cSpectrum::WaveLength(int i) { 
	return spectrum.x(i);
}

double cSpectrum::Data(double wl) { 
	return spectrum.y(wl);
}

cSpectrum& cSpectrum::AddData(double wl,double value){
	spectrum.AddData(wl,value);
	return *this;
}

cSpectrum& cSpectrum::AddConstantData(double wl_start,double wl_end,double wl_step,double value){
	double wl=wl_start;
	while( wl <= wl_end+wl_step*0.001 ){
		spectrum.AddData(wl,value);
		wl+=wl_step;
	}
	return *this;
}

cSpectrum& cSpectrum::RemoveAllData() {
	spectrum.RemoveAll();
	return *this;
}

cSpectrum cSpectrum::Copy() const {
	cSpectrum x=*this;
	return x;
}

cSpectrum& cSpectrum::Apply(cXYList &x) {
	spectrum.Apply(x);
	return *this;
}

int cSpectrum::Apply(const std::string& filename){
	cSpectrum x;
	if( x.Open(filename) ){
		Apply(x.spectrum);
		return 1;
	}
	else{
		return 0;
	}
}

cSpectrum& cSpectrum::ApplyBlackBodyRadiation(double temp) {
	int i;
	cXYData x;
	if(temp>0){
		for(i=1; i<=spectrum.GetSize(); ++i){
			spectrum.GetData(x,i);
			x.y*=BlackBodyRadiation(temp,x.x,560);
			spectrum.SetData(x,i);
		}
	}
	return *this;
}

cSpectrum& cSpectrum::ApplyLuminousEfficiency() {
	Apply(y_bar);
	return *this;
}

cSpectrum& cSpectrum::ApplyColorGlass(std::string filename,double thickness) {
	int i;
	cXYData x;
	cXYList li;
	double t0;
	std::ifstream from( (cMaterial::ColorGlassDataFolder+'/'+filename).c_str() );
	if(from){
		from>>t0;
		while(from>>x) li.AddTail(x);
		for(i=1; i<=li.GetSize(); ++i) {
			li.GetData(x,i);
			x.y=pow( x.y, thickness/t0 );
			li.SetData(x,i);
		}
		spectrum.Apply(li);
	}
	return *this;
}

cSpectrum& cSpectrum::Whiten() {
	int i;
	cXYData x;
	for(i=1; i<=spectrum.GetSize(); ++i){
		spectrum.GetData(x,i);
		x.y=1;
		spectrum.SetData(x,i);
	}
	return *this;
}

cSpectrum& cSpectrum::Normalize(double new_max) {
	spectrum.Normalize(new_max);
	return *this;
}
cSpectrum& cSpectrum::NormalizeByLocal(double wl_start,double wl_end,double new_max) {
	spectrum.NormalizeByLocal(wl_start,wl_end,new_max);
	return *this;
}

cSpectrum& cSpectrum::Shift(double multiplier){
	// multiplier倍シフトさせる．
	// 光学薄膜で膜厚をmultiplier倍にした場合に対応する．
	int i;
	cXYList li;
	double wl,val;
	if(multiplier==0) return *this;
	li=spectrum;
	spectrum.RemoveAll();
	li.Sort();
	for(i=1; i<=li.GetSize(); i++){
		wl=li.x(i);
		val=li.y(wl/multiplier);
		spectrum.AddData(wl,val);
	}
	return *this;
}

cSpectrum& cSpectrum::ReverseTR(){
	int i;
	cXYData x;
	for(i=1; i<=spectrum.GetSize(); i++){
		spectrum.GetData(x,i);
		x.y=1-x.y;
		spectrum.SetData(x,i);
	}
	return *this;
}

double cSpectrum::x() {
	// CIE 1931 xy色度図
	int i;
	cXYData s;
	cXYList x,y,z;
	double X=0,Y=0,Z=0;
	x=y=z=spectrum;
	x.Apply(x_bar);
	y.Apply(y_bar);
	z.Apply(z_bar);
	for(i=1; i<=spectrum.GetSize(); ++i){
		x.GetData(s,i); X+=s.y;
		y.GetData(s,i); Y+=s.y;
		z.GetData(s,i); Z+=s.y;
	}
	return X+Y+Z==0 ? 0 : X/(X+Y+Z);
}

double cSpectrum::y() {
	// CIE 1931 xy色度図
	int i;
	cXYData s;
	cXYList x,y,z;
	double X=0,Y=0,Z=0;
	x=y=z=spectrum;
	x.Apply(x_bar);
	y.Apply(y_bar);
	z.Apply(z_bar);
	for(i=1; i<=spectrum.GetSize(); ++i){
		x.GetData(s,i); X+=s.y;
		y.GetData(s,i); Y+=s.y;
		z.GetData(s,i); Z+=s.y;
	}
	return X+Y+Z==0 ? 0 : Y/(X+Y+Z);
}


double cSpectrum::u() {
	// CIE 1960 UCS色度図
	double x,y;

	x=this->x();
	y=this->y();
	return 4*x/(-2*x+12*y+3);
}

double cSpectrum::v() {
	// CIE 1960 UCS色度図
	double x,y;

	x=this->x();
	y=this->y();
	return 6*y/(-2*x+12*y+3);
}

double cSpectrum::Tcp(){
	// 相関色温度 (JIS Z 8725:2015 による)
	double dummy;

	return Tcp(dummy);
}

double cSpectrum::duv(){
	// 黒体からの偏差 (JIS Z 8725:2015 による)
	double x;

	Tcp(x);
	return x;
}

double cSpectrum::LuminousFlux() {
	int i;
	cXYData s;
	cXYList y;
	double Y=0;
	y=spectrum;
	y.Apply(y_bar);
	for(i=1; i<=y.GetSize(); ++i){
		y.GetData(s,i); Y+=s.y;
	}
	return Y;
}

double cSpectrum::x(std::string IlluminantName,double temp){
	if(IlluminantName!=""){
		if(IlluminantName=="A")     return Copy().Apply(A).x();
		if(IlluminantName=="C")     return Copy().Apply(C).x();
		if(IlluminantName=="D65")   return Copy().Apply(D65).x();
		if(IlluminantName=="WHITE") return Copy().x();
		return x();
	}
	else if(temp!=0){
		return Copy().ApplyBlackBodyRadiation(temp).x();
	}
	else{
		return x();
	}
}

double cSpectrum::y(std::string IlluminantName,double temp){
	if(IlluminantName!=""){
		if(IlluminantName=="A")     return Copy().Apply(A).y();
		if(IlluminantName=="C")     return Copy().Apply(C).y();
		if(IlluminantName=="D65")   return Copy().Apply(D65).y();
		if(IlluminantName=="WHITE") return Copy().y();
		return y();
	}
	else if(temp!=0){
		return Copy().ApplyBlackBodyRadiation(temp).y();
	}
	else{
		return y();
	}
}

double cSpectrum::u(std::string IlluminantName,double temp){
	if(IlluminantName!=""){
		if(IlluminantName=="A")     return Copy().Apply(A).u();
		if(IlluminantName=="C")     return Copy().Apply(C).u();
		if(IlluminantName=="D65")   return Copy().Apply(D65).u();
		if(IlluminantName=="WHITE") return Copy().u();
		return u();
	}
	else if(temp!=0){
		return Copy().ApplyBlackBodyRadiation(temp).u();
	}
	else{
		return u();
	}
}

double cSpectrum::v(std::string IlluminantName,double temp){
	if(IlluminantName!=""){
		if(IlluminantName=="A")     return Copy().Apply(A).v();
		if(IlluminantName=="C")     return Copy().Apply(C).v();
		if(IlluminantName=="D65")   return Copy().Apply(D65).v();
		if(IlluminantName=="WHITE") return Copy().v();
		return v();
	}
	else if(temp!=0){
		return Copy().ApplyBlackBodyRadiation(temp).v();
	}
	else{
		return v();
	}
}

double cSpectrum::Tcp(std::string IlluminantName,double temp){
	if(IlluminantName!=""){
		if(IlluminantName=="A")     return Copy().Apply(A).Tcp();
		if(IlluminantName=="C")     return Copy().Apply(C).Tcp();
		if(IlluminantName=="D65")   return Copy().Apply(D65).Tcp();
		if(IlluminantName=="WHITE") return Copy().Tcp();
		return v();
	}
	else if(temp!=0){
		return Copy().ApplyBlackBodyRadiation(temp).Tcp();
	}
	else{
		return Tcp();
	}
}

double cSpectrum::duv(std::string IlluminantName,double temp){
	if(IlluminantName!=""){
		if(IlluminantName=="A")     return Copy().Apply(A).duv();
		if(IlluminantName=="C")     return Copy().Apply(C).duv();
		if(IlluminantName=="D65")   return Copy().Apply(D65).duv();
		if(IlluminantName=="WHITE") return Copy().duv();
		return v();
	}
	else if(temp!=0){
		return Copy().ApplyBlackBodyRadiation(temp).duv();
	}
	else{
		return duv();
	}
}

double cSpectrum::LuminousFluxRatio(std::string IlluminantName,double temp){
	if(IlluminantName!=""){
		if(IlluminantName=="A"){
			return Copy().Apply(A).LuminousFlux()/Copy().Whiten().Apply(A).LuminousFlux();
		}
		if(IlluminantName=="C"){
			return Copy().Apply(C).LuminousFlux()/Copy().Whiten().Apply(C).LuminousFlux();
		}
		if(IlluminantName=="D65"){
			return Copy().Apply(D65).LuminousFlux()/Copy().Whiten().Apply(D65).LuminousFlux();
		}
		return LuminousFlux()/Copy().Whiten().LuminousFlux();
	}
	else if(temp!=0){
		return Copy().ApplyBlackBodyRadiation(temp).LuminousFlux()
		      /Copy().Whiten().ApplyBlackBodyRadiation(temp).LuminousFlux();
	}
	else{
		return LuminousFlux()/Copy().Whiten().LuminousFlux();
	}
}

double cSpectrum::Integral(double wl_start,double wl_end) {
	return spectrum.Integral(wl_start,wl_end);
}

double cSpectrum::Average(double wl_start,double wl_end) {
	return spectrum.Average(wl_start,wl_end);
}

double cSpectrum::Peak(){
	return spectrum.Peak();
}

double cSpectrum::PeakLocal(double wl_start,double wl_end){
	return spectrum.PeakLocal(wl_start,wl_end);
}

double cSpectrum::FW(double threshold){
	return spectrum.FW(threshold);
}

double cSpectrum::FWStart(double threshold){
	return spectrum.FWStart(threshold);
}

double cSpectrum::FWEnd(double threshold){
	return spectrum.FWEnd(threshold);
}

double cSpectrum::FWCenter(double threshold){
	return spectrum.FWCenter(threshold);
}

double cSpectrum::FWLocal(double wl_start,double wl_end,double threshold){
	return spectrum.FWLocal(wl_start,wl_end,threshold);
}

double cSpectrum::FWStartLocal(double wl_start,double wl_end,double threshold){
	return spectrum.FWStartLocal(wl_start,wl_end,threshold);
}

double cSpectrum::FWEndLocal(double wl_start,double wl_end,double threshold){
	return spectrum.FWEndLocal(wl_start,wl_end,threshold);
}

double cSpectrum::FWCenterLocal(double wl_start,double wl_end,double threshold){
	return spectrum.FWCenterLocal(wl_start,wl_end,threshold);
}

double cSpectrum::FWHM(){
	return spectrum.FWHM();
}

double cSpectrum::FWHMLocal(double wl_start,double wl_end){
	return spectrum.FWHMLocal(wl_start,wl_end);
}

double cSpectrum::TransitionPoint(double val){
	return spectrum.TransitionPoint(val);
}

double cSpectrum::TransitionPointLocal(double wl_start,double wl_end,double val){
	return spectrum.TransitionPointLocal(wl_start,wl_end,val);
}

double cSpectrum::TransitionInterval(double val1,double val2){
	return spectrum.TransitionInterval(val1,val2);
}

double cSpectrum::TransitionIntervalLocal(double wl_start,double wl_end,double val1,double val2){
	return spectrum.TransitionIntervalLocal(wl_start,wl_end,val1,val2);
}

std::string cSpectrum::Str() {
	int i;
	cXYData x;
	std::string s;
	char buf[1000];
	for(i=1; i<=spectrum.GetSize(); ++i) {
		spectrum.GetData(x,i);
		sprintf(buf,"%g\t%g\n", x.x,x.y); s+=buf;
	}
	return s;
}

std::istream& operator>>(std::istream& from, cSpectrum& x){
	from >> x.spectrum;
	return from;
}
std::ostream& operator<<(std::ostream& to, const cSpectrum& x){
	to << x.spectrum;
	return to;
}

int cSpectrum::Open(const std::string& filename){
	std::ifstream from(filename.c_str());
	if(from){
		from >> *this;
		return 1;
	}
	else{
		return 0;
	}
}
int cSpectrum::Save(const std::string& filename){
	std::ofstream to(filename.c_str());
	if(to){
		to << *this;
		return 1;
	}
	else{
		return 0;
	}
}

////////////////////////////////////////////////////////////////////////////////
void cSpectrum::make_constants() {
	// x_bar,y_bar,z_barは http://cvrl.ioo.ucl.ac.uk/ より
	// (CIE 1931 2-deg, XYZ Color Matching Functionsとうたっている)
	// JOEMテキストより桁数が多く，スペクトル軌跡長波長側の(x,y)を計算できる．

	x_bar.AddData( 380, 0.00136800000 );
	x_bar.AddData( 390, 0.00424300000 );
	x_bar.AddData( 400, 0.01431000000 );
	x_bar.AddData( 410, 0.04351000000 );
	x_bar.AddData( 420, 0.13438000000 );
	x_bar.AddData( 430, 0.28390000000 );
	x_bar.AddData( 440, 0.34828000000 );
	x_bar.AddData( 450, 0.33620000000 );
	x_bar.AddData( 460, 0.29080000000 );
	x_bar.AddData( 470, 0.19536000000 );
	x_bar.AddData( 480, 0.09564000000 );
	x_bar.AddData( 490, 0.03201000000 );
	x_bar.AddData( 500, 0.00490000000 );
	x_bar.AddData( 510, 0.00930000000 );
	x_bar.AddData( 520, 0.06327000000 );
	x_bar.AddData( 530, 0.16550000000 );
	x_bar.AddData( 540, 0.29040000000 );
	x_bar.AddData( 550, 0.43344990000 );
	x_bar.AddData( 560, 0.59450000000 );
	x_bar.AddData( 570, 0.76210000000 );
	x_bar.AddData( 580, 0.91630000000 );
	x_bar.AddData( 590, 1.02630000000 );
	x_bar.AddData( 600, 1.06220000000 );
	x_bar.AddData( 610, 1.00260000000 );
	x_bar.AddData( 620, 0.85444990000 );
	x_bar.AddData( 630, 0.64240000000 );
	x_bar.AddData( 640, 0.44790000000 );
	x_bar.AddData( 650, 0.28350000000 );
	x_bar.AddData( 660, 0.16490000000 );
	x_bar.AddData( 670, 0.08740000000 );
	x_bar.AddData( 680, 0.04677000000 );
	x_bar.AddData( 690, 0.02270000000 );
	x_bar.AddData( 700, 0.01135916000 );
	x_bar.AddData( 710, 0.00579034600 );
	x_bar.AddData( 720, 0.00289932700 );
	x_bar.AddData( 730, 0.00143997100 );
	x_bar.AddData( 740, 0.00069007860 );
	x_bar.AddData( 750, 0.00033230110 );
	x_bar.AddData( 760, 0.00016615050 );
	x_bar.AddData( 770, 0.00008307527 );
	x_bar.AddData( 780, 0.00004150994 );

	y_bar.AddData( 380, 0.00003900000 );
	y_bar.AddData( 390, 0.00012000000 );
	y_bar.AddData( 400, 0.00039600000 );
	y_bar.AddData( 410, 0.00121000000 );
	y_bar.AddData( 420, 0.00400000000 );
	y_bar.AddData( 430, 0.01160000000 );
	y_bar.AddData( 440, 0.02300000000 );
	y_bar.AddData( 450, 0.03800000000 );
	y_bar.AddData( 460, 0.06000000000 );
	y_bar.AddData( 470, 0.09098000000 );
	y_bar.AddData( 480, 0.13902000000 );
	y_bar.AddData( 490, 0.20802000000 );
	y_bar.AddData( 500, 0.32300000000 );
	y_bar.AddData( 510, 0.50300000000 );
	y_bar.AddData( 520, 0.71000000000 );
	y_bar.AddData( 530, 0.86200000000 );
	y_bar.AddData( 540, 0.95400000000 );
	y_bar.AddData( 550, 0.99495010000 );
	y_bar.AddData( 560, 0.99500000000 );
	y_bar.AddData( 570, 0.95200000000 );
	y_bar.AddData( 580, 0.87000000000 );
	y_bar.AddData( 590, 0.75700000000 );
	y_bar.AddData( 600, 0.63100000000 );
	y_bar.AddData( 610, 0.50300000000 );
	y_bar.AddData( 620, 0.38100000000 );
	y_bar.AddData( 630, 0.26500000000 );
	y_bar.AddData( 640, 0.17500000000 );
	y_bar.AddData( 650, 0.10700000000 );
	y_bar.AddData( 660, 0.06100000000 );
	y_bar.AddData( 670, 0.03200000000 );
	y_bar.AddData( 680, 0.01700000000 );
	y_bar.AddData( 690, 0.00821000000 );
	y_bar.AddData( 700, 0.00410200000 );
	y_bar.AddData( 710, 0.00209100000 );
	y_bar.AddData( 720, 0.00104700000 );
	y_bar.AddData( 730, 0.00052000000 );
	y_bar.AddData( 740, 0.00024920000 );
	y_bar.AddData( 750, 0.00012000000 );
	y_bar.AddData( 760, 0.00006000000 );
	y_bar.AddData( 770, 0.00003000000 );
	y_bar.AddData( 780, 0.00001499000 );

	z_bar.AddData( 380, 0.00645000100 );
	z_bar.AddData( 390, 0.02005001000 );
	z_bar.AddData( 400, 0.06785001000 );
	z_bar.AddData( 410, 0.20740000000 );
	z_bar.AddData( 420, 0.64560000000 );
	z_bar.AddData( 430, 1.38560000000 );
	z_bar.AddData( 440, 1.74706000000 );
	z_bar.AddData( 450, 1.77211000000 );
	z_bar.AddData( 460, 1.66920000000 );
	z_bar.AddData( 470, 1.28764000000 );
	z_bar.AddData( 480, 0.81295010000 );
	z_bar.AddData( 490, 0.46518000000 );
	z_bar.AddData( 500, 0.27200000000 );
	z_bar.AddData( 510, 0.15820000000 );
	z_bar.AddData( 520, 0.07824999000 );
	z_bar.AddData( 530, 0.04216000000 );
	z_bar.AddData( 540, 0.02030000000 );
	z_bar.AddData( 550, 0.00874999900 );
	z_bar.AddData( 560, 0.00390000000 );
	z_bar.AddData( 570, 0.00210000000 );
	z_bar.AddData( 580, 0.00165000100 );
	z_bar.AddData( 590, 0.00110000000 );
	z_bar.AddData( 600, 0.00080000000 );
	z_bar.AddData( 610, 0.00034000000 );
	z_bar.AddData( 620, 0.00019000000 );
	z_bar.AddData( 630, 0.00004999999 );
	z_bar.AddData( 640, 0.00002000000 );
	z_bar.AddData( 650, 0.00000000000 );
	z_bar.AddData( 660, 0.00000000000 );
	z_bar.AddData( 670, 0.00000000000 );
	z_bar.AddData( 680, 0.00000000000 );
	z_bar.AddData( 690, 0.00000000000 );
	z_bar.AddData( 700, 0.00000000000 );
	z_bar.AddData( 710, 0.00000000000 );
	z_bar.AddData( 720, 0.00000000000 );
	z_bar.AddData( 730, 0.00000000000 );
	z_bar.AddData( 740, 0.00000000000 );
	z_bar.AddData( 750, 0.00000000000 );
	z_bar.AddData( 760, 0.00000000000 );
	z_bar.AddData( 770, 0.00000000000 );
	z_bar.AddData( 780, 0.00000000000 );

	A.AddData( 380,   9.80 );
	A.AddData( 390,  12.09 );
	A.AddData( 400,  14.71 );
	A.AddData( 410,  17.68 );
	A.AddData( 420,  20.99 );
	A.AddData( 430,  24.67 );
	A.AddData( 440,  28.70 );
	A.AddData( 450,  33.09 );
	A.AddData( 460,  37.81 );
	A.AddData( 470,  42.87 );
	A.AddData( 480,  48.24 );
	A.AddData( 490,  53.91 );
	A.AddData( 500,  59.86 );
	A.AddData( 510,  66.06 );
	A.AddData( 520,  72.50 );
	A.AddData( 530,  79.13 );
	A.AddData( 540,  85.95 );
	A.AddData( 550,  92.91 );
	A.AddData( 560, 100.00 );
	A.AddData( 570, 107.18 );
	A.AddData( 580, 114.44 );
	A.AddData( 590, 121.73 );
	A.AddData( 600, 129.04 );
	A.AddData( 610, 136.35 );
	A.AddData( 620, 143.62 );
	A.AddData( 630, 150.84 );
	A.AddData( 640, 157.98 );
	A.AddData( 650, 165.03 );
	A.AddData( 660, 171.96 );
	A.AddData( 670, 178.77 );
	A.AddData( 680, 185.43 );
	A.AddData( 690, 191.93 );
	A.AddData( 700, 198.26 );
	A.AddData( 710, 204.41 );
	A.AddData( 720, 210.36 );
	A.AddData( 730, 216.12 );
	A.AddData( 740, 221.67 );
	A.AddData( 750, 227.00 );
	A.AddData( 760, 232.12 );
	A.AddData( 770, 237.01 );
	A.AddData( 780, 241.68 );

	D65.AddData( 380,  49.98 );
	D65.AddData( 390,  54.65 );
	D65.AddData( 400,  82.75 );
	D65.AddData( 410,  91.49 );
	D65.AddData( 420,  93.43 );
	D65.AddData( 430,  86.68 );
	D65.AddData( 440, 104.86 );
	D65.AddData( 450, 117.01 );
	D65.AddData( 460, 117.81 );
	D65.AddData( 470, 114.86 );
	D65.AddData( 480, 115.92 );
	D65.AddData( 490, 108.81 );
	D65.AddData( 500, 109.35 );
	D65.AddData( 510, 107.80 );
	D65.AddData( 520, 104.79 );
	D65.AddData( 530, 107.69 );
	D65.AddData( 540, 104.41 );
	D65.AddData( 550, 104.05 );
	D65.AddData( 560, 100.00 );
	D65.AddData( 570,  96.33 );
	D65.AddData( 580,  95.79 );
	D65.AddData( 590,  88.69 );
	D65.AddData( 600,  90.01 );
	D65.AddData( 610,  89.60 );
	D65.AddData( 620,  87.70 );
	D65.AddData( 630,  83.29 );
	D65.AddData( 640,  83.70 );
	D65.AddData( 650,  80.03 );
	D65.AddData( 660,  80.21 );
	D65.AddData( 670,  82.28 );
	D65.AddData( 680,  78.28 );
	D65.AddData( 690,  69.72 );
	D65.AddData( 700,  71.61 );
	D65.AddData( 710,  74.35 );
	D65.AddData( 720,  61.60 );
	D65.AddData( 730,  69.89 );
	D65.AddData( 740,  75.09 );
	D65.AddData( 750,  63.59 );
	D65.AddData( 760,  46.42 );
	D65.AddData( 770,  66.81 );
	D65.AddData( 780,  63.38 );

	C.AddData( 380,  33.00 );
	C.AddData( 390,  47.40 );
	C.AddData( 400,  63.30 );
	C.AddData( 410,  80.60 );
	C.AddData( 420,  98.10 );
	C.AddData( 430, 112.40 );
	C.AddData( 440, 121.50 );
	C.AddData( 450, 124.00 );
	C.AddData( 460, 123.10 );
	C.AddData( 470, 123.80 );
	C.AddData( 480, 123.90 );
	C.AddData( 490, 120.70 );
	C.AddData( 500, 112.10 );
	C.AddData( 510, 102.30 );
	C.AddData( 520,  96.90 );
	C.AddData( 530,  98.00 );
	C.AddData( 540, 102.10 );
	C.AddData( 550, 105.20 );
	C.AddData( 560, 105.30 );
	C.AddData( 570, 102.30 );
	C.AddData( 580,  97.80 );
	C.AddData( 590,  93.20 );
	C.AddData( 600,  89.70 );
	C.AddData( 610,  88.40 );
	C.AddData( 620,  88.10 );
	C.AddData( 630,  88.00 );
	C.AddData( 640,  87.80 );
	C.AddData( 650,  88.20 );
	C.AddData( 660,  87.90 );
	C.AddData( 670,  86.30 );
	C.AddData( 680,  84.00 );
	C.AddData( 690,  80.20 );
	C.AddData( 700,  76.30 );
	C.AddData( 710,  72.40 );
	C.AddData( 720,  68.30 );
	C.AddData( 730,  64.40 );
	C.AddData( 740,  61.50 );
	C.AddData( 750,  59.20 );
	C.AddData( 760,  58.10 );
	C.AddData( 770,  58.20 );
	C.AddData( 780,  59.10 );

	for(int i=0; i<=40; i++) WHITE.AddData(380+10*i, 1);
}