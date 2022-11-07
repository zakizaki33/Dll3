#include "stdafx.h"
#include "general_func.h"

int CreateConsole() {
	// プロセスにコンソールを割り当て，標準出力を向ける
	FreeConsole();
	if(AllocConsole()){   // 【64bit環境では実行時になぜか必ず失敗する．32bitでは問題ない】
		freopen("CONOUT$","w",stdout); /* 標準出力(stdout)を新しいコンソールに向ける */
		freopen("CONOUT$","w",stderr); /* 標準エラー出力(stderr)を新しいコンソールに向ける */
		return 1; // 成功
	}
	else{
		return 0; // 失敗
	}
}

std::string RemoveExtension(std::string filename){
	// filenameより拡張子を除いた文字列を返す．
	return remove_extension(filename);
}

std::string Date(){
	std::string s;
	char buf[100];
	time_t now;
	struct tm *pnow;

	now=time(NULL);
	pnow=localtime(&now);
	sprintf(buf, "%4d/%02d/%02d", pnow->tm_year+1900, pnow->tm_mon+1, pnow->tm_mday);  // 1月が0
	return s=buf;
}

std::string Time(){
	std::string s;
	char buf[100];
	time_t now;
	struct tm *pnow;

	now=time(NULL);
	pnow=localtime(&now);
	sprintf(buf, "%02d:%02d:%02d", pnow->tm_hour, pnow->tm_min, pnow->tm_sec);
	return s=buf;
}

std::string LapTime(clock_t start){
	// startからの経過時間を秒単位で返す
	char buf[100];

	sprintf(buf,"%.3fsec",(double)(clock()-start)/CLOCKS_PER_SEC);
	return buf;
}

void dwrite(std::ostream to,double x){ // バイナリモード double型書き出し
	to.write((char*)&x, sizeof(double));
}

double dread(std::istream from){ // バイナリモード double型読み込み
	double x;

	from.read((char*)&x, sizeof(double));
	return x;
}

void iwrite(std::ostream to,int x){ // バイナリモード int型書き出し
	to.write((char*)&x, sizeof(int));
}

int iread(std::istream from){ // バイナリモード int型読み込み
	int x;

	from.read((char*)&x, sizeof(int));
	return x;
}

void swrite(std::ostream to,std::string s){ // バイナリモード std::string型書き出し
	char *buf;
	int n;

	n=(int)strlen(s.c_str());
	buf=new char[n+1];
	strcpy(buf,s.c_str());
	to.write((char*)&n, sizeof(int));
	to.write((char*)buf, n+1);
	delete [] buf;
}

std::string sread(std::istream from){ // バイナリモード std::string型読み込み
	std::string s;
	char *buf;
	int n;

	from.read((char*)&n, sizeof(int)); 
	buf=new char[n+1];
	from.read((char*)buf,n+1);
	s=buf;
	delete [] buf;
	return s;
}

bool is_even(const int& i){
	return (i>0 ? i : -i)%2==0 ? true : false;
}

bool is_odd(const int& i){
	return !(is_even(i));
}

double factorial(const int& i){
	if(i<0){
		return 0;
	}
	else if(i==0){
		return 1;
	}
	else{
		double a;  // "int a" とすると例えば 15! はオーバーフローする
		int j;
		a=1;
		for(j=1; j<=i; j++) a*=(double)j;
		return a;
	}
}

double pw(const double& x,const int& n){
	// xのn乗，n<0では正しくない
	int i;
	double p=1;

	for(i=1; i<=n; i++) p*=x;
	return p;
}

double fceil(const double& x){
	// x以上で最小の
	//    (正負の一桁の整数) * (10の階乗)
	// である数値を返す．
	double xabs= x>=0 ? x : -x;
	return x==0 ? 0 : ceil( x/pow(10,floor(log10(xabs)+1e-15))-1e-15 ) * pow(10,floor(log10(xabs)+1e-15));
}

double ffloor(const double& x){
	// x以下で最大の
	//     (正負の一桁の整数) * (10の階乗)
	// である数値を返す．
	double xabs= x>=0 ? x : -x;
	return x==0 ? 0 : floor( x/pow(10,floor(log10(xabs)+1e-15))+1e-15 ) * pow(10,floor(log10(xabs)+1e-15));
}

double fceil_mantissa(const double& x){
	// x以上で最小の
	//     a * (10の階乗), a=1,2,3,..,9,10
	// である数値を求め，aを返す．
	// グラフの軸の分割数を求める場合等に使用する．
	double xabs= x>=0 ? x : -x;
	return x==0 ? 0 : ceil( x/pow(10,floor(log10(xabs)+1e-15))-1e-15 );
}

double sgn(const double& x) { return x>=0 ? 1 : -1; };

double Round(const double& x,const int& n){
	// n<=0 : 小数点以下|n|桁に丸める
	// n>=1 : 有効桁をn桁に丸める
	double absx;

	if(n<=0){
		double a=1; int i;
		double sign= x>=0 ? 1 : -1;
		absx=fabs(x);
		for(i=1; i<=-n; ++i) a*=10;
		return sign*floor(absx*a+0.5 )/a;
	}
	else{
		double a=1;
		double sign= x>=0 ? 1 : -1;
		absx=fabs(x);
		if(absx>=10){
			while(absx>=10){
				absx/=10; a/=10;
			}
		}
		if(0<absx && absx<1){
			while(absx<1){
				absx*=10; a*=10;
			}
		}
		return sign*Round(absx,-n+1)/a;
	}
}

int ToInt(const double& x){
	// 小数点以下1位を四捨五入した整数を返す
	return x>=0 ? (int)(x+0.5) : (int)(x-0.5);
	// （doubleからintへの変換は，小数点以下を切捨てて0の方向へ向かう
	//   ことがJIS規格で規定されている．）
}

double ArcCos(const double& x){
	// cos(x)よりacos関数を用いて入射角xを求めようとするとき，
	// 数値計算誤差によりxがわずかでも1を超えていると定義域エラーになってしまう．
	// 本関数ではxが1を超えているときは，xの有効桁を減らしてからacosを実行することで
	// エラーの発生を抑える．
	if(x<=1){
		return acos(x);
	}
	else{
		return acos(Round(x,8));  // xの有効桁を8桁に丸める
	}
}

void InverseMatrix(double a[4][4]){
	// 3x3行列の逆行列を計算する
	double det, b[4][4];

	det= a[1][1]*a[2][2]*a[3][3] +a[2][1]*a[3][2]*a[1][3] +a[3][1]*a[1][2]*a[2][3]
		-a[1][1]*a[3][2]*a[2][3] -a[3][1]*a[2][2]*a[1][3] -a[2][1]*a[1][2]*a[3][3];

	if(det==0) return;

	b[1][1]=a[1][1]; b[1][2]=a[1][2]; b[1][3]=a[1][3];
	b[2][1]=a[2][1]; b[2][2]=a[2][2]; b[2][3]=a[2][3];
	b[3][1]=a[3][1]; b[3][2]=a[3][2]; b[3][3]=a[3][3];

	a[1][1]=b[2][2]*b[3][3]-b[2][3]*b[3][2]; a[1][2]=b[1][3]*b[3][2]-b[1][2]*b[3][3]; a[1][3]=b[1][2]*b[2][3]-b[1][3]*b[2][2];
	a[2][1]=b[2][3]*b[3][1]-b[2][1]*b[3][3]; a[2][2]=b[1][1]*b[3][3]-b[1][3]*b[3][1]; a[2][3]=b[1][3]*b[2][1]-b[1][1]*b[2][3];
	a[3][1]=b[2][1]*b[3][2]-b[2][2]*b[3][1]; a[3][2]=b[1][2]*b[3][1]-b[1][1]*b[3][2]; a[3][3]=b[1][1]*b[2][2]-b[1][2]*b[2][1];	

	a[1][1]/=det; a[1][2]/=det; a[1][3]/=det;
	a[2][1]/=det; a[2][2]/=det; a[2][3]/=det;
	a[3][1]/=det; a[3][2]/=det; a[3][3]/=det;
}

void Srand(){
	::srand( (unsigned)time(NULL) );   // 1秒ごとに乱数の並びが変わる
}

double Random(const double& semi_width,const int& IsEndNotUni){
	if(IsEndNotUni) return (double)rand()/(double)RAND_MAX>0.5 ? semi_width : -semi_width;
	else            return -semi_width + 2*semi_width*(double)rand()/(double)RAND_MAX;
}

double Random(const double& limit1,const double& limit2,const int& IsEndNotUni){
	if(IsEndNotUni) return (double)rand()/(double)RAND_MAX>0.5 ? limit1 : limit2;
	else            return limit1 + (limit2-limit1)*(double)rand()/(double)RAND_MAX;
}

double RandomGauss(const double& sigma){
	// |x|<Kσ の範囲の乱数ｘを発生する．
	// xは正規分布 y=exp{-x^2/(2σ^2)} の確率で発生する．
	// この確率は，半径σのガウスビームの強度分布に等しい．
	// |x|<σ である確率は 68.26%, <2σの確率は95.44%(レーザのビーム直径), <3σでは99.74%
	const double K=4;  
	double x,y;

	if( sigma==0 ) return 0;
	do{
		x= -K*sigma +2*K*sigma*(double)rand()/(double)RAND_MAX;
		y= exp( -2*x*x/sigma/sigma );
	} while ( (double)rand()/(double)RAND_MAX > y );
	return x;

	// ガウスビームの式と混同しないこと．
	// ガウスビームのパワー分布は y = exp(-2w/wo) となる（分子に2が付く）．
	// 直径2woの中に入るパワーは全体の86.5%といわれており，
	// これは上記“<2σの確率95.44%”の二乗，91%よりも低いが
	// 91%は一辺が2woの "正方形" に含まれるパワーであるためである．
}

long rgb(const long& r,const long& g,const long& b){
	long rr,gg,bb;

	rr=r; if(rr>255) rr=255; if(rr<0) rr=0;
	gg=g; if(gg>255) gg=255; if(gg<0) gg=0;
	bb=b; if(bb>255) bb=255; if(bb<0) bb=0;
	return bb*256*256+gg*256+rr;
}

long Rrgb(const long& rgb){
	return rgb%256;
}

long Grgb(const long& rgb){
	return (rgb/256)%256;
}

long Brgb(const long& rgb){
	return rgb/256/256;
}

long RGBComplement(const long& rgb){
	// 補色を返す．
	// 参考： http://appakumaturi.hatenablog.com/entry/20120121/1327143125
	long x, r,g,b;

	r=Rrgb(rgb); g=Grgb(rgb); b=Brgb(rgb);
	x=(long)Max(r,g,b)+(long)Min(r,g,b);
	r=x-r; g=x-g; b=x-b;
	return ::rgb(r,g,b);
}


double Max(const double& x1,const double& x2) { return x1>x2 ? x1:x2; }

int Max(const int& i1,const int& i2) { return i1>i2 ? i1:i2; }

double Max(const double& x1,const double& x2,const double& x3){
	double max;
	max=x1;
	if(x2>max) max=x2;
	if(x3>max) max=x3;
	return max;
}

double Max(const double& x1,const double& x2,const double& x3,const double& x4){
	double max;
	max= x1>x2  ? x1:x2;
	max= x3>max ? x3:max;
	max= x4>max ? x4:max;
	return max;
}
	
double Min(const double& x1,const double& x2) { return x1<x2 ? x1:x2; }

int Min(const int& i1,const int& i2) { return i1<i2 ? i1:i2; }

double Min(const double& x1,const double& x2,const double& x3){
	double min;
	min=x1;
	if(x2<min) min=x2;
	if(x3<min) min=x3;
	return min;
}

double Min(const double& x1,const double& x2,const double& x3,const double& x4){
	double min;
	min= x1<x2  ? x1:x2;
	min= x3<min ? x3:min;
	min= x4<min ? x4:min;
	return min;
}

double Median(double *buf,int n){
	// buf[0]からbuf[n-1]の中央値を返す．
	Sort(buf,n);
	return buf[n/2];
}

int QuadraticEq(double& x1,double& x2,double a,double b,double c){
	// 2次方程式 ax^2+bx+c=0 を解く
	if( b*b-4*a*c >= 0 ){
		// -bとsprt(...)との加減算での桁落ちをふせぐため
		//    x1=(-b+sqrt(b*b-4*a*c))/2/a;
		//    x2=(-b-sqrt(b*b-4*a*c))/2/a;
		// とはしない．
		// 参考：“FORTRAN 基本＋応用(刀根薫)”p108
		if(b>=0){
			x1=2*c/(-b-sqrt(b*b-4*a*c));
			x2=(-b-sqrt(b*b-4*a*c))/2/a;
		}
		else{
			x1=(-b+sqrt(b*b-4*a*c))/2/a;
			x2=2*c/(-b+sqrt(b*b-4*a*c));
		}
		return 1;
	}
	else{
		x1=0;
		x2=0;
		return 0;
	}
}
double QuadraticX1(double a,double b,double c){
	double x1,x2;
	if(QuadraticEq(x1,x2,a,b,c)) return x1; else return 0;
}
double QuadraticX2(double a,double b,double c){
	double x1,x2;
	if(QuadraticEq(x1,x2,a,b,c)) return x2; else return 0;
}

void DFT(complex *a,int n,int inv,int optical,double zoom) {
	// *aを含みaからｎ個のデータ(a[0]～a[n-1])に対し離散フーリエ変換を行い，結果を上書きする．
	//     nは任意．nが2の階乗のときはFFTで計算する．
	//     inv==0のとき順変換，inv!=0のとき逆変換．
	//     optical==0のときa[0]が直流成分，optical!=0のとき中央が直流成分．
	//     zoomが0,1以外のとき，周波数軸最大値が1/zoom倍になる．このときFFTは使えない．
	const double PI=3.14159265358979;
	int i,j,k,w,j1,j2;
	int length,ex;
	int numb,lenb,timb;
	complex *w_table;
	complex x,y;
	double nrml;
	bool fft=true;

	length=1; ex=0;
	do{
		ex+=1;
		length*=2;
		if(length>n) { fft=false; break; }
	} while(n!=length);

	if(zoom!=1 && zoom!=0) fft=false;

	if(fft)
	{
		if(optical) for(i=1; i<n; i+=2) a[i]=-a[i];

		w_table=new complex[n+1];
		{
			int i;
			double xx,arg;

			xx=-PI*2.0/(double)n;
			if(inv) xx=-xx;
			for(i=0; i<n; ++i){
				arg=(double)i*xx;
				w_table[i]=complex(cos(arg),sin(arg));
			}
		}

		numb=1;
		lenb=n;
		for(i=0; i<ex; ++i){
			lenb/=2;
			timb=0;
			for(j=0; j<numb; ++j){
				w=0;
				for(k=0; k<lenb; ++k){
					j1=timb+k;
					j2=j1+lenb;
					x=a[j1];
					y=a[j2];
					a[j1]=x+y;
					x=x-y;
					a[j2]=x*w_table[w];
					w+=numb;
				}
				timb+=(2*lenb);
			}
			numb*=2;
		}
		{
			int i,ii,k,bit;
			complex *b=new complex[n+1];

			for(i=0; i<n; ++i){
				for(k=0,ii=i,bit=0; ; bit<<=1,ii>>=1){
					bit= (ii&1) | bit;
					if(++k==ex) break;
				}
				b[i]=a[bit];
			}
			for(i=0; i<n; ++i) a[i]=b[i];
			delete[] b;
		}
		delete[] w_table;

		if(optical) for(i=1; i<n; i+=2) a[i]=-a[i];
	}
	else
	{
/*		{
			// 離散フーリエ変換定義式をそのままプログラムする
			int i,j, i1,j1;
			double w,w0;
			complex *b=new complex[n+1];

			if(inv) w0=-2.0*PI/n; else w0=2.0*PI/n;
			if(zoom!=0) w0/=zoom;

			for(i=0; i<n; i++){
				b[i]=0;
				for(j=0; j<n; j++){
					if(optical){
						i1=i-n/2;
						j1=j-n/2;
					}
					else{
						i1=i;
						j1=j;
					}
					w=w0*i1*j1;
					b[i]=b[i]+a[j]*complex( cos(w),sin(w) );  // b[i]+=a[i]*exp{I*(±2*PI*i1*j1/n)}
				}
			}

			for(i=0; i<n; i++) a[i]=b[i];
			delete[] b;
		}
*/
		{
			// 離散フーリエ変換
			// 定義式どおりにプログラムするのではなく，可読性は落ちるが，
			// ・forループから出せるものは出す．
			// ・三角関数の加法定理を利用して, sin(),cos()の呼び出しを減らす．
			// ことで高速化を図る．
			
			// 順変換         a[j] = Σ(j=0～n-1) a[i]exp( -(2πI/n)ij )
			// 逆変換         a[j] = Σ(j=0～n-1) a[i]exp(  (2πI/n)ij )
			// 順変換(光学的) a[j] = Σ(j=0～n-1) a[i]exp( -(2πI/n)(i-n/2)(j-n/2) )
			// 逆変換(光学的) a[j] = Σ(j=0～n-1) a[i]exp(  (2πI/n)(i-n/2)(j-n/2) )
			int i,j, p;
			double w;
			double dcs,dsn,cs,sn, dcs0,dsn0,cs0,sn0, cs_temp,sn_temp;
			complex *b=new complex[n+1];
			
			w= inv ? 2.0*PI/n : -2.0*PI/n;
			w/=(zoom==0 ? 1 : zoom);
			
			p= optical ? n/2 : 0;

			// w*(-1-p)*(-1-p) は，i=j=-1のときの位相
			cs0=cos(w*(-1-p)*(-1-p));
			sn0=sin(w*(-1-p)*(-1-p));
			// -w*(-1-p)は，j=-1のときの位相(数列)のiに対する公差
			dcs0=cos(w*(-1-p));
			dsn0=sin(w*(-1-p));
			
			for(i=0; i<n; ++i){
				b[i]=0;
				cs_temp=cs0;
				sn_temp=sn0;
				cs0=cs_temp*dcs0-sn_temp*dsn0;  // cos(a+b)=cos(a)cos(b)-sin(a)sin(b)
				sn0=sn_temp*dcs0+cs_temp*dsn0;  // sin(a+b)=sin(a)cos(b)+cos(a)sin(b)
				cs=cs0;
				sn=sn0;
				// w*(i-p) は，位相(数列)のjに対する公差
				dcs=cos(w*(i-p));
				dsn=sin(w*(i-p));
				for(j=0; j<n; ++j){
					cs_temp=cs;
					sn_temp=sn;
					cs=cs_temp*dcs-sn_temp*dsn;
					sn=sn_temp*dcs+cs_temp*dsn;
					b[i]+=a[j]*complex(cs,sn);
				}
			}
			for(i=0; i<n; ++i) a[i]=b[i];
			delete[] b;
		}

	}
	
	nrml=1/sqrt((double)n);
	for(i=0; i<n; ++i) a[i]=a[i]*nrml;
}

void DFT(double *real,double *image,int n,int inv,int optical,double zoom){
	complex* a;
	int i;
	
	a=new complex[n];
	for(i=0; i<n; i++){
		a[i]=complex(real[i],image[i]);
	}
	DFT(a,n,inv,optical,zoom);
	for(i=0; i<n; i++){
		real[i] =Re(a[i]);
		image[i]=Im(a[i]);
	}
	delete [] a;
}

void DFTRow(complex **a,int m,int n,int inv,int optical,double zoom){
	// a[0][0]からa[m-1][n-1]までのm x n個のデータの各行をフーリエ変換する．
	int i;

	for(i=0; i<=m-1; ++i){
		DFT(a[i],n,inv,optical,zoom);
	}
}

void DFTColumn(complex **a,int m,int n,int inv,int optical,double zoom){
	// a[0][0]からa[m-1][n-1]までのm x n個のデータの各列をフーリエ変換する．
	int i,j;
	complex *buf=new complex[m];

	for(j=0; j<=n-1; ++j){
		for(i=0; i<=m-1; ++i) buf[i]=a[i][j];
		DFT(buf,m,inv,optical,zoom);
		for(i=0; i<=m-1; ++i) a[i][j]=buf[i];
	}
	delete [] buf;
}

void DFT(complex **a,int m,int n,int inv,int optical,double zoom){
	// a[0][0]からa[m-1][n-1]までのm x n個のデータをフーリエ変換する．
	DFTRow(a,m,n,inv,optical,zoom);
	DFTColumn(a,m,n,inv,optical,zoom);
}

void HilbertT(complex *a,int n){
	// *aを含みaからｎ個のデータ(a[0]～a[n-1])に対しヒルベルト変換を行い，結果を上書きする．
	int i;
	complex I(0,1);

	DFT(a,n,0,0,0);
	
	for(i=0; i<=n-1; ++i){
		if(i==0)        a[i]=0;
		else if(i<=n/2) a[i]*=I;
		else            a[i]*=-I;  // if(i>n/2)
	}
	
	DFT(a,n,1,0,0);
}

void HilbertT(double *real,double *image,int n){
	complex* a;
	int i;
	
	a=new complex[n];
	for(i=0; i<n; i++){
		a[i]=complex(real[i],image[i]);
	}
	HilbertT(a,n);
	for(i=0; i<n; i++){
		real[i] =Re(a[i]);
		image[i]=Im(a[i]);
	}
	delete [] a;
}

void DispersionCor(double *a,int n,double A2,double A3){
	// *aを含みaからｎ個のデータ(a[0]～a[n-1])に対し分散補正を行い，結果を上書きする．
	// 加える位相量ΔΦは，
	//     ΔΦ = 2π{(A2)x^2 +(A3)x~3} ,  x=i/n (i=0～n-1)
	// とする．
	const double PI=3.141592654;
	const complex I=complex(0,1);

	int i;
	complex *tmp,*c;
	double x,phi;

	if(A2==0 && A3==0) return;

	tmp=new complex[n];
	c=new complex[n];

	for(i=0; i<n; ++i) tmp[i]=c[i]=a[i];

	DFT(c,n,0,0,0);   // フーリエ順変換
	c[0]=0;           // 直流成分を除去 (cが実関数なので，フーリエ変換の定義より，c[0]は実数）
	DFT(c,n,1,0,0);   // フーリエ逆変換 (cから実数成分を除去したことになるから，逆変換後のcは実関数)

	HilbertT(tmp,n);  // ヒルベルト変換（直流成分除去も含まれる）

	for(i=0; i<n; ++i){
		Im(c[i])=Re(tmp[i]);         // 虚数部にヒルベルト変換結果を代入
		x=((double)i-(double)(n-1)/2)/((double)(n-1)/2);  // i=0,n-1でx=±1, i=n/2でx≒0
		phi=2*PI*(A2*x*x+A3*x*x*x);
		c[i]*=exp(I*phi);            // 位相を加える
		a[i]=Re(c[i]);
	}

	delete [] tmp;
	delete [] c;
}

void EnFace(int m,int n,std::string in_filename,std::string out_filename,double gamma/*=1*/){
	// m = en face画像の縦画素数 = Bスキャン画像数
	// n = en face画像の横画素数 = Bスキャン画像の横画素数
	// in_filename  = 入力ファイル名(拡張子を除いたbmpファイル名) = 一連のBスキャン画像のファイル名 
	//               (ex. "oct_a_h#" = "oct_a_h001" ～ "oct_a_h256", "#"が数字で置き換えられる)
	// out_filename = 出力ファイル名(拡張子を除いたbmpファイル名) = en face画像のファイル名
	// gamma = 出力ファイル保存時のガンマ補正値
	cBitmap x;
	int i,j, m_depth,ii;
	matrix<double> A(m,n);
	std::string num;
	char buf[10];
	double max;
	
	for(i=1; i<=m; ++i){
		if     (m<100)   sprintf(buf,"%02d",i);
		else if(m<1000)  sprintf(buf,"%03d",i);
		else if(m<10000) sprintf(buf,"%04d",i);
		else             return;
		num=buf;                                             // 例：m=256のとき，num="001","002",..,"256"
		if( x.Open( replace(in_filename,"#",num)+".bmp" ) ){ // Bスキャン画像を開く
			m_depth=x.GetM();                                // Bスキャン深さ方向画素数
			for(j=1; j<=n; ++j){
				for(ii=1; ii<=m_depth; ++ii) A[i][j]+=x.GetG(ii,j);
			}
		}
	}

	max=0;
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		if(A[i][j]>max) max=A[i][j];
	}
	x.Resize(m,n);
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		A[i][j]*=(255/max);  // 最大値が255となるように規格化
		x.SetRGB(i,j,(int)A[i][j],(int)A[i][j],(int)A[i][j]);	
	}
	x.Save(out_filename+".bmp",gamma);  // en face画像を保存する
}

double SphSag(double r,double h){
	// 球面のサグ量を計算する
	double z;

	if(r==0){
		z=0;
	}
	else{
		z=h*h/r/(1+sqrt( 1-h*h/r/r + 1e-30 ));
	}
	return z;
}

double OverlapAreaCircles(double D,double a){
	// 直径Dの2円の中心間距離がaのときの円の重なり部分の円の面積に対する比を計算する．
	// 無収差円形開口のMTFの計算などに使う．
	// 例えば NA/0.61λ=1.6393NA/λ に対応するMTFは，cut off=2NA/λ だから，
	//   OverlapAreaCircles(2,1.6393) = 0.0894
	// となる．
	double r;
	r=D/2;
	a= a>=0 ? a : -a;
	if(a>=r*2){
		return 0;
	}
	else{
		return ( 2*r*r*acos(a/2/r)-a*sqrt(r*r-a*a/4) )/(PI*r*r);
	}
}

double Distance(double &dx,double &dy,double &dz,double x1,double y1,double z1,double x2,double y2,double z2){
	dx=x2-x1; dy=y2-y1; dz=z2-z1;    // (x1,y1,z1) から見た (x2,y2,z2) の位置
	return sqrt(dx*dx+dy*dy+dz*dz);
}

double DistanceX(double x1,double y1,double z1,double x2,double y2,double z2){
	double dx,dy,dz;
	
	Distance(dx,dy,dz,x1,y1,z1,x2,y2,z2);
	return dx;
}
double DistanceY(double x1,double y1,double z1,double x2,double y2,double z2){
	double dx,dy,dz;
	
	Distance(dx,dy,dz,x1,y1,z1,x2,y2,z2);
	return dy;
}
double DistanceZ(double x1,double y1,double z1,double x2,double y2,double z2){
	double dx,dy,dz;
	
	Distance(dx,dy,dz,x1,y1,z1,x2,y2,z2);
	return dz;
}

void IntersectionsCircles(double &x1,double &y1,double &x2,double &y2,
						         double X1,double Y1,double R1,double X2,double Y2,double R2){
	// 2つの円  (x-X1)^2 +(y-Y1)^2 =R1^2  ... (1)
	//          (x-X2)^2 +(y-Y2)^2 =R2^2  ... (2)
	// の交点 (x1,y1),(x2,y2) を求める．

	// 本関数は結果文字列の出力はしない（ sprintf(...) の処理に時間がかかるため）．
	// 文字列出力には
	// std::string IntersectionsCircles(double X1,double Y1,double R1,double X2,double Y2,double R2);
	// を使う．
	double A,B, a,b,c;
	
	if(Y1==Y2){ // 2円の中心が水平に並んでいる場合
		// (1)(2)より x=-B/A;
		A= 2.0*(X1-X2);
		B= X2*X2-X1*X1 +Y2*Y2-Y1*Y1 -R2*R2+R1*R1;
		// これを(1)に代入し，2次方程式の解の公式でyを求める．
		a= 1;
		b= -2.0*Y1;
		c= Y1*Y1 +(-B/A-X1)*(-B/A-X1) -R1*R1;
		if(b*b-4*a*c){
			x1=y1=x2=y2=0;   // 解なし
		}
		else{
			x1=-B/A; y1=(-b+sqrt(b*b-4*a*c))/2/a;
			x2=-B/A; y2=(-b-sqrt(b*b-4*a*c))/2/a;
		}
	}
	else{
		// 2円の交点を結ぶ直線は y=Ax+B  ((1)(2)から x^2, y^2 を消去し，地道に計算する）
		A= -(X1-X2)/(Y1-Y2);
		B= ( X2*X2-X1*X1 +Y2*Y2-Y1*Y1 -R2*R2+R1*R1 )/2/(Y2-Y1);
		// これを(1)に代入し，2次方程式の解の公式でxを求める．
		a= A*A+1;
		b= -2*X1+2*A*B-2*A*Y1;
		c= X1*X1 +B*B -2*B*Y1 +Y1*Y1 -R1*R1;
		if(b*b-4*a*c <0){
			x1=y1=x2=y2=0;   // 解なし
		}
		else{
			x1=(-b+sqrt(b*b-4*a*c))/2/a;  y1=A*x1+B;
			x2=(-b-sqrt(b*b-4*a*c))/2/a;  y2=A*x2+B;
		}
	}
}

std::string IntersectionsCircles(double X1,double Y1,double R1,double X2,double Y2,double R2){
	double x1,y1,x2,y2;
	char buf[200];

	IntersectionsCircles(x1,y1,x2,y2,X1,Y1,R1,X2,Y2,R2);
	sprintf(buf,"x1=%g y1=%g  x2=%g y2=%g\n", x1,y1,x2,y2);
	return buf;
}

double Circle3Points(double& xc,double& yc,double& r,
				   double x1,double y1,double x2,double y2,double x3,double y3){
	// (x1,y1),(x2,y2),(x3,y3)を通る円の中心(xc,yc)および半径rを計算する．また，rを戻り値とする．
	// 方法：
	//   円の中心は(x1,y1)と(x2,y2)の中点を通り，線分(x1,y1)-(x2,y2)に垂直な直線L1と，
	//   (x1,y1)と(x3,y3)についての同様な直線L2の交点である．
	//   L1の方程式は，((x1+x2)/2,(y1+y2)/2) を通り，(x1,y1)-(x2,y2)に垂直(内積が0)であることから，
	//     (x2-x1)x +(y2-y1)y = (y2^2-y1^2)/2 + (x2^2-x1^2)/2
    //   となる．L2についても同様．

	double delta,b1,b2;

	delta=(x2-x1)*(y3-y1)-(x3-x1)*(y2-y1);
	b1=(y2*y2-y1*y1)/2+(x2*x2-x1*x1)/2;
	b2=(y3*y3-y1*y1)/2+(x3*x3-x1*x1)/2;

	if(delta!=0){
		xc=(1/delta)*( (y3-y1)*b1-(y2-y1)*b2);
		yc=(1/delta)*(-(x3-x1)*b1+(x2-x1)*b2);
		r=sqrt((x1-xc)*(x1-xc)+(y1-yc)*(y1-yc));
	}
	else{
		r=0;
	}

	return r;
}

double TriangleArea(double x1,double y1,double z1,double x2,double y2,double z2,
					double x3,double y3,double z3){
	// (x1,y1,z1),(x2,y2,z2),(x3,y3,z3)を頂点とする三角形の面積
	vector<double> a,b;
	a=vector<double>(x2-x1,y2-y1,z2-z1);
	b=vector<double>(x3-x1,y3-y1,z3-z1);
	return abs(vProduct(a,b))/2;
}

double QuadrangleArea(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4){
	// （x1,y1),(x2,y2),(x3,y3),(x4,y4)を順に結んだ四角形の面積
	return TriangleArea(x1,y1,0,x2,y2,0,x4,y4,0) + TriangleArea(x2,y2,0,x3,y3,0,x4,y4,0);
}

double PolygonArea(int n,double *x,double *y){
	// n個の点 (x[1],y[1]),(x[2],y[2]),  ... (x[n],y[n]) を "この順番で"
	// 直線にて結んだ多角形の面積を求める．
	// 原点と(xi,yi)と(xi+1,yi+1) による三角形の面積を外積により計算し，
	// iについて和をとることによる．
	// 外積は符号付きであるから余分な部分は相殺され必要な面積が求められる．
	int i;
	double S=0;

	for(i=1; i<=n; ++i){
		if(i==1){
			S+= x[n]*y[1]-y[n]*x[1];
		}
		else{
			S+= x[i-1]*y[i]-y[i-1]*x[i];
		}
	}
	return fabs(S)/2;
}

double PolygonPointDistance(int n,double *x,double *y,double X,double Y){
	// n個の点 (x[1],y[1]),(x[2],y[2]),  ... (x[n],y[n]) による多角形と点(X,Y)の距離を求める．
	//  ・点と多角形の各辺までの距離の最小値を返す
	//  ・点が多角形の内側にあるときは正，外側にあるときは負の値を返す
	const double LN=1e30;
	int i,count;
	double min, d;
	vector<double> v,v1,v2;

	min=LN;
	v=vector<double>(X,Y,0);

	for(i=1; i<=n; ++i){
		v1=vector<double>(x[i],y[i],0);
		v2=vector<double>(x[i==n ? 1:i+1],y[i==n ? 1:i+1],0);
		d=DistancePointLinesegment(v,v1,v2);
		if(d<min) min=d;    // 最も近い辺までの距離
	}

	count=0;

	for(i=1; i<=n; ++i){
		if(IsCrossLines(X,Y, X+LN,Y, x[i],y[i], x[i==n ? 1:i+1],y[i==n ? 1:i+1])) count++;
	}

	return count%2==0 ? -min : min;   // Crossing Number Algorithm（交点数が偶数ならば点は多角形の外にある）
}

void ParaxialR(double &r1,double &r2,double &axis_deg,double axx,double ayy,double axy){
	// 曲面z(x,y)の展開の2次の項が，
	//   z=axx*x^2+ayy*y^2+axy*x*y
	// であるとき，主曲率半径r1,r2とr1の方向がx軸となす角axis_degを求める．
	// 参考： 草川徹 “レンズ光学”	(3.4.1～3.4.3式)
	double a1xx,a1yy, th;

	if(axx==ayy && axy>=0){
		axis_deg=45;
	}
	else if(axx==ayy && axy<0){
		axis_deg=-45;
	}
	else{
		axis_deg=atan(axy/(axx-ayy))/2*180/PI;
	}
	th=axis_deg*PI/180;
	a1xx=axx*cos(th)*cos(th)+ayy*sin(th)*sin(th)+axy*sin(th)*cos(th);
	a1yy=axx*sin(th)*sin(th)+ayy*cos(th)*cos(th)-axy*sin(th)*cos(th);
	r1= a1xx==0 ? 0 : 1.0/2.0/a1xx;
	r2= a1yy==0 ? 0 : 1.0/2.0/a1yy;
}

vector<double> NearestPointOnLine(vector<double> P0,vector<double> P1,vector<double> P2){
	// 点P1,P2を通る直線上で点P0に最も近い点（垂線の足）を返す．
	// (P2-P1) と P1+s(P2-P1)-Po は垂直であることより，
	//   (P2-P1)・(P1+s(P2-P1)-Po) =0
	//   s = {(P2-P1)・(P0-P1)} / {(P2-P1)・(P2-P1)}
	double s;
	
	s=sProduct(P2-P1,P0-P1)/sProduct(P2-P1,P2-P1);
	return P1+s*(P2-P1);
}

vector<double> NearestPointLineLine(vector<double> P,vector<double> V,vector<double> Po,vector<double> Vo){
	// 直線 Po+sVo に最も近い直線 P+tV 上の点を返す.
	// ここで，Po,Pは直線上の点，Vo,Vは方向ベクトル(単位ベクトルでなくてもよい)．
	// 最近接点を Po+tVo, P+sV とすると，最近接点同士を結ぶ直線は元の2直線と垂直だから，
	//    (P-Po+sV-tVo)・Vo =0
	//    (P-Po+sV-tVo)・V  =0
	// となる．これらから t を消去し，s について解けば求める点は P+sV となる．
	double s;

	s= ( sProduct(V,Po-P)+sProduct(Vo,P-Po)*sProduct(Vo,V)/sProduct(Vo,Vo) )
	  /( sProduct(V,V)-sProduct(V,Vo)/sProduct(Vo,Vo) );
	return P+s*V;
}

double DistancePointLinesegment(vector<double> P0,vector<double> P1,vector<double> P2){
	// 点P0と線分(P1,P2)の距離
	vector<double> p=NearestPointOnLine(P0,P1,P2);  // P0から直線(P1,P2)に下した垂線の足

	if( sProduct(p-P1,p-P2)<0 ){  // 垂線の足は線分上にある
		return abs(P0-p);
	}
	else{
		return Min( abs(P0-P1), abs(P0-P2) );
	}
}

bool IsCrossLines(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4){
	// 線分(x1,y1)-(x2,y2) と 線分(x3,y3)-(x4,y4) が交わるかどうかを判定する．
	vector<double> v1(x1,y1,0);
	vector<double> v2(x2,y2,0);
	vector<double> v3(x3,y3,0);
	vector<double> v4(x4,y4,0);
	double a,b;

	if(v1==v2 && v3==v4){
		return v1==v3;
	}
	else if(v1==v2){
		return abs(vProduct(v1-v3,v4-v3))==0 && sProduct(v1-v3,v1-v4)<=0;
	}
	else if(v3==v4){
		return abs(vProduct(v3-v1,v3-v2))==0 && sProduct(v3-v1,v3-v2)<=0;
	}
	else{
		a=sProduct(vProduct(v3-v1,v2-v1),vProduct(v4-v1,v2-v1)); // 直線v1-v2と線分v3-v4が交わるか(v3,v4が直線v1-v2の違う側にあるか）
		b=sProduct(vProduct(v1-v3,v4-v3),vProduct(v2-v3,v4-v3)); // 線分v1-v2と直線v3-v4が交わるか
		return a<=0 && b<=0;
	}
}

matrix<double> Tmatrix(double rox,double roy,double roz){
	// 空間ベクトル成分を回転rox,roy,rozを行った座標系に変換する行列（“レンズ光学”3.3.1 )
	matrix<double> Rx(3,3),Ry(3,3),Rz(3,3);
	rox*=PI/180;
	roy*=PI/180;
	roz*=PI/180;
	Rx.a[1][1]=1;         Rx.a[1][2]=0;         Rx.a[1][3]=0;
	Rx.a[2][1]=0;         Rx.a[2][2]=cos(rox);  Rx.a[2][3]=sin(rox);
	Rx.a[3][1]=0;         Rx.a[3][2]=-sin(rox); Rx.a[3][3]=cos(rox);
	Ry.a[1][1]=cos(roy);  Ry.a[1][2]=0;         Ry.a[1][3]=-sin(roy);
	Ry.a[2][1]=0;         Ry.a[2][2]=1;         Ry.a[2][3]=0;
	Ry.a[3][1]=sin(roy);  Ry.a[3][2]=0;         Ry.a[3][3]=cos(roy);
	Rz.a[1][1]=cos(roz);  Rz.a[1][2]=sin(roz);  Rz.a[1][3]=0;
	Rz.a[2][1]=-sin(roz); Rz.a[2][2]=cos(roz);  Rz.a[2][3]=0;
	Rz.a[3][1]=0;         Rz.a[3][2]=0;         Rz.a[3][3]=1;
	return Rz*Ry*Rx;
}

double Spline(double *x,double *y,int N,double xx,int SPDerivativeZero/*=0*/,int EPDerivativeZero/*=0*/){
	// N組のデータ x[0](==*x),x[1], ... x[N-1]
	//             y[0](==*x),y[1], ... y[N-1]
	// の3次スプライン補間により，xxに対応するyの値を求める．
	// SPDerivativeZero : 始点の一次導関数を0とするかどうか（偽のときは二次導関数が0）
	// EPDerivativeZero : 終点の一次導関数を0とするかどうか（偽のときは二次導関数が0）
	int i;
	cSpline spline;

	spline.StartPointDerivativeZero=SPDerivativeZero;
	spline.EndPointDerivativeZero=EPDerivativeZero;

	spline.SetN(N);
	for(i=1; i<=N; ++i){
		spline.SetX(i,x[i-1]);
		spline.SetY(i,y[i-1]);
	}
	return spline.y(xx);
}

double Lagrange(double *x,double *y,int N,double xx){
	// N組のデータ x[0](==*x),x[1], ... x[N-1]
	//             y[0](==*x),y[1], ... y[N-1]
	// のラグランジュ補間により，xxに対応するyの値を求める．
	int i,j;
	double yy,l;

	if(N<2) return 0;

	yy=0;
	for(i=0; i<=N-1; i++){
		l=1;
		for(j=0; j<=N-1; j++){
			if(i!=j){
				l*=((xx-x[j])/(x[i]-x[j]));
			}
		}
		yy+=y[i]*l;
	}

	return yy;
}

void Unwrap(double &ppre,double p0,double &psuc){
	// p0を基点として固定し，ppre,psucを変更して位相接続する．
	// p0とppreの間，またはp0とpsucの間の少なくとも片方には位相飛びがないとする．
	// 引数の単位はdegとする．
	double dp1,dp2;
	const double K=180;

	dp1=p0-ppre;
	dp2=psuc-p0;

	if(dp1*dp2<0){                  // 少なくともp0が尖端であるとする
		if(fabs(dp1)>fabs(dp2)){    // 変化が大きい方が位相飛びを起こしているとする
			if(dp1> K) ppre+=360;
			if(dp1<-K) ppre-=360;
		}
		else{
			if(dp2> K) psuc-=360;
			if(dp2<-K) psuc+=360;
		}
	}
}

void Unwrap(double *th_deg,int N){
	// N個のデータ th_deg[0]....th_deg[N-1]を位相接続し上書きする．th_deg[0]を固定する．
	int i,ii;
	double th0,thsuc0;
	
	th0=th_deg[0];

	for(i=1; i<=N-2; i++){
		thsuc0=th_deg[i+1];
		Unwrap(th_deg[i-1],th_deg[i],th_deg[i+1]);
		for(ii=i+2; ii<=N-1; ii++){
			th_deg[ii]+=th_deg[i+1]-thsuc0;
		}
	}

	for(i=0; i<=N-1; i++){
		th_deg[i]+=th_deg[0]-th0;
	}
}

void AddSCA(double& S,double& C,double& A_deg,
            double S1,double C1,double A1_deg,char plus_minus,double S2,double C2,double A2_deg){
	// (S1,C1,A1)と(S2,C2,A2)を加減算する．
	// 加算か減算かはplus_minusを '+' または '-' として指定する．	        
    double A,A1,A2, sgn,a,b;

    if     (plus_minus=='+') sgn=1;
    else if(plus_minus=='-') sgn=-1;
    else                     return;

    A1=A1_deg*PI/180;
    A2=A2_deg*PI/180;

    a=(C1/2)*sin(A1*2)+sgn*(C2/2)*sin(A2*2);
    b=(C1/2)*cos(A1*2)+sgn*(C2/2)*cos(A2*2);

    A= b==0 ? PI/4 : (1.0/2.0)*atan(a/b);
    C= 2.0*((C1/2)*cos(A1*2)+sgn*(C2/2)*cos(A2*2))/cos(A*2);

    if(C>0){
        C=-C;       // マイナス表示にする
        A+=PI/2;
    }

    if(A<0){
        A+=PI;      // A=0～π にする
    }
    
    S=(S1+C1/2)+sgn*(S2+C2/2)-C/2;
    A_deg=A*180/PI;
}
/*
void AddSCA(double& S,double& C,double& A,
            double S1,double C1,double A1,char plus_minus,double S2,double C2,double A2){
	// (S1,C1,A1)と(S2,C2,A2)を加減算する．
	// 加算か減算かはplus_minusを '+' または '-' として指定する．
	// 結果Cはマイナス表示とするが，引数C1,C2のプラスマイナスは不問．
	//
	// 円柱成分を大きさがCYL, 方位角がAXIS*2 のベクトルで表したとき，
	// 円柱成分の加減は，ベクトルの加減算となることによる．
	double sign, V1x,V1y,V2x,V2y,Vx,Vy;

	if(plus_minus=='+'){
		sign=1;
	}
	else if(plus_minus=='-'){
		sign=-1;
	}
	else{
		S=C=A=0;
		return;
	}

	A1*=PI/180;
	A2*=PI/180;

	V1x=C1*cos(A1*2); V1y=C1*sin(A1*2);
	V2x=C2*cos(A2*2); V2y=C2*sin(A2*2);
	Vx=V1x+sign*V2x;
	Vy=V1y+sign*V2y;
	C=sqrt(Vx*Vx+Vy*Vy);
	A=atan2(Vy,Vx)/2;
	C*=-1; A-=PI/2;                  // Cを-表示にする
	S=(S1+C1/2)+sign*(S2+C2/2)-C/2;  // 等価球面度数は保存されるので，(S1+C1/2)±(S2+C2/2)=S+C/2 より．

	A*=180/PI;
	if(A<0) A+=180;   // Aは0～180degとする．
}
*/
double dBToT(double dB){
	// デシベルで表された減衰量を透過率に換算する
	return pow(10,-dB/10);
}

int CSVReadLine(std::ifstream& from,double* data,int datas,int& valid_line){
	// csvファイル from から1行を読み込み，
	// data[0],data[1], ... ,data[datas-1] にデータを代入する
	int i;
	std::string s;

	if(std::getline(from,s)){
		s=replace(s,",",";");   // csvファイルの区切り "," を，word()での区切り ";" に置き換える
		valid_line=1; for(i=1; i<=datas; ++i) valid_line = valid_line && is_numeric(word(s,i,0));
		if(valid_line){  // 有効な行である（datas個の数値データがある）
			for(i=1; i<=datas; ++i){
				data[i-1]=atof(word(s,i,0).c_str());
			}
		}
		return 1;  // 行を読み込めた
	}
	else{
		return 0;  // 行を読み込めなかった（ファイル終端を超えている）
	}
}

std::string scmd_general_func(const std::string& com,int val){
	// val=tureのとき，一部のコマンドはBasic側にてval関数で処理するために
	// 数値を表す文字列のみ返す．
	//
	// 注意：書式%fは小数点以下の桁数が固定なので，値が小さいと有効桁がなくなる．
	//       %.15gがよい．

	std::string s;
	char buf[1000];
	std::string s0,s1,s2,s3,s4,s5,s6,s7;
	bool b1,b2,b3,b4,b5,b6,b7;

	s0=arg(com,0);

	if(s0=="??"){
		s+="AddSCA Circle3Points dBToT EnFace IntersectionsCircles OverlapAreaCircles SphSag ";
		s+="\n";
		return s;
	}
	if(s0=="AddSCA" || s0=="addsca"){
		double S,C,A,S1,C1,A1,S2,C2,A2;
		char plus_minus;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="AddSCA S1 C1 A1 '+'or'-' S2 C2 A2\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7){
			S1=atof(s1.c_str());
			C1=atof(s2.c_str());
			A1=atof(s3.c_str());
			plus_minus=s4[0];
			S2=atof(s5.c_str());
			C2=atof(s6.c_str());
			A2=atof(s7.c_str());
			AddSCA(S,C,A,S1,C1,A1,plus_minus,S2,C2,A2);
			sprintf(buf,"S=%g C=%g A=%g\n", S,C,A); s=buf;
		}
		return s;
	}
	if(s0=="Circle3Points" || s0=="circle3points"){
		double x1,y1,x2,y2,x3,y3,xc,yc,r;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="Circle3Points x1 y1 x2 y2 x3 y3\n";
		}
		else if(b1 && b2 && b3 && b5 && b6){
			x1=atof(s1.c_str());
			y1=atof(s2.c_str());
			x2=atof(s3.c_str());
			y2=atof(s4.c_str());
			x3=atof(s5.c_str());
			y3=atof(s6.c_str());
			sprintf(buf,"%.15g\n", Circle3Points(xc,yc,r,x1,y1,x2,y2,x3,y3)); s=buf;
		}
		return s;
	}
	if(s0=="dBToT" || s0=="dbtot"){
		double dB;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="dBToT dB_value\n";
		}
		else if(b1){
			dB=atof(s1.c_str());
			sprintf(buf,"%.15g\n",dBToT(dB)); s=buf;
		}
		return s;
	}
	if(s0=="EnFace" || s0=="enface"){
		int m,n;
		std::string in_filename,out_filename;
		double gamma;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3);
		s4=arg(com,4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s= "EnFace m n in_filename,out_filename [gamma=1]\n";
			s+="(ex) EnFace 256 256 oct_a_h# result 1\n";
		}
		else if(b1 && b2 && b5){
			m=atoi(s1.c_str());
			n=atoi(s2.c_str());
			in_filename=s3;
			out_filename=s4;
			gamma=atof(s5.c_str());
			EnFace(m,n,in_filename,out_filename,gamma);
		}
		else if(b1 && b2){
			m=atoi(s1.c_str());
			n=atoi(s2.c_str());
			in_filename=s3;
			out_filename=s4;
			EnFace(m,n,in_filename,out_filename);
		}
		return s;
	}
	if(s0=="IntersectionsCircles" || s0=="intersectionscircles"){
		double X1,Y1,R1,X2,Y2,R2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="IntersectionsCircles X1 Y1 R1 X2 Y2 R2\n";
		}
		else if(b1 && b2 && b3 && b5 && b6){
			X1=atof(s1.c_str());
			Y1=atof(s2.c_str());
			R1=atof(s3.c_str());
			X2=atof(s4.c_str());
			Y2=atof(s5.c_str());
			R2=atof(s6.c_str());
			s+=IntersectionsCircles(X1,Y1,R1,X2,Y2,R2);
		}
		return s;
	}
	if(s0=="OverlapAreaCircles" || s0=="overlapareacircles"){
		double D,a;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="OverlapAreaCircles diameter distance_between_centers\n";
		}
		else if( b1 && b2 ) {
			D=atof(s1.c_str());
			a=atof(s2.c_str());
			sprintf(buf,"%.15g\n",OverlapAreaCircles(D,a)); s=buf;
		}
		return s;
	}
	if(s0=="SphSag" || s0=="sphsag"){
		double r,h;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="SphSag r h\n";
		}
		else if( b1 && b2 ) {
			r=atof(s1.c_str());
			h=atof(s2.c_str());
			sprintf(buf,"%.15g\n",SphSag(r,h)); s=buf;
		}
		return s;
	}

	return s;
}



