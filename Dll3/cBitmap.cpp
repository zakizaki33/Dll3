#include "stdafx.h"
#include "cBitmap.h"
#include "cImage.h"  // これを cBitmap.h に書くと，cBitmap.h と cImage.h が相互参照となり，ビルドエラーが発生

void cBitmap::alloc(){
	int i,j;
	R=new BYTE* [m+1];
	G=new BYTE* [m+1];
	B=new BYTE* [m+1];
	for(i=0; i<=m; ++i){
		R[i]=new BYTE [n+1];
		G[i]=new BYTE [n+1];
		B[i]=new BYTE [n+1];
	}
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		R[i][j]=0;
		G[i][j]=0;
		B[i][j]=0;
	}
}

void cBitmap::free(){
	int i;
	for(i=0; i<=m; ++i){
		delete [] R[i];
		delete [] G[i];
		delete [] B[i];
	}
	delete [] R;
	delete [] G;
	delete [] B;
}

void cBitmap::make_gamma_table(){
	int val;

	if( gamma<=0 || gamma==1 ){  // gamma<=0 または gamma==1であれば何もしない
		for(val=0; val<=255; val++){
			a[val]=val;
		}
	}
	else{
		// 高速化のためのテーブル
		a[0]=0;

		for(val=1; val<=254; ++val){
			a[val]=BYTE( pow( (double)val/255.0, 1/gamma ) * 255 );
		}
		
		a[255]=255;  // 上のループに入れるとなぜか255にならない
	}	
}

cBitmap cBitmap::buf=cBitmap();
cBitmap cBitmap::buf2=cBitmap();

cBitmap::cBitmap(){
	m=n=4;
	alloc();
	SetGamma(1);
}

cBitmap::cBitmap(const cBitmap &x){
	int i,j;

	m=x.m;
	n=x.n;
	alloc();
	SetGamma(x.gamma);
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		B[i][j]=x.B[i][j];
		G[i][j]=x.G[i][j];
		R[i][j]=x.R[i][j];
	}
}

cBitmap::~cBitmap(){
	free();
}

cBitmap& cBitmap::operator=(const cBitmap& x){ 
	int i,j;

	free();
	m=x.m;
	n=x.n;
	alloc();
	SetGamma(x.gamma);
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		B[i][j]=x.B[i][j];
		G[i][j]=x.G[i][j];
		R[i][j]=x.R[i][j];
	}
	return *this;
}

int cBitmap::GetM() const{
	return m;
}
void cBitmap::SetM(int m){
	free();
	this->m=m;
	alloc();
}

int cBitmap::GetN() const{
	return n;
}
void cBitmap::SetN(int n){
	free();
	this->n=n;
	alloc();
}

double cBitmap::GetGamma(){
	return gamma;
}
void cBitmap::SetGamma(double gamma){
	// ガンマテーブル a は変えるが，
	// データ R,G,B 自体は変えない．
	// ガンマ値はGetRGB(),GetR(),... の戻り値に反映される．
	this->gamma=gamma;
	make_gamma_table();
}

long cBitmap::GetRGB(int i,int j){
	// 引数(i,j)=(1,1)を画像の左上とする．
	// 一方，bmpファイルの画像データ部分の先頭は画像の“左下”である．
	long rr,gg,bb;

	if( 1<=i && i<=m && 1<=j && j<=n ){
		rr=(long)a[R[m+1-i][j]];
		gg=(long)a[G[m+1-i][j]];
		bb=(long)a[B[m+1-i][j]];
		return bb*65536+gg*256+rr;
	}
	else{
		return 0;
	}
}

void cBitmap::SetRGB(int i,int j,int R,int G,int B){
	// 引数(i,j)=(1,1)を画像の左上とする．
	// 一方，bmpファイルの画像データ部分の先頭は画像の“左下”である．
	if( 1<=i && i<=m && 1<=j && j<=n ){
		if(R<0)        this->R[m+1-i][j]=0;
		else if(R>255) this->R[m+1-i][j]=255;
		else           this->R[m+1-i][j]=R;

		if(G<0)        this->G[m+1-i][j]=0;
		else if(G>255) this->G[m+1-i][j]=255;
		else           this->G[m+1-i][j]=G;

		if(B<0)        this->B[m+1-i][j]=0;
		else if(B>255) this->B[m+1-i][j]=255;
		else           this->B[m+1-i][j]=B;
	}
}

void cBitmap::SetRGB(int i,int j,long col){
	// col = B*256*256 + G*256 + R
	int B,G,R;

	B=col/65536;
	G=(col/256)%256;
	R=col%256;
	SetRGB(i,j,R,G,B);
}

void cBitmap::SetGrayScale(int i,int j,int val){
	SetRGB(i,j,val,val,val);
}


void cBitmap::SetLine(double x1,double y1,double x2,double y2,int R,int G,int B){
	// 左上を原点として，点(x1画素，y2画素) と 点(x2画素，y2画素) を結ぶ線分を描く
	// x1,y1, x2,y2は double型で指定できる
	int i1,j1,i2,j2, i,j;

	i1=(int)y1; j1=(int)x1;
	i2=(int)y2; j2=(int)x2;

	if(fabs(x1-x2)>=fabs(y1-y2)){  // 全体に横長
		if(j1>j2){
			Swap(i1,i2);
			Swap(j1,j2);
		}
		for(j=j1; j<=j2; ++j){
			i=(int)( (double)i1 +(double)(i2-i1)*(double)(j-j1)/(double)(j2-j1) );
			SetRGB(i,j,R,G,B);
		}
	}
	else{ // 全体に縦長
		if(i1>i2){
			Swap(i1,i2);
			Swap(j1,j2);
		}
		for(i=i1; i<=i2; ++i){
			j=(int)( (double)j1 +(double)(j2-j1)*(double)(i-i1)/(double)(i2-i1) );
			SetRGB(i,j,R,G,B);
		}
	}
}

void cBitmap::SetCircle(double x0,double y0,double radius,int R,int G,int B,int fill/*=0*/){
	// 左上を原点として，点(x0画素,y0画素)を中心とする半径radius画素の円を描く
	// x0,y0,rはdouble型で指定できる
	// fill=1のとき円の内側を塗りつぶす．fill=-1のときは円の外側を塗りつぶす．
	int i,j;
	double xpre,ypre, x,y;
	int th, step=2;
	
	if(fill==0){
		xpre=x0+radius; ypre=y0;
		for(th=step; th<=360; th+=step){
			x=x0+radius*cos((double)th*PI/180);
			y=y0+radius*sin((double)th*PI/180);
			SetLine(xpre,ypre,x,y,R,G,B);
			xpre=x; ypre=y;
		}
	}
	else{
		for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
			x=(double)j;
			y=(double)i;

			if(fill==1){
				if( (x-x0)*(x-x0)+(y-y0)*(y-y0)<=radius*radius ) SetRGB(i,j,R,G,B);
			}
			else if(fill==-1){
				if( (x-x0)*(x-x0)+(y-y0)*(y-y0)>=radius*radius ) SetRGB(i,j,R,G,B);
			}
		}
	}
}

int cBitmap::GetR(int i,int j) const{
	if( 1<=i && i<=m && 1<=j && j<=n ){
		return (int)a[R[m+1-i][j]];
	}
	else{
		return 0;
	}
}

int cBitmap::GetG(int i,int j) const{
	if( 1<=i && i<=m && 1<=j && j<=n ){
		return (int)a[G[m+1-i][j]];
	}
	else{
		return 0;
	}
}

int cBitmap::GetB(int i,int j) const{
	if( 1<=i && i<=m && 1<=j && j<=n ){
		return (int)a[B[m+1-i][j]];
	}
	else{
		return 0;
	}
}

cBitmap& cBitmap::ToR(){
	// R成分による画像に変換する．
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		B[i][j]=G[i][j]=R[i][j];
	}
	return *this;
}

cBitmap& cBitmap::ToG(){
	// G成分による画像に変換する．
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		R[i][j]=B[i][j]=G[i][j];
	}
	return *this;
}

cBitmap& cBitmap::ToB(){
	// B成分による画像に変換する．
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		R[i][j]=G[i][j]=B[i][j];
	}
	return *this;
}

cBitmap& cBitmap::ToWhite(){
	// 白色均一にする
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		B[i][j]=255;
		G[i][j]=255;
		R[i][j]=255;
	}
	return *this;
}

cBitmap& cBitmap::ToNeg(){
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		B[i][j]=255-B[i][j];
		G[i][j]=255-G[i][j];
		R[i][j]=255-R[i][j];
	}
	return *this;
}

cBitmap& cBitmap::Expand(){
	// 最大明るさが255になるようにする
	int i,j;
	BYTE max;
	double c;

	max=0;
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		if(R[i][j]>max) max=R[i][j];
		if(G[i][j]>max) max=G[i][j];
		if(B[i][j]>max) max=B[i][j];
	}
	c=255.0/(double)max;
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		R[i][j]*=(BYTE)c;
		G[i][j]*=(BYTE)c;
		B[i][j]*=(BYTE)c;
	}
	return *this;
}

int cBitmap::Height(){
	return m;
}

int cBitmap::Width(){
	return n;
}

void cBitmap::Resize(int m,int n){
	// 画素数(解像度)を mxn に変更する
	int i,j,ii,jj;
	double **rr,**gg,**bb,**count;
	double ratio_y,ratio_x;

	if(n==0) n=m;
	
	rr=new double* [m+1];
	gg=new double* [m+1];
	bb=new double* [m+1];
	count=new double* [m+1];
	
	for(ii=0; ii<=m; ++ii){
		rr[ii]=new double [n+1];
		gg[ii]=new double [n+1];
		bb[ii]=new double [n+1];
		count[ii]=new double [n+1];
	}

	for(ii=1; ii<=m; ++ii) for(jj=1; jj<=n; ++jj){
		rr[ii][jj]=0;
		gg[ii][jj]=0;
		bb[ii][jj]=0;
		count[ii][jj]=0;
	}

	ratio_y=(double)m/(double)(this->m);
	ratio_x=(double)n/(double)(this->n);

	if(m<=this->m && n<=this->n){
		// x,y方向とも画素数が減少するとき
		for(i=1; i<=this->m; ++i) for(j=1; j<=this->n; ++j){
			ii=(int)((double)i*ratio_y);
			jj=(int)((double)j*ratio_x);
			if(ii<1) ii=1;
			if(ii>m) ii=m;
			if(jj<1) jj=1;
			if(jj>n) jj=n;
			rr[ii][jj]+=(double)(int)R[i][j];
			gg[ii][jj]+=(double)(int)G[i][j];
			bb[ii][jj]+=(double)(int)B[i][j];
			count[ii][jj]+=1;
		}
	}
	else{
		// x,y方向で画素数増加と減少が混ざっているとき，
		// 減少側での平均処理は行われない．
		for(ii=1; ii<=m; ++ii) for(jj=1; jj<=n; ++jj){
			i=(int)((double)ii/ratio_y);
			j=(int)((double)jj/ratio_x);
			if(i<1) i=1;
			if(i>this->m) i=this->m;
			if(j<1) j=1;
			if(j>this->n) j=this->n;
			rr[ii][jj]=(double)(int)R[i][j];
			gg[ii][jj]=(double)(int)G[i][j];
			bb[ii][jj]=(double)(int)B[i][j];
			count[ii][jj]=1;
		}
	}

	SetM(m);
	SetN(n);

	for(ii=1; ii<=m; ++ii) for(jj=1; jj<=n; ++jj){
		R[ii][jj]=(BYTE)(rr[ii][jj]/count[ii][jj]);
		G[ii][jj]=(BYTE)(gg[ii][jj]/count[ii][jj]);
		B[ii][jj]=(BYTE)(bb[ii][jj]/count[ii][jj]);
	}

	for(ii=0; ii<=m; ++ii){
		delete [] rr[ii];
		delete [] gg[ii];
		delete [] bb[ii];
		delete [] count[ii];
	}

	delete [] rr;
	delete [] gg;
	delete [] bb;
	delete [] count;
}

void cBitmap::Resize(double ratio){
	int m,n;

	m=(int)((double)this->m*ratio);
	n=(int)((double)this->n*ratio);

	Resize(m,n);
}

void cBitmap::Trim(int i0,int j0,int height,int width){
	// 左上が(i0,j0)で幅がwidth,高さがheightの矩形部分を抜き出す
	int i,j;
	cBitmap X;

	X.SetM(height);
	X.SetN(width);
	for(i=1; i<=X.m; ++i) for(j=1; j<=X.n; ++j){
		X.R[X.m+1-i][j]=GetR(i0+i-1,j0-1+j);
		X.G[X.m+1-i][j]=GetG(i0+i-1,j0-1+j);
		X.B[X.m+1-i][j]=GetB(i0+i-1,j0-1+j);
	}
	*this=X;
}

void cBitmap::TrimCircle(){
	// 内接円でトリミングする
	int i,j;
	double rx,ry;

	rx=n/2; ry=m/2;
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		if( (i-ry)*(i-ry)/ry/ry + (j-rx)*(j-rx)/rx/rx > 1 ){
			B[i][j]=0;
			G[i][j]=0;
			R[i][j]=0;
		}
	}
}

void cBitmap::GammaCompensation(double gamma){
	// データR,G,B自体をガンマ補正する．gamma=1が無変換．
	// 一方，SetGamma() ではR,G,B自体は変更されない．
	int i,j;
	double gamma0;

	gamma0=this->gamma;
	SetGamma(gamma);

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		B[i][j]=a[B[i][j]];
		G[i][j]=a[G[i][j]];
		R[i][j]=a[R[i][j]];
	}

	SetGamma(gamma0);
}

void cBitmap::Zoom(int i1,int j1,int i2,int j2){
	// 新しい画像範囲を (i1,j1)-(i2,j2) の矩形とする．
	// (i1,j1)-(i2,j2)と元画像の縦横比が異なるとき，同じになるように (i1,j1)-(i2,j2) が調整される．

	int I,J, i,j;
	double x1,y1,x2,y2, h,w, ax,bx,ay,by;
	cBitmap temp;

	if(i1>i2) Swap(i1,i2);  // i1<i2に直す
	if(j1>j2) Swap(j1,j2);

	y1=i1; x1=j1; y2=i2; x2=j2;   // double型へ代入する

	if( (y2-y1)/(x2-x1) > Height()/Width() ){
		// 矩形が元画像より縦長なので幅を広げる
		w=(y2-y1)* (double)n/(double)m;  // 新しい幅
		x1=(x1+x2)/2-w/2;                // 矩形の中央 (x1+x2)/2 は動かさない
		x2=x1+w;
	}
	else{
		// 矩形が元画像より横長なので高さを広げる
		h=(x2-x1)* (double)m/(double)n;  // 新しい高さ
		y1=(y1+y2)/2-h/2;
		y2=y1+h;
	}

	temp.Resize(m,n);
	
	// 新画像の画素(X,Y)に元画像の画素(x,y)の画素値を代入する
	//   y=ay*Y+by
	//   x=ax*X+bx

	ay=(y2-y1)/((double)m-1);  by=y1-ay;
	ax=(x2-x1)/((double)n-1);  bx=x1-ax;
	
	for(I=1; I<=m; ++I){
		i=(int)(ay*(double)I+by);
		
		for(J=1; J<=n; ++J){
			j=(int)(ax*(double)J+bx);

			if(1<=i && i<=m && 1<=j && j<=n){
				temp.B[m+1-I][J]=B[m+1-i][j];
				temp.G[m+1-I][J]=G[m+1-i][j];
				temp.R[m+1-I][J]=R[m+1-i][j];
			}
			else{
				temp.B[m+1-I][J]=0;
				temp.G[m+1-I][J]=0;
				temp.R[m+1-I][J]=0;
			}
		} // next J
	} // next I

	*this=temp;
}

void cBitmap::SeeColorShift(int i,int j){
	// (i,j)を通る十字線を境にして画像をB画像とR画像に分割する．色ずれの評価に使う．
	int p,q;
	
	for(p=1; p<=m; ++p) for(q=1; q<=n; ++q){
		if     (j-n/10.0<q && q<j ) G[m+1-p][q]=R[m+1-p][q]=B[m+1-p][q]; // B画像にする
		else if(j<=q && q<j+n/10.0) B[m+1-p][q]=G[m+1-p][q]=R[m+1-p][q]; // R画像にする

		else if(i-m/10.0<p && p<i ) G[m+1-p][q]=R[m+1-p][q]=B[m+1-p][q];
		else if(i<=p && p<i+m/10.0) B[m+1-p][q]=G[m+1-p][q]=R[m+1-p][q];
	}
}

void cBitmap::ColorShiftCorrect(std::string filename,int inv/*=0*/){
	// 色ずれを補正する．invが真のときは逆変換をする（補正後画像から元画像を推定する）．
	cBitmap x;
	cImage b,g,r;

	x=*this;
	b.FromBitmap(x.ToB());
	b.SetCxCy(filename,"BLUE");
	b.Transform(0,inv);

	x=*this;
	g.FromBitmap(x.ToG());
	g.SetCxCy(filename,"GREEN");
	g.Transform(0,inv);

	x=*this;
	r.FromBitmap(x.ToR());
	r.SetCxCy(filename,"RED");
	r.Transform(0,inv);
	
	*this=cImage::ToBitmap(r,g,b);
}

void cBitmap::DistCorrect(double dist){
	// dist*100% の3次歪曲を補正する
	cBitmap x;
	cImage b,g,r;

	x=*this;
	b.FromBitmap(x.ToB());
	b.DistCorrect(dist);

	x=*this;
	g.FromBitmap(x.ToG());
	g.DistCorrect(dist);

	x=*this;
	r.FromBitmap(x.ToR());
	r.DistCorrect(dist);
	
	*this=cImage::ToBitmap(r,g,b);
}

void cBitmap::Edge(){
	// エッジを強調する．
	cBitmap x;
	cImage b,g,r;

	x=*this;
	b.FromBitmap(x.ToB());
	b.Edge();
	x=*this;
	g.FromBitmap(x.ToG());
	g.Edge();
	x=*this;
	r.FromBitmap(x.ToR());
	r.Edge();

	*this=cImage::ToBitmap(r,g,b);
}

void cBitmap::Median(){
	// メディアンフィルタを作用させる
	cBitmap x;
	cImage b,g,r;

	x=*this;
	b.FromBitmap(x.ToB());
	b.MedianFilter();
	x=*this;
	g.FromBitmap(x.ToG());
	g.MedianFilter();
	x=*this;
	r.FromBitmap(x.ToR());
	r.MedianFilter();

	*this=cImage::ToBitmap(r,g,b);
}

void cBitmap::Binarize(double threshold){
	// 二値化する
	cBitmap x;
	cImage b,g,r;

	if(threshold<0) threshold=0;
	if(threshold>255) threshold=255;

	x=*this;
	b.FromBitmap(x.ToB());
	b.Binarize(threshold);
	x=*this;
	g.FromBitmap(x.ToG());
	g.Binarize(threshold);
	x=*this;
	r.FromBitmap(x.ToR());
	r.Binarize(threshold);

	*this=cImage::ToBitmap(r,g,b);
}

void cBitmap::MarkSaturation(){
	// 飽和画素に印をつける
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		if( B[i][j]==255 || G[i][j]==255 || R[i][j]==255 ){
			B[i][j]=0;
			G[i][j]=0;
			R[i][j]=255;
		}
		else if( B[i][j]==0 && G[i][j]==0 && R[i][j]==0 ){
			B[i][j]=255;
			G[i][j]=0;
			R[i][j]=0;
		}
	}

}

int cBitmap::Save(std::string filename,double gamma/*=1*/){
	// gamma=1が無変換
	FILE* fp;
	BMPFILEHEADER bf;
	BMPINFOHEADER bi;
	int i,j,k,imgsize;
	BYTE *buf,*buf0;
	int pad;
	const int BITCOUNT=24;  // 一色あたり1バイト(0-255)で，1画素あたり3バイト=24ビット とする．

	if ( (fp=fopen((filename).c_str(), "wb"))==NULL ) return 0;
	// <note> "b"(バイナリモード)を指定しないと，デフォルトはテキストモードらしく，
	//        意図しない変換が行なわれ，画像がずれたりする．

	GammaCompensation(gamma);

	// bmpファイルの仕様上，水平方向1ラインのデータのバイト数は4の倍数，
	// 即ちビット数は4x8=32の倍数である必要があり，
	// 各ラインのデータの終端に0のデータをpadバイト加える．
	pad= (n*BITCOUNT)%32 ==0 ? 0 : ( ((n*BITCOUNT)/32+1)*32-(n*BITCOUNT) )/8;
	imgsize=(n*BITCOUNT+pad*8)/8*m;

	// ヘッダ作成
	bf.bfType=*(WORD*)"BM";
	bf.bfSize=14+40+imgsize;
	bf.bfReserved1=0;
	bf.bfReserved2=0;
	bf.bfOffBits=14+40;

	bi.biSize=40;
	bi.biWidth=n;
	bi.biHeight=m;
	bi.biPlanes=1;
	bi.biBitCount=(WORD)BITCOUNT;    
	bi.biCompression=0;
	bi.biSizeImage=imgsize;
	bi.biXPelsPerMeter=0;
	bi.biYPelsPerMeter=0;
	bi.biClrUsed=0;
	bi.biClrImportant=0;
	
	// ヘッダ出力
	//   一気に書き込もうとして，
	//       fwrite((void*)&bf,sizeof(BMPFILEHEADER),1,fp)
	//   のようにするとbmpとして開けないファイルが出来る.
	//   (ファイルダンプするとbf.bfTypeのあとにゴミが2バイト入っていた)
	//   また，sizeof(BMPFILEHEADER)は正しくは14のところ16を返していた．
	//   sizeofはメモリー上のサイズを返すとのことで，
	//   構造体では必ずしも各メンバのサイズの和にはならないようである．
	//   次のようにメンバをひとつずつ書き込むのが確実のようである．

	fwrite((void*)&(bf.bfType),          sizeof(WORD), 1,fp);
	fwrite((void*)&(bf.bfSize),          sizeof(DWORD),1,fp);
	fwrite((void*)&(bf.bfReserved1),     sizeof(WORD), 1,fp);
	fwrite((void*)&(bf.bfReserved2),     sizeof(WORD), 1,fp);
	fwrite((void*)&(bf.bfOffBits),       sizeof(DWORD),1,fp);

	fwrite((void*)&(bi.biSize),          sizeof(DWORD),1,fp);
	fwrite((void*)&(bi.biWidth),         sizeof(LONG), 1,fp);
	fwrite((void*)&(bi.biHeight),        sizeof(LONG), 1,fp);
	fwrite((void*)&(bi.biPlanes),        sizeof(WORD), 1,fp);
	fwrite((void*)&(bi.biBitCount),      sizeof(WORD), 1,fp);
	fwrite((void*)&(bi.biCompression),   sizeof(DWORD),1,fp);
	fwrite((void*)&(bi.biSizeImage),     sizeof(DWORD),1,fp);
	fwrite((void*)&(bi.biXPelsPerMeter), sizeof(LONG), 1,fp);
	fwrite((void*)&(bi.biYPelsPerMeter), sizeof(LONG), 1,fp);
	fwrite((void*)&(bi.biClrUsed),       sizeof(DWORD),1,fp);
	fwrite((void*)&(bi.biClrImportant),  sizeof(DWORD),1,fp);
	
	// 画素データ出力
	/*
	buf0=buf=new BYTE [imgsize+1];
	for(i=1; i<=m; ++i) {
		for(j=1; j<=n; ++j){		
			*(buf++)=B[i][j];
			*(buf++)=G[i][j];
			*(buf++)=R[i][j];	
		}
		for(k=1; k<=pad; ++k){
			*(buf++)=0;
		}
	}
	fwrite((void*)buf0,sizeof(BYTE),imgsize,fp);
	delete [] buf0;  
	// <note> delete [] buf とするとエラー．bufは++により既に配列の先頭ではない．
	*/
	
	buf0=buf=new BYTE [3*n+pad];  // new BYTE[imgsize+1]; として一気に処理することも可能だが，メモリを浪費するので，
	                              // 1ラインずつ処理する．
	for(i=1; i<=m; ++i) {
		for(j=1; j<=n; ++j){		
			*(buf++)=B[i][j];
			*(buf++)=G[i][j];
			*(buf++)=R[i][j];	
		}
		for(k=1; k<=pad; ++k){
			*(buf++)=0;
		}
		fwrite((void*)buf0,sizeof(BYTE),3*n+pad,fp);
		buf=buf0;
	}
	
	delete [] buf0;

	fclose(fp);
	return 1;
}

int cBitmap::Open(std::string filename){
	FILE* fp;
	BMPFILEHEADER bf;
	BMPINFOHEADER bi;
	int i,j,k;
	BYTE buf,buf1;
	int pad,bitcount;
	PALLET pallet[256];
	int pallets;

	if ( (fp=fopen((filename).c_str(), "rb"))==NULL ) return 0;

	// ヘッダ読み込み
	fread((void*)&(bf.bfType),          sizeof(WORD), 1,fp);
	fread((void*)&(bf.bfSize),          sizeof(DWORD),1,fp);
	fread((void*)&(bf.bfReserved1),     sizeof(WORD), 1,fp);
	fread((void*)&(bf.bfReserved2),     sizeof(WORD), 1,fp);
	fread((void*)&(bf.bfOffBits),       sizeof(DWORD),1,fp);
	
	if(bf.bfType != *(WORD*)"BM") return 0;

	fread((void*)&(bi.biSize),          sizeof(DWORD),1,fp);
	fread((void*)&(bi.biWidth),         sizeof(LONG), 1,fp);
	fread((void*)&(bi.biHeight),        sizeof(LONG), 1,fp);
	fread((void*)&(bi.biPlanes),        sizeof(WORD), 1,fp);
	fread((void*)&(bi.biBitCount),      sizeof(WORD), 1,fp);
	fread((void*)&(bi.biCompression),   sizeof(DWORD),1,fp);
	fread((void*)&(bi.biSizeImage),     sizeof(DWORD),1,fp);
	fread((void*)&(bi.biXPelsPerMeter), sizeof(LONG), 1,fp);
	fread((void*)&(bi.biYPelsPerMeter), sizeof(LONG), 1,fp);
	fread((void*)&(bi.biClrUsed),       sizeof(DWORD),1,fp);
	fread((void*)&(bi.biClrImportant),  sizeof(DWORD),1,fp);

	bitcount=(int)bi.biBitCount;

	if     (bitcount==1) pallets=2;
	else if(bitcount==4) pallets=16;
	else if(bitcount==8) pallets=256;
	else if(bitcount==24 || bitcount==32) pallets=0;
	else return 0;

	// パレット読み込み
	for(i=0; i<=pallets-1; ++i){
		fread((void*)&pallet[i].rgbBlue,     sizeof(BYTE),1,fp);
		fread((void*)&pallet[i].rgbGreen,    sizeof(BYTE),1,fp);
		fread((void*)&pallet[i].rgbRed,      sizeof(BYTE),1,fp);
		fread((void*)&pallet[i].rgbReserved, sizeof(BYTE),1,fp);
	}

	m=bi.biHeight;
	n=bi.biWidth;
	alloc();
	pad= (n*bitcount)%32 ==0 ? 0 : ( ((n*bitcount)/32+1)*32-(n*bitcount) )/8; 
	fseek(fp,bf.bfOffBits,SEEK_SET);  // データ部の先頭へ移動

	// 画素データ読み込み
	switch(bitcount){
		case 1:
			for(i=1; i<=m; ++i) {
				for(j=1; j<=n; ++j){
					if(j%8==1) fread((void*)&buf,sizeof(BYTE),1,fp);
					switch(j%8){
						case 1: buf1=buf>>7;      break; // 先頭のビットの値
						case 2: buf1=(buf>>6)&1u; break;
						case 3: buf1=(buf>>5)&1u; break;
						case 4: buf1=(buf>>4)&1u; break;
						case 5: buf1=(buf>>3)&1u; break;
						case 6: buf1=(buf>>2)&1u; break;
						case 7: buf1=(buf>>1)&1u; break;
						case 0: buf1=buf&1u;      break; // 最後のビットの値
					}
					B[i][j]=pallet[buf1].rgbBlue;
					G[i][j]=pallet[buf1].rgbGreen;
					R[i][j]=pallet[buf1].rgbRed;
				}
				for(k=1; k<=pad; ++k){
					fread((void*)&buf,sizeof(BYTE),1,fp);
				}
			}
			break;
		
		case 4:
			for(i=1; i<=m; ++i){
				for(j=1; j<=n; ++j){
					if(j%2==1) fread((void*)&buf,sizeof(BYTE),1,fp);
					switch(j%2){
						case 1: buf1=(buf>>4)&15u; break; // 前半4ビットの値
						case 0: buf1=buf&15u;      break; // 後半4ビットの値
					}
					B[i][j]=pallet[buf1].rgbBlue;
					G[i][j]=pallet[buf1].rgbGreen;
					R[i][j]=pallet[buf1].rgbRed;
				}
				for(k=1; k<=pad; ++k){
					fread((void*)&buf,sizeof(BYTE),1,fp);
				}
			}
			break;
		
		case 8:
			for(i=1; i<=m; ++i){
				for(j=1; j<=n; ++j){
					fread((void*)&buf,sizeof(BYTE),1,fp);
					B[i][j]=pallet[buf].rgbBlue;
					G[i][j]=pallet[buf].rgbGreen;
					R[i][j]=pallet[buf].rgbRed;
				}
				for(k=1; k<=pad; ++k){
					fread((void*)&buf,sizeof(BYTE),1,fp);
				}
			}
			break;
		
		case 24:
			for(i=1; i<=m; ++i){
				for(j=1; j<=n; ++j){
					fread((void*)&B[i][j],sizeof(BYTE),1,fp);
					fread((void*)&G[i][j],sizeof(BYTE),1,fp);
					fread((void*)&R[i][j],sizeof(BYTE),1,fp);
				}
				for(k=1; k<=pad; ++k){
					fread((void*)&buf,sizeof(BYTE),1,fp);
				}
			}
			break;
		
		case 32:
			for(i=1; i<=m; ++i) {
				for(j=1; j<=n; ++j){
					fread((void*)&B[i][j],sizeof(BYTE),1,fp);
					fread((void*)&G[i][j],sizeof(BYTE),1,fp);
					fread((void*)&R[i][j],sizeof(BYTE),1,fp);
					fread((void*)&buf    ,sizeof(BYTE),1,fp);
				}
			}
			break;
	}

	fclose(fp);

	return 1;
}

int cBitmap::Merge(std::string filename_R,std::string filename_G,std::string filename_B){
	// ファイル名が filename_R,filename_G,filename_B のモノクロ画像をそれぞれ R,G,Bに割り当てて合成する．
	cBitmap R,G,B;
	int i,j, m,n;

	m=n=0;
	if(R.Open(filename_R)){ m=R.m;        n=R.n;        }
	if(G.Open(filename_G)){ m=Min(m,G.m); n=Min(n,G.n); }
	if(B.Open(filename_B)){ m=Min(m,B.m); n=Min(n,B.n); }
	
	if(m==0 && n==0){
		return 0;   // 全てのオープンに失敗した場合
	}
	else{
		SetM(m);
		SetN(n);
	}

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		this->R[i][j]=R.G[i][j];  // R,G,Bの内，G[i][j]を代表とする
		this->G[i][j]=G.G[i][j];
		this->B[i][j]=B.G[i][j];
	}

	return 1;
}

void cBitmap::ToBuf(){
	buf=*this;
}

void cBitmap::FromBuf(){
	*this=buf;
}

void cBitmap::ToBuf2(){
	buf2=*this;
}

void cBitmap::FromBuf2(){
	*this=buf2;
}

