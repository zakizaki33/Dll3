#include "stdafx.h"
#include "cBitmap.h"
#include "cImage.h"  // ����� cBitmap.h �ɏ����ƁCcBitmap.h �� cImage.h �����ݎQ�ƂƂȂ�C�r���h�G���[������

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

	if( gamma<=0 || gamma==1 ){  // gamma<=0 �܂��� gamma==1�ł���Ή������Ȃ�
		for(val=0; val<=255; val++){
			a[val]=val;
		}
	}
	else{
		// �������̂��߂̃e�[�u��
		a[0]=0;

		for(val=1; val<=254; ++val){
			a[val]=BYTE( pow( (double)val/255.0, 1/gamma ) * 255 );
		}
		
		a[255]=255;  // ��̃��[�v�ɓ����ƂȂ���255�ɂȂ�Ȃ�
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
	// �K���}�e�[�u�� a �͕ς��邪�C
	// �f�[�^ R,G,B ���͕̂ς��Ȃ��D
	// �K���}�l��GetRGB(),GetR(),... �̖߂�l�ɔ��f�����D
	this->gamma=gamma;
	make_gamma_table();
}

long cBitmap::GetRGB(int i,int j){
	// ����(i,j)=(1,1)���摜�̍���Ƃ���D
	// ����Cbmp�t�@�C���̉摜�f�[�^�����̐擪�͉摜�́g�����h�ł���D
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
	// ����(i,j)=(1,1)���摜�̍���Ƃ���D
	// ����Cbmp�t�@�C���̉摜�f�[�^�����̐擪�͉摜�́g�����h�ł���D
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
	// ��������_�Ƃ��āC�_(x1��f�Cy2��f) �� �_(x2��f�Cy2��f) �����Ԑ�����`��
	// x1,y1, x2,y2�� double�^�Ŏw��ł���
	int i1,j1,i2,j2, i,j;

	i1=(int)y1; j1=(int)x1;
	i2=(int)y2; j2=(int)x2;

	if(fabs(x1-x2)>=fabs(y1-y2)){  // �S�̂ɉ���
		if(j1>j2){
			Swap(i1,i2);
			Swap(j1,j2);
		}
		for(j=j1; j<=j2; ++j){
			i=(int)( (double)i1 +(double)(i2-i1)*(double)(j-j1)/(double)(j2-j1) );
			SetRGB(i,j,R,G,B);
		}
	}
	else{ // �S�̂ɏc��
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
	// ��������_�Ƃ��āC�_(x0��f,y0��f)�𒆐S�Ƃ��锼�aradius��f�̉~��`��
	// x0,y0,r��double�^�Ŏw��ł���
	// fill=1�̂Ƃ��~�̓�����h��Ԃ��Dfill=-1�̂Ƃ��͉~�̊O����h��Ԃ��D
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
	// R�����ɂ��摜�ɕϊ�����D
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		B[i][j]=G[i][j]=R[i][j];
	}
	return *this;
}

cBitmap& cBitmap::ToG(){
	// G�����ɂ��摜�ɕϊ�����D
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		R[i][j]=B[i][j]=G[i][j];
	}
	return *this;
}

cBitmap& cBitmap::ToB(){
	// B�����ɂ��摜�ɕϊ�����D
	int i,j;

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		R[i][j]=G[i][j]=B[i][j];
	}
	return *this;
}

cBitmap& cBitmap::ToWhite(){
	// ���F�ψ�ɂ���
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
	// �ő喾�邳��255�ɂȂ�悤�ɂ���
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
	// ��f��(�𑜓x)�� mxn �ɕύX����
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
		// x,y�����Ƃ���f������������Ƃ�
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
		// x,y�����ŉ�f�������ƌ������������Ă���Ƃ��C
		// �������ł̕��Ϗ����͍s���Ȃ��D
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
	// ���オ(i0,j0)�ŕ���width,������height�̋�`�����𔲂��o��
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
	// ���ډ~�Ńg���~���O����
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
	// �f�[�^R,G,B���̂��K���}�␳����Dgamma=1�����ϊ��D
	// ����CSetGamma() �ł�R,G,B���͕̂ύX����Ȃ��D
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
	// �V�����摜�͈͂� (i1,j1)-(i2,j2) �̋�`�Ƃ���D
	// (i1,j1)-(i2,j2)�ƌ��摜�̏c���䂪�قȂ�Ƃ��C�����ɂȂ�悤�� (i1,j1)-(i2,j2) �����������D

	int I,J, i,j;
	double x1,y1,x2,y2, h,w, ax,bx,ay,by;
	cBitmap temp;

	if(i1>i2) Swap(i1,i2);  // i1<i2�ɒ���
	if(j1>j2) Swap(j1,j2);

	y1=i1; x1=j1; y2=i2; x2=j2;   // double�^�֑������

	if( (y2-y1)/(x2-x1) > Height()/Width() ){
		// ��`�����摜���c���Ȃ̂ŕ����L����
		w=(y2-y1)* (double)n/(double)m;  // �V������
		x1=(x1+x2)/2-w/2;                // ��`�̒��� (x1+x2)/2 �͓������Ȃ�
		x2=x1+w;
	}
	else{
		// ��`�����摜��艡���Ȃ̂ō������L����
		h=(x2-x1)* (double)m/(double)n;  // �V��������
		y1=(y1+y2)/2-h/2;
		y2=y1+h;
	}

	temp.Resize(m,n);
	
	// �V�摜�̉�f(X,Y)�Ɍ��摜�̉�f(x,y)�̉�f�l��������
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
	// (i,j)��ʂ�\���������ɂ��ĉ摜��B�摜��R�摜�ɕ�������D�F����̕]���Ɏg���D
	int p,q;
	
	for(p=1; p<=m; ++p) for(q=1; q<=n; ++q){
		if     (j-n/10.0<q && q<j ) G[m+1-p][q]=R[m+1-p][q]=B[m+1-p][q]; // B�摜�ɂ���
		else if(j<=q && q<j+n/10.0) B[m+1-p][q]=G[m+1-p][q]=R[m+1-p][q]; // R�摜�ɂ���

		else if(i-m/10.0<p && p<i ) G[m+1-p][q]=R[m+1-p][q]=B[m+1-p][q];
		else if(i<=p && p<i+m/10.0) B[m+1-p][q]=G[m+1-p][q]=R[m+1-p][q];
	}
}

void cBitmap::ColorShiftCorrect(std::string filename,int inv/*=0*/){
	// �F�����␳����Dinv���^�̂Ƃ��͋t�ϊ�������i�␳��摜���猳�摜�𐄒肷��j�D
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
	// dist*100% ��3���c�Ȃ�␳����
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
	// �G�b�W����������D
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
	// ���f�B�A���t�B���^����p������
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
	// ��l������
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
	// �O�a��f�Ɉ������
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
	// gamma=1�����ϊ�
	FILE* fp;
	BMPFILEHEADER bf;
	BMPINFOHEADER bi;
	int i,j,k,imgsize;
	BYTE *buf,*buf0;
	int pad;
	const int BITCOUNT=24;  // ��F������1�o�C�g(0-255)�ŁC1��f������3�o�C�g=24�r�b�g �Ƃ���D

	if ( (fp=fopen((filename).c_str(), "wb"))==NULL ) return 0;
	// <note> "b"(�o�C�i�����[�h)���w�肵�Ȃ��ƁC�f�t�H���g�̓e�L�X�g���[�h�炵���C
	//        �Ӑ}���Ȃ��ϊ����s�Ȃ��C�摜�����ꂽ�肷��D

	GammaCompensation(gamma);

	// bmp�t�@�C���̎d�l��C��������1���C���̃f�[�^�̃o�C�g����4�̔{���C
	// �����r�b�g����4x8=32�̔{���ł���K�v������C
	// �e���C���̃f�[�^�̏I�[��0�̃f�[�^��pad�o�C�g������D
	pad= (n*BITCOUNT)%32 ==0 ? 0 : ( ((n*BITCOUNT)/32+1)*32-(n*BITCOUNT) )/8;
	imgsize=(n*BITCOUNT+pad*8)/8*m;

	// �w�b�_�쐬
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
	
	// �w�b�_�o��
	//   ��C�ɏ����������Ƃ��āC
	//       fwrite((void*)&bf,sizeof(BMPFILEHEADER),1,fp)
	//   �̂悤�ɂ����bmp�Ƃ��ĊJ���Ȃ��t�@�C�����o����.
	//   (�t�@�C���_���v�����bf.bfType�̂��ƂɃS�~��2�o�C�g�����Ă���)
	//   �܂��Csizeof(BMPFILEHEADER)�͐�������14�̂Ƃ���16��Ԃ��Ă����D
	//   sizeof�̓������[��̃T�C�Y��Ԃ��Ƃ̂��ƂŁC
	//   �\���̂ł͕K�������e�����o�̃T�C�Y�̘a�ɂ͂Ȃ�Ȃ��悤�ł���D
	//   ���̂悤�Ƀ����o���ЂƂ��������ނ̂��m���̂悤�ł���D

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
	
	// ��f�f�[�^�o��
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
	// <note> delete [] buf �Ƃ���ƃG���[�Dbuf��++�ɂ����ɔz��̐擪�ł͂Ȃ��D
	*/
	
	buf0=buf=new BYTE [3*n+pad];  // new BYTE[imgsize+1]; �Ƃ��Ĉ�C�ɏ������邱�Ƃ��\�����C��������Q���̂ŁC
	                              // 1���C������������D
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

	// �w�b�_�ǂݍ���
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

	// �p���b�g�ǂݍ���
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
	fseek(fp,bf.bfOffBits,SEEK_SET);  // �f�[�^���̐擪�ֈړ�

	// ��f�f�[�^�ǂݍ���
	switch(bitcount){
		case 1:
			for(i=1; i<=m; ++i) {
				for(j=1; j<=n; ++j){
					if(j%8==1) fread((void*)&buf,sizeof(BYTE),1,fp);
					switch(j%8){
						case 1: buf1=buf>>7;      break; // �擪�̃r�b�g�̒l
						case 2: buf1=(buf>>6)&1u; break;
						case 3: buf1=(buf>>5)&1u; break;
						case 4: buf1=(buf>>4)&1u; break;
						case 5: buf1=(buf>>3)&1u; break;
						case 6: buf1=(buf>>2)&1u; break;
						case 7: buf1=(buf>>1)&1u; break;
						case 0: buf1=buf&1u;      break; // �Ō�̃r�b�g�̒l
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
						case 1: buf1=(buf>>4)&15u; break; // �O��4�r�b�g�̒l
						case 0: buf1=buf&15u;      break; // �㔼4�r�b�g�̒l
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
	// �t�@�C������ filename_R,filename_G,filename_B �̃��m�N���摜�����ꂼ�� R,G,B�Ɋ��蓖�Ăč�������D
	cBitmap R,G,B;
	int i,j, m,n;

	m=n=0;
	if(R.Open(filename_R)){ m=R.m;        n=R.n;        }
	if(G.Open(filename_G)){ m=Min(m,G.m); n=Min(n,G.n); }
	if(B.Open(filename_B)){ m=Min(m,B.m); n=Min(n,B.n); }
	
	if(m==0 && n==0){
		return 0;   // �S�ẴI�[�v���Ɏ��s�����ꍇ
	}
	else{
		SetM(m);
		SetN(n);
	}

	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		this->R[i][j]=R.G[i][j];  // R,G,B�̓��CG[i][j]���\�Ƃ���
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

