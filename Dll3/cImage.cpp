#include "stdafx.h"
#include "cImage.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

void cImage::gravity_center(int& ig,int& jg) {
	int i,j;
	double gi=0,gj=0,x=0;

	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		gi+=GetIntensity(i,j)*i;
		gj+=GetIntensity(i,j)*j;
		x+=GetIntensity(i,j);
	}
	if(x==0){
		gi=gj=0;
	}
	else{
		gi/=x;
		gj/=x;
	}
	ig=(int)Round(gi,0); jg=(int)Round(gj,0);
}

cImage cImage::buf=cImage();

// public members /////////////////////////////////////////

cImage::cImage(){
	// �y���Ӂz
	//   ypixels,xpixels��傫�Ȓl(���S�Ȃ�)�ŏ���������ƁC
	//   cImage�I�v�W�F�N�g�������o�Ɏ��N���X����ʂ̃��������g��
	ypixels=1;
	xpixels=1;
	ypitch=xpitch=0;
	A.redim(ypixels,xpixels);

	int i,j;
	for(i=0; i<=MAX_ORDER; ++i) for(j=0; j<=MAX_ORDER; ++j){
		Cx[i][j]=Cy[i][j]=0;
	}
}

cImage::~cImage(){}

cImage& cImage::ToSquare(){
	// A�̗v�f���Βl�̓��ɕϊ�����D�ʑ���0�Ƃ���D
	int i,j;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		A[i][j]=sqabs(A[i][j]);
	}
	return *this;
}

cImage& cImage::ToSqrt(){
	// A�̗v�f���Βl�̕������ɕϊ�����D�ʑ���0�Ƃ���D
	int i,j;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		A[i][j]=sqrt(abs(A[i][j]));
	}
	return *this;
}

cImage& cImage::ToAbs(){
	// A�̗v�f���Βl�ɕϊ�����D�ʑ���0�Ƃ���D
	int i,j;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		A[i][j]=abs(A[i][j]);
	}
	return *this;
}

void cImage::SetData(const matrix<complex>& a) {
	A=a;
	ypixels=a.rows();
	xpixels=a.columns();
}

int cImage::GetYpixels() const{
	return ypixels;
}
void cImage::SetYpixels(int value){
	if(value>0){
		ypixels=value;
		A.redim(ypixels,xpixels);
	}
}

int cImage::GetXpixels() const{
	return xpixels;
}
void cImage::SetXpixels(int value){
	if(value>0){
		xpixels=value;
		A.redim(ypixels,xpixels);
	}
}

int cImage::GetPixels() const{
	return ypixels==xpixels ? ypixels : 0;
}
void cImage::SetPixels(int value){
	if(value>0){
		ypixels=xpixels=value;
		A.redim(ypixels,xpixels);
	}
}

double cImage::GetYpitch() const{
	return ypitch;
}
void cImage::SetYpitch(double value){
	if(value>0){
		ypitch=value;
	}	
}

double cImage::GetXpitch() const{
	return xpitch;
}
void cImage::SetXpitch(double value){
	if(value>0){
		xpitch=value;
	}	
}

double cImage::GetPitch() const{
	return ypitch==xpitch ? ypitch : 0;
}
void cImage::SetPitch(double value){
	if(value>0){
		ypitch=xpitch=value;
	}
}

double cImage::GetHeight() const{
	return ypitch*((double)ypixels-1);
}
void cImage::SetHeight(double value){
	ypitch=value/((double)ypixels-1);
	if(xpitch==0) xpitch=ypitch;
}

double cImage::GetWidth() const{
	return xpitch*((double)xpixels-1);
}
void cImage::SetWidth(double value){
	xpitch=value/((double)xpixels-1);
	if(ypitch==0) ypitch=xpitch;
}

void cImage::Clear(){
	A=zero(A);
}

void cImage::Resize(int ypixels,int xpixels){
	// ��f��(�𑜓x)�� ypixels x xpixels �ɕύX����
	int i,j,ii,jj;
	matrix<complex> AA(ypixels,xpixels);
	matrix<double> count(ypixels,xpixels);
	double ratio_y,ratio_x;
	
	ratio_y=(double)ypixels/(double)(this->ypixels);
	ratio_x=(double)xpixels/(double)(this->xpixels);

	if(ypixels<=this->ypixels && xpixels<=this->xpixels){
		// x,y�����Ƃ���f������������Ƃ�
		for(i=1; i<=this->ypixels; ++i) for(j=1; j<=this->xpixels; ++j){
			ii=(int)((double)i*ratio_y);
			jj=(int)((double)j*ratio_x);
			if(ii<1)       ii=1;
			if(ii>ypixels) ii=ypixels;
			if(jj<1)       jj=1;
			if(jj>xpixels) jj=xpixels;
			AA[ii][jj]+=A[i][j];
			count[ii][jj]+=1;
		}
	}
	else{
		// x,y�����ŉ�f�������ƌ������������Ă���Ƃ��C
		// �������ł̕��Ϗ����͍s���Ȃ��D
		for(ii=1; ii<=ypixels; ++ii) for(jj=1; jj<=xpixels; ++jj){
			i=(int)((double)ii/ratio_y);
			j=(int)((double)jj/ratio_x);
			if(i<1)             i=1;
			if(i>this->ypixels) i=this->ypixels;
			if(j<1)             j=1;
			if(j>this->xpixels) j=this->xpixels;
			AA[ii][jj]=A[i][j];
			count[ii][jj]=1;
		}
	}

	SetYpixels(ypixels);
	SetXpixels(xpixels);

	for(ii=1; ii<=ypixels; ++ii) for(jj=1; jj<=xpixels; ++jj){
		A[ii][jj]=AA[ii][jj]/count[ii][jj];
	}
}

void cImage::ZeroFill(int ypixels,int xpixels){
	// ���͂ɍ���f��ǉ����C�S�̂� ypixels x xpixels �Ƃ���D
	cImage X;

	X.SetYpixels(ypixels); X.SetXpixels(xpixels);
	X.Merge(*this);
	*this=X;
}

void cImage::Trim(int ypixels,int xpixels){
	// �g���~���O���� ypixels x xpixels �ɂ���D
	// ���ʂ�ZeroFill(ypixels,xpixels)�Ɠ����D
	ZeroFill(ypixels,xpixels);
}

void cImage::TrimCircle(int ypixels,int xpixels){
	// �g���~���O���� ypixels x xpixels �ɂ���D
	// ����ɁC���a(ypixels,xpixels)�̑ȉ~�̊O�̉�f�����ɂ���D
	int i,j;
	double cy,cx, x,y;
	
	Trim(ypixels,xpixels);
	cy=(double)ypixels/2;
	cx=(double)xpixels/2;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		y=(double)i-cy;
		x=(double)j-cx;
		if( y*y/cy/cy + x*x/cx/cx >1) A[i][j]=0;
	}
}

complex cImage::GetComplexAmplitude(int i,int j) const{
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ) {
		return A[i][j];
	}
	else {
		return 0;
	}

}
void cImage::SetComplexAmplitude(int i,int j,complex value){
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ) {
		A[i][j]=value;
	}
}

double cImage::GetAmplitude(int i,int j) const{
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ) {
		return abs(A[i][j]);
	}
	else {
		return 0;
	}
}
void cImage::SetAmplitude(int i,int j,double value) {
	const complex I=complex(0,1);
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ) {
		A[i][j]=value*exp(I*arg(A[i][j]));
	}
}

double cImage::GetIntensity(int i,int j) const{
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ){
		return sqabs(A[i][j]);
	}
	else {
		return 0;
	}
}
void cImage::SetIntensity(int i,int j,double value){
	SetAmplitude(i,j, value>0 ? sqrt(value) : 0 );
}

void cImage::AddIntensity(int i,int j,double value){
	value+=GetIntensity(i,j);
	SetAmplitude(i,j, value>0 ? sqrt(value) : 0 );
}

double cImage::GetPhase(int i,int j) const{
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ) {
		return arg(A[i][j]);
	}
	else {
		return 0;
	}
}
void cImage::SetPhase(int i,int j,double value) {
	const complex I=complex(0,1);
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ) {
		A[i][j]=abs(A[i][j])*exp(I*value);
	}
}

int cImage::PixelLocation(int &i,int &j,double x,double y){
	// (x,y)�ɑΉ������f�ʒu(i,j)���Q�ƈ����ɑ������D
	// (x,y)���摜�͈͂ɂ���Ƃ���1�C�Ȃ����0��Ԃ��D
	double xwidth,ywidth;
	
	xwidth=xpitch*xpixels; 
	ywidth=ypitch*ypixels;
	if( x<-xwidth/2 || xwidth/2<x ) return 0;
	if( y<-ywidth/2 || ywidth/2<y ) return 0;
	x= x+xwidth/2; 
	y=-y+ywidth/2;  // y�͔��]����
	i=(int)ceil(y/ypitch);  
	j=(int)ceil(x/xpitch);
	return 1;
}

int cImage::Add(double x,double y,double weight){
	int i,j;
	
	if(PixelLocation(i,j,x,y)==0) return 0;
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ){
		SetIntensity(i,j,GetIntensity(i,j)+weight);
		return 1;
	}	
	else return 0;
}

int cImage::Add(double x,double y){
	return Add(x,y,1);
}

int cImage::AddComplexAmplitude(double x,double y,complex value){
	int i,j;
	
	if(PixelLocation(i,j,x,y)==0) return 0;
	if( 1<=i && i<=ypixels && 1<=j && j<=xpixels ){
		A[i][j]+=value;
		return 1;
	}	
	else return 0;
}

void cImage::SetLine(double x1,double y1,double x2,double y2,double weight){
	double xwidth,ywidth;
	int i1,j1,i2,j2,i,j, iwidth,jwidth;
		
	xwidth=ypitch*xpixels; 
	ywidth=xpitch*ypixels;

	x1= x1+xwidth/2; y1=-y1+ywidth/2;
	x2= x2+xwidth/2; y2=-y2+ywidth/2;
	
	i1=(int)ceil(y1/ypitch); j1=(int)ceil(x1/xpitch);
	i2=(int)ceil(y2/ypitch); j2=(int)ceil(x2/xpitch);

	iwidth=i2-i1>=0 ? i2-i1 : i1-i2;
	jwidth=j2-j1>=0 ? j2-j1 : j1-j2;

	if(iwidth==0 && jwidth==0) return;

	if(iwidth>=jwidth){
		if(i2>i1){
			for(i=i1; i<=i2; ++i){
				j=j1+(j2-j1)*(i-i1)/(i2-i1);
				SetIntensity(i,j,weight);
				// SetIntensity(i,j,GetIntensity(i,j)+weight) �Ƃ����
				// �����̐����������Ƃ���_�������邳������Ă��܂�
			}			
		}
		else{
			for(i=i1; i>=i2; --i){
				j=j1+(j2-j1)*(i-i1)/(i2-i1);
				SetIntensity(i,j,weight);
			}
		}
	}
	else{
		if(j2>j1){
			for(j=j1; j<=j2; ++j){
				i=i1+(i2-i1)*(j-j1)/(j2-j1);
				SetIntensity(i,j,weight);
			}			
		}
		else{
			for(j=j1; j>=j2; --j){
				i=i1+(i2-i1)*(j-j1)/(j2-j1);
				SetIntensity(i,j,weight);
			}
		}
	}
}

void cImage::SetCircle(double xo,double yo,double Dia,double Intensity,
			 double GaussDia/*=0*/,double Tilt_deg/*=0*/,double TiltAzm_deg/*=0*/,double Wl_nm/*=632.8*/){
	// �~��`��
	// xo,yo,Dia : �~�̒��S����ђ��a
	// Intensity : ���x�C�K�E�X���z�ł͒����̋��x
	// GaussDia  : �K�E�X���z�� 1/e^2 ���a�D�ψꕪ�z�ł�0�Ƃ���D
	// Tilt_deg  : ���z�����ʔg�̏Ǝ˂ɂ��Ƃ��C�X�N���[���ƕ��ʔg�i�s�������Ȃ��p�x
	// TiltAzm_deg : Tilt_deg ��0�łȂ��Ƃ��C���̌X���̕��ʊp�Dx���̐��̕�����0�Ƃ���D
	// Wl_nm       : ���ʔg�̔g���DTilt_deg ��ݒ肵���Ƃ��݈̂Ӗ�������D

	int i,j;
	double x,y,rr,I, th,phi,a;

	Clear();
	for(i=1; i<=xpixels; ++i) for(j=1; j<=ypixels; ++j){
		x= (i-xpixels/2)*xpitch;
		y=-(j-ypixels/2)*ypitch;
		rr=x*x+y*y;
		
		if(rr<=Dia*Dia/4){
			I=Intensity;
			if(GaussDia>0) I=I*exp(-2*rr/GaussDia/GaussDia*4);
			SetIntensity(i,j,I);

			if(Tilt_deg!=0){
				th=Tilt_deg*PI/180;
				phi=TiltAzm_deg*PI/180;

				a=cos(phi)*x+sin(phi)*y;   // ���ʊp�����̒P�ʃx�N�g��(cos��,sin��)��(x,y)�̓���
				SetPhase(i,j, 2*PI*a*th/(Wl_nm/1000000));
			}
		}
	} // next i,j
}

void cImage::MaskXCut(){
	// �摜�̋��x��X���ɉ�����0�ɂ���D
	// �t�[���G�ϊ��摜�ɑ΂��ēK�p���邱�ƂŁC�摜�̎΂ߕ����̐����폜����ȂǁD
	int i,j;
	double x,y, th;

	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		x=(double)j-(double)xpixels/2;
		y=-((double)i-(double)ypixels/2);  // (x,y)�͉摜���S�����_�C�E��𐳂Ƃ���
		th= x==0 ? PI/2 : atan(y/x);

		if     ( fabs(x)<2 && fabs(y)<2 );  // �����͎c��(�c���Ȃ��Ƌt�t�[���G�ϊ���ɔ��������]���邱�Ƃ�����D
		                                    //            <1 ���Ƃ��܂������Ȃ����Ƃ�����D)
		else if( (-3.0/8.0)*PI<th && th<(-1.0/8.0)*PI ) SetIntensity(i,j,0);
		else if( (1.0/8.0)*PI<th  && th<(3.0/8.0)*PI  ) SetIntensity(i,j,0);
	}
}

void cImage::Mask(double InnerR,double OuterR){
	// �~��Ƀ}�X�L���O����(�摜�̍����ƕ����قȂ�ꍇ�͑ȉ~��)
	// �~�̓����aInnerR,�O���aOuterR�͉摜�����C�����ꂼ��̔�����1�Ƃ���D
	// OuterR=0�̂Ƃ��́C�}�X�N���O�a�𖳌���Ƃ���D
	// DFT�Ƒg�ݍ��킹�ă��[�p�X�C�n�C�p�X�C�o���h�p�X�t�B���^����邱�Ƃ��ł���D
	int i,j;
	double cy,cx, y,x;
	
	cy=ypixels/2+1;  // �摜������1/2
	cx=xpixels/2+1;  // �摜����1/2
	// �����E������̂Ƃ��́Ccy,cx �͒����̉�f
	//     �V    �����̂Ƃ��́Ccy,cx �͒�����2��f�̉�f�ԍ��̑傫�����i�t�[���G�ϊ��摜�̒萔������f�ƈ�v�j
	
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		y=((double)i-cy)/cy;  // �摜���S����̋K�i������
		x=((double)j-cx)/cx;

		if(OuterR==0){
			if( y*y+x*x > InnerR*InnerR ) A[i][j]=0;
		}
		else{
			if( (y*y+x*x > InnerR*InnerR) && (y*y+x*x < OuterR*OuterR) ) A[i][j]=0;
		}
	}
}

void cImage::MaskXPass(){
	// �摜�̋��x���\���ɉ�����0�ɂ���D
	// �t�[���G�ϊ��摜�ɑ΂��ēK�p���邱�ƂŁC�摜�̏c�������̐����폜����ȂǁD
	int i,j;
	double x,y, th;

	for(i=1; i<=xpixels; ++i) for(j=1; j<=ypixels; ++j){
		x=(double)i-(double)xpixels/2;
		y=-((double)j-(double)ypixels/2);  // (x,y)�͉摜���S�����_�C�E��𐳂Ƃ���
		th= x==0 ? PI/2 : atan(y/x);
		
		if     ( fabs(x)<2 && fabs(y)<2 );  // �����͎c��(�c���Ȃ��Ƌt�t�[���G�ϊ���ɔ��������]���邱�Ƃ�����D
		                                    //            <1 ���Ƃ��܂������Ȃ����Ƃ�����D)
		else if( th<(-3.0/8.0)*PI ) SetIntensity(i,j,0);
		else if( (-1.0/8.0)*PI<th && th<(1.0/8.0)*PI ) SetIntensity(i,j,0);
		else if( (3.0/8.0)*PI<th ) SetIntensity(i,j,0);
	}
}

void cImage::NormalizeAmplitude(double new_max){
	// �U���̍ő傪new_max�ɂȂ�悤�ɂ���D
	double max=0;
	int i,j;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		if( GetAmplitude(i,j)>max ) max=GetAmplitude(i,j);
	}
	if( max>0 ){
		for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
			A[i][j]=A[i][j]*new_max/max;
		}
	}
}

void cImage::Normalize(double new_max){
	// ���x�̍ő傪new_max�ɂȂ�悤�ɂ���D
	double max;
	int i,j;
	max=MaxIntensity();
	if( max>0 ){
		for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
			A[i][j]=A[i][j]*sqrt(new_max/max);
		}
	}
}

void cImage::NormalizeTotalIntensity(double total/*=1*/){
	// ���x�̘a�� total �ɂȂ�悤�ɂ���D
	// (ypitch,xpitch�͍l�����Ȃ��D�P���Ɋe��f�̋��x�����v����D�j
	int i,j;
	double a;

	a=this->TotalIntensity()/total;
	if(a==0) return;
	a=sqrt(a);

	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		A[i][j]/=a;
	}
}

void cImage::GCenterToCenter(){
	// ���x���z�̏d�S���摜�̒��S(ypixels/2,xpixels/2)�Ƃ���D
	int i,j,i0,j0,ig,jg;
	cImage buf=*this;
	
	gravity_center(ig,jg);

	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		i0=i-ypixels/2+ig;   // i=ypixels/2 �̂Ƃ� i0=ig
		j0=j-xpixels/2+jg;
		if(1<=i0 && i0<=ypixels && 1<=j0 && j0<=xpixels){
			A[i][j]=buf.A[i0][j0];
		}
		else{
			A[i][j]=0;
		}
	}
}

void cImage::ToNeg(double new_max){
	// intensity�̕����𔽓]���CIntensity���0��new_max�ɂȂ�悤�ɕ��s�ړ�����D
	int i,j;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		SetIntensity(i,j,new_max-GetIntensity(i,j));
	}
}

void cImage::ReverseXY(){
	// �㉺���E�𔽓]����i�|�����𐳑��ɒ����Ȃǁj
	int i,j;
	cImage tmp=*this;

	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		A[i][j]=tmp.A[ypixels+1-i][xpixels+1-j];
	}
}

void cImage::Zoom(double my,double mx,double io,double jo){
	// �ʒu(io,jo)(������)�𒆐S�Ƃ��ďc����my�{,������mx�{�g�傷��D
	int i,j, i1,j1;
	double ii,jj, p;
	matrix<complex> a(ypixels,xpixels);
	complex z1,z2;
	if(mx==0 || my==0) return;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		ii=io+((double)i-io)/my;
		jj=jo+((double)j-jo)/mx;
		if(1<=ii && ii<=ypixels && 1<=jj && jj<=xpixels){
			i1=(int)ii; if(i1==ypixels) i1=ypixels-1;
			j1=(int)jj; if(j1==xpixels) j1=xpixels-1;
			p=ii-i1;
			z1=A[i1][j1]  *(1-p)+A[i1+1][j1]  *p;
			z2=A[i1][j1+1]*(1-p)+A[i1+1][j1+1]*p;
			p=jj-j1;
			a[i][j]=z1*(1-p)+z2*p;
		}
		else{
			a[i][j]=0;
		}
	}
	A=a;
}

void cImage::Zoom(double m,double io,double jo){
	Zoom(m,m,io,jo);
}

void cImage::Zoom(double m){
	// �摜�����𒆐S�Ƃ���m�{�g�傷��D
	Zoom(m,(double)ypixels/2,(double)xpixels/2);
}

double cImage::Transform(int OriginIsCenter,int inv/*=0*/){
	// �摜��ό`����i�{���F�����̕␳���j�D
	//  x0 = Cx(0,0) +Cx(1,0)x +Cx(0,1)y +Cx(2,0)x^2 +Cx(1,1)xy +Cx(0,2)y^2 + ...
	//  y0 = Cy(0,0) +Cy(1,0)x +Cy(0,1)y +Cy(2,0)x^2 +Cy(1,1)xy +Cy(0,2)y^2 + ...
	//  �����ŁC(x,y)�͕ό`��摜�̍��W�C(x0,y0)��(x,y)�ɑΉ�����ό`�O�摜�̍��W�D
	//  OriginIsCeter=true  �̂Ƃ��C�摜������(x,y)=(-1,-1) �E���(1,1)�Ƃ���D
	//  OriginIsCeter=fasle �̂Ƃ��C�摜�����(x,y)=(0,0)   �E����(1,1)�Ƃ���D

	// �ł������̖�����f�i�k�������ɂ��L���łȂ��Ȃ�����f�j�̈ʒu���摜���ɑ΂����ŕԂ��D

	// inv==true �̂Ƃ��͋t�ϊ�������. 

	int i,j, i0,j0, p,q, max_order;
	matrix<complex> a(ypixels,xpixels);
	double x,y, x0,y0;
	double ratio,min=1;
	double Cx1[MAX_ORDER+1][MAX_ORDER+1],Cy1[MAX_ORDER+1][MAX_ORDER+1], cx,cy;

	if(inv!=0){
		// �ŏ����@�ŋt�ϊ��W�� Cx1[p][q],Cy1[p][q] �����߂�
		cFitting Fx,Fy;

		for(p=0; p<=MAX_ORDER; ++p) for(q=0; q<=MAX_ORDER; ++q){ Cx1[p][q]=0; Cy1[p][q]=0; }

		max_order=0;
		for(p=0; p<=MAX_ORDER; ++p) for(q=0; q<=MAX_ORDER; ++q){
			// Cx,Cy�̒l�������Ă���ő原�� max_order �����߂�
			if( (Cx[p][q]!=0 || Cy[p][q]!=0) && p+q>max_order ) max_order=p+q;
		}

		Fx.SetOrder(max_order);  // �ŏ����@�� max_order���Ƃ���D MAX_ORDER���Ƃ���Ƌߎ��덷���傫���Ȃ邩������Ȃ��D
		Fy.SetOrder(max_order);

		for(i=0; i<=10; ++i) for(j=0; j<=10; ++j){ // �T���v���_���� 11x11=121
			if(OriginIsCenter!=0){
				x=-1+0.2*i; y=-1+0.2*j;  // x,y= -1�`+1 �X�e�b�v 0.2
			}
			else{
				x=i*0.1; y=j*0.1;        // x,y = 0�`+1 �X�e�b�v 0.1
			}

			x0=y0=0;
			for(p=0; p<=MAX_ORDER; ++p) for(q=0; q<=MAX_ORDER; ++q){
				x0+=Cx[p][q]*pw(x,p)*pw(y,q);
				y0+=Cy[p][q]*pw(x,p)*pw(y,q);
			}

			Fx.AddData(x0,y0,x);
			Fy.AddData(x0,y0,y);
		}

		Fx.CalcCoefficients();
		Fy.CalcCoefficients();

		for(p=0; p<=max_order; ++p) for(q=0; q<=max_order; ++q){
			Cx1[p][q]=Fx.coefficientGet(p,q,0);
			Cy1[p][q]=Fy.coefficientGet(p,q,0);
		}
	}

	for(j=1; j<=xpixels; ++j) for(i=1; i<=ypixels; ++i){  // (i,j)�͉摜���オ(1,1)
		if(OriginIsCenter){
			x=-1+2*(double)(j-1)/(double)(xpixels-1);  // j=1��x=-1, j=xpixels��x=1
			y= 1-2*(double)(i-1)/(double)(ypixels-1);  // i=1��y=1,  i=ypixels��y=-1
		}
		else{
			x=(double)(j-1)/(double)(xpixels-1);  // j=1��x=0, j=xpixels��x=1
			y=(double)(i-1)/(double)(ypixels-1);  // i=1��y=0, i=ypixels��y=1
		}

		x0=y0=0;
		for(p=0; p<=MAX_ORDER; ++p) for(q=0; q<=MAX_ORDER; ++q){
			cx= inv==0 ? Cx[p][q] : Cx1[p][q];
			cy= inv==0 ? Cy[p][q] : Cy1[p][q];
			x0+=cx*pw(x,p)*pw(y,q);
			y0+=cy*pw(x,p)*pw(y,q);
		}
		
		if(OriginIsCenter){
			j0=ToInt(1+(x0+1)*(double)(xpixels-1)/2);
			i0=ToInt(1-(y0-1)*(double)(ypixels-1)/2);
		}
		else{
			j0=ToInt(x0*xpixels);
			i0=ToInt(y0*ypixels);
		}
		
		if(1<=j0 && j0<=xpixels && 1<=i0 && i0<=ypixels){
			a[i][j]=A[i0][j0];
		}
		else{
			a[i][j]=0;

			// �����\������ōł������̖�����f�ʒu�𒲂ׂ�
			if(i==ypixels/2){     // ��f(i,j)�������̉�����ɂ���Ƃ�
				ratio=fabs((double)(j-xpixels/2)/((double)xpixels/2));
				if(ratio<min) min=ratio;
			}
			else if(j==xpixels/2){ // ��f(i,j)�������̏c����ɂ���Ƃ�
				ratio=fabs((double)(i-ypixels/2)/((double)ypixels/2));
				if(ratio<min) min=ratio;
			}
		}
	}
	A=a;
	return min;
}

double cImage::GetCx(int i,int j){
	if(i<=MAX_ORDER && j<=MAX_ORDER) return Cx[i][j]; else return 0;
}
void cImage::SetCx(int i,int j,double val){
	if(i<=MAX_ORDER && j<=MAX_ORDER) Cx[i][j]=val;
}

double cImage::GetCy(int i,int j){
	if(i<=MAX_ORDER && j<=MAX_ORDER) return Cy[i][j]; else return 0;
}
void cImage::SetCy(int i,int j,double val){
	if(i<=MAX_ORDER && j<=MAX_ORDER) Cy[i][j]=val;
}

int cImage::SetCxCy(std::string filename,std::string tagname){
	int i,j;
	std::string s,tag, tmp;
	// �ݒ�t�@�C���̗�
	// [BLUE]
	// A00=.01047
	// A10=.9787
	// A01=.0004607
	// ..
	// B00=.01287
	// B10=-.002801
	// B01=.9703
	// ..
	// [GREEN]
	// ..

	std::ifstream from(filename.c_str());

	if(from){
		for(i=0; i<=MAX_ORDER; ++i) for(j=0; j<=MAX_ORDER; ++j){ Cx[i][j]=0; Cy[i][j]=0; }
		do{
			from >> s;
			if(s[0]=='[') tag=trim(s,1);  // s="[BLUE]" �̂Ƃ��C tag="BLUE"
			if(tag==tagname){
				if(s[0]=='A' || s[0]=='a'){
					i=atoi((tmp=s[1]).c_str());          // s="A10=.9787" �̂Ƃ� i=1
					j=atoi((tmp=s[2]).c_str());          // s="A10=.9787" �̂Ƃ� j=0
					s=replace(s,"="," ");                // s="A10=.9787" �̂Ƃ� s="A10 .9787"
					if(0<=i && i<=MAX_ORDER && 0<=j && j<=MAX_ORDER){
						Cx[i][j]=atof(word(s,2,0).c_str());  // s="A10=.9787" �̂Ƃ� Cx[1][0]=.9787
					}
				}	
				else if(s[0]=='B' || s[0]=='b'){
					i=atoi((tmp=s[1]).c_str());
					j=atoi((tmp=s[2]).c_str());
					s=replace(s,"="," ");
					if(0<=i && i<=MAX_ORDER && 0<=j && MAX_ORDER){
						Cy[i][j]=atof(word(s,2,0).c_str());
					}
				}
			}
		} while(!from.eof());
		return 1;
	}
	else{
		return 0;
	}
}

void cImage::DistCorrect(double dist){
	// dist*100% ��3���c�Ȃ�␳����
	// ���Ȃ݂ɁCz=x^3+y^3 �͉�]��Ώ̂Ȃ̂ŁCCx[3][0],Cy[0][3]��ݒ肵��Transform()�����Ă��C
	// �␳�ł��Ȃ��D
	int i,j, ii,jj;
	double y,x,oy,ox,r,th;
	cImage image=*this;

	oy=(double)ypixels/2; ox=(double)xpixels/2;

	dist/=((double)ypixels*(double)ypixels+(double)xpixels*(double)xpixels)/4;  // 4���� dist*100% �̂���ƂȂ�

	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		y=double(i)-oy;
		x=double(j)-ox;

		r=sqrt(y*y+x*x);
		th=atan2(y,x);

		r+=r*r*r*dist;
		y=r*sin(th);
		x=r*cos(th);
		
		ii=(int)(y+oy);
		jj=(int)(x+ox);

		if(1<=ii && ii<=ypixels && 1<=jj && jj<=xpixels){
			this->A[i][j]=image.A[ii][jj];
		}
		else{
			this->A[i][j]=0;
		}
	}
}

cImage cImage::DFTAmplitude(int optical,double zoom) {
	// �U�����t�[���G�ϊ����CA[i][j]�Ƃ���D
	int i;
	for(i=1; i<=ypixels; ++i) ::DFT(A[i]+1,xpixels,0,optical,zoom);
	A=t(A);
	for(i=1; i<=xpixels; ++i) ::DFT(A[i]+1,ypixels,0,optical,zoom);
	A=t(A);
	return *this;
}

cImage cImage::InvDFTAmplitude(int optical,double zoom) {
	// �U�����t�t�[���G�ϊ����CA[i][j]�Ƃ���D
	int i;
	for(i=1; i<=ypixels; ++i) ::DFT(A[i]+1,xpixels,1,optical,zoom);
	A=t(A);
	for(i=1; i<=xpixels; ++i) ::DFT(A[i]+1,ypixels,1,optical,zoom);
	A=t(A);
	return *this;
}

cImage cImage::DFTIntensity(int optical,double zoom) {
	// ���x���t�[���G�ϊ����CA[i][j]�Ƃ���D
	int i;
	ToSquare();
	for(i=1; i<=ypixels; ++i) ::DFT(A[i]+1,xpixels,0,optical,zoom);
	A=t(A);
	for(i=1; i<=xpixels; ++i) ::DFT(A[i]+1,ypixels,0,optical,zoom);
	A=t(A);
	return *this;
}

cImage cImage::InvDFTIntensity(int optical,double zoom) {
	// ���x���t�t�[���G�ϊ����CA[i][j]�Ƃ���D
	int i;
	ToSquare();
	for(i=1; i<=ypixels; ++i) ::DFT(A[i]+1,xpixels,1,optical,zoom);
	A=t(A);
	for(i=1; i<=xpixels; ++i) ::DFT(A[i]+1,ypixels,1,optical,zoom);
	A=t(A);
	return *this;
}

double cImage::Diffraction(double width,double width1,double wl,double z0,double z,double f,int IsFraunhofer)
{
	// A��(������)���f�U�����ߗ����z�����ܑ��̕��f�U�����z�ɕϊ�����D
	//   width  : A(���U�����ߗ����z)�̈�ӂ̎�����
	//   width1 : A(��ܑ��U�����z)�̈�ӂ̎�����(0�̏ꍇ��FFT�̃f�t�H���g�l)
	//   wl : �g��(nm)
	//   z0 : ������_�����܂ł̋����D�����͌�����Ƃ���D
	//   z : ������ώ@�ʂ܂ł̋���
	//   f : ���ɖ������������Y�̏œ_����
	//   IsFraunhofer : 0�ȊO�̂Ƃ�Fraunhofer�ߎ��Ōv�Z����D���Ȃ킿�œ_����z�̃����Y�̏œ_�̉�ܑ����v�Z����D
	//                  ���̂Ƃ�z0,f�͖����ł���D
	//   �߂�l : A(��ܑ��U�����z)�̈�ӂ̎�����
	//
	// �����f�U�����ߗ����zP��z,z0,f�̌��ʂ���荞��P'�����DP'���t�[���G�ϊ���P'��Fraunhofer��ܑ��𓾂�D
	// �����P��Fresnel��ܑ��Ƃ�������D
	//

	const double PI=3.14159265358979;
	double width_im, c, arg;
	complex x;
	int n,i,j;
	double zoom;

	if( ypixels==xpixels ){
		n=ypixels;
	}
	else{
		return 0;
	}

	c=width/sqrt( (wl/1000000)*fabs(z)*n );
	width_im=(wl/1000000)*fabs(z)*n/width;

	for(i=1; i<=n; ++i)
	{
		for(j=1; j<=n; ++j)
		{
			if(IsFraunhofer)
			{
				x=1;
			}
			else
			{
				arg=-PI/n*c*c*((i-n/2.0)*(i-n/2.0)+(j-n/2.0)*(j-n/2.0))
				   *( 1-( z0==0 ? 0 : z/z0 )-(f==0 ? 0 : z/f) )*(z>=0 ? 1 : -1);
				x=complex(cos(arg),sin(arg));
			}
			A[i][j]=A[i][j]*x;
		}
	}

	if(width1==0){
		zoom=1;
	}
	else{
		zoom=width_im/width1;
		width_im=width1;
	}
	
	DFTAmplitude(1,zoom);

	SetWidth(width_im);
	return width_im;
}

double cImage::DiffractionCircleAperture
       (double diameter,int IsGauss,double width1,double wl,double z0,double z,double f,int IsFraunhofer){
	//   diameter : �J���̒��a�D�U�����ߗ��͈��Ƃ���D
	//   IsGauss  : true�̂Ƃ��U�����ߗ���diameter��1/e^2���a�Ƃ���K�E�X���z�̕�����
	//   width1 : A(��ܑ��U�����z)�̈�ӂ̎�����(0�̏ꍇ��FFT�̃f�t�H���g�l)
	//   wl : �g��(nm)
	//   z0 : ������_�����܂ł̋����D�����͌�����Ƃ���D
	//   z : ������ώ@�ʂ܂ł̋���
	//   f : ���ɖ������������Y�̏œ_����
	//   IsFraunhofer : 0�ȊO�̂Ƃ�Fraunhofer�ߎ��Ōv�Z����D���Ȃ킿�œ_����z�̃����Y�̏œ_�̉�ܑ����v�Z����D
	//                  ���̂Ƃ�z0,f�͖����ł���D
	//   �߂�l : A(��ܑ��U�����z)�̈�ӂ̎�����

	double width;
	int n,i,j;
	double ii,jj,nn;
	double rr,rr0;

	if( ypixels==xpixels ){
		n=ypixels;
	}
	else{
		return 0;
	}

	// ��ܑ��摜�ɑ΂���Fraunhofer��ܑ��̑傫�����J���摜�ɑ΂���J���̑傫���Ƃقړ����ƂȂ�悤
	// �J���摜�̈�ӂ̒���width��ݒ肷��D����̂��FWx,FWy�֐����g�����v�Z�̐��x��ǂ�����D
	width=sqrt((double)n)*diameter;

	for(i=1; i<=n; ++i) for(j=1; j<=n; ++j){
		nn=n;
		ii=i;
		jj=j;
		rr0=nn*nn*diameter*diameter/width/width/4;
		rr =(ii-nn/2)*(ii-nn/2)+(jj-nn/2)*(jj-nn/2);
		if(IsGauss){
			A[i][j]=exp( -rr/rr0 );  // A[i][j]�͐U���摜
		}
		else{
			if( rr <= rr0 ){
				A[i][j]=1;
			}
			else{
				A[i][j]=0;
			}
		}
	}

	return Diffraction(width,width1,wl,z0,z,f,IsFraunhofer);
}

double cImage::FWDiffractionCircleAperture
       (int size,double diameter,int IsGauss,double wl,double z0,double z,double f,int IsFraunhofer,double threshold){
	double width1;
	cImage x;
	x.SetYpixels(size);
	x.SetXpixels(size);
	width1=x.DiffractionCircleAperture(diameter,IsGauss,0,wl,z0,z,f,IsFraunhofer);
	return x.FWxIntensity(threshold,width1);
}

double cImage::FWxIntensity(double threshold,double width)
{
	// �d�S���܂ސ���ł�thereshold�ɑΉ�����S����Ԃ��D
	// 1�`n�ɑΉ����钷����width�ŗ^����D
	int ig,jg, j;
	cXYList li;

	gravity_center(ig,jg);
	for(j=1; j<=xpixels; ++j){
		li.AddData(j,GetIntensity(ig,j));
	}
	return width*li.FW(threshold)/(xpixels-1); // width��(xpixels-1)����
}

double cImage::FWyIntensity(double threshold,double width)
{
	// �d�S���܂ސ���ł�thereshold�ɑΉ�����S����Ԃ��D
	// 1�`n�ɑΉ����钷����width�ŗ^����D
	int ig,jg, i;
	cXYList li;

	gravity_center(ig,jg);
	for(i=1; i<=ypixels; ++i){
		li.AddData(i,GetIntensity(i,jg));
	}
	return width*li.FW(threshold)/(ypixels-1); // width��(ypixels-1)����
}

complex cImage::Total(){
	int i,j;
	complex s;
	s=0;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		s+=A[i][j];
	}
	return s;
}

double cImage::TotalAmplitude(){
	int i,j;
	double s;
	s=0;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		s+=GetAmplitude(i,j);
	}
	return s;
}

double cImage::TotalIntensity(){
	int i,j;
	double s;
	s=0;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		s+=GetIntensity(i,j);
	}
	return s;
}

double cImage::MaxIntensity(){
	// ���x�̍ő��Ԃ�
	double max=0;
	int i,j;

	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		if( GetIntensity(i,j)>max ) max=GetIntensity(i,j);
	}
	return max;
}

double cImage::GetBrightness(){
	// Intensity�̕��ς�Brightness�Ƃ���D
	int i,j,n;
	double z;

	n=0;
	z=0;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		z+=GetIntensity(i,j);
		n++;
	}
	return z/n;
}

void cImage::SetBrightness(double value){
	int i,j;
	double dz;

	dz=value-GetBrightness();
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		SetIntensity(i,j,GetIntensity(i,j)+dz);
		// GetIntensity(i,j)+dz<0 ���ƕs�t
	}
}

double cImage::GetContrast(){
	// A[i][j]��RMS���R���g���X�g�Ƃ���D
	int i,j,n;
	double z,ave;

	ave=GetBrightness();
	n=0;
	z=0;
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		z+=(GetIntensity(i,j)-ave)*(GetIntensity(i,j)-ave);
		n++;
	}
	return sqrt(z/n);
}

void cImage::SetContrast(double value){
	int i,j;
	double ave,m;

	ave=GetBrightness();
	m=value/GetContrast();
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		SetIntensity(i,j,ave+(GetIntensity(i,j)-ave)*m); 
		// ave+(GetIntensity(i,j)-ave)*m<0 ���ƕs�t
	}
}

void cImage::SetBrightnessContrast(double brightness,double contrast){
	// SetBrightness() �� SetContrast() �̌��ʂ�0�ȉ��̋��x��񂪗����Ă���̂ŁC
	// �����̘A�p�̍ۂ̏�񗎂����Ȃ�ׂ����Ȃ����邽�߂ɖ{�֐����쐬�����D
	int i,j;
	double m;

	m=contrast/GetContrast();
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		SetIntensity(i,j,brightness+(GetIntensity(i,j)-brightness)*m); 
		// brightness+(GetIntensity(i,j)-brightness)*m<0 ���ƕs�t
	}
}

void cImage::Intensify(double multiplier){
	int i,j;
	cImage a=*this;  // ���[�J���R�s�[�𑀍삵���ق����킸���ɑ����D
	
	if(multiplier>0){
		for(i=1; i<=a.ypixels; ++i) for(j=1; j<=a.xpixels; ++j){
			a.SetIntensity(i,j,a.GetIntensity(i,j)*multiplier);
		}
	}
	*this=a;
}

void cImage::BiasIntensity(double offset){
	int i,j;
	cImage a=*this;  // ���[�J���R�s�[�𑀍삵���ق����킸���ɑ����D
	
	for(i=1; i<=a.ypixels; ++i) for(j=1; j<=a.xpixels; ++j){
		a.SetIntensity(i,j,a.GetIntensity(i,j)+offset);  // offset�ɂ�蕉�ɂȂ�Ƃ���0����������(SetIntensity�̎����ɂ��)
	}
	*this=a;
}

void cImage::BiasAmplitude(double amp,double phase_deg){
	int i,j;
	cImage a=*this;  // ���[�J���R�s�[�𑀍삵���ق����킸���ɑ����D
	complex z;
	
	for(i=1; i<=a.ypixels; ++i) for(j=1; j<=a.xpixels; ++j){
		z=complex(amp*cos(phase_deg*PI/180),amp*sin(phase_deg*PI/180));
		a.SetComplexAmplitude(i,j,a.GetComplexAmplitude(i,j)+z);
	}
	*this=a;
}

void cImage::Merge(const cImage &X){
	// *this��X�𒆉����킹�ŋ��x�Ń}�[�W����D
	int i_offset,j_offset, i1,j1,i2,j2, i,j;
	cImage a=*this;
	
	// X�����͈� i1�`i2, j1�`j2 ��*this�̍��W�ŕ\���D
	i1= X.ypixels>=a.ypixels ? 1 : (a.ypixels-X.ypixels)/2;
	i2= X.ypixels>=a.ypixels ? a.ypixels : i1+X.ypixels-1;
	j1= X.xpixels>=a.xpixels ? 1 : (a.xpixels-X.xpixels)/2;
	j2= X.xpixels>=a.xpixels ? a.xpixels : i2+X.xpixels-1;

	i_offset=-(a.ypixels-X.ypixels)/2;
	j_offset=-(a.xpixels-X.xpixels)/2;

	for(i=i1; i<=i2; ++i) for(j=j1; j<=j2; ++j){
		a.SetIntensity(i,j, a.GetIntensity(i,j)+X.GetIntensity(i+i_offset,j+j_offset));
	}

	*this=a;
}

int cImage::MergeBmp(std::string filename){
	cImage X;

	if(X.OpenFromBmp(filename)){
		Merge(X);
		return 1;
	}
	else{
		return 0;
	}
}

cImage cImage::MedianFilter(){
	int i,j;
	double buf[10];
	cImage a=*this;

	for(i=2; i<=ypixels-1; i++) for(j=2; j<=xpixels-1; j++){  // ����1��C1�s�͏������Ȃ�
		buf[0]=GetIntensity(i-1,j-1);
		buf[1]=GetIntensity(i-1,j  );
		buf[2]=GetIntensity(i-1,j+1);
		buf[3]=GetIntensity(i  ,j-1);
		buf[4]=GetIntensity(i  ,j  );
		buf[5]=GetIntensity(i  ,j+1);
		buf[6]=GetIntensity(i+1,j-1);
		buf[7]=GetIntensity(i+1,j  );
		buf[8]=GetIntensity(i+1,j+1);
		a.SetIntensity(i,j,Median(buf,9));
	}
	return *this=a;
}

void cImage::LowPassFilter(double RelativeR/*=1*/){
	// �t�[���G�ϊ��摜�̑��Δ��a RelativeR ���O���̕���(�摜�̍����ƕ����قȂ�ꍇ�͑ȉ~�̊O���ƂȂ�)
	// ��0�Ƃ��邱�Ƃɂ�郍�[�p�X�t�B���^�D
	// RelativeR�͉摜�����C�����ꂼ��̔�����1�Ƃ���D
	int i,j;
	double cy,cx, y,x;

	DFTAmplitude(1,0);
	
	cy=(double)ypixels/2;  // �摜������1/2
	cx=(double)xpixels/2;  // �摜����1/2
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		y=(double)i-cy;  // �摜���S����̋���
		x=(double)j-cx;
		if( y*y/cy/cy + x*x/cx/cx > RelativeR*RelativeR ) A[i][j]=0;
	}

	InvDFTAmplitude(1,0);
}

void cImage::SetGauss(double dia){
	// ���x���z�𒼌adia�̃K�E�X���z�Ƃ���D�������͏�������D
	int i,j, m,n;
	double ry,rx,rr,r0;

	Clear();
	m=ypixels;
	n=xpixels;
	r0=dia/2;
	for(i=1; i<=m; ++i) for(j=1; j<=n; ++j){
		ry=(i-m/2)*ypitch;
		rx=(j-n/2)*xpitch;
		rr=rx*rx+ry*ry;
		SetIntensity(i,j,exp(-2*rr/r0/r0));
	}
}

cImage cImage::Edge(){
	int i,j;
	double fx,fy;
	cImage x=*this; x.Clear();  // x�̓T�C�Y����A�Ɠ����ŉ�f�l��0�ɂ�������

	for(i=2; i<=ypixels-1; ++i) for(j=2; j<=xpixels-1; ++j){ // �I�y���[�^�̐����ɂ�艏�͏���
		// Sobel�I�y���[�^����p������
		fx= -GetIntensity(i-1,j-1) +GetIntensity(i-1,j+1)
			-2*GetIntensity(i,j-1) +2*GetIntensity(i,j+1)
			-GetIntensity(i+1,j-1) +GetIntensity(i+1,j+1);
		fy= -GetIntensity(i-1,j-1) -2*GetIntensity(i-1,j) -GetIntensity(i-1,j+1)
			+GetIntensity(i+1,j-1) +2*GetIntensity(i+1,j) +GetIntensity(i+1,j+1);

		x.SetIntensity(i,j, sqrt(fx*fx+fy*fy));
	}
	return *this=x;;
}

cImage cImage::Binarize(double threshold){
	// 2�l������
	int i,j;
	
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		if(GetIntensity(i,j)<threshold){
			SetIntensity(i,j,0);
		}
		else{
			SetIntensity(i,j,255);  // cBitmap�ł̗��p��z������255�Ƃ���
		}
	}
	return *this;
}

int cImage::SaveIntensity(std::string filename) const{
	int i,j;
	matrix<double> a(ypixels,xpixels);
	
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		a[i][j]=GetIntensity(i,j);
	}
	return a.save(filename);
}

std::ostream& operator<<(std::ostream& to,const cImage& x){
	to<<x.ypixels<<' '<<x.xpixels<<std::endl;
	to<<x.ypitch<<' '<<x.xpitch<<std::endl;
	to<<x.A<<std::endl;
	return to;
}

std::istream& operator>>(std::istream& from,cImage& x){
	from>>x.ypixels; x.SetYpixels(x.ypixels);
	from>>x.xpixels; x.SetXpixels(x.xpixels);
	from>>x.ypitch>>x.xpitch;
	from>>x.A;
	return from;
}

int cImage::Open(std::string filename){
	std::ifstream from(filename.c_str());
	if(from) {
		from >> *this;
		return 1;
	}
	else return 0;
}

int cImage::Save(std::string filename) const{
	std::ofstream to(filename.c_str());
	if(to) {
		to << *this;
		return 1;
	}
	else return 0;
}

cImage& cImage::FromBitmap(const cBitmap &x,int ByAmplitude/*=0*/){
	// �Ⴆ�΁C
	// ByAmplitude==1�ł���΁C
	// DFTAmplitude()�̌�CSaveAsBmpByAmplitude()��DFT�摜(��Βl�̂�)�̕ۑ����ł���D
	// ���̑��ɂ��CA[][]��bmp�t�@�C�������ڑΉ����Ă���ƕ֗��Ȃ��Ƃ�����D
	cBitmap b;
	int i,j;

	SetYpixels(x.GetM());
	SetXpixels(x.GetN());
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		if(ByAmplitude==0){
			SetIntensity(i,j,x.GetG(i,j));  // �Ƃ肠����G���g���D���x�ɂ͉�f�l(0�`255)���ݒ肳���D
		}
		else{
			SetAmplitude(i,j,x.GetG(i,j));
		}
	}
	return *this;
}

int cImage::OpenFromBmp(std::string filename,int ByAmplitude){
	cBitmap b;

	if( b.Open(filename)==0 ) return 0;
	FromBitmap(b,ByAmplitude);
	return 1;
}

int cImage::OpenFromBmp(std::string filename){
	return OpenFromBmp(filename,0);
}

int cImage::OpenFromBmpByAmplitude(std::string filename){
	return OpenFromBmp(filename,1);
}

cBitmap& cImage::ToBitmap(cBitmap& bitmap,int ByAmplitude/*=0*/) const{
	// ToBitmap(int ByAmplitude) ������������ߖ�ł���i���[�J���� cBitmap x ���s�v�j
	int i,j;

	bitmap.SetM(ypixels);
	bitmap.SetN(xpixels);
	for(i=1; i<=ypixels; ++i) for(j=1; j<=xpixels; ++j){
		if(ByAmplitude==0){
			bitmap.SetGrayScale(i,j,(int)(GetIntensity(i,j)));
		}
		else{
			bitmap.SetGrayScale(i,j,(int)(GetAmplitude(i,j)));
		}
	}
	return bitmap;
}

cBitmap cImage::ToBitmap(int ByAmplitude/*=0*/) const{
	cBitmap x;
	return ToBitmap(x,ByAmplitude);
}

cBitmap cImage::ToBitmap(const cImage &r,const cImage &g,const cImage &b){
	int i,j;
	cBitmap x;

	x.SetM(r.ypixels);
	x.SetN(r.xpixels);
	for(i=1; i<=r.ypixels; ++i) for(j=1; j<=r.xpixels; ++j){
		x.SetRGB(i,j,(int)(r.GetIntensity(i,j)),(int)(g.GetIntensity(i,j)),(int)(b.GetIntensity(i,j)));
	}
	return x;
}

int cImage::SaveAsBmp(std::string filename,double gamma,int ByAmplitude) const{
	cBitmap bitmap;
	return ToBitmap(bitmap,ByAmplitude).Save(filename,gamma);
}

int cImage::SaveAsBmp(std::string filename,double gamma/*=1*/) const{
	return SaveAsBmp(filename,gamma,0);
}

int cImage::SaveAsBmpByAmplitude(std::string filename,double gamma/*=1*/) const{
	return SaveAsBmp(filename,gamma,1);
}

cImage operator+(const cImage& A,const cImage& B){
	int i,j;
	cImage X;
	if(A.ypixels==B.ypixels && A.xpixels==B.xpixels){
		X.SetYpixels(A.ypixels);
		X.SetXpixels(A.xpixels);
		for(i=1; i<=A.ypixels; ++i) for(j=1; j<=A.xpixels; ++j){
			X.A[i][j]=A.A[i][j]+B.A[i][j];
		}
		X.xpitch=A.xpitch; X.ypitch=A.ypitch;   // A�̃s�b�`���p������
		return X;
	}
	else{
		return cImage();
	}
}

cImage SumOnIntensity(const cImage& A,const cImage& B){
	int i,j;
	cImage X;
	if(A.ypixels==B.ypixels && A.xpixels==B.xpixels){
		X.SetYpixels(A.ypixels);
		X.SetXpixels(A.xpixels);
		for(i=1; i<=A.ypixels; ++i) for(j=1; j<=A.xpixels; ++j){
			X.A[i][j]=A.GetIntensity(i,j)+B.GetIntensity(i,j);
		}
		X.ToSqrt();
		X.xpitch=A.xpitch; X.ypitch=A.ypitch;    // A�̃s�b�`���p������
		return X;
	}
	else{
		return cImage();
	}	
}

cImage operator*(const cImage& A,const cImage& B){
	int i,j;
	cImage X;
	if(A.ypixels==B.ypixels && A.xpixels==B.xpixels){
		X.SetYpixels(A.ypixels);
		X.SetXpixels(A.xpixels);
		for(i=1; i<=A.ypixels; ++i) for(j=1; j<=A.xpixels; ++j){
			X.A[i][j]=A.A[i][j]*B.A[i][j];
		}
		X.xpitch=A.xpitch; X.ypitch=A.ypitch;    // A�̃s�b�`���p������
		return X;
	}
	else{
		return cImage();
	}
}

cImage operator*(double x,const cImage& A){
	int i,j;
	cImage X;
	X.SetYpixels(A.ypixels);
	X.SetXpixels(A.xpixels);
	for(i=1; i<=A.ypixels; ++i) for(j=1; j<=A.xpixels; ++j){
		X.A[i][j]=x*A.A[i][j];
	}
	X.xpitch=A.xpitch; X.ypitch=A.ypitch;
	return X;
}

cImage ProductOnIntensity(const cImage& A,const cImage& B){
	int i,j;
	cImage X;
	if(A.ypixels==B.ypixels && A.xpixels==B.xpixels){
		X.SetYpixels(A.ypixels);
		X.SetXpixels(A.xpixels);
		for(i=1; i<=A.ypixels; ++i) for(j=1; j<=A.xpixels; ++j){
			X.A[i][j]=A.GetIntensity(i,j)*B.GetIntensity(i,j);
		}
		X.ToSqrt();
		X.xpitch=A.xpitch; X.ypitch=A.ypitch;    // A�̃s�b�`���p������
		return X;
	}
	else{
		return cImage();
	}	
}

cImage Convolution(const cImage& A,const cImage& B){
	cImage a,b;
	a=A;
	b=B;
	return ( a.DFTAmplitude(1,1) * b.DFTAmplitude(1,1) ).InvDFTAmplitude(1,1);
}

cImage ConvolutionOnIntensity(const cImage& A,const cImage& B){
	cImage a,b;
	a=A;
	b=B;
	return ( a.DFTIntensity(1,1) * b.DFTIntensity(1,1) ).InvDFTAmplitude(1,1).ToSqrt();
}

double cImage::NCC(const cImage& A,const cImage& B,int i1/*=0*/,int j1/*=0*/,int i2/*=0*/,int j2/*=0*/){
	// (i1,j1)-(i2,j2)�͈̔͂�NCC�i���K�����ݑ��ցCnormalized cross correlation) ���v�Z����D
	// NCC��1�ɋ߂��قǉ摜�͑��ւ������D

	// NCC = ��{ (f-f_ave)(g-g_ave)/��((f-f_ave)^2)/��((g-g_ave)^2) }
	//     = { n(��fg) -(��f)(��g) } / ��{ n(��ff) -(��f)(��f) } / ��{ n(��gg) -(��g)(��g) }
	//         f,g �͊e��f�̋��x
	//         f_ave, g_ave �͏����͈͂ɂ����镽�ϒl
	//         ���͏����͈͂ł̘a
	//         n �͏����͈͓��̉�f��
	//    �� ��2���͉�f�Q�Ƃ��ꏄ�ōς܂��悤�ɑ�1����ό`��������
	int i,j;
	double SumA,SumB,SumAA,SumBB,SumAB, a,b,n;

	if(i1==0 && j1==0 && i2==0 && j2==0){
		// �f�t�H���g�̏����͈͉͂摜�S��
		i1=1; j1=1;
		i2=A.ypixels; j2=A.xpixels;
	}

	if(i2>Min(A.ypixels,B.ypixels)) i2=Min(A.ypixels,B.ypixels);  // �����͈͂��摜�O�ɋy�Ԃ��Ƃ�h��
	if(j2>Min(A.xpixels,B.xpixels)) j2=Min(A.xpixels,B.xpixels);

	SumA=SumB=SumAA=SumBB=SumAB=0;

	for(i=i1; i<=i2; ++i) for(j=j1; j<=j2; ++j){
		a=A.GetIntensity(i,j);
		b=B.GetIntensity(i,j);

		SumA+=a;
		SumB+=b;
		SumAA+=a*a;
		SumBB+=b*b;
		SumAB+=a*b;
	}

	n=(i2-i1+1)*(j2-j1+1);  // ��f��
	return (n*SumAB-SumA*SumB)/sqrt(n*SumAA-SumA*SumA)/sqrt(n*SumBB-SumB*SumB);
}

cImage cImage::POC(const cImage &A,const cImage &B){
	// �ʑ����葊�֖@�DA��B�̈ʑ����葊�։摜��Ԃ�
	// A��B�̃T�C�Y�͓����ł��邱�ƁD
	int i,j;
	cImage a=A;
	cImage b=B;
	cImage c;

	a.DFTAmplitude(0,1);
	b.DFTAmplitude(0,1);

	for(i=1; i<=a.ypixels; ++i) for(j=1; j<=a.xpixels; ++j){
		a.A[i][j]/=abs(a.A[i][j]);	
	}

	for(i=1; i<=b.ypixels; ++i) for(j=1; j<=b.xpixels; ++j){
		b.A[i][j]/=abs(b.A[i][j]);
		b.A[i][j]=conj(b.A[i][j]);
	}

	c=a*b;
	c.InvDFTAmplitude(0,1);
	return c;
}

int cImage::POC(std::string bmp_filename1,std::string bmp_filename2){
	// �����Ŏw�肷��2�̃r�b�g�}�b�v�摜�̈ʑ����葊�։摜�𐶐����C�r�b�g�}�b�v�Ƃ��ĕۑ�����
	cImage a,b,c;
	std::string fname;

	if(a.OpenFromBmpByAmplitude(bmp_filename1) && b.OpenFromBmpByAmplitude(bmp_filename2)){
		c=POC(a,b);
		fname="poc_"+date()+"_"+time()+".bmp";
		c.SaveAsBmpByAmplitude(fname,1);
		return 1;
	}
	else{
		return 0;
	}
}


void cImage::ToBuf(){
	buf=*this;
}

void cImage::FromBuf(){
	*this=buf;
}

void cImage::ToBufPlus(){
	buf=buf+(*this);
}

void cImage::ToBufPlusIntensity(){
	buf=SumOnIntensity(buf,*this);
}

void cImage::ToBufMultiply(){
	buf=buf*(*this);
}

void cImage::ToBufConvolute(){
	buf=Convolution(buf,*this);
}

void cImage::ToBufConvoluteIntensity(){
	buf=ConvolutionOnIntensity(buf,*this);
}
