#include "stdafx.h"
#include "cZernike.h"

void cZernike::ConvertR0(){
	double *a;
	int j;

	if(r0fix!=0 && r0!=r0fix){
		// �W�J�a�̕ύX
		//   �W�J�ar0fix�̑������ɉ��߂čŏ����t�B�b�e�B���O�����ɁC
		//   �ϊ������ɂ��W�J�ar0fix�̑������̌W���ɕϊ�����D
		//   �W�J���̕\���֐��͕ϊ��O���exact�Ɉ�v���邪�C
		//   r0fix�̑������ɉ��߂ăt�B�b�e�B���O�����Ƃ��ƈ�v����Ƃ͌���Ȃ��D
		//   ���̕ϊ�����]�����铙�Ɏg����D(10.04.30)
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
	// order���܂ł�Zernike�������̑�������Ԃ��D
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
	// (l,n) �̍����ʂ��ԍ��ŉ��Ԗڂ̍�����Ԃ��D
	int j;

	if(IsFringeOrder){
		//   Fringe order�̋K��
		//     �܂� n+m=n+|l| �̏�������
		//     ���� n�̏�������
		//     ���� m>0, m<0 �̏�
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
		//  ���͂��̕\�����Ȃ���l���o�����D
		//   ( (n+m)/2+1 )*( (n+m)/2+1 ) ��n+m���������e�O���[�v�̍Ō�̍���j�ł���D

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
	// Finge no. �� Standard no. �ɕϊ�����D
	//
	//   Fringe order�̋K��
	//     �܂� n+m=n+|l| �̏�������
	//     ���� n�̏�������
	//     ���� m>0, m<0 �̏�
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
	// �ʂ��ԍ�j�Ԗڂ̍���l��Ԃ��D
	if(jBase==0) j+=1;
	if(IsFringeOrder) j=StandardNo(j,1);
	return -nNumber(j,0,1) + ( j-TotalTerms(nNumber(j,0,1)-1,0)-1 )*2;
}

int cZernike::nNumber(int j,int IsFringeOrder,int jBase){
	// j�Ԗڂ̍���n(����)��Ԃ��D
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
		// �t�����W�I�[�_�[�ł� order=n+m
		return nNumber(j,IsFringeOrder,jBase)+mNumber(j,IsFringeOrder,jBase);
	}
	else{
		// �X�^���_�[�h�I�[�_�[�ł� order=n
		return nNumber(j,IsFringeOrder,jBase);
	}
}

void cZernike::ConvertR0(double *a,int terms,double r1,double r2,int Normalized,int IsFringeOrder){
	// �W�J���ar1�ɂ�����Zernike�W�J�W����W�J���ar2�ɂ����邻��ɕϊ�����D
	// ���� a[1],a[2], .... a[terms], terms=����
	// �ɂ��s�X�g������菇�ԂɌ��̌W����^���C�ϊ���̌W����a[]�ɏ㏑�������D
	// Normalized!=0 ���ς�ς̒P�ʉ~�����ςƂ���Ƃ����ȓ���(2�敽��)��1�ɂȂ鐳�K��
	//           ==0 �P�ʉ~���̍ő�ŏ����}1�ƂȂ鐳�K��
	// IsFringeOrder!=0 �t�����W�I�[�_�[
	//              ==0 �X�^���_�[�h�I�[�_�[
	//
	// �ϊ����͕����FScaling Zernike expansion coefficients to different pupil sizes
	//              (Jim Schwiegerling, Vol.19,No.10/October 2002/J.Opt.Soc.Am.A)
	// ��(A8)���ɂ��D
	// (A8)���͕��G�ł��邽�߁C�����ł́C(A8)������K�i���Ɋւ���v�f�i�����̕t�����v�f�j
	// �͎�菜��������p���C�K�i���Ɋւ��鏈���͕ϊ����K�p�O��ɍs�Ȃ��Ă���D
	// (�������邱�ƂŋK�i������Ă��Ȃ��W���ɂ��ȒP�ɑΉ��ł���D�j
	// �܂��C�g���w�̌����U(M.Born,E.Wolf)�h�̕\�L�Ƃ̐����̂��߁C
	//   �����ł�l=������m
	//   �����ł�m=������|m|
	//   �W��a,b�̓Y���͕����ł� ���a����,�p�x���g�� �̏������C�����ł� �p�x���g��,���a����
	// �Ƃ��Ă���D

	int n,nmax, l,lmin,lmax, m,i2,j,k;
	double *b;
	double sum1,sum2;

	b=new double [terms+1];

	// �K�i���Ɋւ��鏈��
	if(Normalized){
		for(j=1; j<=terms; j++){
			a[j]*=sqrt( 2*(nNumber(j,IsFringeOrder,1)+1)/(1+( lNumber(j,IsFringeOrder,1)==0 ? 1:0 )) );
		}
	}

	// l�̍ŏ��C�ő�����߂�
	lmin=lmax=0;
	for(j=1; j<=terms; j++){
		l=lNumber(j,IsFringeOrder,1);
		if(l<lmin) lmin=l;
		if(l>lmax) lmax=l;
	}

	// �ϊ�����K�p����
	for(l=lmin; l<=lmax; l++){
		m= l>=0 ? l : -l;

		// l�̗�̍ō�����nmax�����߂�
		nmax=0;
		for(j=1; j<=terms; j++){
			if( l==lNumber(j,IsFringeOrder,1) ){
				n=nNumber(j,IsFringeOrder,1);
				if(n>nmax) nmax=n;
			}
		}
		
		for(i2=nmax-m-(is_even(nmax-m) ? 0 :1); i2>=0; i2-=2){
			// �s���~�b�h�}��l�̗�̎����̍���������ϊ���̌W��b(l,n)=b(l,m+i2)�����߂Ă����D

			// (A8)��{}���̑�1��
			sum1=0;
			for(j=0; j<=(nmax-m-i2)/2; j++){
				sum1+=(is_even(j) ? 1 : -1)*a[jNumber(l,m+i2+j*2,IsFringeOrder,1)]
				      *factorial(m+i2+j)/factorial(j)/factorial(m+i2);
			}
			
			// (A8)��{}���̑�2��
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

	// �K�i���Ɋւ��鏈��
	if(Normalized){
		for(j=1; j<=terms; j++){
			b[j]/=sqrt( 2*(nNumber(j,IsFringeOrder,1)+1)/(1+( lNumber(j,IsFringeOrder,1)==0 ? 1:0 )) );
		}
	}

	// ���ʂ������ɏ㏑������
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
	// Zernike���a������
	// �Q�l��: �g���w�̌����U(M.Born,E.Wolf)�h 9.2.2 Zernike��circle polynomial
	// ��(5)��
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
	// Zernike���a��������rho�Δ���
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
	// Zernike���a��������rho2�K�Δ���
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
	// Zernike��������
	// �Q�l��: �g���w�̌����U(M.Born,E.Wolf)�h 9.2.2 Zernike��circle polynomial
	
	// n=0  l=0             piston
	// n=1  l=-1,1          tilt
	// n=2  l=-2,0,2        defocus,astigmatism
	// n=3  l=-3,-1,1,3     coma aberration etc.
	// n=4  l=-4,-2,0,2,4   spherical aberration etc.

	// Zernike�֐��̂悤�ȁCR(r)f(��)�̊֐��`�́C
	// R(r)��r=0��0�łȂ��Cf(��)�����łȂ����
	// xy���W�̌��_(r=0)�ŃƂɂ���Ċ֐��l���ω�����(�܂�xy���W�ł͒�`�ł��Ȃ�)�D
	// ���̂��Ƃ́Cxy���W�n�̌��_�͂��ƍ��W��r=0�����ɑΉ����C
	// �����W�n�͈�Έ�̑Ή��ɂȂ��Ă��Ȃ����ƂɊ֌W����D
	// �������C�S�Ă�Zernike�֐��͂����炭x=y=0�ŘA���ł���D
	// ���������āCx=y=0�ŃƂɂ���Ēl�͕ω����Ȃ��D
	// �����ł́Cx=y=0�Ń�=atan(0,0)��0 ���\�l�Ƃ��ė^���Ă���D
	// �������C�Ⴆ�ΕΔ���Ur�ɂ��Ă͖��炩�Ƀe�B���g���̂����xy���W�̌��_�ŕs�A���ł��邩��C
	// ���ӂ��Ȃ���΂Ȃ�Ȃ��D

	int m;
	double rho,th, R;
	
	if(r0==0) return 0;
	m= l>=0 ? l : -l;
	if( n>=0 && n>=m && is_even(n-m) ){
		rho=sqrt(x*x+y*y)/r0;
		th =atan2(y,x); 
		R=this->R(l,n,rho);
		// normalize!=0 ���ς�ς̒P�ʉ~�����ςƂ���Ƃ����ȓ���(2�敽��)��1�ɂȂ鐳�K��
		// normalize==0 �P�ʉ~���̍ő�ŏ����}1�ƂȂ鐳�K��
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
	// U��r�Δ���
	int m;
	double rho,th, Rr;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // ���_x=y=0�ł͑����֐��ł��邽�ߎg���Ȃ��DUr_pol()���g�����ƁD
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
	// U��r�Δ���
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
	// U�̃ƕΔ���
	int m;
	double rho,th, R;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // x=y=0�ł�0�ł���DUr()�ƈ���Ă���͐������D
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
	// U��r2�K�Δ���
	int m;
	double rho,th, Rrr;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // ���_x=y=0�ł͑����֐��ł��邽�ߎg���Ȃ��DUr_pol()���g�����ƁD
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
	// U��r2�K�Δ���
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
	// U�̃�2�K�Δ���
	int m;
	double rho,th, R;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // x=y=0�ł�0�ł���DUrr()�ƈ���Ă���͐������D
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
	// U��r,th�Δ���
	int m;
	double rho,th, Rr;

	if(r0==0) return 0;
	if(x==0 && y==0) return 0;  // ���_x=y=0�ł͑����֐��ł��邽�ߎg���Ȃ��DUrth_pol()���g�����ƁD
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
	// U��r,th�Δ���
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
	// U��x�Δ���
	double r;

	r=sqrt(x*x+y*y);
	if(r==0){
		// dr/dx,dth/dx�͌��_�Œ�`����Ȃ��D
		// �������CUx�͂����炭���_�Ŋ��炩�Ȋ֐��ł���D
		// �����ŁCdr/dx,dth/dx��x>0�ł̕Б��������Ƃ��āC
		// dr/dx=1, dth/dx=0�Ƃ���D
		// Ur(x,y)�����_�Œ�`����Ȃ��DUr_pol(r,th)���g���D
		return Ur_pol(l,n,0,0);
	}
	else{		
		// Ux=Ur(dr/dx)+Uth(dth/dx)
		// dr/dx=x/r
		// �܂��Ctan(th)=y/x ���C d(th)/cos(th)^2 = (-y/x^2)dx+(1/x)dy
		// ������Cdth/dx=-y/r^2
		// �y���z x=rcos(th)���dx=-rsin(th)dth -> dth/dx=-1/rsin(th)=-1/y �͊ԈႢ
		//        �Δ����ł�y���萔�ł���Cr���萔�ł͂Ȃ��D
		return Ur(l,n,x,y)*(x/r)+Uth(l,n,x,y)*(-y/r/r);
	}
}

double cZernike::Ux(int j,double x,double y) const{
	return Ux(lNumber(j),nNumber(j),x,y);
}

double cZernike::Uy(int l,int n,double x,double y) const{
	// U��y�Δ���
	double r;

	r=sqrt(x*x+y*y);
	if(r==0){
		// dr/dy,dth/dy�͌��_�Œ�`����Ȃ��D
		// �������CUy�͂����炭���_�Ŋ��炩�Ȋ֐��ł���D
		// �����ŁCdr/dy,dth/dy��y>0�ł̕Б��������Ƃ��āC
		// dr/dy=1, dth/dy=0�Ƃ���D
		// Ur(x,y)�����_�Œ�`����Ȃ��DUr_pol(r,th)���g���D
		return Ur_pol(l,n,0,PI/2);
	}
	else{
		// Uy=Ur(dr/dy)+Uth(dth/dy)
		// dr/dy=d(xx+yy)/dy=2y
		// �܂��Ctan(th)=y/x ���C d(th)/cos(th)^2 = (-y/x^2)dx+(1/x)dy
		// ������Cdth/dy=x/r^2
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
		// ���_�ł̃e�[���[�W�J��2���̍��́C
		//     z = (1/2)Uxx(dxdx)+Uxy(dxdy)+(1/2)Uyy(dydy)
		// �ƂȂ�Cdx=dy=1�̂Ƃ��C
		//     z(1,1) = (1/2)Uxx+Uxy+(1/2)Uyy
		// ����CUrr���g���ƁCdr^2=dx^2+dy^2=2 ������C
		//     z(1,1) = (1/2)Urr(��=PI/4)*2=Urr(��=PI/4)
		// ���������āC
		//     (1/2)Uxx+Uxy+(1/2)Uyy = Urr(��=PI/4)
		//     Uxy = Urr(��=PI/4)-(1/2)Uxx-(1/2)Uyy
		//         = Urr(��=PI/4)-(1/2)Urr(��=0)-(1/2)Urr(��=PI/2)
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
		X.datacopy(temp);   // X�̃f�[�^�͉\�Ȕ͈͂ňێ�����
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
	// this->x,y,z�ɌW��X[j][1]���t�B�b�e�B���O����Dj��jBase�ɂ�����炸1����terms�Ƃ���D

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
	// ��F com="Normalize 0; r0 1.5; C 4 0.1; C 12 0.1"

	std::string s,s0;
	int i,j, j_max;
	std::string s1,s2;
	bool b1,b2;

	// �ŏ���com�Ɋ܂܂��ő�̍��ԍ��𒲂ׂāC���������m�ۂ��Ă����K�v������D
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
	// z��Zernike�W���ŕ\�����`��Ƃ���
	int i;

	for(i=1; i<=NumberOfData(); ++i){
		this->z[i]=Z(GetXData(i),GetYData(i));
	}
	fitted=true;
}

double cZernike::Z(double x,double y,double Rref) const{
	// Rref = ����ʂ̋ȗ����a
	//        �Ⴆ�΁C���v��Zernike�W������΍��ł͂Ȃ���Ό`����v�Z����̂Ɏg���D
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
	// ����������`����ʂ�Axial�ȗ����aRa���v�Z����D
	// ��ʋȖʂ��ɍ��W�ɂ��C
	//     Z=(r,��)
	// �ƕ\���Ƃ��C
	//     Ra=r/[ sin{ atan(dZ/dr) } ]
	// �Ƃ���D
	double zr;

	zr=Zr(x,y,Rref,RemoveTilt);
	return zr==0 ? 0 : sqrt(x*x+y*y)/sin(atan(zr));
}

void cZernike::ParaxialExpansion(double &cx,double &cxy,double &cy){
	//  z�̋ߎ��W�J�C
	//      z = (cx/2)*x^2 +cxy*x*y +(cy/2)*y^2
	//  �̌W�������߂�D
	//  �Q�l���F�g���w�̌����U(W.Born E.Wolf ����O��)�h
	
	// �����ł̓e�[���[�W�J�����Ƃ���2���̌W�������߂Ă���D
	// �ł��t�B�b�g�����i��������������Ƃ��̎c�]rms���ŏ���)
	// 2���ȖʂƂ͈Ⴄ���̂ł���D
	// ���������āC��̎������ł��⏞����ዾ��S,C,A�����߂�C
	// �Ƃ������ł���Αf����Zernike�W�J��2���̌W��������
	// ���������悢�D

	int s,m,l, j,n;
	double b;

	if(r0==0) return;
	cx=cy=cxy=0;
	
	for(j=1; j<=terms; ++j){
		l=lNumber(j);
		n=nNumber(j);
		m= l>=0 ? l : -l;

		for(s=0; s<=(n-m)/2; ++s){       // 9.2.1 (5)��
			if( n-2*s==2 ){              // n-2*s : r�̎���
				b=(is_even(s) ? 1: -1)
				  *factorial(n-s)/factorial(s)/factorial((n+m)/2-s)/factorial((n-m)/2-s);
				if(normalize) b*= m==0 ? sqrt((double)(n+1)) : sqrt((double)(n+1)*2);
				switch(l){
				case 0:  // r^2=x^2+y^2, C*b*(x^2+y^2) = (cx/2)*x^2+(cy/2)*y^2
					cx+=b*2*GetC(j);
					cy+=b*2*GetC(j);
					break;
				case -2: // r^2*sin(2��)=2*x*y, C*b*2*x*y = cxy*x*y
					cxy+=b*2*GetC(j);
					break;
				case 2:  // r^2*cos(2��)=x^2-y^2, C*b*(x^2-y^2)=(cx/2)*x^2+(cy/2)*y^2
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
	// �t�B�b�e�B���O��Cthis->z���piston,tilt,sph���̕�����������D
	//     piston<>0  piston����������
	//     tilt  <>0  tilt����������
	//     sph   <>0  sph����������
	// �W�J����Orders��1���ȉ��̂Ƃ��͉������Ȃ��D
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
	// �t�B�b�e�B���O��C
	// piston,tilt,sph�����𒲐����āCthis->z�𕪎U��x,y���z�̈���ōŏ��ɂȂ�悤�ɂ���D
	//     tilt=ture�Ȃ��piston���C
	//     sph=true�Ȃ��piston,tilt��
	// ��������D
	// �Ⴆ�΁Cx,y�̕��z�̈悪�~���炸��Ă���(�ȉ~�Ȃ�)�Ƃ��C
	// Zernike�������͂��̗̈�ł͒������Ȃ��D
	// ���������āCRemoveTerms()�֐��ŁC
	//     z����piston���̕����������Ă�z�̕��z�������Ώ̂ɂȂ�Ȃ��D
	//     z����sph���̕����������Ă��p����rms��������D
	// ���Ƃ����������ߖ{�֐����쐬�����D (080910)
	//
	// �y���Ӂz
	//    z���g�ʎ����Cx,y�������W�ł���Ƃ��C
	//    ���ʂƎ�����������łȂ��Ƃ�(�Ⴆ��FT�̃E�G�n�ʂ𓵂Ƃ����Ƃ�),
	//    sph�����͌��̐i�s�������猩��Δ�_�������܂ށD
	//    ���������āCsph�̒����ɂ���_��������������C
	//    ���������ۂ��ǂ������邱�Ƃ�����D

	int i,j, m,n;
	matrix<double> A1,X1,F1;
	double x,y,z;

	if(Fit() && terms>=4){
		// �f�[�^��m, ����n�̐ݒ�D�s��p�������̊m�ہD
		m=this->x.GetSize();	
		n=0;
		if(piston) n=1;
		if(tilt)   n=3;
		if(sph)    n=4;
		if(n==0) return 0;
		A1.redim(m,n);
		F1.redim(m,1);
		// �ŏ����@�̕Δ����W���s��̐ݒ�
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
		// �ŏ����@�̖ڕW�l��z=0;
		for(i=1; i<=m; ++i){
			F1[i][1]=-(this->z[i]);
		}
		// �ŏ����@�̎��s
		if( rank(t(A1)*A1)<n ) return 0; // �񐳑��Ȃ�Ή������Ȃ��ŏI��
		X1=inv( t(A1)*A1 )*t(A1)*F1;
		// �W�����C������D(X[]�̓Y���̓s�X�g����1�ł��邱�Ƃɒ���)
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
		// z��V�����W�����g�����W�J���̒l�ɂ���D
		for(i=1; i<=m; ++i){
			x=this->x[i];
			y=this->y[i];
			// �����ŁCz=zApproximate(x,y) �Ƃ����̂��l�����邪�C
			// ���������RMSError(),PVError()��0�ɂȂ��Ă��܂��D�i���Ȃ킿�����̎�����������D�j
			// �t�B�b�e�B���O�G���[�ʂ�ۑ����邽�ߎ��̂悤�ɂ���D
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
	// tilt !=0       �̂Ƃ���tilt��(n==1)��������
	// sph  !=0       �̂Ƃ���sph��(n==2;l==0)��������
	// cyl  !=0       �̂Ƃ���cyl��(n==2;l==-2,2)��������
	// high_order !=0 �̂Ƃ��͍�����(n>=3)��������
	// (���𐫂��g���Ă��Ȃ��̂ŁCx,y�̕��z�̈悪�~�łȂ��Ă��L��)
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
	// *this�̂��f�[�^��dx,dy�������炵���ꏊ��z�f�[�^�Ƃ̍����Ƃ�C*this�֑������D
	// x,y�f�[�^��*this�Ɠ����ɂȂ�D
    // ���炵��z�f�[�^�ɂ�Zernike�W�J��p����D���̓W�J��pv�G���[��Ԃ��D
	// �W�J�̃I�[�_�[���w�肷�邽�߁C�R�s�[���g���Čv�Z����D
	// z���g�ʎ����ł����Shearing���v�̃V�~�����[�V�����ƂȂ�D
	cZernike a;
	int i;
	double pv;

	a=*this;
	a.SetMaxOrder(order);
	pv=a.PVError();

	for(i=1; i<=a.NumberOfData(); ++i){
		this->z[i]=a.z[i]-a.ZApproximate(a.x[i]+dx,a.y[i]+dy);
	}
	
	AdjustTerms(1,1,0);   // piston,tilt�����͏�������
	this->fitted=false;
	return pv;
}




