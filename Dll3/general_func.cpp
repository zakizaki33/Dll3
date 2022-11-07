#include "stdafx.h"
#include "general_func.h"

int CreateConsole() {
	// �v���Z�X�ɃR���\�[�������蓖�āC�W���o�͂�������
	FreeConsole();
	if(AllocConsole()){   // �y64bit���ł͎��s���ɂȂ����K�����s����D32bit�ł͖��Ȃ��z
		freopen("CONOUT$","w",stdout); /* �W���o��(stdout)��V�����R���\�[���Ɍ����� */
		freopen("CONOUT$","w",stderr); /* �W���G���[�o��(stderr)��V�����R���\�[���Ɍ����� */
		return 1; // ����
	}
	else{
		return 0; // ���s
	}
}

std::string RemoveExtension(std::string filename){
	// filename���g���q���������������Ԃ��D
	return remove_extension(filename);
}

std::string Date(){
	std::string s;
	char buf[100];
	time_t now;
	struct tm *pnow;

	now=time(NULL);
	pnow=localtime(&now);
	sprintf(buf, "%4d/%02d/%02d", pnow->tm_year+1900, pnow->tm_mon+1, pnow->tm_mday);  // 1����0
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
	// start����̌o�ߎ��Ԃ�b�P�ʂŕԂ�
	char buf[100];

	sprintf(buf,"%.3fsec",(double)(clock()-start)/CLOCKS_PER_SEC);
	return buf;
}

void dwrite(std::ostream to,double x){ // �o�C�i�����[�h double�^�����o��
	to.write((char*)&x, sizeof(double));
}

double dread(std::istream from){ // �o�C�i�����[�h double�^�ǂݍ���
	double x;

	from.read((char*)&x, sizeof(double));
	return x;
}

void iwrite(std::ostream to,int x){ // �o�C�i�����[�h int�^�����o��
	to.write((char*)&x, sizeof(int));
}

int iread(std::istream from){ // �o�C�i�����[�h int�^�ǂݍ���
	int x;

	from.read((char*)&x, sizeof(int));
	return x;
}

void swrite(std::ostream to,std::string s){ // �o�C�i�����[�h std::string�^�����o��
	char *buf;
	int n;

	n=(int)strlen(s.c_str());
	buf=new char[n+1];
	strcpy(buf,s.c_str());
	to.write((char*)&n, sizeof(int));
	to.write((char*)buf, n+1);
	delete [] buf;
}

std::string sread(std::istream from){ // �o�C�i�����[�h std::string�^�ǂݍ���
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
		double a;  // "int a" �Ƃ���ƗႦ�� 15! �̓I�[�o�[�t���[����
		int j;
		a=1;
		for(j=1; j<=i; j++) a*=(double)j;
		return a;
	}
}

double pw(const double& x,const int& n){
	// x��n��Cn<0�ł͐������Ȃ�
	int i;
	double p=1;

	for(i=1; i<=n; i++) p*=x;
	return p;
}

double fceil(const double& x){
	// x�ȏ�ōŏ���
	//    (�����̈ꌅ�̐���) * (10�̊K��)
	// �ł��鐔�l��Ԃ��D
	double xabs= x>=0 ? x : -x;
	return x==0 ? 0 : ceil( x/pow(10,floor(log10(xabs)+1e-15))-1e-15 ) * pow(10,floor(log10(xabs)+1e-15));
}

double ffloor(const double& x){
	// x�ȉ��ōő��
	//     (�����̈ꌅ�̐���) * (10�̊K��)
	// �ł��鐔�l��Ԃ��D
	double xabs= x>=0 ? x : -x;
	return x==0 ? 0 : floor( x/pow(10,floor(log10(xabs)+1e-15))+1e-15 ) * pow(10,floor(log10(xabs)+1e-15));
}

double fceil_mantissa(const double& x){
	// x�ȏ�ōŏ���
	//     a * (10�̊K��), a=1,2,3,..,9,10
	// �ł��鐔�l�����߁Ca��Ԃ��D
	// �O���t�̎��̕����������߂�ꍇ���Ɏg�p����D
	double xabs= x>=0 ? x : -x;
	return x==0 ? 0 : ceil( x/pow(10,floor(log10(xabs)+1e-15))-1e-15 );
}

double sgn(const double& x) { return x>=0 ? 1 : -1; };

double Round(const double& x,const int& n){
	// n<=0 : �����_�ȉ�|n|���Ɋۂ߂�
	// n>=1 : �L������n���Ɋۂ߂�
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
	// �����_�ȉ�1�ʂ��l�̌ܓ�����������Ԃ�
	return x>=0 ? (int)(x+0.5) : (int)(x-0.5);
	// �idouble����int�ւ̕ϊ��́C�����_�ȉ���؎̂Ă�0�̕����֌�����
	//   ���Ƃ�JIS�K�i�ŋK�肳��Ă���D�j
}

double ArcCos(const double& x){
	// cos(x)���acos�֐���p���ē��ˊpx�����߂悤�Ƃ���Ƃ��C
	// ���l�v�Z�덷�ɂ��x���킸���ł�1�𒴂��Ă���ƒ�`��G���[�ɂȂ��Ă��܂��D
	// �{�֐��ł�x��1�𒴂��Ă���Ƃ��́Cx�̗L���������炵�Ă���acos�����s���邱�Ƃ�
	// �G���[�̔�����}����D
	if(x<=1){
		return acos(x);
	}
	else{
		return acos(Round(x,8));  // x�̗L������8���Ɋۂ߂�
	}
}

void InverseMatrix(double a[4][4]){
	// 3x3�s��̋t�s����v�Z����
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
	::srand( (unsigned)time(NULL) );   // 1�b���Ƃɗ����̕��т��ς��
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
	// |x|<K�� �͈̗̔͂������𔭐�����D
	// x�͐��K���z y=exp{-x^2/(2��^2)} �̊m���Ŕ�������D
	// ���̊m���́C���a�Ђ̃K�E�X�r�[���̋��x���z�ɓ������D
	// |x|<�� �ł���m���� 68.26%, <2�Ђ̊m����95.44%(���[�U�̃r�[�����a), <3�Ђł�99.74%
	const double K=4;  
	double x,y;

	if( sigma==0 ) return 0;
	do{
		x= -K*sigma +2*K*sigma*(double)rand()/(double)RAND_MAX;
		y= exp( -2*x*x/sigma/sigma );
	} while ( (double)rand()/(double)RAND_MAX > y );
	return x;

	// �K�E�X�r�[���̎��ƍ������Ȃ����ƁD
	// �K�E�X�r�[���̃p���[���z�� y = exp(-2w/wo) �ƂȂ�i���q��2���t���j�D
	// ���a2wo�̒��ɓ���p���[�͑S�̂�86.5%�Ƃ����Ă���C
	// ����͏�L�g<2�Ђ̊m��95.44%�h�̓��C91%�����Ⴂ��
	// 91%�͈�ӂ�2wo�� "�����`" �Ɋ܂܂��p���[�ł��邽�߂ł���D
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
	// ��F��Ԃ��D
	// �Q�l�F http://appakumaturi.hatenablog.com/entry/20120121/1327143125
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
	// buf[0]����buf[n-1]�̒����l��Ԃ��D
	Sort(buf,n);
	return buf[n/2];
}

int QuadraticEq(double& x1,double& x2,double a,double b,double c){
	// 2�������� ax^2+bx+c=0 ������
	if( b*b-4*a*c >= 0 ){
		// -b��sprt(...)�Ƃ̉����Z�ł̌��������ӂ�������
		//    x1=(-b+sqrt(b*b-4*a*c))/2/a;
		//    x2=(-b-sqrt(b*b-4*a*c))/2/a;
		// �Ƃ͂��Ȃ��D
		// �Q�l�F�gFORTRAN ��{�{���p(�����O)�hp108
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
	// *a���܂�a���炎�̃f�[�^(a[0]�`a[n-1])�ɑ΂����U�t�[���G�ϊ����s���C���ʂ��㏑������D
	//     n�͔C�ӁDn��2�̊K��̂Ƃ���FFT�Ōv�Z����D
	//     inv==0�̂Ƃ����ϊ��Cinv!=0�̂Ƃ��t�ϊ��D
	//     optical==0�̂Ƃ�a[0]�����������Coptical!=0�̂Ƃ����������������D
	//     zoom��0,1�ȊO�̂Ƃ��C���g�����ő�l��1/zoom�{�ɂȂ�D���̂Ƃ�FFT�͎g���Ȃ��D
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
			// ���U�t�[���G�ϊ���`�������̂܂܃v���O��������
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
					b[i]=b[i]+a[j]*complex( cos(w),sin(w) );  // b[i]+=a[i]*exp{I*(�}2*PI*i1*j1/n)}
				}
			}

			for(i=0; i<n; i++) a[i]=b[i];
			delete[] b;
		}
*/
		{
			// ���U�t�[���G�ϊ�
			// ��`���ǂ���Ƀv���O��������̂ł͂Ȃ��C�ǐ��͗����邪�C
			// �Efor���[�v����o������̂͏o���D
			// �E�O�p�֐��̉��@�藝�𗘗p����, sin(),cos()�̌Ăяo�������炷�D
			// ���Ƃō�������}��D
			
			// ���ϊ�         a[j] = ��(j=0�`n-1) a[i]exp( -(2��I/n)ij )
			// �t�ϊ�         a[j] = ��(j=0�`n-1) a[i]exp(  (2��I/n)ij )
			// ���ϊ�(���w�I) a[j] = ��(j=0�`n-1) a[i]exp( -(2��I/n)(i-n/2)(j-n/2) )
			// �t�ϊ�(���w�I) a[j] = ��(j=0�`n-1) a[i]exp(  (2��I/n)(i-n/2)(j-n/2) )
			int i,j, p;
			double w;
			double dcs,dsn,cs,sn, dcs0,dsn0,cs0,sn0, cs_temp,sn_temp;
			complex *b=new complex[n+1];
			
			w= inv ? 2.0*PI/n : -2.0*PI/n;
			w/=(zoom==0 ? 1 : zoom);
			
			p= optical ? n/2 : 0;

			// w*(-1-p)*(-1-p) �́Ci=j=-1�̂Ƃ��̈ʑ�
			cs0=cos(w*(-1-p)*(-1-p));
			sn0=sin(w*(-1-p)*(-1-p));
			// -w*(-1-p)�́Cj=-1�̂Ƃ��̈ʑ�(����)��i�ɑ΂������
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
				// w*(i-p) �́C�ʑ�(����)��j�ɑ΂������
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
	// a[0][0]����a[m-1][n-1]�܂ł�m x n�̃f�[�^�̊e�s���t�[���G�ϊ�����D
	int i;

	for(i=0; i<=m-1; ++i){
		DFT(a[i],n,inv,optical,zoom);
	}
}

void DFTColumn(complex **a,int m,int n,int inv,int optical,double zoom){
	// a[0][0]����a[m-1][n-1]�܂ł�m x n�̃f�[�^�̊e����t�[���G�ϊ�����D
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
	// a[0][0]����a[m-1][n-1]�܂ł�m x n�̃f�[�^���t�[���G�ϊ�����D
	DFTRow(a,m,n,inv,optical,zoom);
	DFTColumn(a,m,n,inv,optical,zoom);
}

void HilbertT(complex *a,int n){
	// *a���܂�a���炎�̃f�[�^(a[0]�`a[n-1])�ɑ΂��q���x���g�ϊ����s���C���ʂ��㏑������D
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
	// *a���܂�a���炎�̃f�[�^(a[0]�`a[n-1])�ɑ΂����U�␳���s���C���ʂ��㏑������D
	// ������ʑ��ʃ����́C
	//     ���� = 2��{(A2)x^2 +(A3)x~3} ,  x=i/n (i=0�`n-1)
	// �Ƃ���D
	const double PI=3.141592654;
	const complex I=complex(0,1);

	int i;
	complex *tmp,*c;
	double x,phi;

	if(A2==0 && A3==0) return;

	tmp=new complex[n];
	c=new complex[n];

	for(i=0; i<n; ++i) tmp[i]=c[i]=a[i];

	DFT(c,n,0,0,0);   // �t�[���G���ϊ�
	c[0]=0;           // �������������� (c�����֐��Ȃ̂ŁC�t�[���G�ϊ��̒�`���Cc[0]�͎����j
	DFT(c,n,1,0,0);   // �t�[���G�t�ϊ� (c������������������������ƂɂȂ邩��C�t�ϊ����c�͎��֐�)

	HilbertT(tmp,n);  // �q���x���g�ϊ��i���������������܂܂��j

	for(i=0; i<n; ++i){
		Im(c[i])=Re(tmp[i]);         // �������Ƀq���x���g�ϊ����ʂ���
		x=((double)i-(double)(n-1)/2)/((double)(n-1)/2);  // i=0,n-1��x=�}1, i=n/2��x��0
		phi=2*PI*(A2*x*x+A3*x*x*x);
		c[i]*=exp(I*phi);            // �ʑ���������
		a[i]=Re(c[i]);
	}

	delete [] tmp;
	delete [] c;
}

void EnFace(int m,int n,std::string in_filename,std::string out_filename,double gamma/*=1*/){
	// m = en face�摜�̏c��f�� = B�X�L�����摜��
	// n = en face�摜�̉���f�� = B�X�L�����摜�̉���f��
	// in_filename  = ���̓t�@�C����(�g���q��������bmp�t�@�C����) = ��A��B�X�L�����摜�̃t�@�C���� 
	//               (ex. "oct_a_h#" = "oct_a_h001" �` "oct_a_h256", "#"�������Œu����������)
	// out_filename = �o�̓t�@�C����(�g���q��������bmp�t�@�C����) = en face�摜�̃t�@�C����
	// gamma = �o�̓t�@�C���ۑ����̃K���}�␳�l
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
		num=buf;                                             // ��Fm=256�̂Ƃ��Cnum="001","002",..,"256"
		if( x.Open( replace(in_filename,"#",num)+".bmp" ) ){ // B�X�L�����摜���J��
			m_depth=x.GetM();                                // B�X�L�����[��������f��
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
		A[i][j]*=(255/max);  // �ő�l��255�ƂȂ�悤�ɋK�i��
		x.SetRGB(i,j,(int)A[i][j],(int)A[i][j],(int)A[i][j]);	
	}
	x.Save(out_filename+".bmp",gamma);  // en face�摜��ۑ�����
}

double SphSag(double r,double h){
	// ���ʂ̃T�O�ʂ��v�Z����
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
	// ���aD��2�~�̒��S�ԋ�����a�̂Ƃ��̉~�̏d�Ȃ蕔���̉~�̖ʐςɑ΂������v�Z����D
	// �������~�`�J����MTF�̌v�Z�ȂǂɎg���D
	// �Ⴆ�� NA/0.61��=1.6393NA/�� �ɑΉ�����MTF�́Ccut off=2NA/�� ������C
	//   OverlapAreaCircles(2,1.6393) = 0.0894
	// �ƂȂ�D
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
	dx=x2-x1; dy=y2-y1; dz=z2-z1;    // (x1,y1,z1) ���猩�� (x2,y2,z2) �̈ʒu
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
	// 2�̉~  (x-X1)^2 +(y-Y1)^2 =R1^2  ... (1)
	//          (x-X2)^2 +(y-Y2)^2 =R2^2  ... (2)
	// �̌�_ (x1,y1),(x2,y2) �����߂�D

	// �{�֐��͌��ʕ�����̏o�͂͂��Ȃ��i sprintf(...) �̏����Ɏ��Ԃ������邽�߁j�D
	// ������o�͂ɂ�
	// std::string IntersectionsCircles(double X1,double Y1,double R1,double X2,double Y2,double R2);
	// ���g���D
	double A,B, a,b,c;
	
	if(Y1==Y2){ // 2�~�̒��S�������ɕ���ł���ꍇ
		// (1)(2)��� x=-B/A;
		A= 2.0*(X1-X2);
		B= X2*X2-X1*X1 +Y2*Y2-Y1*Y1 -R2*R2+R1*R1;
		// �����(1)�ɑ�����C2���������̉��̌�����y�����߂�D
		a= 1;
		b= -2.0*Y1;
		c= Y1*Y1 +(-B/A-X1)*(-B/A-X1) -R1*R1;
		if(b*b-4*a*c){
			x1=y1=x2=y2=0;   // ���Ȃ�
		}
		else{
			x1=-B/A; y1=(-b+sqrt(b*b-4*a*c))/2/a;
			x2=-B/A; y2=(-b-sqrt(b*b-4*a*c))/2/a;
		}
	}
	else{
		// 2�~�̌�_�����Ԓ����� y=Ax+B  ((1)(2)���� x^2, y^2 ���������C�n���Ɍv�Z����j
		A= -(X1-X2)/(Y1-Y2);
		B= ( X2*X2-X1*X1 +Y2*Y2-Y1*Y1 -R2*R2+R1*R1 )/2/(Y2-Y1);
		// �����(1)�ɑ�����C2���������̉��̌�����x�����߂�D
		a= A*A+1;
		b= -2*X1+2*A*B-2*A*Y1;
		c= X1*X1 +B*B -2*B*Y1 +Y1*Y1 -R1*R1;
		if(b*b-4*a*c <0){
			x1=y1=x2=y2=0;   // ���Ȃ�
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
	// (x1,y1),(x2,y2),(x3,y3)��ʂ�~�̒��S(xc,yc)����є��ar���v�Z����D�܂��Cr��߂�l�Ƃ���D
	// ���@�F
	//   �~�̒��S��(x1,y1)��(x2,y2)�̒��_��ʂ�C����(x1,y1)-(x2,y2)�ɐ����Ȓ���L1�ƁC
	//   (x1,y1)��(x3,y3)�ɂ��Ă̓��l�Ȓ���L2�̌�_�ł���D
	//   L1�̕������́C((x1+x2)/2,(y1+y2)/2) ��ʂ�C(x1,y1)-(x2,y2)�ɐ���(���ς�0)�ł��邱�Ƃ���C
	//     (x2-x1)x +(y2-y1)y = (y2^2-y1^2)/2 + (x2^2-x1^2)/2
    //   �ƂȂ�DL2�ɂ��Ă����l�D

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
	// (x1,y1,z1),(x2,y2,z2),(x3,y3,z3)�𒸓_�Ƃ���O�p�`�̖ʐ�
	vector<double> a,b;
	a=vector<double>(x2-x1,y2-y1,z2-z1);
	b=vector<double>(x3-x1,y3-y1,z3-z1);
	return abs(vProduct(a,b))/2;
}

double QuadrangleArea(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4){
	// �ix1,y1),(x2,y2),(x3,y3),(x4,y4)�����Ɍ��񂾎l�p�`�̖ʐ�
	return TriangleArea(x1,y1,0,x2,y2,0,x4,y4,0) + TriangleArea(x2,y2,0,x3,y3,0,x4,y4,0);
}

double PolygonArea(int n,double *x,double *y){
	// n�̓_ (x[1],y[1]),(x[2],y[2]),  ... (x[n],y[n]) �� "���̏��Ԃ�"
	// �����ɂČ��񂾑��p�`�̖ʐς����߂�D
	// ���_��(xi,yi)��(xi+1,yi+1) �ɂ��O�p�`�̖ʐς��O�ςɂ��v�Z���C
	// i�ɂ��Ęa���Ƃ邱�Ƃɂ��D
	// �O�ς͕����t���ł��邩��]���ȕ����͑��E����K�v�Ȗʐς����߂���D
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
	// n�̓_ (x[1],y[1]),(x[2],y[2]),  ... (x[n],y[n]) �ɂ�鑽�p�`�Ɠ_(X,Y)�̋��������߂�D
	//  �E�_�Ƒ��p�`�̊e�ӂ܂ł̋����̍ŏ��l��Ԃ�
	//  �E�_�����p�`�̓����ɂ���Ƃ��͐��C�O���ɂ���Ƃ��͕��̒l��Ԃ�
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
		if(d<min) min=d;    // �ł��߂��ӂ܂ł̋���
	}

	count=0;

	for(i=1; i<=n; ++i){
		if(IsCrossLines(X,Y, X+LN,Y, x[i],y[i], x[i==n ? 1:i+1],y[i==n ? 1:i+1])) count++;
	}

	return count%2==0 ? -min : min;   // Crossing Number Algorithm�i��_���������Ȃ�Γ_�͑��p�`�̊O�ɂ���j
}

void ParaxialR(double &r1,double &r2,double &axis_deg,double axx,double ayy,double axy){
	// �Ȗ�z(x,y)�̓W�J��2���̍����C
	//   z=axx*x^2+ayy*y^2+axy*x*y
	// �ł���Ƃ��C��ȗ����ar1,r2��r1�̕�����x���ƂȂ��paxis_deg�����߂�D
	// �Q�l�F ����O �g�����Y���w�h	(3.4.1�`3.4.3��)
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
	// �_P1,P2��ʂ钼����œ_P0�ɍł��߂��_�i�����̑��j��Ԃ��D
	// (P2-P1) �� P1+s(P2-P1)-Po �͐����ł��邱�Ƃ��C
	//   (P2-P1)�E(P1+s(P2-P1)-Po) =0
	//   s = {(P2-P1)�E(P0-P1)} / {(P2-P1)�E(P2-P1)}
	double s;
	
	s=sProduct(P2-P1,P0-P1)/sProduct(P2-P1,P2-P1);
	return P1+s*(P2-P1);
}

vector<double> NearestPointLineLine(vector<double> P,vector<double> V,vector<double> Po,vector<double> Vo){
	// ���� Po+sVo �ɍł��߂����� P+tV ��̓_��Ԃ�.
	// �����ŁCPo,P�͒�����̓_�CVo,V�͕����x�N�g��(�P�ʃx�N�g���łȂ��Ă��悢)�D
	// �ŋߐړ_�� Po+tVo, P+sV �Ƃ���ƁC�ŋߐړ_���m�����Ԓ����͌���2�����Ɛ���������C
	//    (P-Po+sV-tVo)�EVo =0
	//    (P-Po+sV-tVo)�EV  =0
	// �ƂȂ�D����炩�� t ���������Cs �ɂ��ĉ����΋��߂�_�� P+sV �ƂȂ�D
	double s;

	s= ( sProduct(V,Po-P)+sProduct(Vo,P-Po)*sProduct(Vo,V)/sProduct(Vo,Vo) )
	  /( sProduct(V,V)-sProduct(V,Vo)/sProduct(Vo,Vo) );
	return P+s*V;
}

double DistancePointLinesegment(vector<double> P0,vector<double> P1,vector<double> P2){
	// �_P0�Ɛ���(P1,P2)�̋���
	vector<double> p=NearestPointOnLine(P0,P1,P2);  // P0���璼��(P1,P2)�ɉ����������̑�

	if( sProduct(p-P1,p-P2)<0 ){  // �����̑��͐�����ɂ���
		return abs(P0-p);
	}
	else{
		return Min( abs(P0-P1), abs(P0-P2) );
	}
}

bool IsCrossLines(double x1,double y1,double x2,double y2,double x3,double y3,double x4,double y4){
	// ����(x1,y1)-(x2,y2) �� ����(x3,y3)-(x4,y4) ������邩�ǂ����𔻒肷��D
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
		a=sProduct(vProduct(v3-v1,v2-v1),vProduct(v4-v1,v2-v1)); // ����v1-v2�Ɛ���v3-v4������邩(v3,v4������v1-v2�̈Ⴄ���ɂ��邩�j
		b=sProduct(vProduct(v1-v3,v4-v3),vProduct(v2-v3,v4-v3)); // ����v1-v2�ƒ���v3-v4������邩
		return a<=0 && b<=0;
	}
}

matrix<double> Tmatrix(double rox,double roy,double roz){
	// ��ԃx�N�g����������]rox,roy,roz���s�������W�n�ɕϊ�����s��i�g�����Y���w�h3.3.1 )
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
	// N�g�̃f�[�^ x[0](==*x),x[1], ... x[N-1]
	//             y[0](==*x),y[1], ... y[N-1]
	// ��3���X�v���C����Ԃɂ��Cxx�ɑΉ�����y�̒l�����߂�D
	// SPDerivativeZero : �n�_�̈ꎟ���֐���0�Ƃ��邩�ǂ����i�U�̂Ƃ��͓񎟓��֐���0�j
	// EPDerivativeZero : �I�_�̈ꎟ���֐���0�Ƃ��邩�ǂ����i�U�̂Ƃ��͓񎟓��֐���0�j
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
	// N�g�̃f�[�^ x[0](==*x),x[1], ... x[N-1]
	//             y[0](==*x),y[1], ... y[N-1]
	// �̃��O�����W����Ԃɂ��Cxx�ɑΉ�����y�̒l�����߂�D
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
	// p0����_�Ƃ��ČŒ肵�Cppre,psuc��ύX���Ĉʑ��ڑ�����D
	// p0��ppre�̊ԁC�܂���p0��psuc�̊Ԃ̏��Ȃ��Ƃ��Е��ɂ͈ʑ���т��Ȃ��Ƃ���D
	// �����̒P�ʂ�deg�Ƃ���D
	double dp1,dp2;
	const double K=180;

	dp1=p0-ppre;
	dp2=psuc-p0;

	if(dp1*dp2<0){                  // ���Ȃ��Ƃ�p0����[�ł���Ƃ���
		if(fabs(dp1)>fabs(dp2)){    // �ω����傫�������ʑ���т��N�����Ă���Ƃ���
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
	// N�̃f�[�^ th_deg[0]....th_deg[N-1]���ʑ��ڑ����㏑������Dth_deg[0]���Œ肷��D
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
	// (S1,C1,A1)��(S2,C2,A2)�������Z����D
	// ���Z�����Z����plus_minus�� '+' �܂��� '-' �Ƃ��Ďw�肷��D	        
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
        C=-C;       // �}�C�i�X�\���ɂ���
        A+=PI/2;
    }

    if(A<0){
        A+=PI;      // A=0�`�� �ɂ���
    }
    
    S=(S1+C1/2)+sgn*(S2+C2/2)-C/2;
    A_deg=A*180/PI;
}
/*
void AddSCA(double& S,double& C,double& A,
            double S1,double C1,double A1,char plus_minus,double S2,double C2,double A2){
	// (S1,C1,A1)��(S2,C2,A2)�������Z����D
	// ���Z�����Z����plus_minus�� '+' �܂��� '-' �Ƃ��Ďw�肷��D
	// ����C�̓}�C�i�X�\���Ƃ��邪�C����C1,C2�̃v���X�}�C�i�X�͕s��D
	//
	// �~��������傫����CYL, ���ʊp��AXIS*2 �̃x�N�g���ŕ\�����Ƃ��C
	// �~�������̉����́C�x�N�g���̉����Z�ƂȂ邱�Ƃɂ��D
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
	C*=-1; A-=PI/2;                  // C��-�\���ɂ���
	S=(S1+C1/2)+sign*(S2+C2/2)-C/2;  // �������ʓx���͕ۑ������̂ŁC(S1+C1/2)�}(S2+C2/2)=S+C/2 ���D

	A*=180/PI;
	if(A<0) A+=180;   // A��0�`180deg�Ƃ���D
}
*/
double dBToT(double dB){
	// �f�V�x���ŕ\���ꂽ�����ʂ𓧉ߗ��Ɋ��Z����
	return pow(10,-dB/10);
}

int CSVReadLine(std::ifstream& from,double* data,int datas,int& valid_line){
	// csv�t�@�C�� from ����1�s��ǂݍ��݁C
	// data[0],data[1], ... ,data[datas-1] �Ƀf�[�^��������
	int i;
	std::string s;

	if(std::getline(from,s)){
		s=replace(s,",",";");   // csv�t�@�C���̋�؂� "," ���Cword()�ł̋�؂� ";" �ɒu��������
		valid_line=1; for(i=1; i<=datas; ++i) valid_line = valid_line && is_numeric(word(s,i,0));
		if(valid_line){  // �L���ȍs�ł���idatas�̐��l�f�[�^������j
			for(i=1; i<=datas; ++i){
				data[i-1]=atof(word(s,i,0).c_str());
			}
		}
		return 1;  // �s��ǂݍ��߂�
	}
	else{
		return 0;  // �s��ǂݍ��߂Ȃ������i�t�@�C���I�[�𒴂��Ă���j
	}
}

std::string scmd_general_func(const std::string& com,int val){
	// val=ture�̂Ƃ��C�ꕔ�̃R�}���h��Basic���ɂ�val�֐��ŏ������邽�߂�
	// ���l��\��������̂ݕԂ��D
	//
	// ���ӁF����%f�͏����_�ȉ��̌������Œ�Ȃ̂ŁC�l���������ƗL�������Ȃ��Ȃ�D
	//       %.15g���悢�D

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



