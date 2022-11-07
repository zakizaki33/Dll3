#include "stdafx.h"
#include "cSelfoc.h"

double cSelfoc::N(vector<double> r){
	return N(no,rA,r.x,r.y);					
}

vector<double> cSelfoc::gradN(vector<double> r){
	return vector<double>(-no*rA*rA*r.x, -no*rA*rA*r.y, 0);
}

vector<double> cSelfoc::NgradN(vector<double> r){
	return N(r)*gradN(r);
}

double cSelfoc::N(double no,double rA,double x,double y){
	return no*(1-rA*rA/2*(x*x+y*y));		
}

cSelfoc::cSelfoc(){}

cSelfoc::cSelfoc(double no,double rA,double phi){
	Set_no(no);
	Set_rA(rA);
	Set_phi(phi);
}

cSelfoc::~cSelfoc(){}

double cSelfoc::Get_no(){ return no; }
void   cSelfoc::Set_no(double value){ no=value; }

double cSelfoc::Get_rA(){ return rA; }
void   cSelfoc::Set_rA(double value){ rA=value; }

double cSelfoc::Get_phi(){ return phi; }
void   cSelfoc::Set_phi(double value){ phi=value; }

double cSelfoc::Get_L(){ return L; }
void   cSelfoc::Set_L(double value){
	L=value;
	Pn=L/(2*PI/rA);
}

double cSelfoc::Get_Pn(){ return Pn; }
void   cSelfoc::Set_Pn(double value){
	Pn=value;
	L=(2*PI/rA)*Pn;
}

int cSelfoc::trace(vector<double> &r,vector<double> &Q,double ds){
	//   r(x,y,z) : �ʒu���W
	//   Q(X,Y,Z) : ���������]��
	// ��^���C�����ɉ�����r����ds�������ꂽ�_�ɂ�����,r,Q���v�Z����D

	// ����������  (d/ds){n(r)(d/ds)r(s)} = grad(n(r)) �ɂ��D
	// ����ɒu�� dt=ds/n ���g���ƁC(d/dt)(d/dt)r = n*grad(n)
	// �ƂȂ�D����ɏ����lr,dr/dt��^�������Q�N�b�^�@�ŉ����D

	vector<double> v;
	double dt;

	v=Q/abs(Q)*N(r);   // v=dr/dt=N(dr/ds)=NQ �̏����l
	dt=ds/N(r);

	{
		// 2�K����������� y''(x) =f(x,y(x),y'(x)) �������Q�N�b�^�@�ŉ����D
		// �����ł́C y''(x) = (d/dt)(d/dt)r
		//                x  = t
		//              y(x) = r
		//             y'(x) = dr/dt = v
		//           f(x,y(x),y'(x))=NgradN(r) (���Ȃ킿��2�����݂̂Ɉˑ�)
		// �Q�l���F�gFORTRAN�ɂ�鉉�K���l�v�Z �F�V�����j��(�T�C�G���X��)�h 9.4��
		double dx;
		vector<double> y,y1;
		vector<double> k11,k21,k12,k22,k13,k23,k14,k24;

		y=r;
		y1=v;
		dx=dt;
		
		k11=y1;
		k21=NgradN(y);
		
		k12=y1+(dx/2)*k21;
		k22=NgradN(y+(dx/2)*k11);

		k13=y1+(dx/2)*k22;
		k23=NgradN(y+(dx/2)*k12);

		k14=y1+dx*k23;
		k24=NgradN(y+dx*k13);

		y +=(dx/6)*(k11+2.0*k12+2.0*k13+k14);
		y1+=(dx/6)*(k21+2.0*k22+2.0*k23+k24);

		r=y;
		v=y1;
	}
	
	Q=v/abs(v);

	return r.x*r.x+r.y*r.y<=phi*phi/4 ? 1 : 0;  // �߂�l�Ń����Y�a�𒴂������ǂ�����Ԃ��D
}

int cSelfoc::Trace(double &x,double &y,double &Qx,double &Qy,double &Qz,double N,double N1,double ds){
	// ���˖ʂɂ����� x,y,Q(���ˑ��}����)�ɂ������ǐՂ����C�o�˖ʂ�x,y,Q�ŏ㏑������D
	//   N  : ���ˑ��}�����ܗ�
	//   N1 : �o�ˑ��}�����ܗ�
	//   ds : �����Q�N�b�^�@�̃X�e�b�v(�����ɉ�����������)
	double xi,xi1,n,G;
	vector<double> r,r0,Q;
	
	Q=vector<double>(Qx,Qy,Qz);
	xi=Q.z/abs(Q);                                        // (3.20) "���� �����Y�݌v�@"
	r=vector<double>(x,y,0);
	n=this->N(r);
	xi1=N*n*xi/fabs(N*n*xi)*sqrt(1-(N*N/n/n)*(1-xi*xi));  // (3.21)
	G=xi1-fabs(N/n)*xi;                                   // (3.22)
	Q=fabs(N/n)*Q+vector<double>(0,0,G);                  // (3.23)

	do{
		r0=r;
		trace(r,Q,ds);
		if (r.x*r.x+r.y*r.y>phi*phi/4) return 0;
	} while( fabs(r.z)<L );

	r=r0+(r-r0)*(L-fabs(r0.z))/fabs(r.z-r0.z);
	x=r.x;
	y=r.y;

	xi=Q.z/abs(Q);                                            // (3.20)
	r=vector<double>(r.x,r.y,0);
	n=this->N(r);
	xi1=n*N1*xi/fabs(n*N1*xi)*sqrt(1-(n*n/N1/N1)*(1-xi*xi));  // (3.21)
	G=xi1-fabs(n/N1)*xi;                                      // (3.22)
	Q=fabs(n/N1)*Q+vector<double>(0,0,G);                     // (3.23)
	Qx=Q.x; Qy=Q.y; Qz=Q.z;

	return 1;
}

