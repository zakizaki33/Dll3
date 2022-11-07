#include "stdafx.h"
#include "cLeneq.h"

cLeneq::cLeneq() {
	int i;

	NumberOfEq=2;
	NumberOfVar=2;
	A.redim(NumberOfEq,NumberOfVar);
	W.redim(NumberOfEq,1); for(i=1; i<=NumberOfEq; ++i) W[i][1]=1;
	C.redim(NumberOfEq,1); for(i=1; i<=NumberOfEq; ++i) C[i][1]=DLS;
	F.redim(NumberOfEq,1);
	X.redim(NumberOfVar,1);
	calced=false;
}

cLeneq::~cLeneq() {
}

void cLeneq::SetNumberOfEq(int m) {
	int i;

	if(m>0){
		NumberOfEq=m;
		A.redim(NumberOfEq,NumberOfVar);
		W.redim(NumberOfEq,1); for(i=1; i<=NumberOfEq; ++i) W[i][1]=1;
		C.redim(NumberOfEq,1); for(i=1; i<=NumberOfEq; ++i) C[i][1]=DLS;
		F.redim(NumberOfEq,1);
		calced=false;
	}
}
int cLeneq::GetNumberOfEq(){
	return NumberOfEq;
}

void cLeneq::SetNumberOfVar(int n){
	int i;

	if(n>0){
		NumberOfVar=n;
		A.redim(NumberOfEq,NumberOfVar);
		W.redim(NumberOfEq,1); for(i=1; i<=NumberOfEq; ++i) W[i][1]=1;
		C.redim(NumberOfEq,1); for(i=1; i<=NumberOfEq; ++i) C[i][1]=DLS;
		X.redim(NumberOfVar,1);
		calced=false;
	}
}
int cLeneq::GetNumberOfVar(){
	return NumberOfVar;
}

void cLeneq::SetA(int i,int j,double value){
	if( 1<=i && i<=NumberOfEq && 1<=j && j<=NumberOfVar){
		A.a[i][j]=value;
		calced=false;
	}
}
double cLeneq::GetA(int i,int j) const{
	return A[i][j];
}

matrix<double> cLeneq::GetAMatrix() const{
	return A;
}

void cLeneq::SetB(int i,double value){
	if(1<=i && i<=NumberOfEq){
		F[i][1]=value;
		calced=false;
	}
}
double cLeneq::GetB(int i){
	return F[i][1];
}

void cLeneq::SetWeight(int i,double value){
	if(1<=i && i<=NumberOfEq){
		W[i][1]=value;
		calced=false;
	}
}
double cLeneq::GetWeight(int i){
	if(1<=i && i<=NumberOfEq){
		if(C[i][1]==DLS) return W[i][1]; else return 0;
		// W���C��D�悷��D
		// �������Ă����΁C�Ⴆ��DLS�@�Ń����b�g�֐����v�Z����Ƃ��C�R���X�g���C���c�ɐݒ肵��
		// �]���֐���Wi ���[���ɂȂ邽�߁CDLS�̊֐����R���X�g���C���c�̊֐����ꍇ�������Ȃ��čς�.
	}
	else{
		return 0;
	}
}

void cLeneq::SetConstrKind(int i,int value){
	if(1<=i && i<=NumberOfEq){
		C[i][1]=value;
		calced=false;
	}
}
int cLeneq::GetConstrKind(int i){
	if(1<=i && i<=NumberOfEq){
		return C[i][1];
	}
	else{
		return 0;
	}
}

double cLeneq::GetDampedX(int j,double rho,int GradientDescent/*=0*/){
	// GradientDescent���^�̂Ƃ���DLS�@�̑���ɍŋ}�~���@��p����
	const int MAXTIMES=100;    // �L������@�̌J��Ԃ��񐔏��
	int ii,jj,kk, m,ma,mb1,mb2,mb2_, n;
	bool e1;
	matrix<double> A,B1,B2,B2_, W, Fa,Fb1,Fb2,Fb2_, A1,Fa1, AA,FF,XX, BX;
	matrix<int> C;
	// Ci=DLS && Wi>0 : �ŏ����@�ɂ��œK��
	// Ci=EQ          : �ꎟ�����S��(Lagrange����搔�@�jAx=F
	// Ci=LT,GT       : �ꎟ�s�����S��(�L������@�j    Ax<F(Ci=LT), Ax>F(Ci=GT)
	// �����ȊO     : �������Ȃ�

	if( !(calced && calced_rho==rho) ){
		
		A=this->A; B1=this->A; B2=this->A;
		C=this->C;
		W=this->W;
		Fa=this->F; Fb1=this->F; Fb2=this->F;

		m=NumberOfEq;
		n=NumberOfVar;
		ma=mb1=mb2=0;

		for(ii=m; ii>=1; --ii){
			if(C[ii][1]==DLS && W[ii][1]>0){ // �ŏ����@�ōœK��
				B1.delrow(ii);
				B2.delrow(ii);
				Fb1.delrow(ii);
				Fb2.delrow(ii);
				ma++;               // �ŏ����@�ōœK������]���֐��̐�
			}
			else if(C[ii][1]==EQ){  // 1�������ōS��
				A.delrow(ii);
				B2.delrow(ii);
				W.delrow(ii);
				Fa.delrow(ii);
				Fb2.delrow(ii);
				mb1++;              // 1�������ōS������]���֐��̐�
			}
			else if(C[ii][1]==LT || C[ii][1]==GT){    // �ꎟ�s�����ōS��
				if(C[ii][1]==GT){
					// CX>Fc�̗��ӂ�-1���|���ĕs�����̌������t�ɂ���
					for(jj=1; jj<=n; ++jj) B2[ii][jj]=-B2[ii][jj];
					Fb2[ii][1]=-Fb2[ii][1];
				}
				A.delrow(ii);
				B1.delrow(ii);
				W.delrow(ii);
				Fa.delrow(ii);
				Fb1.delrow(ii);
				mb2++;
			}
			else{                   // �����ȊO�̕]���֐�(C=0����W=0�Ȃ�)�͏�������D
				A.delrow(ii);
				B1.delrow(ii);
				B2.delrow(ii);
				W.delrow(ii);
				Fa.delrow(ii);
				Fb1.delrow(ii);
				Fb2.delrow(ii);
			}
		}
			
		for(kk=1; kk<=MAXTIMES; ++kk){
			//  ======== �L������@�i�s�����S���̏����j ===============================
			//    �Ⴆ�΁C2����(x,y)�� ��=(x-a)^2+(y-b)^2 ((x,y)=(a,b)�Ńӂ��ŏ�)
			//    �� y<a1*x+b1, y<a2*x+b2 �̏����̂��Ƃōŏ����������x,y���ʂŐ}�ɕ`���Ă݂��
			//    �����I�ɔc�����₷��. (�Q�l�F�g��b�����ϕ�(�F�V�����j)�h�̐}4.9�j
			if(kk==1){      // kk���[�v1��ڂ�
				mb2_=0;     // �s������S�Ė�������
			}
			else if(kk==2){ // kk���[�v2���
				B2_=B2;
				Fb2_=Fb2;
				mb2_=mb2;
				BX=B2*X;
				for(ii=mb2; ii>=1; --ii){
					if(Round(BX[ii][1],5)<=Round(Fb2[ii][1],5)){   // �s�����𖞑��ł����
						B2_.delrow(ii);         // ���̎����폜
						Fb2_.delrow(ii);
						--mb2_;
					}
				}
				// ���̎��_�ŁC��������Ȃ��s�����̂ݎc���Ă���D
				if(mb2_==0) break;  // DLS�@�Ɠ����ōœK����C�s�������������Ă��邱�ƂɂȂ�̂ŏI��
			}
			else{ // kk���[�v3��ڈȍ~
				e1=true;                  // �I���t���O
				for(ii=mb2_; ii>=1; --ii){
					BX=B2_*X;
					if(XX[n+mb1+ii][1]<=0 && Round(BX[ii][1],5)<=Round(Fb2_[ii][1],5)){ // Lagrange�搔�����ŁC
						                                                                // ���s�����𖞑��ł���΁C
						    // Round()���Ȃ��ƁCdouble�̗L�����܂ł͐��l�����肵�Ȃ����߁C
						    // �R���p�C���I�v�V�����Ȃǂɂ��s�����̕]�����قȂ邱�Ƃ�����D
						e1=false;         // �I���t���O��|���C
						B2_.delrow(ii);   // ���̎����폜����D
						Fb2_.delrow(ii);
						--mb2_;
					}
				}
				if(e1) break;
			}
			//  ========= �L������@ �I��� ==========================

			if(ma>0){
				matrix<double> MU(n,n);
				matrix<double> WA(ma,n);
				matrix<double> WFa(ma,1);

				for(ii=1; ii<=ma; ++ii) for(jj=1; jj<=n; ++jj) WA[ii][jj]=W[ii][1]*A[ii][jj];
				// ���F
				//   W��Ίp�s��Ƃ��C�Ίp�����ɃE�G�C�g��ݒ肷��΁C
				//   WA�͍s��̐�W*A �ƂȂ�C�R�[�h�͊Ȍ��ɂȂ�D
				//   �������CW�̑Ίp�����ȊO��0�Ȃ̂�NumberOfEq���傫���Ƃ�
				//   �������̖��ʁC�����W*A�ɂ����鑽�����0�Ƃ̏�Z�ɂ�閳�ʂ��傫���Ȃ�D

				for(jj=1; jj<=n; ++jj){
					for(ii=1; ii<=ma; ++ii){
						MU[jj][jj]+=WA[ii][jj]*WA[ii][jj];  // ����O�g�����Y���w�h A.6 DLS�@ �� ��j
					}
					if(MU[jj][jj]==0) MU[jj][jj]=1;  // ���֌W(���ׂĂ̕]���֐��̔����W����0)�ȕϐ�Xj
													 // ���܂ނƂ��ł��CA1�𐳑��Ƃ���
				}

				if(GradientDescent==0){
					A1=t(WA)*WA+rho*MU;
				}
				else{
					// �ŋ}�~���@�ł�A1���قȂ�̂݁i�g�����Y�݌v (�����F���C���C��w�o�ŉ�j�h�\8.2 ���ސ������j
					// DLS�@�Ɠ��l�Ƀς��傫���ق�X���������������̂ŁC�t�����|����D
					A1.redim(n,n);
					A1=unit(A1)*(rho==0 ? 1e30 : 1.0/rho);
				}

				for(ii=1; ii<=ma; ++ii) WFa[ii][1]=W[ii][1]*Fa[ii][1];
				Fa1=t(WA)*WFa;	
			}
			else{
				A1.redim(n,n);
				Fa1.redim(n,1);
			}

			// �g�����Y�݌v (�����F���C���C��w�o�ŉ�j�h(8.45)���ƑΉ�������ƁC
			//     A1      = t(A)*A+pE  (�A���C�ŋ}�~���@�ł�(8.10)��� 2pE )
			//     B1,B2   = B
			//     Fa1     = t(A)*O
			//     Fb1,Fb2 = F-F(0)
			//  �ƂȂ�D
			AA.redim(n+mb1+mb2_,n+mb1+mb2_);   // (8.45)�����ӂ̍ŏ��̍s��
			for(ii=1; ii<=n;    ++ii)     for(jj=1; jj<=n;    ++jj)     AA[ii][jj]=A1[ii][jj];
			for(ii=1; ii<=mb1;  ++ii)     for(jj=1; jj<=n;    ++jj)     AA[n+ii][jj]=B1[ii][jj];
			for(ii=1; ii<=mb2_; ++ii)     for(jj=1; jj<=n;    ++jj)     AA[n+mb1+ii][jj]=B2_[ii][jj];
			for(ii=1; ii<=n;    ++ii)     for(jj=1; jj<=mb1;  ++jj)     AA[ii][n+jj]=t(B1)[ii][jj];
			for(ii=1; ii<=n;    ++ii)     for(jj=1; jj<=mb2_; ++jj)     AA[ii][n+mb1+jj]=t(B2_)[ii][jj];
			for(ii=1; ii<=mb1+mb2_; ++ii) for(jj=1; jj<=mb1+mb2_; ++jj) AA[n+ii][n+jj]=0;

			FF.redim(n+mb1+mb2_,1);     // (8.45)���E��
			for(ii=1; ii<=n;    ++ii) FF[ii      ][1]=Fa1[ii][1];
			for(ii=1; ii<=mb1;  ++ii) FF[n+ii    ][1]=Fb1[ii][1];
			for(ii=1; ii<=mb2_; ++ii) FF[n+mb1+ii][1]=Fb2_[ii][1];

			XX=inv(AA)*FF;        // (8.46)��
			for(jj=1; jj<=n; ++jj) X[jj][1]=XX[jj][1];

			calced=true;
			calced_rho=rho; 

			if(mb2==0) break;  // ���Ƃ��ƕs�����������Ȃ��ꍇ��kk���[�v���񂳂��ɏI��
		} // next kk
	}

	return X[j][1];
}

double cLeneq::GetX(int j) {
	return GetDampedX(j,0);
}

double cLeneq::pvError(double rho){
	int i;
	double peak,valley;
	GetDampedX(1,rho);
	matrix<double> F1;
	
	F1=A*X;
	peak=-1e30;
	valley=1e30;
	for(i=1; i<=GetNumberOfEq(); i++){
		if( F1[i][1]-F[i][1]>peak   ) peak  =F1[i][1]-F[i][1];
		if( F1[i][1]-F[i][1]<valley ) valley=F1[i][1]-F[i][1];
	}
	return peak-valley;	
}

double cLeneq::pvError(){
	return pvError(0);
}