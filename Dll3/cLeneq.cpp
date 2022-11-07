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
		// WよりCを優先する．
		// こうしておけば，例えばDLS法でメリット関数を計算するとき，コンストレインツに設定した
		// 評価関数のWi がゼロになるため，DLSの関数かコンストレインツの関数か場合分けしなくて済む.
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
	// GradientDescentが真のときはDLS法の代わりに最急降下法を用いる
	const int MAXTIMES=100;    // 有効制約法の繰り返し回数上限
	int ii,jj,kk, m,ma,mb1,mb2,mb2_, n;
	bool e1;
	matrix<double> A,B1,B2,B2_, W, Fa,Fb1,Fb2,Fb2_, A1,Fa1, AA,FF,XX, BX;
	matrix<int> C;
	// Ci=DLS && Wi>0 : 最小二乗法による最適化
	// Ci=EQ          : 一次等式拘束(Lagrange未定乗数法）Ax=F
	// Ci=LT,GT       : 一次不等式拘束(有効制約法）    Ax<F(Ci=LT), Ax>F(Ci=GT)
	// これら以外     : 何もしない

	if( !(calced && calced_rho==rho) ){
		
		A=this->A; B1=this->A; B2=this->A;
		C=this->C;
		W=this->W;
		Fa=this->F; Fb1=this->F; Fb2=this->F;

		m=NumberOfEq;
		n=NumberOfVar;
		ma=mb1=mb2=0;

		for(ii=m; ii>=1; --ii){
			if(C[ii][1]==DLS && W[ii][1]>0){ // 最小二乗法で最適化
				B1.delrow(ii);
				B2.delrow(ii);
				Fb1.delrow(ii);
				Fb2.delrow(ii);
				ma++;               // 最小二乗法で最適化する評価関数の数
			}
			else if(C[ii][1]==EQ){  // 1次等式で拘束
				A.delrow(ii);
				B2.delrow(ii);
				W.delrow(ii);
				Fa.delrow(ii);
				Fb2.delrow(ii);
				mb1++;              // 1次等式で拘束する評価関数の数
			}
			else if(C[ii][1]==LT || C[ii][1]==GT){    // 一次不等式で拘束
				if(C[ii][1]==GT){
					// CX>Fcの両辺に-1を掛けて不等号の向きを逆にする
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
			else{                   // これら以外の評価関数(C=0かつW=0など)は消去する．
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
			//  ======== 有効制約法（不等式拘束の処理） ===============================
			//    例えば，2次元(x,y)で φ=(x-a)^2+(y-b)^2 ((x,y)=(a,b)でφが最小)
			//    を y<a1*x+b1, y<a2*x+b2 の条件のもとで最小化する問題をx,y平面で図に描いてみると
			//    直感的に把握しやすい. (参考：“基礎微分積分(洲之内治男)”の図4.9）
			if(kk==1){      // kkループ1回目は
				mb2_=0;     // 不等式を全て無視する
			}
			else if(kk==2){ // kkループ2回目
				B2_=B2;
				Fb2_=Fb2;
				mb2_=mb2;
				BX=B2*X;
				for(ii=mb2; ii>=1; --ii){
					if(Round(BX[ii][1],5)<=Round(Fb2[ii][1],5)){   // 不等式を満足であれば
						B2_.delrow(ii);         // その式を削除
						Fb2_.delrow(ii);
						--mb2_;
					}
				}
				// この時点で，満足されない不等式のみ残っている．
				if(mb2_==0) break;  // DLS法と等式で最適化後，不等式も満たしていることになるので終了
			}
			else{ // kkループ3回目以降
				e1=true;                  // 終了フラグ
				for(ii=mb2_; ii>=1; --ii){
					BX=B2_*X;
					if(XX[n+mb1+ii][1]<=0 && Round(BX[ii][1],5)<=Round(Fb2_[ii][1],5)){ // Lagrange乗数が負で，
						                                                                // かつ不等式を満足であれば，
						    // Round()がないと，doubleの有効桁までは数値が安定しないため，
						    // コンパイラオプションなどにより不等式の評価が異なることがある．
						e1=false;         // 終了フラグを倒し，
						B2_.delrow(ii);   // その式を削除する．
						Fb2_.delrow(ii);
						--mb2_;
					}
				}
				if(e1) break;
			}
			//  ========= 有効制約法 終わり ==========================

			if(ma>0){
				matrix<double> MU(n,n);
				matrix<double> WA(ma,n);
				matrix<double> WFa(ma,1);

				for(ii=1; ii<=ma; ++ii) for(jj=1; jj<=n; ++jj) WA[ii][jj]=W[ii][1]*A[ii][jj];
				// 注：
				//   Wを対角行列とし，対角成分にウエイトを設定すれば，
				//   WAは行列の積W*A となり，コードは簡潔になる．
				//   しかし，Wの対角成分以外は0なのでNumberOfEqが大きいとき
				//   メモリの無駄，およびW*Aにおける多数回の0との乗算による無駄が大きくなる．

				for(jj=1; jj<=n; ++jj){
					for(ii=1; ii<=ma; ++ii){
						MU[jj][jj]+=WA[ii][jj]*WA[ii][jj];  // 草川徹“レンズ光学” A.6 DLS法 の μj
					}
					if(MU[jj][jj]==0) MU[jj][jj]=1;  // 無関係(すべての評価関数の微分係数が0)な変数Xj
													 // を含むときでも，A1を正則とする
				}

				if(GradientDescent==0){
					A1=t(WA)*WA+rho*MU;
				}
				else{
					// 最急降下法ではA1が異なるのみ（“レンズ設計 (高橋友刀，東海大学出版会）”表8.2 より類推した）
					// DLS法と同様にρが大きいほどXを小さくしたいので，逆数を掛ける．
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

			// “レンズ設計 (高橋友刀，東海大学出版会）”(8.45)式と対応させると，
			//     A1      = t(A)*A+pE  (但し，最急降下法では(8.10)より 2pE )
			//     B1,B2   = B
			//     Fa1     = t(A)*O
			//     Fb1,Fb2 = F-F(0)
			//  となる．
			AA.redim(n+mb1+mb2_,n+mb1+mb2_);   // (8.45)式左辺の最初の行列
			for(ii=1; ii<=n;    ++ii)     for(jj=1; jj<=n;    ++jj)     AA[ii][jj]=A1[ii][jj];
			for(ii=1; ii<=mb1;  ++ii)     for(jj=1; jj<=n;    ++jj)     AA[n+ii][jj]=B1[ii][jj];
			for(ii=1; ii<=mb2_; ++ii)     for(jj=1; jj<=n;    ++jj)     AA[n+mb1+ii][jj]=B2_[ii][jj];
			for(ii=1; ii<=n;    ++ii)     for(jj=1; jj<=mb1;  ++jj)     AA[ii][n+jj]=t(B1)[ii][jj];
			for(ii=1; ii<=n;    ++ii)     for(jj=1; jj<=mb2_; ++jj)     AA[ii][n+mb1+jj]=t(B2_)[ii][jj];
			for(ii=1; ii<=mb1+mb2_; ++ii) for(jj=1; jj<=mb1+mb2_; ++jj) AA[n+ii][n+jj]=0;

			FF.redim(n+mb1+mb2_,1);     // (8.45)式右辺
			for(ii=1; ii<=n;    ++ii) FF[ii      ][1]=Fa1[ii][1];
			for(ii=1; ii<=mb1;  ++ii) FF[n+ii    ][1]=Fb1[ii][1];
			for(ii=1; ii<=mb2_; ++ii) FF[n+mb1+ii][1]=Fb2_[ii][1];

			XX=inv(AA)*FF;        // (8.46)式
			for(jj=1; jj<=n; ++jj) X[jj][1]=XX[jj][1];

			calced=true;
			calced_rho=rho; 

			if(mb2==0) break;  // もともと不等式条件がない場合はkkループを回さずに終了
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