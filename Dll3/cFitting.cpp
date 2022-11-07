#include "stdafx.h"
#include "cFitting.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#include <math.h>


//////////// cFitting members /////////////////////////////////////////////////////////////////////

///// private ///////////////////////////

int cFitting::N(){
	return NumberOfTerms+2;
}

void cFitting::alloc_terms() {
	int j;
	int m=GetNumberOfData();
	int n=NumberOfTerms;
	x1dim=new int [n+1];
	x2dim=new int [n+1];
	for(j=0; j<=n; j++) x1dim[j]=x2dim[j]=0;
	Digits=new int [n+1];
	for(j=0; j<=n; j++) Digits[j]=0;
	C.redim(n,1);
}

void cFitting::free_terms() {
	delete [] x1dim; delete [] x2dim;
	delete [] Digits;
}

double cFitting::f(double x1,double x2,double y,double dummy){
	int j;
	double a;

	a=0;
	for(j=1; j<=NumberOfTerms; j++){
		a+=C[j][1]*pow( (x1-x10),x1dim[j] )*pow( (x2-x20),x2dim[j] );
	}
	return a-y;
}

void cFitting::SetA(){
	int i,j,jj;

	for(i=1; i<=GetNumberOfData(); i++) for(j=1; j<=N(); j++) {
		if(j==N()-1){
			A[i][N()-1]=0;
			if(x10Variable!=0){
				for(jj=1; jj<=NumberOfTerms; jj++){
					if(x1dim[jj]!=0){
						A[i][N()-1]-=C[jj][1]*x1dim[jj]
						            *pow(x1[i]-x10,x1dim[jj]-1)*pow(x2[i]-x20,x2dim[jj]);
					}
				}
			}
		}
		else if(j==N()){
			A[i][N()]=0;
			if(x20Variable!=0){
				for(jj=1; jj<=NumberOfTerms; jj++){
					if(x2dim[jj]!=0){
						A[i][N()]-=C[jj][1]
						          *pow(x1[i]-x10,x1dim[jj])*x2dim[jj]*pow(x2[i]-x20,x2dim[jj]-1);
					}
				}
			}
		}
		else{
			A[i][j]=pow(x1[i]-x10,x1dim[j])*pow(x2[i]-x20,x2dim[j]);
		}
	}

}

void cFitting::SetX(){
	int j;

	for(j=1; j<=N(); j++) {
		if(j==N()-1){
			X[j]=&x10;
		}
		else if(j==N()){
			X[j]=&x20;
		}
		else{
			X[j]=&C[j][1];
		}
	}
}

int cFitting::NumberOfCoefficients(){
	return N();
}

///// public ///////////////////////////

cFitting::cFitting() : y(x3){
	SetNumberOfData(0);
	NumberOfTerms=3;
	alloc_terms();
	x10=x20=0;
	x10Variable=x20Variable=0;
	CoefficientsIsNew=false;
}

cFitting::~cFitting() {
	free_terms();
}

void cFitting::SetNumberOfTerms(int n) {
	NumberOfTerms=n;
	free_terms();
	alloc_terms();
	CoefficientsIsNew=false;
}

int cFitting::GetNumberOfTerms() const {
	return NumberOfTerms;
}

void cFitting::SetOrder(int order){
	// orderを指定して自動的に SetNumberOfTerms(),dimensionSet()
	// を実行する．
	int j,jj,k,n, x1_order,x2_order;
	
	n=(order+1)*(order+2)/2;
	SetNumberOfTerms(n);

	j=1;
	for(jj=0; jj<=order; jj++){
		x1_order=jj;
		x2_order=0;
		for(k=1; k<=jj+1; k++){
			dimensionSet(j,x1_order,x2_order);
			x1_order--;
			x2_order++;
			j++;
		}
	}
}

void cFitting::dimensionSet(int j,int x1dimension,int x2dimension) {
	if(1<=j && j<=NumberOfTerms) {
		x1dim[j]=x1dimension;
		x2dim[j]=x2dimension;
	}
	CoefficientsIsNew=false;
}

int cFitting::x1dimensionGet(int j) const {
	if(1<=j && j<=NumberOfTerms){
		return x1dim[j];
	}
	else return 0;
}

int cFitting::x2dimensionGet(int j) const {
	if(1<=j && j<=NumberOfTerms){
		return x2dim[j];
	}
	else return 0;
}

void cFitting::RemoveX1Odd(){
	// x1の奇数次項を削除する
	int j, counter;
	int *buf1,*buf2;

	buf1=new int[NumberOfTerms+1];
	buf2=new int[NumberOfTerms+1];
	
	counter=0;
	for(j=1; j<=NumberOfTerms; ++j){
		if(is_even(x1dim[j])){  // x1の偶数次項を残す
			counter++;
			buf1[counter]=x1dim[j];
			buf2[counter]=x2dim[j];
		}
	}

	SetNumberOfTerms(counter);

	for(j=1; j<=NumberOfTerms; ++j){
		x1dim[j]=buf1[j];
		x2dim[j]=buf2[j];
	}

	CoefficientsIsNew=false;
	delete buf1;
	delete buf2;
}

void cFitting::RemoveX1Even(){
	// x1の偶数次項を削除する．ただし定数項は残す．
	int j, counter;
	int *buf1,*buf2;

	buf1=new int[NumberOfTerms+1];
	buf2=new int[NumberOfTerms+1];
	
	counter=0;
	for(j=1; j<=NumberOfTerms; ++j){
		if( is_odd(x1dim[j]) || (x1dim[j]==0 && x2dim[j]==0) ){  // x1の奇数次項と定数項を残す
			counter++;
			buf1[counter]=x1dim[j];
			buf2[counter]=x2dim[j];
		}
	}

	SetNumberOfTerms(counter);

	for(j=1; j<=NumberOfTerms; ++j){
		x1dim[j]=buf1[j];
		x2dim[j]=buf2[j];
	}

	CoefficientsIsNew=false;
	delete buf1;
	delete buf2;
}

void cFitting::RemoveX2Odd(){
	// x2の奇数次項を削除する
	int j, counter;
	int *buf1,*buf2;

	buf1=new int[NumberOfTerms+1];
	buf2=new int[NumberOfTerms+1];
	
	counter=0;
	for(j=1; j<=NumberOfTerms; ++j){
		if(is_even(x2dim[j])){  // x2の偶数次項を残す
			counter++;
			buf1[counter]=x1dim[j];
			buf2[counter]=x2dim[j];
		}
	}

	SetNumberOfTerms(counter);

	for(j=1; j<=NumberOfTerms; ++j){
		x1dim[j]=buf1[j];
		x2dim[j]=buf2[j];
	}

	CoefficientsIsNew=false;
	delete buf1;
	delete buf2;
}

void cFitting::RemoveX2Even(){
	// x2の偶数次項を削除する．ただし定数項は残す．
	int j, counter;
	int *buf1,*buf2;

	buf1=new int[NumberOfTerms+1];
	buf2=new int[NumberOfTerms+1];
	
	counter=0;
	for(j=1; j<=NumberOfTerms; ++j){
		if( is_odd(x2dim[j]) || (x1dim[j]==0 && x2dim[j]==0) ){  // x2の奇数次項と定数項を残す
			counter++;
			buf1[counter]=x1dim[j];
			buf2[counter]=x2dim[j];
		}
	}

	SetNumberOfTerms(counter);

	for(j=1; j<=NumberOfTerms; ++j){
		x1dim[j]=buf1[j];
		x2dim[j]=buf2[j];
	}

	CoefficientsIsNew=false;
	delete buf1;
	delete buf2;
}

void cFitting::dataSet(int i,double x1_value,double x2_value,double y_value) {
	if(1<=i && i<=GetNumberOfData()){
		x1[i]=x1_value;
		x2[i]=x2_value;
		y[i]=y_value;
	}
	CoefficientsIsNew=false;
}

void cFitting::AddData(double x1_value,double x2_value,double y_value){
	cLeastSquares::AddData(x1_value,x2_value,y_value);
	CoefficientsIsNew=false;
}

void cFitting::DataClear(){
	cLeastSquares::DataClear();
	CoefficientsIsNew=false;
}

double cFitting::x1dataGet(int i) const{
	if(1<=i && i<=GetNumberOfData()){
		return x1[i];
	}
	else return 0;
}

double cFitting::x2dataGet(int i) const{
	if(1<=i && i<=GetNumberOfData()){
		return x2[i];
	}
	else return 0;
}

double cFitting::ydataGet(int i) const{
	if(1<=i && i<=GetNumberOfData()){
		return y[i];
	}
	else return 0;
}

void cFitting::DigitsSet(int j,int value) {
	if(1<=j && j<=NumberOfTerms){
		Digits[j]=value;
	}
	CoefficientsIsNew=false;
}

int cFitting::CalcCoefficients() {
	const double RHO=0.1;
	// x0,y0固定の場合は最小二乗法は線形で逐次近似は不要だが，
	// RHO=0とすると解けない．
	// ここではΔx0,Δy0が解けない方程式となるため，t(A)Aが正則でなくなるためか?
	// x0,y0固定のときはAの列数を減らせばいいのかもしないが，場合分けが面倒である．

	int j,improved;
	
	improved=0;
	if(CoefficientsIsNew==false){
		improved=ModifyCoefficients(RHO);
		for(j=1; j<=NumberOfTerms; j++){
			if(Digits[j]!=0) C[j][1]=Round(C[j][1],Digits[j]);
			// 小数点以下四捨五入Round(C[j][1],0)は不可
			// 現状，x10,x20はまるめていない
		}
		CoefficientsIsNew=true;
		return improved;
	}
	else {
		return improved;
	}
}

double cFitting::coefficientGet(int j,int CoefficientCalc) {
	// CoefficientCalc==0 のとき係数を計算しない．
	// C[j][1]を再フィッティングしないで今の値を取得するときなどに CoefficientCalc==0 とする．
	// <参考>
	//    CoefficientCalc!=0の場合，CalcCoefficients()が実行されるが，
	//    CalcCoefficients()は係数が既にフィッティングされていれば無駄に再計算しない．

	if(CoefficientCalc){
		CalcCoefficients();  
	}

	if(1<=j && j<=NumberOfTerms) {
		return C[j][1];
	}
	else{
		return 0;
	}
}

double cFitting::coefficientGet(int x1_order,int x2_order,int CoefficientCalc){
	int j;

	for(j=1; j<=NumberOfTerms; j++){
		if( x1_order==x1dim[j] && x2_order==x2dim[j] ){
			return coefficientGet(j,CoefficientCalc);
		}
	}
	return 0;
}

void cFitting::coefficientSet(int j,double value) {
	if(1<=j && j<=NumberOfTerms) {
		C[j][1]=value;
	}
}

double cFitting::x10Get() const{
	return x10;
}

void cFitting::x10Set(double value){
	x10=value;
	CoefficientsIsNew=false;
}

double cFitting::x20Get() const{
	return x20;
}

void cFitting::x20Set(double value){
	x20=value;
	CoefficientsIsNew=false;
}

int cFitting::x10VariableGet() const{
	return x10Variable;
}

void cFitting::x10VariableSet(int value){
	x10Variable=value;
	CoefficientsIsNew=false;
}

int cFitting::x20VariableGet() const{
	return x20Variable;
}

void cFitting::x20VariableSet(int value){
	x20Variable=value;
	CoefficientsIsNew=false;
}

double cFitting::yApproximate(double x1value,double x2value,int EliminateConst,int CoefficientCalc){
	// EliminateConst!=0 のときは定数項を除去する．
	// CoefficientCalc==0 とすれば係数をフィッティングしない．
	int j;
	double ya;
	if(CoefficientCalc){
		CalcCoefficients();
	}
	ya=0;
	for(j=1; j<=NumberOfTerms; j++){
		if( !( EliminateConst!=0 && x1dim[j]==0 && x2dim[j]==0) ){
			ya+=C[j][1]*pow( x1value-x10,x1dim[j] )*pow( x2value-x20,x2dim[j] );
		}
	}
	return ya;
}

double cFitting::yApproximate(double x1value,double x2value,int CoefficientCalc){
	return yApproximate(x1value,x2value,0,CoefficientCalc);
}

double cFitting::dyApproximate(double x1value,double x2value,int CoefficientCalc) {
	// 0,1次の係数で決まる線形部分との差（結像における収差）を計算する．
	// CoefficientCalc==0 とすれば係数を再計算しない．
	int j;
	double ya0;

	if(CoefficientCalc){
		CalcCoefficients();
	}

	ya0=0;
	for(j=1; j<=NumberOfTerms; j++){
		int p,q;
		p=x1dim[j];
		q=x2dim[j];
		if( (p==0 && q==0) || (p==1 && q==0) || (p==0 && q==1) ){
			ya0+=C[j][1]*pow( x1value-x10,p )*pow( x2value-x20,q );
		}
	}

	return yApproximate(x1value,x2value,CoefficientCalc)-ya0;
}


double cFitting::Error(int i,int CoefficientCalc) {
	// CoefficientCalc==0 とすれば係数を再計算しない．
	// coefficientSet()で指定した係数を評価をしたいときなどに使用する．
	int j;
	double ya=0;

	if(CoefficientCalc){
		CalcCoefficients();
	}
	if(1<=i && i<=GetNumberOfData()) {
		for(j=1; j<=NumberOfTerms; j++){
			ya+=C[j][1]*pow( x1[i]-x10,x1dim[j] )*pow( x2[i]-x20,x2dim[j] );
		}
		return ya-y[i];
	}
	else return 0;
}

double cFitting::rmsError(int CoefficientCalc) {
	// CoefficientCalc==0 とすれば係数を再計算しない．
	// coefficientSet()で指定した係数を評価をしたいときなどに使用する．

	if(CoefficientCalc){
		CalcCoefficients();
	}	
	return cLeastSquares::rmsError();
}

double cFitting::pvError(int CoefficientCalc) {
	// CoefficientCalc==0 とすれば係数を再計算しない．
	// coefficientSet()で指定した係数を評価をしたいときなどに使用する．

	if(CoefficientCalc){
		CalcCoefficients();
	}
	return cLeastSquares::pvError();
}

double cFitting::maxError(int CoefficientCalc) {
	// CoefficientCalc==0 とすれば係数を再計算しない．
	// coefficientSet()で指定した係数を評価をしたいときなどに使用する．

	if(CoefficientCalc){
		CalcCoefficients();
	}
	return cLeastSquares::maxError();
}

double cFitting::R2(int CoefficientCalc) {
	// CoefficientCalc==0 とすれば係数を再計算しない．
	// coefficientSet()で指定した係数を評価をしたいときなどに使用する．
	int i;
	cXYList X;

	if(CoefficientCalc){
		CalcCoefficients();
	}
	for(i=1; i<=GetNumberOfData(); i++){
		X.AddData(y[i],yApproximate(x1[i],x2[i],0));
	}
	return X.R2();
}

void cFitting::RemoveTerms() {
	int i;
	double *temp;
	temp=new double[GetNumberOfData()+1];
	for(i=1; i<=GetNumberOfData(); i++){
		temp[i]=-Error(i,1);
	}
	for(i=1; i<=GetNumberOfData(); i++){
		y[i]=temp[i];
	}
	CoefficientsIsNew=false; CalcCoefficients();
	delete [] temp;
}

std::ostream& operator<<(std::ostream& to,const cFitting& x){
	// データ一覧を出力する
	int i;
	
	for(i=1; i<=x.GetNumberOfData(); i++){
		to << x.x1dataGet(i) << '\t';
		to << x.x2dataGet(i) << '\t';
		to << x.ydataGet(i)  << std::endl;
	}
	return to;
}

