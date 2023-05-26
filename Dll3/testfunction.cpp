#include "testfunction.h"
#include "testfunction.h"
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "complex.h"
#include "Matrix.h"
#include "cLens1.h"

complex_wrap* static_c;
complex* p_complex;

// ↓最終的には、これは使わないで動作出来るようにする（沢山扱える）
//  グローバル変数を使ってのやり取りだと1個だけになる
// 
// もしかするとダブルポインタ**でのやり取りになるかも
//  実体をやり取りさせたいためのポインタ(初期化したい)
// 
//　T **a　のダブルポインタの使い方(多次元配列)とは異なる

matrix_wrap* static_m;
matrix<double>* p_matrix;

void test_01()
{
	printf("\nHello Dll3!!! 2022-11-07\n\n");
}

double Return123() 
{
	return 123;
}

complex_wrap* test_001() {
	static_c = new complex_wrap();
	static_c->p1 = 0.123;
	static_c->p2 = 0.456;
	return static_c;
}

void New_complex()
{
	p_complex = new complex();
}

void Delete_complex()
{
	delete p_complex;
}

void ABS(complex_wrap* p){
	// New_complex();
	p_complex = new complex(p->p1, p->p2);
	// p_complex->x = p->p1;
	// p_complex->y = p->p2;
	double test;
	// complex テスト(1,1);
	// test = abs(テスト);
	test = abs(*p_complex);
	std::cout << "ABSの関数を動かす" << test << std::endl;
	Delete_complex();
}

double ARG(complex_wrap* p){
	p_complex = new complex(p->p1, p->p2);
	double test;
	test = arg(*p_complex);
	std::cout << "argの関数を動かす" << test << std::endl;
	Delete_complex();
	return test;
}

void Conjugate(complex_wrap* comp1, complex_wrap* comp2)
{
	// インプット、アウトプット用の複素数を二つ用意する
	complex* p_c1 = new complex(comp1->p1, comp1->p2);
	complex* p_c2 = new complex(comp2->p1, comp2->p2);

	*p_c2 = conj(*p_c1);
	comp2->p1 = Re(*p_c2);
	comp2->p2 = Im(*p_c2);

	std::cout << "動作確認用メッセージ 2022-11-14 " << std::endl;

	delete p_c1;
	delete p_c2;

}

void New_matrix(matrix_wrap* p)
{
	// matrix<double> A1();
	//matrix<double>* A0 = new matrix<double>();
	// p_matrix = new matrix<double>();
}

void New_matrixMN(matrix_wrap* p, int m, int n)
{
	p_matrix = new matrix<double>(m, n);
	p->m = m;
	p->n = n;
	p->dpmatrix_wrap = p_matrix->a;

}
void New_matrixCopy(matrix_wrap* p1, matrix_wrap* p2)
{
}
void Delete_matrix(matrix_wrap* p)
{
	delete p_matrix;
}

/*
void matrix_init(matrix_wrap* p, int m, int n)
{
	p_matrix = new matrix<double> (m, n);
	std::cout << "matrix_init呼び出し成功"<< std::endl;
	p->m = m;
	p->n = n;

}
*/

// 行列の初期化
void matrix_init(matrix<double>* p, int m, int n) {
	//p_matrix = new matrix<double>(m, n);
	// p_matrix->a[1][1] = 1.1;
	// p_matrix->a[0][1] = 2.2;
	// p_matrix->a[1][0] = 3.3;
	// p_matrix->a[1][1] = 4.4;
	// p = p_matrix;

	// p = new matrix<double>(m, n);

	// 要素m,nで書き換え
	p->redim(m, n);

	p->a[1][1] = 1.0;
    p->a[1][2] = 3.0;
    p->a[2][1] = 2.0;
    p->a[2][2] = 5.0;
	
}

// 逆行列の関数
void matrix_inv(matrix<double>* p1, matrix<double>* p2) 
{
	*p2 = inv(*p1);
}



matrix<double>* CreateMatrix()
{
	return new matrix<double>();
}

matrix<double>* CreateMatrixMN(int m, int n)
{
	return new matrix<double>(m, n);
}

void SetMatrixIJ(matrix<double>* pMatrix, int m, int n, double value)
{
	pMatrix->a[m][n] = value;
}

double GetMatrixIJ(matrix<double>* pMatrix, int m, int n)
{
	return pMatrix->a[m][n];
}


void DeleteMatrix(matrix<double>* pMatrix)
{
	delete pMatrix;
}

// cLen1の生成
cLens1* Create_cLens1()
{
	return new cLens1();
}

void Delete_cLens1(cLens1* p_cLens1)
{
	delete p_cLens1;
}

void SetRadius(cLens1* p_cLens1, int surf_i, double value) 
{
	p_cLens1->Set_r(surf_i, value);
}

double GetRadius(cLens1* p_cLens1, int surf_i)
{
	return p_cLens1->Get_r(surf_i);
}

void SetDistance(cLens1* p_cLens1, int surf_i, double value)
{
	p_cLens1->Set_d(surf_i, value);
}

double GetDistance(cLens1* p_cLens1, int surf_i)
{
	return p_cLens1->Get_d(surf_i);
}

void SetGlassName(cLens1* p_cLens1, int surf_i, const char* pText)
{
	std::string str = pText;
	std::cout << "strの中身を確認する" << std::endl << str << std::endl;
	p_cLens1->Set_gname(surf_i, str);
}

BSTR GetGlassName(cLens1* p_cLens1, int surf_i)
{
	// return p_cLens1->Get_gname(surf_i);
	std::string str;
	str = p_cLens1->Get_gname(surf_i);

	// この100という適当な設定で良いのか？
	wchar_t ws[100];
	mbstowcs(ws, str.c_str(), 100);
	return SysAllocString(ws);
	// return SysAllocString(L"GetGlassNameだよ");
}

int GetK(cLens1* p_cLens1) 
{
	return p_cLens1->Get_k();
}

double focallength(cLens1* p_cLens1)
{
	return p_cLens1->f();
}

double backf(cLens1* p_cLens1) 
{
	return p_cLens1->bf();
}

void SetStop(cLens1* p_cLens1, int surf_i) 
{
	p_cLens1->Set_stop(surf_i);
}

int GetStop(cLens1* p_cLens1)
{
	return p_cLens1->Get_stop();
}

void SetEAy(cLens1* p_cLens1, int surf_i, double value)
{
	p_cLens1->Set_EAy(surf_i, value);
}

double GetEAy(cLens1* p_cLens1, int surf_i)
{
	return p_cLens1->Get_EAy(surf_i);
}

void EPCalculation(cLens1* p_cLens1)
{
	p_cLens1->EPCalc();
}

void Set_s(cLens1* p_cLens1, double value)
{
	p_cLens1->Set_s(value);
}

double Get_s(cLens1* p_cLens1)
{
	return p_cLens1->Get_s();
}

void Set_t(cLens1* p_cLens1, double value)
{
	p_cLens1->Set_t(value);
}

double Get_t(cLens1* p_cLens1)
{
	return p_cLens1->Get_t();
}

void Set_EPD(cLens1* p_cLens1, double value)
{
	p_cLens1->Set_EPD(value);
}

double Get_EPD(cLens1* p_cLens1)
{
	return p_cLens1->Get_EPD();
}

void SetColorN(cLens1* p_cLens1,int num)
{
	p_cLens1->Set_cn(num);
}

int GetColorN(cLens1* p_cLens1)
{
	return p_cLens1->Get_cn();
}

void SetColorWeight(cLens1* p_cLens1, int num, double value)
{
	p_cLens1->Set_colorweight(num, value);
}

double GetColorWeight(cLens1* p_cLens1, int num)
{
	return p_cLens1->Get_colorweight(num);
}

void SetColor(cLens1* p_cLens1, int num, const char* pText)
{
	// std::string str = pText;
	// std::cout << "strの中身を確認する" << std::endl << str << std::endl;
	// p_cLens1->Set_color(num, str);
	p_cLens1->Set_color(num, pText);
}

BSTR GetColor(cLens1* p_cLens1, int num)
{
	std::string str;
	str = p_cLens1->Get_color(num);
	// 100は適当
	wchar_t ws[100]; 
	mbstowcs(ws, str.c_str(), 100);
	return SysAllocString(ws);
}

void MakeSAGraph(cLens1* p_cLens1)
{
	// 一旦は、固定値とする
	int FindPupil = 1;
	double yPupilMax = 1;
	double FullScale = 0.02;
	p_cLens1->MakeSAGraph(FindPupil, yPupilMax, FullScale);
}

void SaveAsBmp(cLens1* p_cLens1)
{
	// CurrentShapeの0は球面収差を意味する
	p_cLens1->CurrentShapes(0)->SaveAsBmp("SAGraph_test_dll3.bmp", 400);
}

void MakeLensView(cLens1* p_cLens1)
{
	// p_cLens1->MakeSAGraph(FindPupil, yPupilMax, FullScale);
}

void SaveAsBmpLensView(cLens1* p_cLens1)
{
	// CurrentShapeの???はLensViewを意味する
	// p_cLens1->CurrentShapes(??)->SaveAsBmp("SAGraph_test_dll3.bmp", 400);
}