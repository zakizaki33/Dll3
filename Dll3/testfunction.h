#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include "complex.h"
#include "Matrix.h"
#include "cLens1.h"
#include <comdef.h>


#ifdef DLL3_EXPORTS
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API __declspec(dllimport)
#endif
// ここ良くわからないので、DLL3_EXPORTSは使わないでおこう

struct complex_wrap {
	double p1;
	double p2;
};

struct matrix_wrap {
	int m = 0;
	int n = 0;
	double** dpmatrix_wrap = nullptr;
};

// 単なるテスト関数
extern "C" MATHLIBRARY_API void test_01();

extern "C" MATHLIBRARY_API double Return123();

// complex用のテスト関数
extern "C" MATHLIBRARY_API complex_wrap* test_001();

extern "C" MATHLIBRARY_API void New_complex();

extern "C" MATHLIBRARY_API void Delete_complex();

extern "C" MATHLIBRARY_API void ABS(complex_wrap*);

extern "C" MATHLIBRARY_API double ARG(complex_wrap*);

extern "C" MATHLIBRARY_API void Conjugate(complex_wrap*, complex_wrap*);

// ここからはmatrixに関する関数
// newの時に、ポインタが必要。それぞれ引数が(python で言うself)1つ増える
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p);
extern "C" MATHLIBRARY_API void New_matrixMN(matrix_wrap* p, int m, int n);
extern "C" MATHLIBRARY_API void New_matrixCopy(matrix_wrap* p1, const matrix_wrap * p2);

extern "C" MATHLIBRARY_API void Delete_matrix(matrix_wrap * p);

extern "C" MATHLIBRARY_API void Matrix_Set_Value(matrix_wrap * p, int m, int n, double num);

/*
// matrix(const matrix<T>&); オリジナル
// 候補1
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p1, const matrix_wrap* p2);
// 候補2
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p1, const matrix_wrap& );
*/

// ↓これは要らない
/*
extern "C" MATHLIBRARY_API void matrix_init(matrix_wrap * p, int m, int n);
*/

// Sammple_26_1_copy\ConsoleApp1の書き方を参考にする
extern "C" MATHLIBRARY_API void matrix_init(matrix<double> * p, int m, int n);

extern "C" MATHLIBRARY_API void matrix_inv(matrix<double> * p1, matrix<double> * p2);

// https://www.youtube.com/watch?v=xkfQFYUuE4o&list=LL&index=22
// ここの動画の真似
extern "C" MATHLIBRARY_API matrix<double>* CreateMatrix();

extern "C" MATHLIBRARY_API matrix<double> * CreateMatrixMN(int m, int n);

extern "C" MATHLIBRARY_API void SetMatrixIJ(matrix<double>* pMatrix,int m, int n, double value);

extern "C" MATHLIBRARY_API double GetMatrixIJ(matrix<double> * pMatrix, int m, int n);

extern "C" MATHLIBRARY_API void DeleteMatrix(matrix<double>* pMatrix);

// ここからはいよいよcLens1を利用できるようにする
extern "C" MATHLIBRARY_API cLens1* Create_cLens1();

extern "C" MATHLIBRARY_API void Delete_cLens1(cLens1* p_cLens1);

extern "C" MATHLIBRARY_API void SetRadius(cLens1 * p_cLens1, int surf_i, double value);

extern "C" MATHLIBRARY_API double GetRadius(cLens1 * p_cLens1, int surf_i);

extern "C" MATHLIBRARY_API void SetDistance(cLens1 * p_cLens1, int surf_i, double value);

extern "C" MATHLIBRARY_API double GetDistance(cLens1 * p_cLens1, int surf_i);

extern "C" MATHLIBRARY_API void SetGlassName(cLens1* p_cLens1, int surf_i, const char* pText);

extern "C" MATHLIBRARY_API BSTR GetGlassName(cLens1 * p_cLens1, int surf_i);

extern "C" MATHLIBRARY_API int GetK(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API double focallength(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API double backf(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void SetStop(cLens1 * p_cLens1, int surf_i);

extern "C" MATHLIBRARY_API int GetStop(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void SetEAy(cLens1 * p_cLens1, int surf_i, double value);

extern "C" MATHLIBRARY_API double GetEAy(cLens1 * p_cLens1, int surf_i);

extern "C" MATHLIBRARY_API void EPCalculation(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void Set_s(cLens1 * p_cLens1,double value);

extern "C" MATHLIBRARY_API double Get_s(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void Set_t(cLens1 * p_cLens1, double value);

extern "C" MATHLIBRARY_API double Get_t(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void Set_EPD(cLens1 * p_cLens1, double value);

extern "C" MATHLIBRARY_API double Get_EPD(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void SetColorN(cLens1 * p_cLens1, int num);

extern "C" MATHLIBRARY_API int GetColorN(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void SetColorWeight(cLens1 * p_cLens1, int num, double value);

extern "C" MATHLIBRARY_API double GetColorWeight(cLens1 * p_cLens1, int num);

extern "C" MATHLIBRARY_API void SetColor(cLens1 * p_cLens1, int num, const char* pText);

extern "C" MATHLIBRARY_API BSTR GetColor(cLens1 * p_cLens1, int num);

// 画像の保存を実装トライする
extern "C" MATHLIBRARY_API void MakeSAGraph(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void SaveAsBmp(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void MakeLensView(cLens1 * p_cLens1);

extern "C" MATHLIBRARY_API void SaveAsBmpLensView(cLens1 * p_cLens1);

// TODO vector.h のreturnを修正した　2023-05-26