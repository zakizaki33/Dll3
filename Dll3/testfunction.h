#pragma once
#include "complex.h"
#include "Matrix.h"

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