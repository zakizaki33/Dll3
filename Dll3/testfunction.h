#pragma once
#include "complex.h"
#include "Matrix.h"
#include "cLens1.h"

#ifdef DLL3_EXPORTS
#define MATHLIBRARY_API __declspec(dllexport)
#else
#define MATHLIBRARY_API __declspec(dllimport)
#endif
// �����ǂ��킩��Ȃ��̂ŁADLL3_EXPORTS�͎g��Ȃ��ł�����

struct complex_wrap {
	double p1;
	double p2;
};

struct matrix_wrap {
	int m = 0;
	int n = 0;
	double** dpmatrix_wrap = nullptr;
};

// �P�Ȃ�e�X�g�֐�
extern "C" MATHLIBRARY_API void test_01();
// complex�p�̃e�X�g�֐�
extern "C" MATHLIBRARY_API complex_wrap* test_001();

extern "C" MATHLIBRARY_API void New_complex();

extern "C" MATHLIBRARY_API void Delete_complex();

extern "C" MATHLIBRARY_API void ABS(complex_wrap*);

extern "C" MATHLIBRARY_API double ARG(complex_wrap*);

extern "C" MATHLIBRARY_API void Conjugate(complex_wrap*, complex_wrap*);

// ���������matrix�Ɋւ���֐�
// new�̎��ɁA�|�C���^���K�v�B���ꂼ�������(python �Ō���self)1������
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p);
extern "C" MATHLIBRARY_API void New_matrixMN(matrix_wrap* p, int m, int n);
extern "C" MATHLIBRARY_API void New_matrixCopy(matrix_wrap* p1, const matrix_wrap * p2);

extern "C" MATHLIBRARY_API void Delete_matrix(matrix_wrap * p);

extern "C" MATHLIBRARY_API void Matrix_Set_Value(matrix_wrap * p, int m, int n, double num);

/*
// matrix(const matrix<T>&); �I���W�i��
// ���1
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p1, const matrix_wrap* p2);
// ���2
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p1, const matrix_wrap& );
*/

// ������͗v��Ȃ�
/*
extern "C" MATHLIBRARY_API void matrix_init(matrix_wrap * p, int m, int n);
*/

// Sammple_26_1_copy\ConsoleApp1�̏��������Q�l�ɂ���
extern "C" MATHLIBRARY_API void matrix_init(matrix<double> * p, int m, int n);

extern "C" MATHLIBRARY_API void matrix_inv(matrix<double> * p1, matrix<double> * p2);

// https://www.youtube.com/watch?v=xkfQFYUuE4o&list=LL&index=22
// �����̓���̐^��
extern "C" MATHLIBRARY_API matrix<double>* CreateMatrix();

extern "C" MATHLIBRARY_API matrix<double> * CreateMatrixMN(int m, int n);

extern "C" MATHLIBRARY_API void SetMatrixIJ(matrix<double>* pMatrix,int m, int n, double value);

extern "C" MATHLIBRARY_API double GetMatrixIJ(matrix<double> * pMatrix, int m, int n);

extern "C" MATHLIBRARY_API void DeleteMatrix(matrix<double>* pMatrix);

// ��������͂��悢��cLens1�𗘗p�ł���悤�ɂ���
extern "C" MATHLIBRARY_API cLens1* Create_cLens1();

extern "C" MATHLIBRARY_API void Delete_cLens1(cLens1* p_cLens1);

extern "C" MATHLIBRARY_API void SetRadius(cLens1 * p_cLens1, int surf_i, double value);

extern "C" MATHLIBRARY_API double GetRadius(cLens1 * p_cLens1, int surf_i);

extern "C" MATHLIBRARY_API void SetDistance(cLens1 * p_cLens1, int surf_i, double value);

extern "C" MATHLIBRARY_API void SetGlassName(cLens1* p_cLens1, int surf_i, const char* pText);

extern "C" MATHLIBRARY_API std::string GetGlassName(cLens1 * p_cLens1, int surf_i);

extern "C" MATHLIBRARY_API double focallength(cLens1 * p_cLens1);