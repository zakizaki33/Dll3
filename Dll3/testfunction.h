#pragma once
#include "complex.h"
#include"Matrix.h"

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
	double** dpmatrix_wrap;
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
// new�̎��ɁA�|�C���^���K�v�B�Ȃ̂ň�����1������
// 
extern "C" MATHLIBRARY_API void New_matrix();
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p);
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p, int m, int n);

// matrix(const matrix<T>&); ??�I���W�i��
// ���2
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p1, const matrix_wrap& );
// ���1
extern "C" MATHLIBRARY_API void New_matrix(matrix_wrap* p1, const matrix_wrap* p2);



extern "C" MATHLIBRARY_API void Delete_matrix();

// ������͗v��Ȃ�
extern "C" MATHLIBRARY_API void matrix_init(matrix_wrap * p, int m, int n);
