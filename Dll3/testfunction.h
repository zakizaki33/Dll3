#pragma once
#include "complex.h"

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

// �P�Ȃ�e�X�g�֐�
extern "C" MATHLIBRARY_API void test_01();
// complex�p�̃e�X�g�֐�
extern "C" MATHLIBRARY_API complex_wrap* test_001();

extern "C" MATHLIBRARY_API void New_complex();

extern "C" MATHLIBRARY_API void Delete_complex();

extern "C" MATHLIBRARY_API void ABS(complex_wrap*);

extern "C" MATHLIBRARY_API double ARG(complex_wrap*);

extern "C" MATHLIBRARY_API void Conjugate(complex_wrap*, complex_wrap*);
