#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "testfunction.h"
#include "complex.h"
#include "Matrix.h"

complex_wrap* static_c;
complex* p_complex;

// ���ŏI�I�ɂ́A����͎g��Ȃ��œ���o����悤�ɂ���i��R������j
//  �O���[�o���ϐ����g���Ă̂���肾��1�����ɂȂ�
// 
// ����������ƃ_�u���|�C���^**�ł̂����ɂȂ邩��
//  ���̂�����肳���������߂̃|�C���^(������������)
// 
//�@T **a�@�̃_�u���|�C���^�̎g����(�������z��)�Ƃ͈قȂ�

matrix_wrap* static_m;
matrix<double>* p_matrix;

void test_01()
{
	printf("\nHello Dll3!!! 2022-11-07\n\n");
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
	// complex �e�X�g(1,1);
	// test = abs(�e�X�g);
	test = abs(*p_complex);
	std::cout << "ABS�̊֐��𓮂���" << test << std::endl;
	Delete_complex();
}

double ARG(complex_wrap* p){
	p_complex = new complex(p->p1, p->p2);
	double test;
	test = arg(*p_complex);
	std::cout << "arg�̊֐��𓮂���" << test << std::endl;
	Delete_complex();
	return test;
}

void Conjugate(complex_wrap* comp1, complex_wrap* comp2)
{
	// �C���v�b�g�A�A�E�g�v�b�g�p�̕��f�����p�ӂ���
	complex* p_c1 = new complex(comp1->p1, comp1->p2);
	complex* p_c2 = new complex(comp2->p1, comp2->p2);

	*p_c2 = conj(*p_c1);
	comp2->p1 = Re(*p_c2);
	comp2->p2 = Im(*p_c2);

	std::cout << "����m�F�p���b�Z�[�W 2022-11-14 " << std::endl;

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
	std::cout << "matrix_init�Ăяo������"<< std::endl;
	p->m = m;
	p->n = n;

}
*/

// �s��̏�����
void matrix_init(matrix<double>* p, int m, int n) {
	//p_matrix = new matrix<double>(m, n);
	// p_matrix->a[1][1] = 1.1;
	// p_matrix->a[0][1] = 2.2;
	// p_matrix->a[1][0] = 3.3;
	// p_matrix->a[1][1] = 4.4;
	// p = p_matrix;

	// p = new matrix<double>(m, n);

	// �v�fm,n�ŏ�������
	p->redim(m, n);

	p->a[1][1] = 1.1;
    p->a[1][2] = 2.2;
    p->a[2][1] = 3.3;
    p->a[2][2] = 4.4;
	
}