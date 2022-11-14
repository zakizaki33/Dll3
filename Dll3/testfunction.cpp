#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "testfunction.h"
#include "complex.h"

complex_wrap* static_c;
complex* p_complex;

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