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
	printf("\nHello Dll3!!! 2022-10-24\n\n");
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