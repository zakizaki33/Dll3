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