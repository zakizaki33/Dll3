#include <iostream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "testfunction.h"
#include "complex.h"
#include "Matrix.h"

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

void New_matrix()
{
	// matrix<double> A1();
	//matrix<double>* A0 = new matrix<double>();
	p_matrix = new matrix<double>();
}

void Delete_matrix()
{
	delete p_matrix;
}

void matrix_init(matrix_wrap* p, int m, int n)
{
	p_matrix = new matrix<double> (m, n);
	std::cout << "matrix_init呼び出し成功"<< std::endl;
	p->m = m;
	p->n = n;

}