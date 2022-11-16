﻿// USING_Dll3.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
// fileの読み込み
#include "../Dll3/testfunction.h"
// #pragma comment(lib, "../x64/Debug/Dll3.lib")

int main()
{
    test_01();

    // testfunction.hをインクルードしているのでcomplex_wrapが使える
    complex_wrap c1;
    c1.p1 = 10;
    c1.p2 = 20;
    std::cout << "c1のp1の値は=" << c1.p1 << std::endl;
    std::cout << "c1のp2の値は=" << c1.p2 << std::endl;

    complex_wrap* c2;
    c2 = test_001();
    std::cout << "c2のp1の値は=" << c2->p1 << std::endl;
    std::cout << "c2のp2の値は=" << c2->p2 << std::endl;

    // complexのABSを計算する

    ABS(c2);
    
    // complexのargを計算する(戻り値がdoubleの場合)
    complex_wrap c3 ;
    c3.p1 = 1;
    c3.p2 = 1;
    double arg_c3 = ARG(&c3);
    std::cout << "c3の角度は　=" << arg_c3 << "rad" << std::endl;

    // complex_wrapのconjugateをGetする
    complex_wrap c4;
    c4.p1 = 0;
    c4.p2 = 0;
    Conjugate(&c3, &c4);
    std::cout << "c3の共役複素数(実数)は　=" << c4.p1 << "です" << std::endl;
    std::cout << "c3の共役複素数(虚数)は　=" << c4.p2 << "です" << std::endl;

    // matrixの練習はここから
    matrix_wrap* mat1 = new matrix_wrap();
    matrix_init(mat1, 2, 2);
    
    std::cout << "mat1の中身" << mat1->m << std::endl;

}
