// USING_Dll3.cpp : このファイルには 'main' 関数が含まれています。プログラム実行の開始と終了がそこで行われます。
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

    // dllの中のcomplex_wrap が使えるようになった
    // そういう物なの？
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
    
    /*
    complex c3 = complex(1,2);
    std::cout << "c3の共役複素数の絶対値は=" << abs(conj(c3)) << std::endl;
    */

}
