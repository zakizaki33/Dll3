#include <iostream>
#include <stdio.h>
#include <stdlib.h>
// HeaderFileの読み込み
#include "../Dll3/testfunction.h"
// libFileの読み込み
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
    std::cout << "c3の共役複素数(虚数)は　=" << c4.p2 << "です" << std::endl << std::endl;

    // matrixの練習はここから
    std::cout << "matrixの練習はここから"<< std::endl;
    matrix_wrap* mat1 = new matrix_wrap();
    New_matrixMN(mat1, 2, 2);
            
    std::cout << "mat1のmの中身" << mat1->m << std::endl;
    std::cout << "mat1のnの中身" << mat1->n << std::endl;

    // matrixの練習(書き方変更：直接matrix<double>を使う）
    std::cout << "matrixの練習その②" << std::endl;
    matrix<double>*  pmat2 =new matrix<double>();
    matrix_init(pmat2, 2, 2);

    std::cout << "mat2のmの行の数" << pmat2->rows() << std::endl;
    std::cout << "mat2のmの中身 a[1][1]" << pmat2->a[1][1] << std::endl;
    std::cout << "mat2のmの中身 a[1][2]" << pmat2->a[1][2] << std::endl;
    std::cout << "mat2のmの中身 a[2][1]" << pmat2->a[2][1] << std::endl;
    std::cout << "mat2のmの中身 a[2][2]" << pmat2->a[2][2] << std::endl;

    std::cout << "逆行列を求める" << std::endl;
    matrix<double>* pmat2_inv = new matrix<double>();
    matrix_init(pmat2_inv, 2, 2);
    matrix_inv(pmat2, pmat2_inv);
    std::cout << "mat2_invのmの中身 a[1][1]" << pmat2_inv->a[1][1] << std::endl;
    std::cout << "mat2_invのmの中身 a[1][2]" << pmat2_inv->a[1][2] << std::endl;
    std::cout << "mat2_invのmの中身 a[2][1]" << pmat2_inv->a[2][1] << std::endl;
    std::cout << "mat2_invのmの中身 a[2][2]" << pmat2_inv->a[2][2] << std::endl;

    matrix<double>* pmat3 = CreateMatrixMN(2, 2);
    // ここでrowsが使えるのは、ヘッダーファイルに関数が実装されているから？
    // C#ではそのまま使えないはず。。。
    std::cout << "mat3のmの行の数" << pmat3->rows() << std::endl;
    SetMatrixIJ(pmat3, 1, 1, 1);
    SetMatrixIJ(pmat3, 1, 2, 2);
    SetMatrixIJ(pmat3, 2, 1, 3);
    SetMatrixIJ(pmat3, 2, 2, 7);
    std::cout << "mat3のmの中身 a[1][1]" << pmat3->a[1][1] << std::endl;
    std::cout << "mat3のmの中身 a[1][2]" << pmat3->a[1][2] << std::endl;
    std::cout << "mat3のmの中身 a[2][1]" << pmat3->a[2][1] << std::endl;
    std::cout << "mat3のmの中身 a[2][2]" << pmat3->a[2][2] << std::endl;

    matrix<double>* pmat3_inv = CreateMatrixMN(2, 2);
    matrix_inv(pmat3, pmat3_inv);
    std::cout << "mat3_invのmの中身 a[1][1]" << pmat3_inv->a[1][1] << std::endl;
    std::cout << "mat3_invのmの中身 a[1][2]" << pmat3_inv->a[1][2] << std::endl;
    std::cout << "mat3_invのmの中身 a[2][1]" << pmat3_inv->a[2][1] << std::endl;
    std::cout << "mat3_invのmの中身 a[2][2]" << pmat3_inv->a[2][2] << std::endl;

    // メモリーリークのテスト
    // http://c-lang.sevendays-study.com/column-28.html
    int i_;
    char* p_;
    for (i_ = 0; i_ < 1000; i_++) {
        //  メモリを生成（開放しない）
        p_ = (char*)malloc(sizeof(char) * 1000);
    }
    std::cout << "メモリの様子チェック" << std::endl;


    /*
    const int j = 10000;
    matrix<double>* array_pmat[j];

    std::cout << ((size_t)sizeof(j)) << std::endl;
    std::cout << ((size_t)sizeof(matrix<double>)) << std::endl;
        
    // メモリーリークのテスト
    for (int i=0 ; i<j ; i++) {
    matrix<double>* pmat = new matrix<double>();
    array_pmat[i] = pmat;
    }

    // delete[] array_pmat;

    std::cout << "メモリの様子チェック" << std::endl;

    // https://brain.cc.kogakuin.ac.jp/~kanamaru/lecture/prog1/11-02.html
    matrix<double>* p = new matrix<double>[j];

    std::cout << "メモリの様子チェック" << std::endl;
    delete[] array_pmat;
    delete[] p;
    std::cout << "メモリの様子チェック" << std::endl;
    */

    // ファイル読み書き
    //  https://stackoverflow.com/questions/6404586/double-to-const-char
    std::ofstream file("output.txt");

    complex_wrap c11;
    c11.p1 = 12345;
    c11.p2 = 24678;
    
    // coutと併記してみる
    std::cout << "c11のp1の値は=" << c11.p1 << std::endl;
    file << "c11のp1の値は=" << c11.p1 << std::endl;
    std::cout << "c11のp2の値は=" << c11.p2 << std::endl;
    file << "c11のp2の値は=" << c11.p2 << std::endl;

    // お片付けを忘れない
    file.close(); 

}
