#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream>
#include "../S04/s04.h"
//#include "../S04/Matrix.h"
//#include "../S04/list.h"
//#include "../S04/cLeastSquares.h"
//#include "../S04/cLeneq.h"
#include "cLens1.h"

// <Tips1>
// Visual Studio で作成する場合、DLLの同じソリューションにもう一つのテスト用のプロジェクトを設けておくとお手軽である
// まず、CUIで良いのでプロジェクトを作成し、ソリューションエクスプローラーから、このプロジェクトのツリーにある
// 参照を右クリックすると「参照の追加」が表示される。
// この参照の追加に、利用するDLLのプロジェクト（今回はS04）を追加しておけば、ややこしいディレクトリとの関連付けはおおよそ不要になる

int main()
{

	// cLens1 lens(2, 3);

	cLens1 lens;

	lens.Set_r(1,100);
	lens.Set_r(2, -100);
	lens.Set_d(1, 10);
	// lens.Set_gname(1, "S-BSL7"); // ガラス面はi-1の値を入れる
	// e-line での設定になっている
	lens.Set_gname(1, "518640"); // ガラスコードでも入れられる
	
	// レンズの有効径の設定　EAy(i)を使う
	// Excel上では　8列目EAx(i), 9列目EAType(i)が指定できる仕様？
	lens.Set_EAy(1, 1.999);
	lens.Set_EAy(2, 1.999);

	lens.Set_stop(1);   // 絞り面の設定
	lens.EPCalc();
	lens.Set_EPD(1.0);// EPD　の設定

	std::cout << "GetEPD" << std::endl;
	std::cout << lens.Get_EPD() << std::endl;


	// 物体距離を指定しているつもり
	lens.Set_s(1000000);

	// 波長(e-line)を指定しているつもり
	lens.Set_cn(1); //波長の数
	lens.Set_color(1, "e");
	lens.Set_colorweight(1, 1);


	double r1;
	r1=lens.Get_r(1);
	std::cout << "Show r1" << std::endl;
	std::cout << r1 << std::endl;

	double focal_length;
	focal_length = lens.f();
	std::cout << focal_length << std::endl;
	
	double bf;
	bf = lens.bf();
	std::cout << bf << std::endl;
	

	
	// ExcelVBAからRを代入したときにこのGet_rにつながっているかを調べる


	// MTFsmaveを使ってみる
	//　構文：
	//  double MTFsmave(double nu,double yObj,double xObj);
	// とりあえず、適当に入れてみる

	double MTF_value;
	// 10LPMM, Y=0,X=0のつもり
	MTF_value = lens.MTFsmave(10,0,0);
    // 値が出てこない、ただしく計算させられていない
	// 何がただしくセット出来ていないのか？
	std::cout << "MTFの値は　＝　" << MTF_value << "です　" << std::endl;

	// こちらで入れてみる
	// double MTFs(double nu, double yObj, double xObj, double defocus,
	// 	int IsGeometrical, int AberrationFree, int ColorStart, int ColorEnd, int FindPupil = 0);
	MTF_value = lens.MTFsmave(10, 0, 0, 0, 0, 0, 1, 1, 1);
	// 値が出てこない、ただしく計算させられていない
	// 何がただしくセット出来ていないのか？
	std::cout << "MTFの値は　＝　" << MTF_value << "です　" << std::endl;

	// ダメもとでMTFyでも計算させてみる
	MTF_value = lens.MTFy(10.0, 0.0, 0.0);
	std::cout << "MTFの値は　＝　" << MTF_value << "です　" << std::endl;
	// やっぱりだめだ

	// 気分を変えて、RayPosYの練習
	double pos_y;
	pos_y = lens.RayPosY(2,1,0,"xxx",1,0,1);
	std::cout << "yの値は　＝　" << pos_y << "です　" << std::endl;

	// 球面収差を計算させる。LSA()を使っていく
	double val_LSA;
	val_LSA = lens.LSA();
	std::cout << "球面収差10割の値は　＝　" << val_LSA << "です　" << std::endl;
	std::cout << "球面収差 9割の値は　＝　" << lens.LSA(0.9, 1, 1) << "です　" << std::endl;
	std::cout << "球面収差 7割の値は　＝　" << lens.LSA70() << "です　" << std::endl;
	std::cout << "球面収差 5割の値は　＝　" << lens.LSA50() << "です　" << std::endl;
	std::cout << "球面収差 3割の値は　＝　" << lens.LSA(0.3, 1, 1) << "です　" << std::endl;
	std::cout << "球面収差 1割の値は　＝　" << lens.LSA(0.1, 1, 1) << "です　" << std::endl;
	//次回聞くこと
	// 球面収差10,7,5割のあたいをZEMAXとくらべてどうなるか？
	//　グラフの書き方を理解するにはどういうことを理解すればよいか？

	// 図形を描く練習をしていく　2022-08-08
	// 利用するclassはcShapes
	// とりあえず、インスタンスを生成出来るか確認
	cShapes shape1;
	cShapes* shape2=new cShapes();

	shape1.AddLine(50, 100, 200, 400, rgb(255, 0, 255), DASH );
	shape1.AddCircle(50, 100, 100, 180, 0L, 0L);  // なぜか出てこない
	// shape1.AddText("text output", 300, 300, 5, 0); //　なぜか出てこない
	shape1.AddLine(25, 200, 100, 100, 1L, 1L);
	shape1.AddLine(25, 200, 250, 200, 0L, 1L);

	// shape1.AllView(); // 何使う？

	shape1.SaveAsBmp("test0808.bmp", 600);
	std::cout << "一つ目のLineのX座標は　＝　" << shape1.X1(1) << "です　" << std::endl;
	std::cout << "2つ目のLineのX座標は　＝　" << shape1.X1(2) << "です　" << std::endl;

	lens.MakeLensView("",0);  // これは一応動く。どうこれを利用する？


	lens.CurrentShapes(2)->SaveAsBmp("test0829.bmp", 600);

	// Excelプログラムと同じように入れてみる
	lens.MakeSAGraph(1, 1, 0.02);
	// SA_GRAPHはゼロと定義されている様子
	// なので、
		lens.CurrentShapes(lens.SAGraph())->SaveAsBmp("test0926.bmp", 800);
	
}
