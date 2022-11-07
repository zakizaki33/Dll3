#if !defined(AFX_CGLASS_H__6CC95E2E_A786_42F4_9CA0_D729E28AD4C9__INCLUDED_)
#define AFX_CGLASS_H__6CC95E2E_A786_42F4_9CA0_D729E28AD4C9__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <string>
#include "cMaterial.h"
#include "cOptics.h"

class cGlass  
{
	int cn;
	double *wl;  // 波長(nm)

	std::string name,name_ini;  // 例："-BK7 50" -> d線屈折率=-1.51680+0.00050
	double nd,nd_ini;           // d線屈折率の絶対値
	double nud,nud_ini;         // アッベ数νd
	double dn,dn_ini;           // 屈折率の誤差
	
	double *n;
	int is_grin;
	double *ra;
	double grin_phi;

	void alloc();
	void initialize();
	void free();
	void assignment(const cGlass& x);
	void set_index();
	void set_index(int j);
	static std::string values_to_name(double nd,double nud,double dn);
	void coordinate();
	bool coordinated;  // 高速化のため

public:
	cGlass();
	cGlass(int cn);
	cGlass(const cGlass& x);
	virtual ~cGlass();
	cGlass& operator=(const cGlass& x);
	friend bool operator==(cGlass a,cGlass b);

	void set_wl(int j,double val);
	void set_cn_wl(int cn,double *pwl);

	// インターフェース
	static bool UsePointer;  // メンバを以下のインターフェースを介さずに
	                         // ポインタなどで操作するとき（自動設計など）はtrueにする．
	                         // そうしないとメンバ間の整合が取れなくなる．
	                         // (trueにすると実行速度が遅くなる）
	static int Digits;   // 屈折率の小数点以下有効桁数．
	                     // アッベ数の小数点以下有効桁数は，Digits-3 とする．
	                     // 初期値は5だが，自動設計時には増やさないと微分係数の精度が出ない．
	std::string& Name();
	double& Nd();
	double& Nud();
	double& dN();
	double N(int j);
	int IsGlass();
	int IsGrin();
	double rA(int j);
	double GrinPhi();

	static double NNuActual(double Nd,double Nud);
	double NNuActual();
};

#endif // !defined(AFX_CGLASS_H__6CC95E2E_A786_42F4_9CA0_D729E28AD4C9__INCLUDED_)
