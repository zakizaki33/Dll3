#pragma once
#include "stdafx.h"
#include "AUserDef.h"

// ここでcLens1のユーザー定義面(成分)のz(x,y)およびその1次,2次導関数を定義する

double UserDefSurfZ(double x,double y){

	return 0;
}

void UserDefSurfDerivative(double &zx,double &zy,double x,double y){
	/*
	int i;
	double r0[10],r,dr,Zr;

	r0[1]=0.5; r0[2]=1.5; r0[3]=3;
	dr=0.02;
	Zr=0.15;	// =3μm/20μm

	r=sqrt(x*x+y*y);

	zx=zy=0;

	for(i=1; i<=3; i++){
		if(r0[i]<r && r<r0[i]+dr){
			zx=Zr*x/r;
			zy=Zr*y/r;
		}
	}
	*/

	zx=zy=0;
	return;
}

void UserDefSurf2ndDerivative(double &zxx,double &zxy,double &zyy,double x,double y){

	zxx=zxy=zyy=0;
	return;
}