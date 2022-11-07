	#if !defined(AFX_CLEASTSQUARES_H__F0CE8FC0_EE37_11D7_BE64_9FFDC1C8AE50__INCLUDED_)
#define AFX_CLEASTSQUARES_H__F0CE8FC0_EE37_11D7_BE64_9FFDC1C8AE50__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "matrix.h"
#include "list.h"

class cLeastSquares
{

public:
	list<double> x1,x2,x3,x4,w;
	matrix<double> A,dX,F;
	double **X;       // X�͔h���N���X�̕ϐ����Q�Ƃ��邽��matrix�ł͂Ȃ��|�C���^�̔z��Ƃ���
	                  // X[]��ModifyCoefficients()�̒��ł̂݃��������m�ۂ��Ă���̂ŁC
	                  // ���̑��̕�������Q�Ƃ��Ȃ��悤���ӂ��邱��
	cLeastSquares();
	virtual ~cLeastSquares();
	int  GetNumberOfData() const;
	void SetNumberOfData(int m);
	double Get_x1(int i) const; void Set_x1(int i,double value);
	double Get_x2(int i) const; void Set_x2(int i,double value);
	double Get_x3(int i) const; void Set_x3(int i,double value);
	double Get_x4(int i) const; void Set_x4(int i,double value);
	double GetWeight(int i) const;
	void SetWeight(int i,double value);
	void SetWeight(double value);
	void AddData(double x1,double x2,double x3,double x4);
	void AddData(double x1,double x2,double x3);
	void AddData(double x1,double x2);
	void DataClear();
	double Error(int i);
	double rmsError();
	double pvError();
	double maxError();
	int ModifyCoefficients(double rho);

protected:
	virtual double f(double x1,double x2,double x3,double x4)=0;
	virtual void SetA()=0;
	virtual void SetX()=0;
public:
	virtual int NumberOfCoefficients()=0;
};


// ��F�Q���֐����Ă͂߂̏ꍇ( y=ax^2+bx+c, z�͕s�v)

class cFittingQuadratic : public cLeastSquares
{
	list<double> &x,&y;  // �R�[�h�����₷�����邽�߂� x=x1, y=x2 �Ƃ���D
	double a,b,c;

	cFittingQuadratic() : x(x1),y(x2) {};
	// �y���Ӂz�R���X�g���N�^�����̂悤�ɏ����ƃG���[ ////////////////////////////
	/*
	cFittingQuadratic()
	{
		x=x1;   // �Q�ƃ����o�̏�������{}���ɏ����ƃG���[�ɂȂ�D
		y=x2;
	};
	*/

	double f(double x,double y,double z,double dummy){
		return a*x*x +b*x +c -y;   // f(x,y,z)=0 �ɂ��Ă͂߂�
	};

	void SetA(){
		for(int i=1; i<=GetNumberOfData(); ++i){
			A[i][1]=x[i]*x[i];
			A[i][2]=x[i];
			A[i][3]=1;
		}
	};

	void SetX(){
		X[1]=&a;
		X[2]=&b;
		X[3]=&c;
	};

	int NumberOfCoefficients(){
		return 3;
	};
};


#endif // !defined(AFX_CLEASTSQUARES_H__F0CE8FC0_EE37_11D7_BE64_9FFDC1C8AE50__INCLUDED_)
