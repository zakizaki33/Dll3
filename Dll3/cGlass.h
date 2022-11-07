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
	double *wl;  // �g��(nm)

	std::string name,name_ini;  // ��F"-BK7 50" -> d�����ܗ�=-1.51680+0.00050
	double nd,nd_ini;           // d�����ܗ��̐�Βl
	double nud,nud_ini;         // �A�b�x����d
	double dn,dn_ini;           // ���ܗ��̌덷
	
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
	bool coordinated;  // �������̂���

public:
	cGlass();
	cGlass(int cn);
	cGlass(const cGlass& x);
	virtual ~cGlass();
	cGlass& operator=(const cGlass& x);
	friend bool operator==(cGlass a,cGlass b);

	void set_wl(int j,double val);
	void set_cn_wl(int cn,double *pwl);

	// �C���^�[�t�F�[�X
	static bool UsePointer;  // �����o���ȉ��̃C���^�[�t�F�[�X�������
	                         // �|�C���^�Ȃǂő��삷��Ƃ��i�����݌v�Ȃǁj��true�ɂ���D
	                         // �������Ȃ��ƃ����o�Ԃ̐��������Ȃ��Ȃ�D
	                         // (true�ɂ���Ǝ��s���x���x���Ȃ�j
	static int Digits;   // ���ܗ��̏����_�ȉ��L�������D
	                     // �A�b�x���̏����_�ȉ��L�������́CDigits-3 �Ƃ���D
	                     // �����l��5�����C�����݌v���ɂ͑��₳�Ȃ��Ɣ����W���̐��x���o�Ȃ��D
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
