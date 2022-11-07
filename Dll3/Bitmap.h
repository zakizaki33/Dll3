#if !defined(AFX_BITMAP_H__E6EE433A_CFC4_4222_813E_0FF580444148__INCLUDED_)
#define AFX_BITMAP_H__E6EE433A_CFC4_4222_813E_0FF580444148__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
// Bitmap.h : �w�b�_�[ �t�@�C��
//
#include "cBitmap.h"

/////////////////////////////////////////////////////////////////////////////
// Bitmap �R�}���h �^�[�Q�b�g

class Bitmap : public CCmdTarget
{
	DECLARE_DYNCREATE(Bitmap)

	Bitmap();           // ���I�����Ɏg�p�����v���e�N�g �R���X�g���N�^

	cBitmap X;

// �A�g���r���[�g
public:

// �I�y���[�V����
public:

// �I�[�o�[���C�h
	// ClassWizard �͉��z�֐��̃I�[�o�[���C�h�𐶐����܂��B
	//{{AFX_VIRTUAL(Bitmap)
	public:
	virtual void OnFinalRelease();
	//}}AFX_VIRTUAL

// �C���v�������e�[�V����
protected:
	virtual ~Bitmap();

	// �������ꂽ���b�Z�[�W �}�b�v�֐�
	//{{AFX_MSG(Bitmap)
		// ���� - ClassWizard �͂��̈ʒu�Ƀ����o�֐���ǉ��܂��͍폜���܂��B
	//}}AFX_MSG

	DECLARE_MESSAGE_MAP()
	DECLARE_OLECREATE(Bitmap)

	// �������ꂽ OLE �f�B�X�p�b�` �}�b�v�֐�
	//{{AFX_DISPATCH(Bitmap)
	afx_msg long GetM();
	afx_msg void SetM(long nNewValue);
	afx_msg long GetN();
	afx_msg void SetN(long nNewValue);
	afx_msg double GetGamma();
	afx_msg void SetGamma(double newValue);
	afx_msg void GammaCompensation(double gamma);
	afx_msg void ToNeg();
	afx_msg long Save(LPCTSTR filename, double gamma);
	afx_msg long Open(LPCTSTR filename);
	afx_msg long GetR(long i, long j);
	afx_msg long GetG(long i, long j);
	afx_msg long GetB(long i, long j);
	afx_msg void Expand();
	afx_msg void Resize(long m, long n);
	afx_msg void ToR();
	afx_msg void ToG();
	afx_msg void ToB();
	afx_msg long Merge(LPCTSTR filename_R, LPCTSTR filename_G, LPCTSTR filename_B);
	afx_msg void SetLine(double x1, double y1, double x2, double y2, long R, long G, long B);
	afx_msg void ToWhite();
	afx_msg void SetCircle(double x0, double y0, double radius, long R, long G, long B, long fill);
	afx_msg long Height();
	afx_msg long Width();
	afx_msg void Zoom(long i1, long j1, long i2, long j2);
	afx_msg void ToBuf();
	afx_msg void FromBuf();
	afx_msg void ToBuf2();
	afx_msg void FromBuf2();
	afx_msg void SeeColorShift(long i, long j);
	afx_msg void Median();
	afx_msg void Edge();
	afx_msg void Binarize(double threshold);
	afx_msg void Trim(long i0, long j0, long height, long width);
	afx_msg void TrimCircle();
	afx_msg void MarkSaturation();
	afx_msg void DistCorrect(double dist);
	afx_msg void ColorShiftCorrect(LPCTSTR filename, long inv);
	afx_msg long GetRGB(long i, long j);
	afx_msg void SetRGB(long i, long j, long nNewValue);
	//}}AFX_DISPATCH
	DECLARE_DISPATCH_MAP()
	DECLARE_INTERFACE_MAP()
};

/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ �͑O�s�̒��O�ɒǉ��̐錾��}�����܂��B

#endif // !defined(AFX_BITMAP_H__E6EE433A_CFC4_4222_813E_0FF580444148__INCLUDED_)
