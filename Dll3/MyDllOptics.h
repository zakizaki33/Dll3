// MyDllOptics.h : MYDLLOPTICS �A�v���P�[�V�����̃��C�� �w�b�_�[ �t�@�C���ł��B
//

#if !defined(AFX_MYDLLOPTICS_H__AD966568_51EE_11D7_BE64_C06041BB0250__INCLUDED_)
#define AFX_MYDLLOPTICS_H__AD966568_51EE_11D7_BE64_C06041BB0250__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifndef __AFXWIN_H__
	#error include 'stdafx.h' before including this file for PCH
#endif

#include "resource.h"		// ���C�� �V���{��

/////////////////////////////////////////////////////////////////////////////
// CMyDllOpticsApp
// ���̃N���X�̓���̒�`�Ɋւ��Ă� MyDllOptics.cpp �t�@�C�����Q�Ƃ��Ă��������B
//

class CMyDllOpticsApp : public CWinApp
{
public:
	CMyDllOpticsApp();

// �I�[�o�[���C�h
	// ClassWizard �͉��z�֐��̃I�[�o�[���C�h�𐶐����܂��B
	//{{AFX_VIRTUAL(CMyDllOpticsApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

	//{{AFX_MSG(CMyDllOpticsApp)
		// ���� -  ClassWizard �͂��̈ʒu�Ƀ����o�֐���ǉ��܂��͍폜���܂��B
		//         ���̈ʒu�ɐ��������R�[�h��ҏW���Ȃ��ł��������B
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ �͑O�s�̒��O�ɒǉ��̐錾��}�����܂��B

#endif // !defined(AFX_MYDLLOPTICS_H__AD966568_51EE_11D7_BE64_C06041BB0250__INCLUDED_)
