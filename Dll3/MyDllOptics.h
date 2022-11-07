// MyDllOptics.h : MYDLLOPTICS アプリケーションのメイン ヘッダー ファイルです。
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

#include "resource.h"		// メイン シンボル

/////////////////////////////////////////////////////////////////////////////
// CMyDllOpticsApp
// このクラスの動作の定義に関しては MyDllOptics.cpp ファイルを参照してください。
//

class CMyDllOpticsApp : public CWinApp
{
public:
	CMyDllOpticsApp();

// オーバーライド
	// ClassWizard は仮想関数のオーバーライドを生成します。
	//{{AFX_VIRTUAL(CMyDllOpticsApp)
	public:
	virtual BOOL InitInstance();
	//}}AFX_VIRTUAL

	//{{AFX_MSG(CMyDllOpticsApp)
		// メモ -  ClassWizard はこの位置にメンバ関数を追加または削除します。
		//         この位置に生成されるコードを編集しないでください。
	//}}AFX_MSG
	DECLARE_MESSAGE_MAP()
};


/////////////////////////////////////////////////////////////////////////////

//{{AFX_INSERT_LOCATION}}
// Microsoft Visual C++ は前行の直前に追加の宣言を挿入します。

#endif // !defined(AFX_MYDLLOPTICS_H__AD966568_51EE_11D7_BE64_C06041BB0250__INCLUDED_)
