// Bitmap.cpp : インプリメンテーション ファイル
//

#include "stdafx.h"
#include "MyDllOptics.h"
#include "Bitmap.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// Bitmap

IMPLEMENT_DYNCREATE(Bitmap, CCmdTarget)

Bitmap::Bitmap()
{
	EnableAutomation();
	
	// OLE オートメーション オブジェクトがアクティブな間はアプリケーションの実行を
	// 終了させないために、コンストラクタから AfxOleLockApp を呼び出します。
	
	AfxOleLockApp();
}

Bitmap::~Bitmap()
{
	// すべてのオブジェクトが OLE オートメーションで作成されている時に
	// アプリケーションを終了するために、デストラクタは AfxOleUnlockApp を呼び出します。
	
	AfxOleUnlockApp();
}


void Bitmap::OnFinalRelease()
{
	// オートメーション オブジェクトの最後の参照が解放された時に
	// OnFinalRelease が呼び出されます。基本クラスは自動的にオブジェクト 
	// を破棄します。オブジェクトをメモリから破棄する前に必要な
	// 終了処理があれば追加して下さい。

	CCmdTarget::OnFinalRelease();
}


BEGIN_MESSAGE_MAP(Bitmap, CCmdTarget)
	//{{AFX_MSG_MAP(Bitmap)
		// メモ - ClassWizard はこの位置にマッピング用のマクロを追加または削除します。
	//}}AFX_MSG_MAP
END_MESSAGE_MAP()

BEGIN_DISPATCH_MAP(Bitmap, CCmdTarget)
	//{{AFX_DISPATCH_MAP(Bitmap)
	DISP_PROPERTY_EX(Bitmap, "M", GetM, SetM, VT_I4)
	DISP_PROPERTY_EX(Bitmap, "N", GetN, SetN, VT_I4)
	DISP_PROPERTY_EX(Bitmap, "Gamma", GetGamma, SetGamma, VT_R8)
	DISP_FUNCTION(Bitmap, "GammaCompensation", GammaCompensation, VT_EMPTY, VTS_R8)
	DISP_FUNCTION(Bitmap, "ToNeg", ToNeg, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "Save", Save, VT_I4, VTS_BSTR VTS_R8)
	DISP_FUNCTION(Bitmap, "Open", Open, VT_I4, VTS_BSTR)
	DISP_FUNCTION(Bitmap, "GetR", GetR, VT_I4, VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "GetG", GetG, VT_I4, VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "GetB", GetB, VT_I4, VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "Expand", Expand, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "Resize", Resize, VT_EMPTY, VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "ToR", ToR, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "ToG", ToG, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "ToB", ToB, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "Merge", Merge, VT_I4, VTS_BSTR VTS_BSTR VTS_BSTR)
	DISP_FUNCTION(Bitmap, "SetLine", SetLine, VT_EMPTY, VTS_R8 VTS_R8 VTS_R8 VTS_R8 VTS_I4 VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "ToWhite", ToWhite, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "SetCircle", SetCircle, VT_EMPTY, VTS_R8 VTS_R8 VTS_R8 VTS_I4 VTS_I4 VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "Height", Height, VT_I4, VTS_NONE)
	DISP_FUNCTION(Bitmap, "Width", Width, VT_I4, VTS_NONE)
	DISP_FUNCTION(Bitmap, "Zoom", Zoom, VT_EMPTY, VTS_I4 VTS_I4 VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "ToBuf", ToBuf, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "FromBuf", FromBuf, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "ToBuf2", ToBuf2, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "FromBuf2", FromBuf2, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "SeeColorShift", SeeColorShift, VT_EMPTY, VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "Median", Median, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "Edge", Edge, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "Binarize", Binarize, VT_EMPTY, VTS_R8)
	DISP_FUNCTION(Bitmap, "Trim", Trim, VT_EMPTY, VTS_I4 VTS_I4 VTS_I4 VTS_I4)
	DISP_FUNCTION(Bitmap, "TrimCircle", TrimCircle, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "MarkSaturation", MarkSaturation, VT_EMPTY, VTS_NONE)
	DISP_FUNCTION(Bitmap, "DistCorrect", DistCorrect, VT_EMPTY, VTS_R8)
	DISP_FUNCTION(Bitmap, "ColorShiftCorrect", ColorShiftCorrect, VT_EMPTY, VTS_BSTR VTS_I4)
	DISP_PROPERTY_PARAM(Bitmap, "RGB", GetRGB, SetRGB, VT_I4, VTS_I4 VTS_I4)
	//}}AFX_DISPATCH_MAP
END_DISPATCH_MAP()

// メモ: IID_IBitmap に対して VBA からタイプセーフ バインディングをサポートするた
// めのサポートを追加します。この IID は .ODL ファイルのディスプインターフェイ
// スにアタッチされている GUID とマッチしなければなりません。

// {0B408802-9A11-4F0F-A428-A7CB6B1E67F9}
static const IID IID_IBitmap =
{ 0xb408802, 0x9a11, 0x4f0f, { 0xa4, 0x28, 0xa7, 0xcb, 0x6b, 0x1e, 0x67, 0xf9 } };

BEGIN_INTERFACE_MAP(Bitmap, CCmdTarget)
	INTERFACE_PART(Bitmap, IID_IBitmap, Dispatch)
END_INTERFACE_MAP()

// {18E55E11-BBE5-4BD2-BADD-13FCDC5FAB38}
IMPLEMENT_OLECREATE(Bitmap, "MyDllOptics.Bitmap", 0x18e55e11, 0xbbe5, 0x4bd2, 0xba, 0xdd, 0x13, 0xfc, 0xdc, 0x5f, 0xab, 0x38)

/////////////////////////////////////////////////////////////////////////////
// Bitmap メッセージ ハンドラ

long Bitmap::GetM() { return X.GetM(); }
void Bitmap::SetM(long nNewValue) { X.SetM(nNewValue); }

long Bitmap::GetN() { return X.GetN(); }
void Bitmap::SetN(long nNewValue) { X.SetN(nNewValue); }

double Bitmap::GetGamma() { return X.GetGamma(); }
void Bitmap::SetGamma(double newValue) { X.SetGamma(newValue); }

long Bitmap::GetRGB(long i, long j) { return X.GetRGB(i,j); }
void Bitmap::SetRGB(long i, long j, long nNewValue) { X.SetRGB(i,j,nNewValue); }

long Bitmap::GetR(long i, long j) { return X.GetR(i,j); }

long Bitmap::GetG(long i, long j) {	return X.GetG(i,j); }

long Bitmap::GetB(long i, long j) {	return X.GetB(i,j); }

void Bitmap::SetLine(double x1, double y1, double x2, double y2, long R, long G, long B){
	X.SetLine(x1,y1,x2,y2,R,G,B);
}

void Bitmap::SetCircle(double x0, double y0, double radius, long R, long G, long B, long fill) {
	X.SetCircle(x0,y0,radius,R,G,B,fill);
}

void Bitmap::ToR() { X.ToR(); }

void Bitmap::ToG() { X.ToG(); }

void Bitmap::ToB() { X.ToB(); }

void Bitmap::ToWhite() { X.ToWhite(); }

void Bitmap::ToNeg() { X.ToNeg(); }

void Bitmap::Expand() { X.Expand(); }

long Bitmap::Height() { return X.Height(); }
long Bitmap::Width()  { return X.Width(); }

void Bitmap::Resize(long m, long n){ X.Resize(m,n); }

void Bitmap::Trim(long i0, long j0, long height, long width) { X.Trim(i0,j0,height,width); }

void Bitmap::TrimCircle() { X.TrimCircle(); }

void Bitmap::GammaCompensation(double gamma){ X.GammaCompensation(gamma); }

void Bitmap::Zoom(long i1, long j1, long i2, long j2) { X.Zoom(i1,j1,i2,j2); }

void Bitmap::SeeColorShift(long i, long j) { X.SeeColorShift(i,j); }

void Bitmap::ColorShiftCorrect(LPCTSTR filename, long inv) { X.ColorShiftCorrect(filename,inv); }

void Bitmap::DistCorrect(double dist) { X.DistCorrect(dist); }

void Bitmap::Edge() { X.Edge(); }

void Bitmap::Median() { X.Median(); }

void Bitmap::Binarize(double threshold) { X.Binarize(threshold); }

void Bitmap::MarkSaturation() { X.MarkSaturation(); } 

long Bitmap::Save(LPCTSTR filename, double gamma){ return X.Save(filename,gamma); }

long Bitmap::Open(LPCTSTR filename){ return X.Open(filename); }

long Bitmap::Merge(LPCTSTR filename_R, LPCTSTR filename_G, LPCTSTR filename_B) {
	return X.Merge(filename_R,filename_G,filename_B);
}

void Bitmap::ToBuf() { X.ToBuf(); }
void Bitmap::FromBuf() { X.FromBuf(); }

void Bitmap::ToBuf2() { X.ToBuf2(); }
void Bitmap::FromBuf2() { X.FromBuf2(); }


