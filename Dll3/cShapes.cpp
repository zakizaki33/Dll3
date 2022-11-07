#include "stdafx.h"
#include "cShapes.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

void cShapes::initialize() {
	OffsetX=OffsetY=0;
	window=view=rect(0,0,640,480);
	T=matrix<double>(3,3);
	T=unit(T);
}

void cShapes::Add(shape newdata) {
	newdata.offset(OffsetX,OffsetY);
	shapes.AddTail(newdata);
	if( newdata.trim(window) ) shapes_in_window.AddTail(newdata);
}

void cShapes::ReSize(double m){
	// 大きさをm倍にする．
	// Zoom()とは異なる．Zoomではshapeの大きさは変わらず，windowの大きさが変わるのみ．
	int i;
	point p1,p2;

	p1=window.TopLeft()*m;
	p2=window.BottomRight()*m;

	for(i=1; i<=shapes.GetSize(); i++){
		(shapes[i].x1)*=m;
		(shapes[i].x2)*=m;
		(shapes[i].y1)*=m;
		(shapes[i].y2)*=m;
		(shapes[i].fontsize)*=m;
	}

	SetWindow(p1.x,p1.y,p2.x,p2.y);
}

// public members /////////////////////////////////

cShapes::cShapes() {
	initialize();
}

void cShapes::SetWindow(double x1,double y1,double x2,double y2){
	int i; shape s;
	point p1=point(x1,y1); point p2=point(x2,y2);
	if(p1!=p2) {
		window=rect(p1,p2);
		shapes_in_window.RemoveAll();
		for(i=1; i<=shapes.GetSize(); ++i){
			shapes.GetData(s,i);
			if(s.trim(window)) shapes_in_window.AddTail(s);
		}
	}
}

void cShapes::SetWindow(double aspectratio){
	// aspectratio=w/h
	int i; double xmax,xmin,ymax,ymin;
	shape s;
	point p1,p2, p;
	xmax=ymax=-1e30;
	xmin=ymin= 1e30;
	for(i=1; i<=shapes.GetSize(); ++i){
		shapes.GetData(s,i);
		switch(s.kind){
		case LINE:
			if(s.x1>xmax) xmax=s.x1;
			if(s.x2>xmax) xmax=s.x2;
			if(s.x1<xmin) xmin=s.x1;
			if(s.x2<xmin) xmin=s.x2;
			if(s.y1>ymax) ymax=s.y1;
			if(s.y2>ymax) ymax=s.y2;
			if(s.y1<ymin) ymin=s.y1;
			if(s.y2<ymin) ymin=s.y2;
		case TEXT:
			if(s.x1>xmax) xmax=s.x1;
			if(s.x1<xmin) xmin=s.x1;
			if(s.y1>ymax) ymax=s.y1;
			if(s.y1<ymin) ymin=s.y1;
		}
		p1=point(xmax,ymax); p2=point(xmin,ymin);
		if(aspectratio!=0){
			if( (p1.x-p2.x)/aspectratio > p1.y-p2.y ){
				p1.y+=((p1.x-p2.x)/aspectratio-(p1.y-p2.y))/2;
				p2.y-=((p1.x-p2.x)/aspectratio-(p1.y-p2.y))/2;
			}
			else{
				p1.x+=((p1.y-p2.y)*aspectratio-(p1.x-p2.x))/2;
				p2.x-=((p1.y-p2.y)*aspectratio-(p1.x-p2.x))/2;
			}
		}
		p=p1-p2;
		p1=p1+p*0.05;
		p2=p2-p*0.05;
		SetWindow(p1.x,p1.y, p2.x,p2.y);
	}
}

double cShapes::WindowAspectRatio(){
	// 横幅／高さ を返す
	return (window.TopRight().x-window.TopLeft().x)/(window.TopRight().y-window.BottomRight().y);
}

void cShapes::SetViewTopLeft(double x,double y){
	// ビューの左上を設定するが，
	// ビューの座標では，yの正方向が下なので，
	// この関数はview.BottomLeft()を設定する．
	point p=view.TopRight()-view.BottomLeft();
	view.p1=point(x,y);
	view.p2=view.p1+p;
}

void cShapes::SetViewWidth(double xw,double yw){
	// xw,ywの一方が正でないときは，ウインドウの縦横比と同じになる．
	point p;
	p=window.TopRight()-window.BottomLeft();
	if(xw>0 && yw<=0){
		p*=xw/p.x;
	}
	else if(xw<=0 && yw>0){
		p*=yw/p.y;
	}
	else if(xw>0 && yw>0){
		p.x*=xw/p.x;
		p.y*=yw/p.y;
	}
	else{
		return;
	}
	view.p1=view.BottomLeft();
	view.p2=view.p1+p;
	return;
}

void cShapes::RemoveTail(){
	shapes.RemoveTail();
	shapes_in_window.RemoveTail();
}

void cShapes::AddLine(double x1,double y1,double x2,double y2,long color,long linestyle){
	Add(shape(x1,y1,x2,y2,color,linestyle));
}

void cShapes::AddBox(double x1,double y1,double x2,double y2,long color/*=0*/,long linestyle/*=SOLID*/,
					 double ChamferX/*=0*/,double ChamferY/*=0=*/){
	// ChamferX,ChamferYは四隅の面取り
	if(x1<x2) Swap(x1,x2);
	if(y1<y2) Swap(y1,y2);

	AddLine(x1,y1-ChamferY, x1,y2+ChamferY, color,linestyle);
	if(ChamferX>0 || ChamferY>0) AddLine(x1,y2+ChamferY, x1-ChamferX,y2, color,linestyle);
	AddLine(x1-ChamferX,y2, x2+ChamferX,y2, color,linestyle);
	if(ChamferX>0 || ChamferY>0) AddLine(x2+ChamferX,y2, x2,y2+ChamferY, color,linestyle);
	AddLine(x2,y2+ChamferY, x2,y1-ChamferY, color,linestyle);
	if(ChamferX>0 || ChamferY>0) AddLine(x2,y1-ChamferY, x2+ChamferX,y1, color,linestyle);
	AddLine(x2+ChamferX,y1, x1-ChamferX,y1, color,linestyle);
	if(ChamferX>0 || ChamferY>0) AddLine(x1-ChamferX,y1, x1,y1-ChamferY, color,linestyle);
}

void cShapes::AddEllipse(double x,double y,double rx,double ry,int div,long color,long linestyle){
	int i;
	double x1,y1,x2,y2;
	for(i=0; i<=div-1; i++){
		x1=rx*cos(2*PI*i/div);
		y1=ry*sin(2*PI*i/div);
		x2=rx*cos(2*PI*(i+1)/div);
		y2=ry*sin(2*PI*(i+1)/div);
		AddLine(x1,y1,x2,y2,color,linestyle);
	}
}

void cShapes::AddCircle(double x,double y,double r,int div,long color,long linestyle){
	AddEllipse(x,y,r,r,div,color,linestyle);
}

void cShapes::AddText(std::string text,double x,double y,double fontsize,long color){
	Add(shape(text,x,y,fontsize,color));
}

void cShapes::Add(cShapes newdata) {
	int i;
	shape s;
	for(i=1; i<=newdata.shapes.GetSize(); i++){
		newdata.shapes.GetData(s,i);
		Add(s);
	}
}

cShapes cShapes::AddInLine(cShapes newdata) const{
	// ・*thisを cShapes s にコピーする．
	// ・newdataのshapesとwindowの大きさを調節して，newdataおよびsのwindow幅が同じになるようにする．
	// ・newdataのshapesをs.window.TopRight()とnesdata.window.TopLeft()がそろうように右横に加える．
	// ・全体を囲むにようにwindowを再設定する．sを返す．
	double xw,yw,xw1,yw1, m;
	cShapes s=*this;

	xw=s.window.TopRight().x - s.window.TopLeft().x;
	yw=s.window.TopRight().y - s.window.BottomRight().y;
	xw1=newdata.window.TopRight().x - newdata.window.TopLeft().x;
	yw1=newdata.window.TopRight().y - newdata.window.BottomRight().y;

	m=xw/xw1;
	newdata.ReSize(m);

	s.OffsetX=s.window.TopRight().x - newdata.window.TopLeft().x;
	s.OffsetY=s.window.TopRight().y - newdata.window.TopLeft().y;
	
	s.Add(newdata);
	s.SetWindow(window.BottomLeft().x,window.BottomLeft().y, window.TopRight().x+xw1*m,window.TopRight().y);

	return s;
}

shape cShapes::Add3DLine(vector<double> v1,vector<double> v2,long color,long linestyle){
	// 3D座標上の線v1-v2を画像座標xy面に投影した後，リストに加える．
	// 投影後の線分データを返す．
	shape s;

	v1=T*v1;
	v2=T*v2;
	s=shape(v1.x,v1.y,v2.x,v2.y,color,linestyle);
	Add(s);
	return s;
}

void cShapes::ViewAngle(double rox,double roy,double roz){
	// 画像座標は画面上右がx方向，上がy方向，手前がz方向となる．
	// 3D座標上の要素は画像座標xy面に投影されて描画される．
	// 3D座標は動かさず画像座標にrox,roy,rozの回転を加えて視線を移動させる．
	T=::Tmatrix(rox,roy,roz)*T;
}

void cShapes::TopView(){
	// 現在の画像を上から見る
	ViewAngle(-90,0,0);
}
void cShapes::BottomView(){
	// 現在の画像を下から見る
	ViewAngle(90,0,0);
}
void cShapes::RightView(){
	// 現在の画像を右から見る
	ViewAngle(0,90,0);
}
void cShapes::LeftView(){
	// 現在の画像を左から見る
	ViewAngle(0,-90,0);
}
void cShapes::BackView(){
	// 現在の画像を裏から見る
	ViewAngle(0,180,0);
}

void cShapes::xView(){
	// 現在の視線方向にかかわらず3D座標のx方向から見る
	//（上はy方向，右は-z方向，視線は-x方向）
	T=unit(T);
	ViewAngle(0,90,0);
}
void cShapes::yView(){
	// 現在の視線方向にかかわらず3D座標のy方向から見る（
	// (上は-z方向，右はx方向，視線は-y方向)
	T=unit(T);
	ViewAngle(-90,0,0);
}
void cShapes::zView(){
	// 現在の視線方向にかかわらず3D座標のz方向から見る
	//（上はy方向，右はx方向，視線は-z方向）
	T=unit(T);
}

matrix<double> cShapes::Tmatrix(){
	return T;
}

void cShapes::Reset() {
	shapes.RemoveAll();
	shapes_in_window.RemoveAll();
	initialize();
}

void cShapes::Clear() {
	shapes.RemoveAll();
	shapes_in_window.RemoveAll();
}

int cShapes::Size() const {
	return shapes_in_window.GetSize();
}

int cShapes::Kind(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) return s.kind;
	else                                return 0;
}

double cShapes::X1(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) {
		return (s.x1-window.BottomLeft().x)/(window.TopRight().x-window.BottomLeft().x)
		      *(view.BottomRight().x-view.TopLeft().x) +view.TopLeft().x;
		// 注：(view.BottomRight().x-view.TopLeft().x) > 0
	}
	else return 0;
}

double cShapes::Y1(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) {
		return (s.y1-window.BottomLeft().y)/(window.TopRight().y-window.BottomLeft().y)
		      *(view.BottomRight().y-view.TopLeft().y) +view.TopLeft().y;			
		// 注: (view.BottomRight().y-view.TopLeft().y) < 0
	}
	else return 0;
}

double cShapes::X2(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) {
		switch(s.kind){
		case  LINE: 		
			return (s.x2-window.BottomLeft().x)/(window.TopRight().x-window.BottomLeft().x)
		          *(view.BottomRight().x-view.TopLeft().x) +view.TopLeft().x;
		default: return 0;
		}
	}
	else return 0;
}

double cShapes::Y2(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) {
		switch(s.kind){
		case  LINE: 		
			return (s.y2-window.BottomLeft().y)/(window.TopRight().y-window.BottomLeft().y)
		          *(view.BottomRight().y-view.TopLeft().y) +view.TopLeft().y;
		default: return 0;
		}
	}
	else return 0;
}

int cShapes::LineStyle(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) {
		switch(s.kind){
		case LINE: return s.linestyle;
		default:   return 0;
		}
	}
	else return 0;
}

std::string cShapes::Text(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) {
		switch(s.kind) {
		case TEXT: return s.text;
		default:   return "";
		}
	}
	else return "";
}

double cShapes::FontSize(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) {
		switch(s.kind) {
		case TEXT: 
			return s.fontsize
			      *(view.TopLeft().y-view.BottomLeft().y)/(window.TopLeft().y-window.BottomLeft().y);
		default: return 0;
		}
	}
	else return 0;
}

long cShapes::Color(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) return s.color;
	else return 0;
}

void cShapes::Zoom(double m) {
	point pc, p1,p2;
	pc=(window.p1+window.p2)/2;
	if(m==0) m=2;
	p1=pc+(window.p1-pc)/m;
	p2=pc+(window.p2-pc)/m;
	SetWindow(p1.x,p1.y,p2.x,p2.y);
}

void cShapes::Translate(double dx,double dy){
	// window幅，高さのdx,dy倍移動する
	point p1,p2;
	dx=(window.TopRight().x-window.TopLeft().x)*dx;
	dy=(window.TopLeft().y-window.BottomLeft().y)*dy;
	p1=window.p1-point(dx,dy);
	p2=window.p2-point(dx,dy);
	SetWindow(p1.x,p1.y,p2.x,p2.y);
}

void cShapes::AllView(){
	SetWindow(4.0/3.0);
}

int cShapes::SaveAsBmp(std::string filename,int yPixels){
	// 全体を縦画素数がyPixelsのビットマップに保存する．
	// windowは既に設定されてある必要がある．
	int k, i,j,x1,y1,x2,y2,xw,yw;
	double a;
	cShapes s;
	cBitmap b;
	long backcolor=rgb(255,255,255); // 背景白

	s=*this;
	s.SetViewTopLeft(0,0);
	s.SetViewWidth(0,yPixels);
	b.SetM(yPixels);
	b.SetN((int)(yPixels*s.WindowAspectRatio()));

	for(i=1; i<=b.GetM(); ++i) for(j=1; j<=b.GetN(); ++j) b.SetRGB(i,j,backcolor);
	
	for(k=1; k<=s.Size(); ++k){
		switch(s.Kind(k)){
			case LINE:  // LineStyleは無視する．
				xw=(int)fabs(s.X1(k)-s.X2(k));
				yw=(int)fabs(s.Y1(k)-s.Y2(k));
				if( (xw>=yw && s.X1(k)<=s.X2(k)) || xw<yw && s.Y1(k)<=s.Y2(k) ){
					x1=(int)s.X1(k); y1=(int)s.Y1(k);
					x2=(int)s.X2(k); y2=(int)s.Y2(k);
				}
				else{
					x1=(int)s.X2(k); y1=(int)s.Y2(k);
					x2=(int)s.X1(k); y2=(int)s.Y1(k);
				}
				if(xw>=yw){
					a=(double)(y2-y1)/(double)(x2-x1);
					for(j=x1; j<=x2; ++j){
						i=(int)(y1+(double)(j-x1)*a);
						b.SetRGB(i,j,s.Color(k));
					}
				}
				else{
					a=(double)(x2-x1)/(double)(y2-y1);
					for(i=y1; i<=y2; ++i){
						j=(int)(x1+(double)(i-y1)*a);
						b.SetRGB(i,j,s.Color(k));
					}
				}
				break;
			default:
				// Kind(i)==TEXTでは何もしない．
				break;
		}
	}

	return b.Save(filename,1);
}