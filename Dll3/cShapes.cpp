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
	// �傫����m�{�ɂ���D
	// Zoom()�Ƃ͈قȂ�DZoom�ł�shape�̑傫���͕ς�炸�Cwindow�̑傫�����ς��̂݁D
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
	// �����^���� ��Ԃ�
	return (window.TopRight().x-window.TopLeft().x)/(window.TopRight().y-window.BottomRight().y);
}

void cShapes::SetViewTopLeft(double x,double y){
	// �r���[�̍����ݒ肷�邪�C
	// �r���[�̍��W�ł́Cy�̐����������Ȃ̂ŁC
	// ���̊֐���view.BottomLeft()��ݒ肷��D
	point p=view.TopRight()-view.BottomLeft();
	view.p1=point(x,y);
	view.p2=view.p1+p;
}

void cShapes::SetViewWidth(double xw,double yw){
	// xw,yw�̈�������łȂ��Ƃ��́C�E�C���h�E�̏c����Ɠ����ɂȂ�D
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
	// ChamferX,ChamferY�͎l���̖ʎ��
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
	// �E*this�� cShapes s �ɃR�s�[����D
	// �Enewdata��shapes��window�̑傫���𒲐߂��āCnewdata�����s��window���������ɂȂ�悤�ɂ���D
	// �Enewdata��shapes��s.window.TopRight()��nesdata.window.TopLeft()�����낤�悤�ɉE���ɉ�����D
	// �E�S�̂��͂ނɂ悤��window���Đݒ肷��Ds��Ԃ��D
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
	// 3D���W��̐�v1-v2���摜���Wxy�ʂɓ��e������C���X�g�ɉ�����D
	// ���e��̐����f�[�^��Ԃ��D
	shape s;

	v1=T*v1;
	v2=T*v2;
	s=shape(v1.x,v1.y,v2.x,v2.y,color,linestyle);
	Add(s);
	return s;
}

void cShapes::ViewAngle(double rox,double roy,double roz){
	// �摜���W�͉�ʏ�E��x�����C�オy�����C��O��z�����ƂȂ�D
	// 3D���W��̗v�f�͉摜���Wxy�ʂɓ��e����ĕ`�悳���D
	// 3D���W�͓��������摜���W��rox,roy,roz�̉�]�������Ď������ړ�������D
	T=::Tmatrix(rox,roy,roz)*T;
}

void cShapes::TopView(){
	// ���݂̉摜���ォ�猩��
	ViewAngle(-90,0,0);
}
void cShapes::BottomView(){
	// ���݂̉摜�������猩��
	ViewAngle(90,0,0);
}
void cShapes::RightView(){
	// ���݂̉摜���E���猩��
	ViewAngle(0,90,0);
}
void cShapes::LeftView(){
	// ���݂̉摜�������猩��
	ViewAngle(0,-90,0);
}
void cShapes::BackView(){
	// ���݂̉摜�𗠂��猩��
	ViewAngle(0,180,0);
}

void cShapes::xView(){
	// ���݂̎��������ɂ�����炸3D���W��x�������猩��
	//�i���y�����C�E��-z�����C������-x�����j
	T=unit(T);
	ViewAngle(0,90,0);
}
void cShapes::yView(){
	// ���݂̎��������ɂ�����炸3D���W��y�������猩��i
	// (���-z�����C�E��x�����C������-y����)
	T=unit(T);
	ViewAngle(-90,0,0);
}
void cShapes::zView(){
	// ���݂̎��������ɂ�����炸3D���W��z�������猩��
	//�i���y�����C�E��x�����C������-z�����j
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
		// ���F(view.BottomRight().x-view.TopLeft().x) > 0
	}
	else return 0;
}

double cShapes::Y1(int i) {
	shape s;
	if( shapes_in_window.GetData(s,i) ) {
		return (s.y1-window.BottomLeft().y)/(window.TopRight().y-window.BottomLeft().y)
		      *(view.BottomRight().y-view.TopLeft().y) +view.TopLeft().y;			
		// ��: (view.BottomRight().y-view.TopLeft().y) < 0
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
	// window���C������dx,dy�{�ړ�����
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
	// �S�̂��c��f����yPixels�̃r�b�g�}�b�v�ɕۑ�����D
	// window�͊��ɐݒ肳��Ă���K�v������D
	int k, i,j,x1,y1,x2,y2,xw,yw;
	double a;
	cShapes s;
	cBitmap b;
	long backcolor=rgb(255,255,255); // �w�i��

	s=*this;
	s.SetViewTopLeft(0,0);
	s.SetViewWidth(0,yPixels);
	b.SetM(yPixels);
	b.SetN((int)(yPixels*s.WindowAspectRatio()));

	for(i=1; i<=b.GetM(); ++i) for(j=1; j<=b.GetN(); ++j) b.SetRGB(i,j,backcolor);
	
	for(k=1; k<=s.Size(); ++k){
		switch(s.Kind(k)){
			case LINE:  // LineStyle�͖�������D
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
				// Kind(i)==TEXT�ł͉������Ȃ��D
				break;
		}
	}

	return b.Save(filename,1);
}