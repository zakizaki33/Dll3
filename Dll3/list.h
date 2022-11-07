#if !defined(AFX_LIST_H__6F2A879E_CFD1_4F46_B233_C068800B1A62__INCLUDED_)
#define AFX_LIST_H__6F2A879E_CFD1_4F46_B233_C068800B1A62__INCLUDED_
#define _CRT_SECURE_NO_WARNINGS
#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#if _MSC_VER > 1200   // コンパイラのバージョン VisualC++6.0 = 1200
#define TEMPLATECLASST template<class T>
#else
#define TEMPLATECLASST
#endif

#include <fstream>
#include <string>
#include <math.h>
#include "general_func.h"

template<class T> struct link
{ 
	T data;
	link<T> *pre, *suc;
	link(const T& x){ data=x; pre=suc=0; }
};

template<class T> class list  
{
	int maxsize;
	int size;
	link<T> *top,*tail;
	mutable link<T> *now; mutable int now_n;
	int point(int n) const;
public:
	list();
	list(const list& x);
	virtual ~list();
	list& operator=(const list& x);
	int GetMaxSize() const;
	void SetMaxSize(int n);
	int GetSize() const;
	int GetData(T& returnValue,int n) const;
	int SetData(const T& newValue,int n);
	int GetTailData(T& returnValue) const;
	T& operator[](int n) const;
	void Add(const T& NewData,int n);
	void AddTail(const T& NewData);
	int Remove(int n);
	int RemoveTail();
	void RemoveAll();
	// 【注意】
	//  Visual Studio 2013 でコンパイルするときは，
	//  friend 関数の先頭に "template<class T>" をつけないとエラーとなる．
	//  (Visual C++ 6.0 でコンパイルするときはつけるとエラーになる．)
	TEMPLATECLASST friend std::ostream& operator<<(std::ostream& to, const list<T>& x);
	TEMPLATECLASST friend std::istream& operator>>(std::istream& from, list<T>& x);
	void save(std::string filename) const;
	void open(std::string filename);
	T Sum();
	T Ave();
	T Stdev();
	T RMS();
	T Max();
	T Min();
	T PV();
	T Mode();
	void Swap(int i,int j);
	void Sort();
	void EraseDuplicates();
};


class stringlist : public list<std::string> 
{
public:
	void Add(const std::string& s);
};


template <class T> class stack
{
	list<T> un,re;
public:
	void SetMaxSize(int n);
	void push(const T& NewData);
	int pop(T& returnValue);
	int undo(T& returnValue);
	int redo(T& returnValue);
};


/////  list member  ///////////////////////////////////////////////////////////////////
template<class T> int list<T>::point(int n) const{
	if( n<1 || size<n ) return 0;
	if( now==0 ) { now=top; now_n=1; }
	while( now_n!=n ){
		if( n>now_n ) {
			now=now->suc; ++now_n;
		}
		else {
			now=now->pre; --now_n;
		}
	}
	return 1;
}

template<class T> list<T>::list(){ 
	maxsize=100000000;
	size=0;
	top=tail=now=0;
}

template<class T> list<T>::list(const list<T>& x) {
	link<T> *p;
	maxsize=100000000; size=0; top=tail=now=0;
	p=x.top;
	while(p!=0){
		AddTail(p->data);
		p=p->suc;
	}
}

template<class T> list<T>::~list(){ 
	RemoveAll(); 
}

template<class T> list<T>& list<T>::operator=(const list<T>& x){
	link<T> *p;
	RemoveAll();
	p=x.top;
	while(p!=0){
		AddTail(p->data);
		p=p->suc;
	}
	now=0;
	return *this;
}

template<class T> int list<T>::GetMaxSize() const{
	return maxsize;
}

template<class T> void list<T>::SetMaxSize(int n){
	while(size>n) { Remove(1); now=0; }
	maxsize=n;
}

template<class T> int list<T>::GetSize() const{
	return size;
}

template<class T> int list<T>::GetData(T& returnValue,int n) const{
	if( point(n) ){
		returnValue=now->data;
		return 1;
	}
	else return 0;
}

template<class T> int list<T>::SetData(const T& newValue,int n) {
	if( point(n) ){
		now->data=newValue;
		return 1;
	}
	else return 0;
}

template<class T> int list<T>::GetTailData(T& returnValue) const {
	if(tail!=0){
		returnValue=tail->data;
		return 1;
	}
	else return 0;
}

template<class T> T& list<T>::operator[](int n) const{
	// 06.12.05 作成
	// GetData(),SetData()よりプログラムを簡潔に書けるが，
	// 成功，不成功は返せない．
	static T x;
	if( point(n) ){
		return now->data;
	}
	else return x;
}

template<class T> void list<T>::Add(const T& NewData,int n){
//  Add an element after the n-th element.
	if(size==maxsize) return;
	if(n>size) n=size;
	if(n<0)    n=0;
	link<T> *pnew =new link<T>(NewData);
	if(size==0){ 
		top=tail=pnew;
	}
	else{
		if(n==size){
			tail->suc=pnew;
			pnew->pre=tail;
			tail=pnew;
		}
		else if(n==0){
			pnew->suc=top;
			top->pre=pnew;
			top=pnew;
			now_n++;
		}
		else{
			point(n);
			pnew->suc=now->suc;
			pnew->pre=now;
			if(now->suc==0){
				tail=pnew;
			}
			else{
				now->suc->pre=pnew;
			}
			now->suc=pnew;
		}
	}
	size++;
}

template<class T> void list<T>::AddTail(const T& NewData){
	// maxsizeを超える場合は先頭のデータを捨てる．
	while(size>=maxsize){
		Remove(1);
	}
	Add(NewData,size);
}

template<class T> int list<T>::Remove(int n){
	if( !point(n) ) return 0;
	if(now->pre==0){
		if(now->suc!=0) now->suc->pre=0;
		top=now->suc;
	}
	else{
		now->pre->suc=now->suc;
	}
	if(now->suc==0){
		if(now->pre!=0) now->pre->suc=0;
		tail=now->pre;
	}
	else{
		now->suc->pre=now->pre;
	}
	delete now; now=0;
	size--;
	return 1;
}

template<class T> int list<T>::RemoveTail(){
	if( Remove(size) ){
		now=tail; now_n=size;
		return 1;
	}
	else{
		return 0;
	}
}

template<class T> void list<T>::RemoveAll(){
	link<T> *p;
	while(tail!=0) {
		p=tail;
		tail=tail->pre;
		delete p;
	};
	size=0;
	top=tail=now=0;
}

template<class T> std::ostream& operator<<(std::ostream& to, const list<T>& x){
	link<T> *p;
//	to<<x.size<<std::endl;
	p=x.top;
	while(p!=0){
		to<<p->data<<std::endl;
		p=p->suc;
	}
	return to;
}

template<class T> std::istream& operator>>(std::istream& from, list<T>& x){
	T data;
	x.RemoveAll();
	while(from>>data) x.AddTail(data);
	/* 次のようにするとなぜか最後のデータを2回読み込んでしまう．
	do{
		from>>data;
		x.AddTail(data);
	} while( !from.eof() );
	*/
	return from;
}

template<class T> void list<T>::save(std::string filename) const{
	std::ofstream to(filename.c_str());
	to << *this;
}

template<class T> void list<T>::open(std::string filename){
	std::ifstream from(filename.c_str());
	from >> *this;
}

template<class T> T list<T>::Sum() {
	if(size==0) return 0;
	int i;
	T x,sum;
	sum=0;
	for(i=1; i<=size; ++i){
		GetData(x,i);
		sum+=x;
	}
	return sum;
}

template<class T> T list<T>::Ave() {
	return size==0 ? 0 : Sum()/size;
}

template<class T> T list<T>::Stdev() {
	int i;
	if(size==0) return 0;
	T x,sum,ave;
	sum=0; ave=this->Ave();
	for(i=1; i<=size; ++i){
		GetData(x,i);
		sum+=(x-ave)*(x-ave);
	}
	return sqrt(sum/size);
}

template<class T> T list<T>::RMS() {
	int i;
	if(size==0) return 0;
	T x,sum;
	sum=0;
	for(i=1; i<=size; ++i){
		GetData(x,i);
		sum+=x*x;
	}
	return sqrt(sum/size);
}

template<class T> T list<T>::Max() {
	int i;
	T x1,x2,max;
	switch(size){
	case 0:
		return 0;
	case 1:
		GetData(x1,1);
		return x1;
	default:
		GetData(x1,1);
		GetData(x2,2);
		max= x1>x2 ? x1 : x2;
		for(i=3; i<=size; ++i){
			GetData(x1,i);
			max= max>x1 ? max : x1;
		};
		return max;
	}
}

template<class T> T list<T>::Min() {
	int i;
	T x1,x2,min;
	switch(size){
	case 0:
		return 0;
	case 1:
		GetData(x1,1);
		return x1;
	default:
		GetData(x1,1); 
		GetData(x2,2);
		min= x1<x2 ? x1 : x2;
		for(i=3; i<=size; ++i){
			GetData(x1,i);
			min= min<x1 ? min : x1;
		};
		return min;
	}
}

template<class T> T list<T>::PV() {
	return Max()-Min();
}

template<class T> T list<T>::Mode() {
	// 最頻値を求める．最頻値が複数ある場合はその内のいずれかを返す．
	list<T> li;
	int i,count;
	T x0,x, max, mode;
	li=*this; // ソートするのでコピーをつくる
	li.Sort();
	if(li.GetData(x0,1)){
		count=1;
		mode=x0;
		max=1;
	}
	else{
		return 0;  // データが1つもない
	}
	for(i=2; i<=li.GetSize(); i++){
		li.GetData(x,i);
		if(x0==x){
			count++;
			if(count>max){
				max=count;
				mode=x;
			}
		}
		else{
			x0=x;
			count=1;
		}
	}
	return mode;
}

template<class T> void list<T>::Swap(int i,int j) {
	T xi,xj;
	if( 1<=i && i<=size && 1<=j && j<=size){
		GetData(xi,i); GetData(xj,j);
		SetData(xj,i); SetData(xi,j);
	}
}

template<class T> void list<T>::Sort() {
	T *buf;
	int i;

	buf=new T[size];
	for(i=1; i<=size; i++) GetData(buf[i-1],i);
	::Sort(buf,size);
	for(i=1; i<=size; i++) SetData(buf[i-1],i);
	delete [] buf;

	/*  
	int i,j;
	T xi,xj;
	for( i=1; i<size; i++ ) {
		for( j=i+1; j<=size; j++ ) {
			GetData(xi,i); GetData(xj,j);
			if(xi>xj) Swap(i,j);
		}
	}
	*/
	//// 高速化 ///////////////////////////
	/*
	int i,j;
	T x,*p;

	p=new T[size+1];
	for(i=1; i<=size; i++){
		GetData(p[i],i);
	}
	for( i=1; i<size; i++ ) {
		for( j=i+1; j<=size; j++ ) {
			if(p[i]>p[j]){
				x=p[i];
				p[i]=p[j];
				p[j]=x;
			}
		}
	}
	for(i=1; i<=size; i++){
		SetData(p[i],i);
	}
	delete [] p;
	*/
/*	
	///// 高速化 ヒープソート /////////////////////////////
	
	int i,j,k,n;
	T *heap, temp;

	heap=new T [size+1];     // heap[1]を根とするヒープ用配列

	for(n=1; n<=size; n++){  // ヒープに値を入れていく
		GetData(heap[n],n);
		i=n;
		j=i/2;                           // heap[i]の親はheap[i/2]
		while(j>=1 && heap[j]>heap[i]){  // 親の方が大きければ子と値を交換
			temp=heap[i];
			heap[i]=heap[j];
			heap[j]=temp;
			i=j;
			j=i/2;
		}
	}

	for(n=1; n<=size; n++){  // 根の値をリストに加えていく
		SetData(heap[1],n);  
		k=size-n;
		heap[1]=heap[k+1];   // ヒープの末尾の値を根に代入
		i=1; j=i*2;          // heap[i]の子はheap[i*2],heap[i*2+1]
		while(j<=k){
			if(j+1<=k && heap[j]>heap[j+1]) j++;   // 値が小さい方の子をheap[j]とする
			if(heap[i]>heap[j]){                   // 子の方が小さければ親と値を交換
				temp=heap[i];
				heap[i]=heap[j];
				heap[j]=temp;
			}
			i=j;
			j=i*2;
		}
	}

	delete [] heap;
	
	////////////////////////////////////////////////////// 
*/
}

template<class T> void list<T>::EraseDuplicates() {
	int i,j;
	T xi,xj;
	for( i=1; i<size; ++i ) {
		for( j=i+1; j<=size; ++j ) {
			GetData(xi,i); GetData(xj,j);
			if(xi==xj) { Remove(j); --j; }
		}
	}
	now=0;
}

/////  stack member  ///////////////////////////////////////////////////////////////////

template<class T> void stack<T>::SetMaxSize(int n){
	un.SetMaxSize(n);
	re.SetMaxSize(n);
}

template<class T> void stack<T>::push(const T& NewData){
	un.AddTail(NewData);
}

template<class T> int stack<T>::pop(T& returnValue){
	if( un.GetTailData(returnValue) ){
		un.RemoveTail();
		return 1;
	}
	else {
		return 0;
	}
}

template<class T> int stack<T>::undo(T& returnValue){
	if( un.GetSize()>0 ){
		re.AddTail(returnValue);
		un.GetTailData(returnValue);
		un.RemoveTail();
		return 1;
	}
	else {
		return 0;
	}
}

template<class T> int stack<T>::redo(T& returnValue){
	if( re.GetSize()>0 ){
		un.AddTail(returnValue);
		re.GetTailData(returnValue);
		re.RemoveTail();
		return 1;
	}
	else {
		return 0;
	}
}

#endif // !defined(AFX_LIST_H__6F2A879E_CFD1_4F46_B233_C068800B1A62__INCLUDED_)