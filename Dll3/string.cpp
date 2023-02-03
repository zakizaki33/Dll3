#include "stdafx.h"
#include "string.h"

bool is_kanji(const char *p,int n){
	// p[n] が全角文字の一部のときtrueを返す．
	int i=0;
	while(p[i]!=0){
		unsigned char c=p[i];
		if( (0x81<=c && c<=0x9F) || (0xE0<=c && c<=0xEF) ){
			// シフトJISコードの全角文字1バイト目 0x81～0x9F, 0xE0～0xEF
			if(p[i+1]==0) return false;
			if(i==n || i+1==n) return true;
			i+=2;
		}
		else{
			if(i==n) return false;
			i++;
		}
	}
	return false;
}

bool is_space(const char& c){
	if( c==' ' || c=='\t' || c=='\n' || c=='\r' ||c==';' ) return true;
	return false;
}

bool is_space(const std::string& s,int i){
	if( i<0 || (int)s.length()-1<i ){
		return true;
	}
	else{
		return is_space(s[i]);
	}
}

bool is_topedge(const std::string& s,int i){
	// s[i-1]が空白でs[i]は空白でないとき真．
	// 引用がない場合はここがwordの頭となる．
	if( i==0 ){
		if( !is_space(s[0]) ) return true;
		else                  return false;
	}
	else if( 0<=i && i<=(int)s.length()-1 ){
		if( is_space(s[i-1]) && !is_space(s[i]) ) return true;
		else                                      return false;
	}
	else return false;
}

bool escaped(const std::string& s,int i){
	if( i<=0 || (int)s.length()-1<i ) return false;
	if( 1<=i && s[i-1]=='\\' ) return true; else return false;
}

std::string erase_yen(const std::string& s){
	// エスケープシーケンスの'\'を除く
	int i;
	std::string x=s;
	for(i=0; i<=(int)x.length()-2; i++){ 
		// x.length()の戻り値はunsigned. (int)がないとlength()<=1のとき予期せぬ結果になる．
		if(x[i]=='\\'){
			if(x[i+1]=='"' ) { x.replace(i,2,"\""); i++; }
			if(x[i+1]=='\'') { x.replace(i,2,"'" ); i++; }
			if(x[i+1]=='\\') { x.replace(i,2,"\\"); i++; } 
		}
	}
	return x;	
}

bool is_singlequote(const std::string& s,int i){
	if(i<0 || (int)s.length()-1<i) return false;
	if( s[i]=='\'' && !escaped(s,i) ) {
		return true;
	}
	else{
		return false;
	}
}

bool is_doublequote(const std::string& s,int i){
	if(i<0 || (int)s.length()-1<i) return false;
	if( s[i]=='"' && !escaped(s,i) ) {
		return true;
	}
	else{
		return false;
	}
}

void step_brackted(const std::string& s,int& i){
	// 与えられたs[i]が括弧( " または ' または ( または { または [ )
	// のとき括弧部の頭と認識し，
	// iを括弧部の最後まで進める．
	// 最後とは次の同じ括弧の位置（"(","{","[" のときは逆向きの括弧となる位置)．
	// それがなければsの最後である．
	if( i<0 || (int)s.length()-1<=i ) return;
	if(is_doublequote(s,i)){
		do i++; while( i<(int)s.length()-1 && !(is_doublequote(s,i)) );
	}
	if(is_singlequote(s,i)){
		do i++; while( i<(int)s.length()-1 && !(is_singlequote(s,i)) );
	}
	if(s[i]=='('){
		do i++; while( i<(int)s.length()-1 && !(s[i]==')') );
	}
	if(s[i]=='{'){
		do i++; while( i<(int)s.length()-1 && !(s[i]=='}') );
	}
	if(s[i]=='['){
		do i++; while( i<(int)s.length()-1 && !(s[i]==']') );
	}
}

std::string trim(const std::string& s,int erase_bracket){
	int i;
	std::string s1;

	// まず左右の空白を削除する．
	s1=s;
	while( is_space(s1[0]) ) s1.erase(0,1);
	while( is_space(s1[s1.length()-1]) ) s1.erase(s1.length()-1,1);
	
	if(erase_bracket){
		// さらに残った部分が括弧部,
		// "*" または "* または '*' または '* または 
		// (*) または (* または {*} または {* または [*] または [*
		// の場合は，括弧を削除する．
		while(true){
			if(s1.length()<2) break;
			i=0;
			step_brackted(s1,i);
			if(i==s1.length()-1){
				if( is_singlequote(s1,0) ){
					s1.erase(0,1);
					if( is_singlequote(s1,(int)s1.length()-1) ){
						s1.erase(s1.length()-1,1);	
					}
				}
				else if( is_doublequote(s1,0) ){
					s1.erase(0,1);
					if( is_doublequote(s1,(int)s1.length()-1) ){
						s1.erase(s1.length()-1,1);
					}
				}
				else if( s1[0]=='(' ){
					s1.erase(0,1);
					if( s1[s1.length()-1]==')' ){
						s1.erase(s1.length()-1,1);
					}
				}
				else if( s1[0]=='{' ){
					s1.erase(0,1);
					if( s1[s1.length()-1]=='}' ){
						s1.erase(s1.length()-1,1);
					}
				}
				else if( s1[0]=='[' ){
					s1.erase(0,1);
					if( s1[s1.length()-1]==']' ){
						s1.erase(s1.length()-1,1);
					}
				}
			}
			else break;
		}
	}

	return s1;
}

std::string word_(const std::string& s,int n){
	int i, m, i1,i2;
	if(n<=0) return "";
	m=0;
	for(i=0; i<(int)s.length(); ++i){
		if( is_topedge(s,i) ) ++m;
		if(m==n){
			i1=i;
			for(i=i1; i<(int)s.length(); ++i){
				step_brackted(s,i);
				if( is_space(s[i]) ){
					i2=i-1;
					return s.substr(i1,i2-i1+1);
				}
				if( i==s.length()-1 ){
					i2=i;
					return s.substr(i1,i2-i1+1);
				}
			}	
		}
		step_brackted(s,i);
	}
	return "";
}

std::string remove_path(const std::string& filename){
	// パス名を含むファイル名よりパス名を除く．
	const int SIZE=1000;
	char buf[SIZE];
	char* p=buf;
	int i;
	strcpy(buf,filename.c_str());
	for(i=0; i<=SIZE-1; ++i){
		if( (buf[i]=='\\' || buf[i]=='/') && i+1<=SIZE-1 && !is_kanji(buf,i) ) p=&buf[i+1];
		// !is_kanji(buf,i) がないと remive_path("表面鏡") == "面鏡" などとなってしまう．
		if(buf[i]==0){
			if( i-1>=0){
				if(buf[i-1]=='.') buf[i-1]=0; // 最後が'.'の場合これを除く．
			}
			break;
		}
	}
	return p;
}

std::string remove_extension(const std::string& filename){
	// ファイル名より拡張子を除く．
	// 拡張子とは最後のピリオドから後ろの部分とする（拡張子の文字数は不問）．
	int i;
	std::string s=filename;

	for(i=(int)s.length()-1; i>0; --i){
		if(s[i]=='.'){
			return s.substr(0,i);
		}
	}
	return s;
}

bool is_numeric(const std::string& s){
	if('0'<=s[0] && s[0]<='9') return true;
	if(s[0]=='+' || s[0]=='-'){
		if('0'<=s[1] && s[1]<='9') return true;
		if(s[1]=='.'){
			if('0'<=s[2] && s[2]<='9') return true;
		}
	}
	if(s[0]=='.'){
		if('0'<=s[1] && s[1]<='9') return true;
	}
	return false;	
}

bool is_string(const std::string& s){
	return s=="" ? false : true;
}

int words(const std::string& s){
	int i, m=0;
	for(i=0; i<(int)s.length(); ++i){
		if( is_topedge(s,i) ) ++m;
		step_brackted(s,i);
	}
	return m;
}

std::string word(const std::string& s,int n,int erase_bracket/*=1*/){
	// 追記　2023-02-03　ここから
	std::string w = word_(s, n);
	if (w == "") {
		return "";
	}
	// ここまで
	return erase_yen( trim(word_(s,n),erase_bracket) ); 
}

int sentences(const std::string& s){
	// 空白の一種である';'で区切られるwordの集合をsentenceとする．
	// ただし文字列終端は";"がなくてもsentence終端となる．
	int i, m=0;
	std::string s1=trim(s,1);

	for(i=0; i<(int)s1.length()-1; ++i){ // 文字列終端の";"は数えない．s1.length()==0のとき(int)必要．
		step_brackted(s1,i);
		if(s1[i]==';') ++m;
	}
	return s1.length()==0 ? 0 : m+1;  // 文字列終端を加算．
}

std::string sentence(const std::string& s,int n){
	int i,i0=0, m=0;
	std::string s1=trim(s,1);
	for(i=0; i<(int)s1.length(); ++i){
		step_brackted(s1,i);
		if(s1[i]==';'){
			++m;
			if(m==n){
				return trim(s1.substr(i0,i-i0),0);
				// 注意：
				//   trim(s1.substr(..),1) とすると，すなわち括弧も除くと不具合が起こる．
				//   例えば，cLens::optimize()の引数，targetにおいて，
				//   "zoo 3"; と zoo 3; は異なる．
				//   前者は第3ズームポジションを指定しているが，後者は
				//   arg(s,0)="zoo"，arg(s,1)="3" と解釈され，無意味のものとなってしまう．
			}
			i0=i;
		}
		if(i==s1.length()-1){
			++m;
			if(m==n){
				return trim(s1.substr(i0,i-i0+1),0);
			}
			i0=i;
		}
	}
	return "";
}

int args(const std::string& sentence){
	return words(sentence)-1;
}

std::string arg(const std::string& sentence,int n,int erase_bracket/*=1*/){
	return word(sentence,n+1,erase_bracket);
}

std::string remove_arg(const std::string& sentence,int i1,int i2){
	// arg(i1)からarg(i2)を削除したものを返す
	int i;
	std::string s;

	for(i=0; i<=args(sentence); i++){
		if(i<i1 || i2<i) s+=word_(sentence,i+1)+" ";
	}
	return s;
}

std::string replace_arg(const std::string& sentence,int i,const std::string& new_val){
	// arg(i)をnew_valで置き換えたものを返す
	int ii;
	std::string s;

	for(ii=0; ii<=args(sentence); ii++){
		if(ii==i){
			s+=new_val+" ";	
		}
		else{
			s+=word_(sentence,ii+1)+" ";
		}
	}
	return s;
}

std::string str(double x){
	char buf[1000];
	sprintf(buf,"%g",x);
	return buf;
}

std::string blank_filled(const std::string& s){
	// 空白を含む文字列をファイルに書き込むと，
	// 読み出しではその文字列は空白の所で分断されてしまう．
	// これを防ぐため空白を文字C1+C2で置き換える．
	std::string s1,tmp;
	int i;
	const char C1=(char)0x82;   // 0x82(アスキーコード)はunsigned, charはsigned なのでキャストが必要
	const char C2=(char)0x01;   
	                    // 0x8201 はシフトJISコードで空きになっている．
	                    // 以前は使う可能性の低い '~' を用いて C='~'=0x7e で置換していた．
						// 半角文字だけであればこれで問題ない．
						// しかし，例えば "ミ" = 0x837e の後半を置換してしまう不具合が発生した．
	                    // やはり空白は空いている全角文字コードで置換するのが確実と考えられる．(2015.05.05)

	s1=s;               // sがconstなのでコピーに対して操作する
	for(i=0; i<(int)s1.length(); ++i){
		if(is_space(s1,i) && s1[i]!=';'){ // sentenceの区切り";"は保存したい(110601)
			s1[i]=C1; 
			tmp=C2;
			s1.insert(i+1,tmp);
		}
	}
	// 空白文字列だと読み込み時に飛ばされ，以降の読み込みにずれが生じてしまうので，
	// 先頭にCをつけたものをファイルへ格納する．
	s1=C1+(C2+s1);
	return s1;
}

std::string inv_blank_filled(const std::string& s){
	// blank_filled()の変換を元に戻す．
	std::string s1;
	int i;
	const char C1=(char)0x82;
	const char C2=(char)0x01;  

	{
		// 旧ファイルフォーマット('~'で置換したもの)に対応する．//////////
		const char C='~';

		s1=s;
		if(s1[0]==C){
			s1.erase(0,1);  // 最初の'~'を削除
			for(i=0; i<(int)s1.length(); ++i){
				if(s1[i]==C) s1[i]=' ';
			}
			return s1;
		}
		// 旧ファイルフォーマットへの対応ここまで /////////////////////////
	}

	s1=s;
	// 最初のC1,C2を削除
	s1.erase(0,2);
	for(i=0; i<(int)s1.length(); ++i){
		// C1,C2を空白に戻す
		if(s1[i]==C1 && s1[i+1]==C2){
			s1[i]=' ';
			s1[i+1]=' ';
			s1.erase(i,1);
		}
	}
	return s1;
}

std::string to_space(const std::string& s,char c){
	// sの中のcを空白に置き換える．
	// 例: s="0.07_-4.2", c='_'  -> "0.07 -4.2"
	std::string s1;
	int i;

	s1=s;
	for(i=0; i<(int)s1.length(); i++){
		if(s1[i]==c) s1[i]=' ';
	}
	return s1;
}

std::string spc(int n){
	// 長さnの空白文字列
	return std::string(n,' ');
}

std::string left(const std::string& s,int n){
	// 先頭からn文字を抜き出す
	return s.substr(0,n);
}

std::string mid(const std::string& s,int start,int n/*=0*/){
	// start番目の文字からn文字(終端に達した時は終端まで)を抜き出す．
	// nが省略されたときは最後まで抜き出す．
	if(n>0){
		if( start-1+n > (int)s.length() ) n=(int)s.length()-start+1;  // nが大きすぎる場合に対応
		return s.substr(start-1,n);
	}
	else{
		return s.substr(start-1,s.length()-start+1);
	}
}

std::string replace(const std::string& s,const std::string& before,const std::string& after){
	// sに含まれる全てのbeforeをafterに置換する．
	int i;
	std::string buf;

	buf=s;
	for(i=1; i<=(int)buf.length()-(int)before.length()+1; ++i){  // length()の戻り値はunsignedなので (int) 必要
		if( mid(buf,i,(int)before.length())==before ){
			buf= left(buf,i-1) + after + mid(buf,i+(int)before.length());
			i+=(int)after.length()-1;
		}
	}
	return buf;
}

bool forward_match(const std::string& s0,const std::string& s){
	// sの前方部分がs0と一致するかどうかを返す
	// 例：  forward_match("abc","abcd")=true,  forward_match("abcd","abc")=false
	return left(s,(int)s0.length())==s0;
}

std::string date(){
	// 現在の年月日を "yymmdd" の形で返す
	int year,month,date;
	char buf[100];
	time_t timer;
	tm *local;

	timer=time(NULL);
	local=localtime(&timer);

	year=local->tm_year+1900-2000;  // 2000～2999年の下2桁
	month=local->tm_mon+1;
	date=local->tm_mday;

	sprintf(buf,"%2d%2d%2d",year,month,date);
	return buf;
}

std::string time(){
	// 現在の時刻を "hhmmss" の形で返す
	int hour,minute,second;
	char buf[100];
	time_t timer;
	tm *local;

	timer=time(NULL);
	local=localtime(&timer);

	hour=local->tm_hour;
	minute=local->tm_min;
	second=local->tm_sec;

	sprintf(buf,"%2d%2d%2d",hour,minute,second);
	return buf;
}