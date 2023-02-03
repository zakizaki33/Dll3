#include "stdafx.h"
#include "string.h"

bool is_kanji(const char *p,int n){
	// p[n] ���S�p�����̈ꕔ�̂Ƃ�true��Ԃ��D
	int i=0;
	while(p[i]!=0){
		unsigned char c=p[i];
		if( (0x81<=c && c<=0x9F) || (0xE0<=c && c<=0xEF) ){
			// �V�t�gJIS�R�[�h�̑S�p����1�o�C�g�� 0x81�`0x9F, 0xE0�`0xEF
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
	// s[i-1]���󔒂�s[i]�͋󔒂łȂ��Ƃ��^�D
	// ���p���Ȃ��ꍇ�͂�����word�̓��ƂȂ�D
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
	// �G�X�P�[�v�V�[�P���X��'\'������
	int i;
	std::string x=s;
	for(i=0; i<=(int)x.length()-2; i++){ 
		// x.length()�̖߂�l��unsigned. (int)���Ȃ���length()<=1�̂Ƃ��\�����ʌ��ʂɂȂ�D
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
	// �^����ꂽs[i]������( " �܂��� ' �܂��� ( �܂��� { �܂��� [ )
	// �̂Ƃ����ʕ��̓��ƔF�����C
	// i�����ʕ��̍Ō�܂Ői�߂�D
	// �Ō�Ƃ͎��̓������ʂ̈ʒu�i"(","{","[" �̂Ƃ��͋t�����̊��ʂƂȂ�ʒu)�D
	// ���ꂪ�Ȃ����s�̍Ō�ł���D
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

	// �܂����E�̋󔒂��폜����D
	s1=s;
	while( is_space(s1[0]) ) s1.erase(0,1);
	while( is_space(s1[s1.length()-1]) ) s1.erase(s1.length()-1,1);
	
	if(erase_bracket){
		// ����Ɏc�������������ʕ�,
		// "*" �܂��� "* �܂��� '*' �܂��� '* �܂��� 
		// (*) �܂��� (* �܂��� {*} �܂��� {* �܂��� [*] �܂��� [*
		// �̏ꍇ�́C���ʂ��폜����D
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
	// �p�X�����܂ރt�@�C�������p�X���������D
	const int SIZE=1000;
	char buf[SIZE];
	char* p=buf;
	int i;
	strcpy(buf,filename.c_str());
	for(i=0; i<=SIZE-1; ++i){
		if( (buf[i]=='\\' || buf[i]=='/') && i+1<=SIZE-1 && !is_kanji(buf,i) ) p=&buf[i+1];
		// !is_kanji(buf,i) ���Ȃ��� remive_path("�\�ʋ�") == "�ʋ�" �ȂǂƂȂ��Ă��܂��D
		if(buf[i]==0){
			if( i-1>=0){
				if(buf[i-1]=='.') buf[i-1]=0; // �Ōオ'.'�̏ꍇ����������D
			}
			break;
		}
	}
	return p;
}

std::string remove_extension(const std::string& filename){
	// �t�@�C�������g���q�������D
	// �g���q�Ƃ͍Ō�̃s���I�h������̕����Ƃ���i�g���q�̕������͕s��j�D
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
	// �ǋL�@2023-02-03�@��������
	std::string w = word_(s, n);
	if (w == "") {
		return "";
	}
	// �����܂�
	return erase_yen( trim(word_(s,n),erase_bracket) ); 
}

int sentences(const std::string& s){
	// �󔒂̈��ł���';'�ŋ�؂���word�̏W����sentence�Ƃ���D
	// ������������I�[��";"���Ȃ��Ă�sentence�I�[�ƂȂ�D
	int i, m=0;
	std::string s1=trim(s,1);

	for(i=0; i<(int)s1.length()-1; ++i){ // ������I�[��";"�͐����Ȃ��Ds1.length()==0�̂Ƃ�(int)�K�v�D
		step_brackted(s1,i);
		if(s1[i]==';') ++m;
	}
	return s1.length()==0 ? 0 : m+1;  // ������I�[�����Z�D
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
				// ���ӁF
				//   trim(s1.substr(..),1) �Ƃ���ƁC���Ȃ킿���ʂ������ƕs����N����D
				//   �Ⴆ�΁CcLens::optimize()�̈����Ctarget�ɂ����āC
				//   "zoo 3"; �� zoo 3; �͈قȂ�D
				//   �O�҂͑�3�Y�[���|�W�V�������w�肵�Ă��邪�C��҂�
				//   arg(s,0)="zoo"�Carg(s,1)="3" �Ɖ��߂���C���Ӗ��̂��̂ƂȂ��Ă��܂��D
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
	// arg(i1)����arg(i2)���폜�������̂�Ԃ�
	int i;
	std::string s;

	for(i=0; i<=args(sentence); i++){
		if(i<i1 || i2<i) s+=word_(sentence,i+1)+" ";
	}
	return s;
}

std::string replace_arg(const std::string& sentence,int i,const std::string& new_val){
	// arg(i)��new_val�Œu�����������̂�Ԃ�
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
	// �󔒂��܂ޕ�������t�@�C���ɏ������ނƁC
	// �ǂݏo���ł͂��̕�����͋󔒂̏��ŕ��f����Ă��܂��D
	// �����h�����ߋ󔒂𕶎�C1+C2�Œu��������D
	std::string s1,tmp;
	int i;
	const char C1=(char)0x82;   // 0x82(�A�X�L�[�R�[�h)��unsigned, char��signed �Ȃ̂ŃL���X�g���K�v
	const char C2=(char)0x01;   
	                    // 0x8201 �̓V�t�gJIS�R�[�h�ŋ󂫂ɂȂ��Ă���D
	                    // �ȑO�͎g���\���̒Ⴂ '~' ��p���� C='~'=0x7e �Œu�����Ă����D
						// ���p���������ł���΂���Ŗ��Ȃ��D
						// �������C�Ⴆ�� "�~" = 0x837e �̌㔼��u�����Ă��܂��s������������D
	                    // ��͂�󔒂͋󂢂Ă���S�p�����R�[�h�Œu������̂��m���ƍl������D(2015.05.05)

	s1=s;               // s��const�Ȃ̂ŃR�s�[�ɑ΂��đ��삷��
	for(i=0; i<(int)s1.length(); ++i){
		if(is_space(s1,i) && s1[i]!=';'){ // sentence�̋�؂�";"�͕ۑ�������(110601)
			s1[i]=C1; 
			tmp=C2;
			s1.insert(i+1,tmp);
		}
	}
	// �󔒕����񂾂Ɠǂݍ��ݎ��ɔ�΂���C�ȍ~�̓ǂݍ��݂ɂ��ꂪ�����Ă��܂��̂ŁC
	// �擪��C���������̂��t�@�C���֊i�[����D
	s1=C1+(C2+s1);
	return s1;
}

std::string inv_blank_filled(const std::string& s){
	// blank_filled()�̕ϊ������ɖ߂��D
	std::string s1;
	int i;
	const char C1=(char)0x82;
	const char C2=(char)0x01;  

	{
		// ���t�@�C���t�H�[�}�b�g('~'�Œu����������)�ɑΉ�����D//////////
		const char C='~';

		s1=s;
		if(s1[0]==C){
			s1.erase(0,1);  // �ŏ���'~'���폜
			for(i=0; i<(int)s1.length(); ++i){
				if(s1[i]==C) s1[i]=' ';
			}
			return s1;
		}
		// ���t�@�C���t�H�[�}�b�g�ւ̑Ή������܂� /////////////////////////
	}

	s1=s;
	// �ŏ���C1,C2���폜
	s1.erase(0,2);
	for(i=0; i<(int)s1.length(); ++i){
		// C1,C2���󔒂ɖ߂�
		if(s1[i]==C1 && s1[i+1]==C2){
			s1[i]=' ';
			s1[i+1]=' ';
			s1.erase(i,1);
		}
	}
	return s1;
}

std::string to_space(const std::string& s,char c){
	// s�̒���c���󔒂ɒu��������D
	// ��: s="0.07_-4.2", c='_'  -> "0.07 -4.2"
	std::string s1;
	int i;

	s1=s;
	for(i=0; i<(int)s1.length(); i++){
		if(s1[i]==c) s1[i]=' ';
	}
	return s1;
}

std::string spc(int n){
	// ����n�̋󔒕�����
	return std::string(n,' ');
}

std::string left(const std::string& s,int n){
	// �擪����n�����𔲂��o��
	return s.substr(0,n);
}

std::string mid(const std::string& s,int start,int n/*=0*/){
	// start�Ԗڂ̕�������n����(�I�[�ɒB�������͏I�[�܂�)�𔲂��o���D
	// n���ȗ����ꂽ�Ƃ��͍Ō�܂Ŕ����o���D
	if(n>0){
		if( start-1+n > (int)s.length() ) n=(int)s.length()-start+1;  // n���傫������ꍇ�ɑΉ�
		return s.substr(start-1,n);
	}
	else{
		return s.substr(start-1,s.length()-start+1);
	}
}

std::string replace(const std::string& s,const std::string& before,const std::string& after){
	// s�Ɋ܂܂��S�Ă�before��after�ɒu������D
	int i;
	std::string buf;

	buf=s;
	for(i=1; i<=(int)buf.length()-(int)before.length()+1; ++i){  // length()�̖߂�l��unsigned�Ȃ̂� (int) �K�v
		if( mid(buf,i,(int)before.length())==before ){
			buf= left(buf,i-1) + after + mid(buf,i+(int)before.length());
			i+=(int)after.length()-1;
		}
	}
	return buf;
}

bool forward_match(const std::string& s0,const std::string& s){
	// s�̑O��������s0�ƈ�v���邩�ǂ�����Ԃ�
	// ��F  forward_match("abc","abcd")=true,  forward_match("abcd","abc")=false
	return left(s,(int)s0.length())==s0;
}

std::string date(){
	// ���݂̔N������ "yymmdd" �̌`�ŕԂ�
	int year,month,date;
	char buf[100];
	time_t timer;
	tm *local;

	timer=time(NULL);
	local=localtime(&timer);

	year=local->tm_year+1900-2000;  // 2000�`2999�N�̉�2��
	month=local->tm_mon+1;
	date=local->tm_mday;

	sprintf(buf,"%2d%2d%2d",year,month,date);
	return buf;
}

std::string time(){
	// ���݂̎����� "hhmmss" �̌`�ŕԂ�
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