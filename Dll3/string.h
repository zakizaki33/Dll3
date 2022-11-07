#ifndef STRING_H_INCLUDED
#define STRING_H_INCLUDED
#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <time.h>
#include "complex.h"

bool is_kanji(const char *p,int n);
bool is_space(const char& c);
bool is_space(const std::string& s,int i);
bool is_topedge(const std::string& s,int i);
bool escaped(const std::string& s,int i);
std::string erase_yen(const std::string& s);
bool is_singlequote(const std::string& s,int i);
bool is_doublequote(const std::string& s,int i);
void step_brackted(const std::string& s,int& i);
std::string trim(const std::string& s,int erase_bracket);
std::string word_(const std::string& s,int n);
std::string remove_path(const std::string& filename);
std::string remove_extension(const std::string& filename);
bool is_numeric(const std::string& s);
bool is_string(const std::string& s);
int words(const std::string& s);
std::string word(const std::string& s,int n,int erase_bracket=1);
int sentences(const std::string& s);
std::string sentence(const std::string& s,int n);
int args(const std::string& sentence);
std::string arg(const std::string& sentence,int n,int erase_bracket=1);
std::string remove_arg(const std::string& sentence,int i1,int i2);
std::string replace_arg(const std::string& sentence,int i,const std::string& new_val);
std::string str(double x);
std::string blank_filled(const std::string& s);
std::string inv_blank_filled(const std::string& s);
std::string to_space(const std::string& s,char c);
std::string spc(int n);
std::string left(const std::string& s,int n);
std::string mid(const std::string& s,int start,int n=0);
std::string replace(const std::string& s,const std::string& before,const std::string& after);
bool forward_match(const std::string& s0,const std::string& s);
std::string time();
std::string date();

#endif // #ifndef STRING_H_INCLUDED