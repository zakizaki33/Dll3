#include "stdafx.h"
#include "list.h"

/////////  stringlist members  ///////////////////////////////////////////////////////////////
void stringlist::Add(const std::string& s){
	int i; std::string buf;
	for(i=0; i<=(int)(s.length())-1; i++){
		if(s[i]!='\n') { buf+=s[i]; }
		else           { AddTail(buf); buf=""; } 
	}		
}


