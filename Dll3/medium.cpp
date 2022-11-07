#include "stdafx.h"
#include "medium.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

int medium::file_ver=0;

medium::medium() { 
	d=delta_d=dVariable=0;
}

medium::medium(int cn){
	d=delta_d=dVariable=0;
	g=cGlass(cn); gVariable=0;
}

bool operator==(const medium& a,const medium& b) {
	bool x;
	x=(a.d==b.d);
	x=x&&(a.delta_d==b.delta_d);
	x=x&&(a.dVariable==b.dVariable);
	x=x&&(a.g==b.g);
	x=x&&(a.gVariable==b.gVariable);
	return x;
}

medium medium::reversed() {
	medium x=*this;
	if( x.g.Name()[0]=='-' ) x.g.Name().erase(0,1); else x.g.Name()='-'+x.g.Name();
	x.d=-d;
	return x;
}

std::ostream& operator<<(std::ostream& to,medium x) {
	to<<x.d<<' '<<x.delta_d<<' '<<x.dVariable<<std::endl;
	to<<blank_filled(x.g.Name())<<' '<<x.gVariable<<std::endl;
	return to;
}

std::istream& operator>>(std::istream& from,medium& x) {
	switch(medium::file_ver){
		case 138:
			from>>x.d>>x.delta_d>>x.dVariable;
			from>>x.g.Name(); x.g.Name()=inv_blank_filled(x.g.Name());
			from>>x.gVariable;
			return from;
			break;
		case 137:
		case 136:
		case 135:
		case 134:	
		case 133:
		case 132:
		case 131:
		case 130:
		case 129:
		case 128:
		case 127:
		case 126:
		case 125:
		case 124:
		case 123:
		case 122:
		case 121:
		case 120:
		case 119:
		case 118:
		case 117:
		case 116:
		case 115:
		case 114:
		case 113:
		case 112:
		case 111:
		case 110:
		case 109:
			from>>x.d>>x.delta_d>>x.dVariable;
			from>>x.g.Name(); x.g.Name()=inv_blank_filled(x.g.Name());
			return from;
			break;
		case 108:
		case 107:
		case 106:
		case 105:
		case 104:
		case 103:
		case 102:
		case 101:
		case 100:
		default:
			from>>x.d>>x.delta_d>>x.dVariable;
			from>>x.g.Name();
			return from;
			break;
	}
}

void medium::scale(double m) {
	if(m==0) return;
	d*=m;
}

