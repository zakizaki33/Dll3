#include "stdafx.h"
#include "cLens1.h"

bool cLens1::is_numeric(std::string &s){
	// ◆sが数を表す文字列とき，trueを返す．
	// ◆次に，いずれかのrem(i)と一致するときは，s=str(i)としtrueを返す．
	//     scmd()のコマンドの面番号を表す引数に使用すると，
	//     面の追加や削除をしてもコマンドを修正しないで済む．
	//     "obj","img",いずれかのrem(i) と一致すると無条件に面番号に変換するので，
	//     不具合の原因にはなりうる（本関数を削除するだけで元に戻る）.
	//     したがって，面の指定以外コマンドの引数にこれらを使わないようにする必要がある．
	// ◆次に，cOptics::is_numeric(s)により，cmd(s)の結果が数である場合，
	//   s=cmd(s)としてtrueを返す．
	int i;

	if(s==""){
		return false;
	}
	if(::is_numeric(s)){
		return true;
	}
	{
		for(i=1; i<=k; i++){
			if(s==rem(i)){
				s=str(i);  // 同じrem(i)が複数の場合，より若い番号となる．
				return true;
			}
		}
		if(s=="obj"){
			s="0";
			return true;
		}
		if(s=="img"){
			s=str(k+1);
			return true;
		}
	}	
	if(cOptics::is_numeric(s)){
		return true;
	}
	return false;
}

int cLens1::IsCoefficient(const std::string& com){
	// comの第0引数(先頭のword)が収差係数の種類を表すものかどうか
	// 戻り値 ： 真の場合は，wordの数 (本来は1でよいが，Len.xlsで収差係数の表を作るときの都合による）
	//           偽の場合は0
	std::string s0;
	s0=arg(com,0);

	if(s0=="hQ"      || s0=="hq"      ||
	   s0=="hQp"     || s0=="hqp"     ||
	   s0=="SA"      || s0=="sa"      ||
	   s0=="CM"      || s0=="cm"      ||
	   s0=="AS"      || s0=="as"      ||
	   s0=="DS"      || s0=="ds"      ||
	   s0=="PT"      || s0=="pt"      ||
	   s0=="LC"      || s0=="lc"      ||
	   s0=="TC"      || s0=="tc"      ||
	   s0=="LC2"     || s0=="lc2"     ||
	   s0=="TC2"     || s0=="tc2"     ||
	   s0=="SAP"     || s0=="sap"     ||
	   s0=="CMP"     || s0=="cmp"     ||
	   s0=="ASP"     || s0=="asp"     ||
	   s0=="DSP"     || s0=="dsp"     ||
	   s0=="LCP"     || s0=="lcp"     ||
	   s0=="SA5"     || s0=="sa5"     ||
	   s0=="CM41"    || s0=="cm41"    ||
	   s0=="CM41Z"   || s0=="cm41z"   ||
	   s0=="CM41ALL" || s0=="cm41all" ||
	   s0=="SA32"    || s0=="sa32"    ||
	   s0=="SA32F"   || s0=="sa32f"   ||
	   s0=="SA32Z"   || s0=="sa32z"   ||
	   s0=="SA32ALL" || s0=="sa32all" ||
	   s0=="CM23"    || s0=="cm23"    ||
	   s0=="CM23P"   || s0=="cm23p"   ||
	   s0=="CM23Z"   || s0=="cm23z"   ||
	   s0=="CM23ALL" || s0=="cm23all" ||
	   s0=="AS5"     || s0=="as5"     ||
	   s0=="SG5"     || s0=="sg5"     ||
	   s0=="DS5"     || s0=="ds5"     ||
	   s0=="PRE"     || s0=="pre"     ||
	   s0=="DSE1"    || s0=="dse1"    ||
	   s0=="DSE2"    || s0=="dse2"    ||
	   s0=="ASE"     || s0=="ase"     ||
	   s0=="PTE"     || s0=="pte"     ||
	   s0=="CME"     || s0=="cme"        )
	{
		return words(com);
	}
	else{
		return 0;
	}
}

std::string cLens1::scmd(std::string com,int val){
	// val=trueのとき，一部のコマンドはBasic側にてval関数で処理するために
	// 数値を表す文字列のみ返す．
	std::string s;
	char buf[1000];
	std::string s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11;
	bool b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11;

	s0=arg(com,0);

	s+=scmd_general_func(com,val);  if(s0!="??" && s!="") return s;
	s+=cOptics::scmd(com,val);      if(s0!="??" && s!="") return s;

	if(s0=="??"){
		s+="SA CM AS DS PT LC TC LC2 TC2 SAP CMP ASP DSP LCP ";
		s+="SA5 CM41 CM41Z CM41ALL SA32 SA32F SA32Z SA32ALL ";
		s+="CM23 CM23P CM23Z CM23ALL AS5 SG5 DS5 ";
		s+="PRE DSE1 DSE2 ASE PTE CME ";
		s+="\n";
		s+="a10 a11 a12 a13 a14 a15 a16 a18 a1 a2 a20 a3 a4 a5 a6 a7 a8 a9 ";
		s+="AbyA0 aCOA Add2ndOrderTerm AdjustFocalLength Afocal Al Alp alpha_h_table ";
		s+="AxialRadiusY Axis ";
		s+="b bCOA bf BFOverF2 bfRatio bTerms ";
		s+="caCOA cbCOA cc cc1 Chamfer ChirpedPulseWidth CHMx CHMy Clearance CM0 CMatrix CMatrixStr ";
		s+="CoatName CoatReverse ";
		s+="Color ConjugateLength ct Cyl Cyl0Deg Cyl45Deg Cyl50 Cyl70 Cyl90 ";
		s+="d da DAbsC dc DconA4 DconA6 DconA8 DconA10 DconA12 DconA14 DconA16 DconA18 DconA20 ";
		s+="DconRn DconSetRn DconToFreeForm DconToPSeries ";
		s+="DeadSpace DefaultBC delta delta1 DeltaH DeltaM DeltaM70 DeltaS DeltaS70 ";
		s+="difforder DiffractObj DispersionTable DispersionToZero Dist Dist70 Distance DistanceX DistanceY DistanceZ dkp ";
		s+="DLSA2to3 dN dNewton dNewtonTol dpower DSagDFringe ";
		s+="DSDMs DSToDM DSToDM70 dx dx1 DXxmax dy dy1 DzToD ";
		s+="DYpr DYpr50 DYpr75 DYunsymmetric DYymax DYymin dz dz1 ";
		s+="e e1 EAdx EAdy EAx EAy EllipseRatio EllipsePhi EncircledEnergy EPCalc EPD EPDx EPx EPy ExitAngle ";
		s+="ExitDirectionCosineX ExitDirectionCosineY ExitDirectionCosineZ ExitPupilDia ExitPupilZ ExitTanX ExitTanY ExpandDist ";
		s+="f ff FFOverF2 ffRatio fideal FNumber FNumberObj FNumberObjXYAve FNumberParaxial FNumberXYAve ";
		s+="FocalLengths Footprint FootprintXdia FootprintYdia FresnelNumber fx fy ";
		s+="g1 g1_hat gamma GaussBeamTruncatedPower GDD gname grating grx gry grz gpitch g_hat ";
		s+="H Hp ";
		s+="ImageInfinity ImagePlaneTiltX IncidentAmplitude ";
		s+="IncidentAngle IncidentAngleMax IncidentAnglePrinMin IncidentAngleX IncidentAngleY IncidentAngleYFanMin ";
		s+="IncidentIntensity IncidentPhase IncidentPointGlobal Index Inflection InflectionPoint ";
		s+="Koba kpx kpy ";
		s+="LCPar LeC LegendreR0 LegendreToFreeForm LensDataZemax LensCenterPowerWorn ";
		s+="LSA LSA2nd LSA3rd LSA50 LSA70 LSAp LSAp2nd LSAp50 LSAp70 LSAps LSAs ";
		s+="M MakeDistChart MakeSpot MartinCyl MartinEqGeneral MartinSph MinMTF ";
		s+="ModelEyeFundusR ModelEyeS1 Mpupil Mpx Mpy Mstop Mstop1 ";
		s+="MTFm MTFs MTFsmave MTFx MTFxyave MTFy Mx My ";
		s+="N NA NAObj NAObjXYAve NAParaxial NAXYAve Nd Newton NewtonTol NewtonToR NNuActual nodal NormH ";
		s+="Note NotThruSurfaces Nud ";
		s+="OffAxialConicAB OffAxialMirrorComa ";
		s+="OPD OPD2 OPDPV OPDRMS OPDRMS0 OPDs OPL OptimizedDefocus OptimizeS1fix OSC OSC70 OSC50 OSCp ";
		s+="ParaxialValues PerturbDxDy pideal power PreformKoba PreformR PSeriesToDcon PSeriesToFreeForm PupilToPupil ";
		s+="qBend qValue ";
		s+="r RayPosX RayPosY RayPosZ ReduceAsphTerms RingsSpecRefSurf RmsPhi RmsPhiOff RmsPhiOn ";
		s+="RotateBlock RotateBlockAroundPupil RotateBlockX RotateBlockY RotateSurface RoundA ";
		s+="rox rox1 roy roy1 roz roz1 RToNewton rx ry ";
		s+="s s1 s1fix s1i S1x S1y SA0 SC ScheimpflugImagePlaneTiltX SchiempflugImagePlaneTiltY ";
		s+="SetHyperboloid SetHyperboloid2 SetM SetParaboloid SetParaboloid2 SetSpheroid SetSpheroid2 SetSpheroid3 ";
		s+="SetSpline ShapeFactor si SingleRayTrace ";
		s+="Sph SphEq SphEq50 SphEq70 SplineDoubleN SplineN SplineH SplineSetHStep SplineZ SPO SPO2 ";
		s+="SpotXGravityCenter SpotYGravityCenter Steepness StopDominate ";
		s+="StrehlDef SurfaceSag SurfaceSagTable SurfaceSagMax SurfaceSagMin ";
		s+="SurfaceSlope SurfaceSlopeMax SurfaceSlopeTable SwapObjPupil ";
		s+="t t1 t1i TangentialCurvature TangentialCurvatureMax tCalc tCOA TelecentricityObj TelephotoF1 TelephotoF2 ";
		s+="ThinT thNewton ti tLack ToAplanaticSurf ToBlock ToConcentricSurf ";
		s+="Toff ToIdealLens ToLMTestLens Ton ToroidRa ToroidZ ToScheimpflugImagePlane ";
		s+="TotalInAirThickness TotalOpticalThickness TotalThickness ToThinLens ";
		s+="TransformACoefficients Transmittance var var1 var2 var3 var4 var5 vdiopter1i vertex vertex1 VertexGlobal Wl ";
		s+="xObjectMax xScan xUsedRange yObjectMax yObjectMaxAng yScan yUsedRange yVignetting ";
		s+="ZC ZernikeC ZernikeCHigh ZernikeMaxOrder ZernikeR0 Zoom ZoomCamTable zScan Zvalue ";
		s+="\n";

		return s;
	}

	if( IsCoefficient(s0) )
	{   // 収差係数　/////////////////////////////////  
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s+=s0+" (no arguments)\n";  // 全系の和
			s+=s0+" i1 [i2]\n";         // i1面の値，またはi1面からi2面までの和
			s+=s0+" e\n";               // 指数形式で表すときの指数
			s+=s0+" rms [i1 i2]\n";     // 全系のrms，またはi1面からi2面までのrms

			// 収差係数を表等で各面の係数を指数一定の仮数にそろえるときは，
			// 例えば"SA e"で指数eを取得し，10^eで割る．
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",Coefficient(s0,1,0,0,i1,i2)); s=buf;
		}
		else if(b1){
			i1=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Coefficient(s0,1,0,0,i1,0)); s=buf;
		}
		else if(s1=="e"){
			sprintf(buf,"%.15g\n",Coefficient(s0,0,1,1,0,0)); s=buf;
		}
		else if(s1=="rms"){
			if(b2 && b3){
				i1=atoi(s2.c_str());
				i2=atoi(s3.c_str());
					sprintf(buf,"%.15g\n",Coefficient(s0,2,0,0,i1,i2)); s=buf;
			}
			else{
				sprintf(buf,"%.15g\n",Coefficient(s0,2,0,1,0,0)); s=buf;
			}
		}
		else{
			sprintf(buf,"%.15g\n",Coefficient(s0,1,0,1,0,0)); s=buf;
		}

		return s;
	}

	{   // double型プロパティ ///////////////////////////
		double *p=0;
		int args;
		std::string s1_0;
		
		s1_0=arg(com,1);       // property関数によってcomは変更されるので，あらかじめ第1引数("?"かどうか)を保存．
		p=property(com,args);  // double型のプロパティを扱う関数propertyは自動設計と共用する．
		                       //  (関数propertyによって，データ以外の部分はcomから削除される．
		                       //   r(1)=10の場合，com="r 1 10" -> com="10")
		if(s1_0=="?"){
			if(args==0){
				s+=s0+" (no arguments)\n";
				s+=s0+" new_value\n";
				return s;
			}
			if(args==1){
				s+=s0+" i\n";
				s+=s0+" i new_value\n";
				return s;
			}
			if(args==2){
				s+=s0+" i j\n";
				s+=s0+" i j new_value\n";
				return s;
			}
			if(args==3){
				s+=s0+" i m n\n";
				s+=s0+" i m n new_value\n";
				return s;
			}
		}
		else if(p!=0){
			s1=arg(com,0); b1=is_numeric(s1);

			if(b1){
				*p=atof(s1.c_str());
				sprintf(buf,"%.15g\n",*p); s=buf;
			}
			else{
				sprintf(buf,"%.15g\n",*p); s=buf;
			}
			return s;
		}
	}

	{   // 0変数プロパティ /////////////////////////////
		int *ip=0;
		std::string *sp=0;
		int ival;
		std::string sval;

		if     (s0=="Afocal" || s0=="afocal") ip=&Afocal;
		else if(s0=="IgnoreTC" || s0=="ignoretc") ip=&IgnoreTC;
		else if(s0=="Note" || s0=="note") sp=&Note;
		else if(s0=="nSpot" || s0=="nspot") ip=&nSpot;
		else if(s0=="stop") ip=&stop;
		else if(s0=="StopDominate" || s0=="stopdominate") ip=&StopDominate;

		if(ip!=0 || sp!=0){
			s1=arg(com,1); b1=is_numeric(s1);

			if(s1=="?"){
				s+=s0+" (no arguments)\n";
				s+=s0+" new_value\n";
			}
			else{
				if(sp!=0){
					if(args(com)==1){
						sval=s1;
						*sp=sval;
						sprintf(buf,"%s\n",(*sp).c_str()); s=buf;
					}
					else{
						sprintf(buf,"%s\n",(*sp).c_str()); s=buf;
					}
				}
				else if( b1 ) {
					if(ip!=0){
						ival=atoi(s1.c_str());
						*ip=ival;
						sprintf(buf,"%d\n",*ip); s=buf;
					}
				}
				else{
					if(ip!=0){
						sprintf(buf,"%d\n",*ip); s=buf;
					}
				}
			}

			return s;
		}
	}

	{   // 1変数プロパティ //////////////////////////////
		int& (cLens1::*ifp)(int)=0;
		std::string& (cLens1::*sfp)(int)=0;
		int i;
		int ival;
		std::string sval;
		
		if     (s0=="asph_type") ifp=&cLens1::asph_type;   // 古いコンパイラ(VC6など)であれば，ifp=asph_type; でも通る
		else if(s0=="CoatName"    || s0=="coatname")    sfp=&cLens1::CoatName;
		else if(s0=="CoatReverse" || s0=="coatreverse") ifp=&cLens1::CoatReverse;		
		else if(s0=="decenter_type")   ifp=&cLens1::decenter_type; 
		else if(s0=="grating")   ifp=&cLens1::grating;
		else if(s0=="difforder") ifp=&cLens1::difforder;

		if(ifp!=0 || sfp!=0){
			s1=arg(com,1); b1=is_numeric(s1);
			s2=arg(com,2); b2=is_numeric(s2);

			if(s1=="?"){
				s+=s0+" i\n";
				s+=s0+" i new_value\n";
			}
			else if(b1){
				i=atoi(s1.c_str());

				if(sfp!=0){
					if(args(com)==2){
						sval=s2;
						(this->*sfp)(i)=sval;
						sprintf(buf,"%s\n",(this->*sfp)(i).c_str()); s=buf;
					}
					else{
						sprintf(buf,"%s\n",(this->*sfp)(i).c_str()); s=buf;
					}
				}
				else if(b2){	
					if(ifp!=0 ){
						ival=atoi(s2.c_str());
						// ifpは int& (cLens1::*)(int) 型であり，
						// メンバ関数をifpで呼ぶときはインスタンス(ここでは*this)の指定を省略できない．
						// 演算子 ->* を使う．
						(this->*ifp)(i)=ival;
						sprintf(buf,"%d\n",(this->*ifp)(i)); s=buf;
					}
				}
				else{
					i=atoi(s1.c_str());
					if(ifp!=0){
						sprintf(buf,"%d\n",(this->*ifp)(i)); s=buf;
					}
				}
			}

			return s;
		}
	}

	if(s0=="AbyA0" || s0=="abya0"){
		s1=arg(com,1);
		if(s1=="?"){
			s="AbyA0 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",AbyA0()); s=buf;
		}
		return s;
	}
	if(s0=="Add2ndOrderTerm" || s0=="add2ndorderterm"){
		int i;
		double a2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Add2ndOrderTerm i a2\n";
		}
		else if( b1 && b2 ){
			i=atoi(s1.c_str());
			a2=atof(s2.c_str());
			Add2ndOrderTerm(i,a2);
		}
		return s;
	}
	if(s0=="AdjustFocalLength" || s0=="adjustfocallength"){
		double fl;
		int i1,i2, KeepD,ByAllSystem,WithEA;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="AdjustFocalLength fl i1 i2 [KeepD=1 ByAllSystem=0 WithEA=0]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6){
			fl=atof(s1.c_str());
			i1=atoi(s2.c_str());
			i2=atoi(s3.c_str());
			KeepD=atoi(s4.c_str());
			ByAllSystem=atoi(s5.c_str());
			WithEA=atoi(s6.c_str());
			AdjustFocalLength(fl,i1,i2,KeepD,ByAllSystem,WithEA);
		}
		else if(b1 && b2 && b3){
			fl=atof(s1.c_str());
			i1=atoi(s2.c_str());
			i2=atoi(s3.c_str());
			AdjustFocalLength(fl,i1,i2);
		}
		return s;
	}
	if(s0=="Al" || s0=="al"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Al i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Al(i)); s=buf;
		}
		return s;
	}
	if(s0=="Alp" || s0=="alp"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Alp i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Alp(i)); s=buf;
		}
		return s;
	}
	if(s0=="alpha_h_table"){
		s1=arg(com,1);
		if(s1=="?"){
			s="alpha_h_table (no argements)\n";
		}
		else {
			s=alpha_h_table();
		}
		return s;
	}
	if(s0=="AxialRadiusY" || s0=="axialradiusy"){
		int i;
		double y;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="AxialRadiusY i y\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			y=atof(s2.c_str());
			sprintf(buf,"%.15g\n",AxialRadiusY(i,y)); s=buf;
		}
		return s;
	}
	if(s0=="Axis" || s0=="axis"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="Axis [yObj=yObjectMax xObj=xObjectMax [findpupil=0]]\n";
		}
		else if(b1 && b2 && b3) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",Axis(yObj,xObj,findpupil)); s=buf;
		}
		else if(b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Axis(yObj,xObj)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Axis()); s=buf;
		}
		return s;
	}
	if(s0=="bf"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="bf [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",bf(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",bf(1,k)); s=buf;
		}
		return s;
	}
	if(s0=="BFOverF2" || s0=="bfoverf2"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="BFOverF2 [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",BFOverF2(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",BFOverF2()); s=buf;
		}
		return s;
	}
	if(s0=="bfRatio" || s0=="bfratio"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="bfRatio [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",bfRatio(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",bfRatio()); s=buf;
		}
		return s;
	}
	if(s0=="bTerms" || s0=="bterms"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="bTerms i\n";
		}
		else if( b1 ){
			i=atoi(s1.c_str());
			s=Get_bTerms(i)+'\n';
		}
		return s;
	}
	if(s0=="cc"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="cc i\n";
		}
		else if( b1 ){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",cc(i)); s=buf;
		}
		return s;
	}
	if(s0=="cc1"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="cc1 i\n";
		}
		else if( b1 ){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",cc1(i)); s=buf;
		}
		return s;
	}
	if(s0=="Chamfer" || s0=="chamfer"){
		int i;
		double Wx,Wy,Cx,Cy,Offset;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="Chamfer i Wx Wy Cx Cy Offset\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6) {
			i=atoi(s1.c_str());
			Wx=atof(s2.c_str());
			Wy=atof(s3.c_str());
			Cx=atof(s4.c_str());
			Cy=atof(s5.c_str());
			Offset=atof(s6.c_str());
			Chamfer(i,Wx,Wy,Cx,Cy,Offset);
		}
		return s;
	}
	if(s0=="ChirpedPulseWidth" || s0=="chirpedpulsewidth"){
		int i1,i2,j;
		double dt0;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="ChirpedPulseWidth i1 i2 j dt0_fs\n";
		}
		else if( b1 && b2 && b3 && b4 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			j =atoi(s3.c_str());
			dt0=atof(s4.c_str());
			sprintf(buf,"%.15g\n",ChirpedPulseWidth(i1,i2,j,dt0)); s=buf;
		}
		return s;
	}
	if(s0=="Clearance" || s0=="clearance"){
		int iPoint,iLine,FindPupil;
		double yObjPoint,xObjPoint,yObjLine,xObjLine;
		std::string SetRayPoint,SetRayLine;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s="Clearance iPoint yObjPoint xObjPoint SetRayPoint iLine yObjLine xObjLine SetRayLine FindPupil\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7 && b9){
			iPoint=atoi(s1.c_str());
			yObjPoint=atof(s2.c_str());
			xObjPoint=atof(s3.c_str());
			SetRayPoint=s4;
			iLine=atoi(s5.c_str());
			yObjLine=atof(s6.c_str());
			xObjLine=atof(s7.c_str());
			SetRayLine=s8;
			FindPupil=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",Clearance(iPoint,yObjPoint,xObjPoint,SetRayPoint,iLine,yObjLine,xObjLine,SetRayLine,FindPupil));
			s=buf;
		}
		return s;
	}
	if(s0=="CM0" || s0=="cm0"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="CM0 i1 i2 (calculated value)\n";
			s+="CM0 i     (surface property)\n";
			s+="CM0       (calculated value)\n";
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",CM0(i1,i2)); s=buf;
		}
		else if( b1 ){
			i1=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",CM0(i1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",CM0()); s=buf;
		}
		return s;
	}
	if(s0=="CMatrix" || s0=="cmatrix"){
		int i1,i2,i,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="CMatrix i1 i2 i j\n";
		}
		else if( b1 && b2 && b3 && b4 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			i =atoi(s3.c_str());
			j =atoi(s4.c_str());
			sprintf(buf,"%.15g\n",CMatrix(i1,i2,i,j)); s=buf;
		}
		return s;
	}
	if(s0=="CMatrixStr" || s0=="cmatrixstr"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="CMatrixStr i1 i2\n";
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			s=CMatrixStr(i1,i2);
		}
		return s;
	}
	if(s0=="Color" || s0=="color"){
		int j;
		std::string new_val;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2);
		if(s1=="?"){
			s="Color j [new_val]\n";
		}
		else if( b1 ){
			j=atoi(s1.c_str());
			new_val=s2;
			if(new_val!="") Set_color(j,new_val);
			sprintf(buf,(Get_color(j)+"\n").c_str()); s=buf;
		}
		return s;
	}
	if(s0=="ConjugateLength" || s0=="conjugatelength"){
		s1=arg(com,1);
		if(s1=="?"){
			s="ConjugateLength (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",ConjugateLength()); s=buf;
		}
		return s;
	}
	if(s0=="ct"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="ct i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",ct(i));
			s=buf;
		}
		return s;
	}
	if(s0=="Cyl" || s0=="cyl"){
		double yObj,xObj;
		int findpupil,i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Cyl [yObj=yObjectMax xObj=xObjectMax [findpupil=0 [i=0]]]\n";
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			i=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",Cyl(yObj,xObj,findpupil,i)); s=buf;
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",Cyl(yObj,xObj,findpupil)); s=buf;
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Cyl(yObj,xObj)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Cyl()); s=buf;
		}
		return s;
	}
	if(s0=="Cyl0Deg" || s0=="cyl0deg"){
		double yObj,xObj;
		int findpupil,i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Cyl0Deg [yObj=yObjectMax xObj=xObjectMax [findpupil=0 [i=0]]]\n";
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			i=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",Cyl0Deg(yObj,xObj,findpupil,i)); s=buf;
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",Cyl0Deg(yObj,xObj,findpupil)); s=buf;
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Cyl0Deg(yObj,xObj)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Cyl0Deg(yObjectMax,xObjectMax)); s=buf;
		}
		return s;
	}
	if(s0=="Cyl45Deg" || s0=="cyl45deg"){
		double yObj,xObj;
		int findpupil,i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Cyl45Deg [yObj=yObjectMax xObj=xObjectMax [findpupil=0 [i=0]]]\n";
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			i=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",Cyl45Deg(yObj,xObj,findpupil,i)); s=buf;
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",Cyl45Deg(yObj,xObj,findpupil)); s=buf;
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Cyl45Deg(yObj,xObj)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Cyl45Deg(yObjectMax,xObjectMax)); s=buf;
		}
		return s;
	}
	if(s0=="Cyl50" || s0=="cyl50"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Cyl50 (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",Cyl50()); s=buf;
		}
		return s;
	}
	if(s0=="Cyl70" || s0=="cyl70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Cyl70 (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",Cyl70()); s=buf;
		}
		return s;
	}
	if(s0=="Cyl90" || s0=="cyl90"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Cyl90 (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",Cyl90()); s=buf;
		}
		return s;
	}
	if(s0=="Cylsgn" || s0=="cylsgn"){
		double yObj,xObj;
		int findpupil,i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="CylSgn [yObj=yObjectMax xObj=xObjectMax [findpupil=0 [i=0]]]\n";
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			i=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",CylSgn(yObj,xObj,findpupil,i)); s=buf;
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",CylSgn(yObj,xObj,findpupil)); s=buf;
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",CylSgn(yObj,xObj)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",CylSgn(yObjectMax,xObjectMax)); s=buf;
		}
		return s;
	}
	if(s0=="da" ){
		int i,n;
		double dz;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="da i n dz\n";
		}
		else if(b1 && b2 && b3) {
			i=atoi(s1.c_str());
			n=atoi(s2.c_str());
			dz=atof(s3.c_str());
			sprintf(buf,"%.15g\n",da(i,n,dz)); s=buf;
		}
		return s;
	}
	if(s0=="DAbsC" || s0=="dabsc"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DAbsC i1 i2\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",DAbsC(i1,i2)); s=buf;
		}
		return s;
	}
	if(s0=="dc" ){
		int i;
		double dz;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="dc i dz\n";
		}
		else if(b1 && b2) {
			i=atoi(s1.c_str());
			dz=atof(s2.c_str());
			sprintf(buf,"%.15g\n",dc(i,dz)); s=buf;
		}
		return s;
	}
	if(s0=="DconSetRn" || s0=="dconsetrn"){
		int i;
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DconSetRn i [val]\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			val=atof(s2.c_str());
			sprintf(buf,"%.15g\n",DconSetRn(i,val)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",DconSetRn(i)); s=buf;
		}
		return s;
	}
	if(s0=="DconToFreeForm" || s0=="dcontofreeform"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="DconToFreeForm i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			DconToFreeForm(i);
		}
		return s;
	}
	if(s0=="DconToPSeries" || s0=="dcontopseries"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="DconToPSeries i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			DconToPSeries(i);
		}
		return s;
	}
	if(s0=="DeadSpace" || s0=="deadspace"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DeadSpace [i1=1 i2=k]\n";
		}
		else if(b1 && b2){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",deadspace(i1,i2,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",deadspace(1,k,1)); s=buf;
		}
		return s;
	}
	if(s0=="DefaultBC" || s0=="defaultbc"){
		int i1,i2;
		double weight,min_koba,min_ct,max_steepness;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="DefaultBC [min_koba=1 min_ct=0.5 max_steepness=0.8 [weight=1 [i1=1 i2=k]]]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6){
			min_koba=atof(s1.c_str());
			min_ct=atof(s2.c_str());
			max_steepness=atof(s3.c_str());
			weight=atof(s4.c_str());
			i1=atoi(s5.c_str());
			i2=atoi(s6.c_str());

			s=DefaultBC(min_koba,min_ct,max_steepness,weight,i1,i2);
		}
		else if(b1 && b2 && b3 && b4){
			min_koba=atof(s1.c_str());
			min_ct=atof(s2.c_str());
			max_steepness=atof(s3.c_str());
			weight=atof(s4.c_str());
			s=DefaultBC(min_koba,min_ct,max_steepness,weight);
		}
		else if(b1 && b2 && b3){
			min_koba=atof(s1.c_str());
			min_ct=atof(s2.c_str());
			max_steepness=atof(s3.c_str());
			s=DefaultBC(min_koba,min_ct,max_steepness);
		}
		else{
			s=DefaultBC();
		}
		return s;
	}
	if(s0=="delta"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="delta [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",delta(i1,i2,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",delta(1,k,1)); s=buf;
		}
		return s;
	}
	if(s0=="delta1"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="delta1 [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",delta1(i1,i2,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",delta1(1,k,1)); s=buf;
		}
		return s;
	}
	if(s0=="DeltaH" || s0=="deltah"){
		int DYIs1DXIs2,FindPupil,j,j0;
		double yObj,xObj,ypNormalized,xpNormalized,defocus;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s="DeltaH DYIs1DXIs2 yObj xObj FindPupil ypNormalized xpNormalized defocus j j0\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9){
			DYIs1DXIs2=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			FindPupil=atoi(s4.c_str());
			ypNormalized=atof(s5.c_str());
			xpNormalized=atof(s6.c_str());
			defocus=atof(s7.c_str());
			j=atoi(s8.c_str());
			j0=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",DeltaH(DYIs1DXIs2,yObj,xObj,FindPupil,ypNormalized,xpNormalized,defocus,j,j0));
			s=buf;
		}
		return s;
	}
	if(s0=="DeltaM" || s0=="deltam"){
		double yObjNormalized;
		int findpupil,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="DeltaM [yObjNomalized=1 findpupil=1 [j=1]]\n";
		}
		else if(b1 && b2 && b3){
			yObjNormalized=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			j=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",DeltaM(yObjNormalized,findpupil,j)); s=buf;
		}
		else if(b1 && b2){
			yObjNormalized=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",DeltaM(yObjNormalized,findpupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",DeltaM()); s=buf;
		}
		return s;
	}
	if(s0=="DeltaM70" || s0=="deltam70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="DeltaM70 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",DeltaM70()); s=buf;
		}
		return s;
	}
	if(s0=="DeltaS" || s0=="deltas"){
		double yObjNormalized;
		int findpupil,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="DeltaS [yObjNomalized=1 findpupil=1 [j=1]]\n";
		}
		else if(b1 && b2 && b3){
			yObjNormalized=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			j=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",DeltaS(yObjNormalized,findpupil,j)); s=buf;
		}
		else if(b1 && b2){
			yObjNormalized=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",DeltaS(yObjNormalized,findpupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",DeltaS()); s=buf;
		}
		return s;
	}
	if(s0=="DeltaS70" || s0=="deltas70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="DeltaS70 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",DeltaS70()); s=buf;
		}
		return s;
	}
	if(s0=="DiffractObj" || s0=="diffractobj"){
		double yObj,xObj;
		int FindPupil,nSpot,DeleteFrom;
		s1=arg(com,1);   b1=is_numeric(s1);
		s2=arg(com,2);   b2=is_numeric(s2);
		s3=arg(com,3);   b3=is_numeric(s3);
		s4=arg(com,4);   b4=is_numeric(s4);
		s5=arg(com,5);   b5=is_numeric(s5);
		if(s1=="?"){
			s="DiffractObj yObj xObj FindPupil nSpot DeleteFrom\n";
		}
		else if(b1 && b2 && b3 && b4 && b5){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			nSpot=atoi(s4.c_str());
			DeleteFrom=atoi(s5.c_str());
			DiffractObj(yObj,xObj,FindPupil,nSpot,DeleteFrom);
		}
		return s;
	}
	if(s0=="DispersionTable" || s0=="dispersiontable"){
		std::string SetRay;
		double yObj,xObj,yPupil,xPupil,yPupil_principal,xPupil_principal,wl_start,wl_end;
		int wl_points;
		s1=arg(com,1);   b1=is_numeric(s1);
		s2=arg(com,2);   b2=is_numeric(s2);
		s3=arg(com,3);
		s4=arg(com,4);   b4=is_numeric(s4);
		s5=arg(com,5);   b5=is_numeric(s5);
		s6=arg(com,6);   b6=is_numeric(s6);
		s7=arg(com,7);   b7=is_numeric(s7);
		s8=arg(com,8);   b8=is_numeric(s8);
		s9=arg(com,9);   b9=is_numeric(s9);
		s10=arg(com,10); b10=is_numeric(s10);
		if(s1=="?"){
			s+="DispersionTable yObj xObj SetRay yPupil xPupil yPupil_principal xPupil_principal\n";
			s+="                     wl_start wl_end wl_points\n";
		}
		else if(b1 && b2 && b4 && b5 && b6 && b7 && b8 && b9 && b10){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3;
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			yPupil_principal=atof(s6.c_str());
			xPupil_principal=atof(s7.c_str());
			wl_start=atof(s8.c_str());
			wl_end=atof(s9.c_str());
			wl_points=atoi(s10.c_str());
			s=DispersionTable(yObj,xObj,SetRay,yPupil,xPupil,yPupil_principal,xPupil_principal,
				               wl_start,wl_end,wl_points);
		}
		return s;
	}
	if(s0=="DispersionToZero" || s0=="dispersiontozero"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="DispersionToZero i1 i2\n";
		}
		else if(b1 && b2){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			DispersionToZero(i1,i2);
		}
		return s;
	}
	if(s0=="Dist" || s0=="dist"){
		double yObjNormalized;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Dist [yObjNomalized=1 findpupil=1]\n";
		}
		else if(b1 && b2){
			yObjNormalized=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",Dist(yObjNormalized,findpupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Dist()); s=buf;
		}
		return s;
	}
	if(s0=="Dist70" || s0=="dist70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Dist70 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",Dist70()); s=buf;
		}return s;
	}
	if(s0=="dkp" ){
		int i;
		double dz;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="dkp i dz\n";
		}
		else if(b1 && b2) {
			i=atoi(s1.c_str());
			dz=atof(s2.c_str());
			sprintf(buf,"%.15g\n",dkp(i,dz)); s=buf;
		}
		return s;
	}
	if(s0=="DLSA2to3" || s0=="dlsa2to3"){
		s1=arg(com,1);
		if(s1=="?"){
			s="DLSA2to3 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",DLSA2to3()); s=buf;
		}
		return s;
	}
	if(s0=="dpower"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="dpower [i1 i2]\n";
		}
		else if(b1 && b2) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",dpower(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",dpower()); s=buf;
		}
		return s;
	}
	if(s0=="DSagDFringe" || s0=="dsagdfringe"){
		int i;
		double yObj,xObj;
		std::string SetRay;
		double yPupil,xPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s="DSagDFringe i [yObj=0 xObj=0 SetRay=\"principal\" yPupil=0 xPupil=0]\n";
		}
		else if( b1 && b2 && b4 && b5 ) {
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			sprintf(buf,"%.15g\n",DSagDFringe(i,yObj,xObj,SetRay,yPupil,xPupil)); s=buf;
		}
		else if( b1 ){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",DSagDFringe(i,0,0,"principal",0,0)); s=buf;
		}
		return s;
	}
	if(s0=="DSDMs" || s0=="dsdms"){
		int findpupil,cols,fine_pitch;
		double s_weight,m_weight;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);

		if(s1=="?"){
			s="DSDMs findpupil [s_weight=1 m_weight=1 [cols=1 [fine_pitch=0]]]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5){
			findpupil=atoi(s1.c_str());
			s_weight=atof(s2.c_str());
			m_weight=atof(s3.c_str());
			cols=atoi(s4.c_str());
			fine_pitch=atoi(s5.c_str());
			s=DSDMs(findpupil,s_weight,m_weight,cols,fine_pitch);
		}
		else if(b1 && b2 && b3 && b4){
			findpupil=atoi(s1.c_str());
			s_weight=atof(s2.c_str());
			m_weight=atof(s3.c_str());
			cols=atoi(s4.c_str());
			s=DSDMs(findpupil,s_weight,m_weight,cols);
		}
		else if(b1 && b2 && b3){
			findpupil=atoi(s1.c_str());
			s_weight=atof(s2.c_str());
			m_weight=atof(s3.c_str());
			s=DSDMs(findpupil,s_weight,m_weight);
		}
		else if(b1){
			findpupil=atoi(s1.c_str());
			s=DSDMs(findpupil);
		}
		return s;
	}
	if(s0=="DSToDM" || s0=="dstodm"){
		double yObjNormalized;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DSToDM [yObjNomalized=1 findpupil=1]\n";
		}
		else if(b1 && b2){
			yObjNormalized=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",DSToDM(yObjNormalized,findpupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",DSToDM()); s=buf;
		}
		return s;
	}
	if(s0=="DSToDM70" || s0=="dstodm70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="DSToDM70 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",DSToDM70()); s=buf;
		}
		return s;
	}
	if(s0=="DSToDM50" || s0=="dstodm50"){
		s1=arg(com,1);
		if(s1=="?"){
			s="DSToDM50 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",DSToDM50()); s=buf;
		}
		return s;
	}
	if(s0=="DXxmax" || s0=="dxxmax"){
		double yObj;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="DXxmax [yObj]\n";
		}
		else if( b1 ){
			yObj=atof(s1.c_str());
			sprintf(buf,"%.15g\n",DXxmax(yObj)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",DXxmax()); s=buf;
		}
		return s;
	}
	if(s0=="DYpr" || s0=="dypr"){
		double yObj;
		int j,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="DYpr [yObj=yObjectMax] j FindPuppil\n";
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			j=atoi(s2.c_str());
			FindPupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",DYpr(yObj,j,FindPupil)); s=buf;
		}
		else if(b1 && b2){
			j=atoi(s1.c_str());
			FindPupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",DYpr(j,FindPupil)); s=buf;
		}
		return s;
	}
	if(s0=="DYpr50" || s0=="dypr50"){
		int j,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DYpr50 j FindPupil\n";
		}
		else{
			j=atoi(s1.c_str());
			FindPupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",DYpr50(j,FindPupil)); s=buf;
		}
		return s;
	}
	if(s0=="DYpr75" || s0=="dypr75"){
		int j,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DYpr75 j FindPupil\n";
		}
		else{
			j=atoi(s1.c_str());
			FindPupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",DYpr75(j,FindPupil)); s=buf;
		}
		return s;
	}
	if(s0=="DYunsymmetric" || s0=="dyunsymmetric"){
		double yObj;
		int FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="DYunsymmetric yObj FindPuppil\n";
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			FindPupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",DYunsymmetric(yObj,FindPupil)); s=buf;
		}
		return s;
	}
	if(s0=="DYymax" || s0=="dyymax"){
		double yObj;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="DYymax [yObj]\n";
		}
		else if( b1 ){
			yObj=atof(s1.c_str());
			sprintf(buf,"%.15g\n",DYymax(yObj)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",DYymax()); s=buf;
		}
		return s;
	}
	if(s0=="DYymin" || s0=="dyymin"){
		double yObj;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="DYymin [yObj]\n";
		}
		else if( b1 ){
			yObj=atof(s1.c_str());
			sprintf(buf,"%.15g\n",DYymin(yObj)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",DYymin()); s=buf;
		}
		return s;
	}
	if(s0=="DzToD" || s0=="dztod"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="DzToD i\n";
		}
		else if( b1 ){
			i=atoi(s1.c_str());
			DzToD(i);
		}
		return s;
	}
	if(s0=="e"){
		int pre_i1,pre_i2,post_i1,post_i2;
		int i1,i2;
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s+="e pre_i1 pre_i2 post_i1 post_i2 [new_value]\n";
			s+="e i1 i2\n";
		}
		else if(b1 && b2 && b3 && b4 && b5){
			pre_i1 =atoi(s1.c_str());
			pre_i2 =atoi(s2.c_str());
			post_i1=atoi(s3.c_str());
			post_i2=atoi(s4.c_str());
			val=atof(s5.c_str());
			Set_e(pre_i1,pre_i2,post_i1,post_i2,val);
			sprintf(buf,"%.15g\n",Get_e(pre_i1,pre_i2,post_i1,post_i2)); s=buf;
		}
		else if(b1 && b2 && b3 && b4){
			pre_i1 =atoi(s1.c_str());
			pre_i2 =atoi(s2.c_str());
			post_i1=atoi(s3.c_str());
			post_i2=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",Get_e(pre_i1,pre_i2,post_i1,post_i2)); s=buf;
		}
		else if(b1 && b2){
			i1 =atoi(s1.c_str());
			i2 =atoi(s2.c_str());
			sprintf(buf,"%.15g\n",e(i1,i2)); s=buf;
		}
		return s;
	}
	if(s0=="e1"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="e1 i1 i2\n";
		}
		else if(b1 && b2){
			i1 =atoi(s1.c_str());
			i2 =atoi(s2.c_str());
			sprintf(buf,"%.15g\n",e1(i1,i2)); s=buf;
		}
		return s;
	}
	if(s0=="EllipseRatio" || s0=="ellipseratio"){
		int findpupil;
		double yObj,xObj,yPupil,xPupil,a,b,phi_deg;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="EllipseRatio [yObj=0 xObj=0 [SetRay=principal findpupil=1 yPupil=0 xPupil=0]]\n";
			s+="             a b phi_deg\n";

		}
		else if( b1 && b2 && b4 && b5 && b6 && b7 && b8 && b9){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3;
			findpupil=atoi(s4.c_str());
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			a=atof(s7.c_str());
			b=atof(s8.c_str());
			phi_deg=atof(s9.c_str());
			sprintf(buf,"%.15g\n",EllipseRatio(yObj,xObj,SetRay,findpupil,yPupil,xPupil,k,1,a,b,phi_deg));
			s=buf;
		}
		else if(b1 && b2 && b3 && b4 && b5){
			xObj=atof(s1.c_str());
			yObj=atof(s2.c_str());
			a=atof(s3.c_str());
			b=atof(s4.c_str());
			phi_deg=atof(s5.c_str());
			sprintf(buf,"%.15g\n",EllipseRatio(yObj,xObj,"principal",1,0,0,k,1,a,b,phi_deg));
			s=buf;
		}
		else if(b1 && b2 && b3){
			a=atof(s1.c_str());
			b=atof(s2.c_str());
			phi_deg=atof(s3.c_str());
			sprintf(buf,"%.15g\n",EllipseRatio(0,0,"principal",1,0,0,k,1,a,b,phi_deg));
			s=buf;
		}
		return s;
	}
	if(s0=="EllipsePhi" || s0=="ellipsephi"){
		int findpupil;
		double yObj,xObj,yPupil,xPupil,a,b,phi_deg;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="EllipsePhi [yObj=0 xObj=0 [SetRay=principal findpupil=1 yPupil=0 xPupil=0]]\n";
			s+="           a b phi_deg\n";

		}
		else if( b1 && b2 && b4 && b5 && b6 && b7 && b8 && b9){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3;
			findpupil=atoi(s4.c_str());
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			a=atof(s7.c_str());
			b=atof(s8.c_str());
			phi_deg=atof(s9.c_str());
			sprintf(buf,"%.15g\n",EllipsePhi(yObj,xObj,SetRay,findpupil,yPupil,xPupil,k,1,a,b,phi_deg));
			s=buf;
		}
		else if(b1 && b2 && b3 && b4 && b5){
			xObj=atof(s1.c_str());
			yObj=atof(s2.c_str());
			a=atof(s3.c_str());
			b=atof(s4.c_str());
			phi_deg=atof(s5.c_str());
			sprintf(buf,"%.15g\n",EllipsePhi(yObj,xObj,"principal",1,0,0,k,1,a,b,phi_deg));
			s=buf;
		}
		else if(b1 && b2 && b3){
			a=atof(s1.c_str());
			b=atof(s2.c_str());
			phi_deg=atof(s3.c_str());
			sprintf(buf,"%.15g\n",EllipsePhi(0,0,"principal",1,0,0,k,1,a,b,phi_deg));
			s=buf;
		}
		return s;
	}	
	if(s0=="EncircledEnergy" || s0=="encircledenergy"){
		double SensorPhi,yObj,xObj,defocus;
		int ColorStart,ColorEnd,FindPupil,OptimizeDefocus;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s+="EncircledEnergy SensorPhi yObj xObj [defocus=0 ColorStart=1 ColorEnd=1 FindPupil=1 OptimizeDefocus=1]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8){
			SensorPhi=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			defocus=atof(s4.c_str());
			ColorStart=atoi(s5.c_str());
			ColorEnd=atoi(s6.c_str());
			FindPupil=atoi(s7.c_str());
			OptimizeDefocus=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",
			        EncircledEnergy(SensorPhi,yObj,xObj,defocus,ColorStart,ColorEnd,FindPupil,OptimizeDefocus));
			s=buf;
		}
		else if(b1 && b2 && b3){
			SensorPhi=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			sprintf(buf,"%.15g\n",EncircledEnergy(SensorPhi,yObj,xObj,0,1,1,1,1));
			s=buf;
		}
		return s;
	}
	if(s0=="EPCalc" || s0=="epcalc"){
		s1=arg(com,1);
		if(s1=="?"){
			s="EPCalc (no argements)\n";
		}
		else{
			EPCalc();
		}
		return s;
	}
	if(s0=="ExitAngle" || s0=="exitangle"){
		double yObj,xObj,yPupil,xPupil;
		int i,j,InAir;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8);
		if(s1=="?"){
			s+="ExitAngle i [yObj=yObjctMax xObj=xObjctMax] SetRay yPupil xPupil j InAir\n";
			s+="ExitAngle i [yObj=yObjctMax xObj=xObjctMax] SetRay yPupil xPupil j (InAir=0)\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			InAir=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",ExitAngle(i,yObj,xObj,SetRay,yPupil,xPupil,j,InAir));
			s=buf;
		}
		else if(b1 && b3 && b4 && b5){
			i=atoi(s1.c_str());
			yObj=this->yObjectMax;
			xObj=this->xObjectMax;
			SetRay=s2;
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			j=atoi(s5.c_str());
			InAir=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",ExitAngle(i,yObj,xObj,SetRay,yPupil,xPupil,j,InAir));
			s=buf;
		}
		return s;
	}
	if(s0=="ExitDirectionCosineX" || s0=="exitdirectioncosinex"){
		int i,j, FindPupil;
		double yObj,xObj,yPupil,xPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="ExitDirectionCosineX i [yObj=yObjectMax xObj=xObjectMax] SetRay FindPupil yPupil xPupil j\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7 && b8){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			FindPupil=atoi(s5.c_str());
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",ExitDirectionCosineX(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		else if(b1 && b3 && b4 && b5 && b6){
			i=atoi(s1.c_str());
			SetRay=s2;
			FindPupil=atoi(s3.c_str());
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			j=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",ExitDirectionCosineX(i,yObjectMax,xObjectMax,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		return s;
	}
	if(s0=="ExitDirectionCosineY" || s0=="exitdirectioncosiney"){
		int i,j, FindPupil;
		double yObj,xObj,yPupil,xPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="ExitDirectionCosineY i [yObj=yObjectMax xObj=xObjectMax] SetRay FindPupil yPupil xPupil j\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7 && b8){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			FindPupil=atoi(s5.c_str());
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",ExitDirectionCosineY(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		else if(b1 && b3 && b4 && b5 && b6){
			i=atoi(s1.c_str());
			SetRay=s2;
			FindPupil=atoi(s3.c_str());
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			j=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",ExitDirectionCosineY(i,yObjectMax,xObjectMax,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		return s;
	}
	if(s0=="ExitDirectionCosineZ" || s0=="exitdirectioncosinez"){
		int i,j, FindPupil;
		double yObj,xObj,yPupil,xPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="ExitDirectionCosineZ i [yObj=yObjectMax xObj=xObjectMax] SetRay FindPupil yPupil xPupil j\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7 && b8){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			FindPupil=atoi(s5.c_str());
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",ExitDirectionCosineZ(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		else if(b1 && b3 && b4 && b5 && b6){
			i=atoi(s1.c_str());
			SetRay=s2;
			FindPupil=atoi(s3.c_str());
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			j=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",ExitDirectionCosineZ(i,yObjectMax,xObjectMax,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		return s;
	}
	if(s0=="ExitPupilDia" || s0=="exitpupildia"){
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="ExitPupilDia findpupil\n";
		}
		else if(b1){
			findpupil=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",ExitPupilDia(findpupil)); s=buf;
		}
		return s;
	}
	if(s0=="ExitPupilZ" || s0=="exitpupilz"){
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="ExitPupilZ findpupil\n";
		}
		else if(b1){
			findpupil=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",ExitPupilZ(findpupil)); s=buf;
		}
		return s;
	}
	if(s0=="ExitTanX" || s0=="exittanx"){
		int i,j, FindPupil;
		double yObj,xObj,yPupil,xPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="ExitTanX i [yObj=yObjectMax xObj=xObjectMax] SetRay FindPupil yPupil xPupil j\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7 && b8){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			FindPupil=atoi(s5.c_str());
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",ExitTanX(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		else if(b1 && b3 && b4 && b5 && b6){
			i=atoi(s1.c_str());
			SetRay=s2;
			FindPupil=atoi(s3.c_str());
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			j=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",ExitTanX(i,yObjectMax,xObjectMax,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		return s;
	}
	if(s0=="ExitTanY" || s0=="exittany"){
		int i,j, FindPupil;
		double yObj,xObj,yPupil,xPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="ExitTanY i [yObj=yObjectMax xObj=xObjectMax] SetRay FindPupil yPupil xPupil j\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7 && b8){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			FindPupil=atoi(s5.c_str());
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",ExitTanY(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		else if(b1 && b3 && b4 && b5 && b6){
			i=atoi(s1.c_str());
			SetRay=s2;
			FindPupil=atoi(s3.c_str());
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			j=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",ExitTanY(i,yObjectMax,xObjectMax,SetRay,FindPupil,yPupil,xPupil,j)); s=buf;
		}
		return s;
	}
	if(s0=="ExpandDist" || s0=="expanddist"){
		double yObjMax,xObjMax,defocus;
		int FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="ExpandDist yObjMax xObjMax FindPupil defocus\n";
		}
		else if( b1 && b2 && b3 && b4 ){
			yObjMax=atof(s1.c_str());
			xObjMax=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			defocus=atof(s4.c_str());
			s=ExpandDist(yObjMax,xObjMax,FindPupil,defocus);
		}
		return s;
	}
	if(s0=="f"){
		int i1,i2,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="f [i1=1 i2=k][j=1]\n";
		}
		else if(b1 && b2 && b3){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			j=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",f(i1,i2,j)); s=buf;
		}
		else if(b1 && b2){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",f(i1,i2)); s=buf;
		}
		else if(b1){
			j=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",f(j)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",f()); s=buf;
		}

		return s;
	}
	if(s0=="ff"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ff [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",ff(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",ff(1,k)); s=buf;
		}
		return s;
	}
	if(s0=="FFOverF2" || s0=="ffoverf2"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="FFOverF2 [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",FFOverF2(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",FFOverF2()); s=buf;
		}
		return s;
	}
	if(s0=="ffRatio" || s0=="ffratio"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ffRatio [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",ffRatio(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",ffRatio()); s=buf;
		}
		return s;
	}
	if(s0=="FNumber" || s0=="fnumber"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="FNumber [yObj=0 xObj=0 [findpupil=1]]\n";
		}
		else if( b1 && b2 && b3 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",FNumber(yObj,xObj,findpupil)); s=buf;
		}
		else if( b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",FNumber(yObj,xObj,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",FNumber(0,0,1)); s=buf;
		}
		return s;
	}
	if(s0=="FNumberObj" || s0=="fnumberobj"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="FNumberObj [yObj=0 xObj=0 [findpupil=1]]\n";
		}
		else if( b1 && b2 && b3 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",FNumberObj(yObj,xObj,findpupil)); s=buf;
		}
		else if( b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",FNumberObj(yObj,xObj,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",FNumberObj(0,0,1)); s=buf;
		}
		return s;
	}
	if(s0=="FNumberObjXYAve" || s0=="fnumberobjxyave"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="FNumberObjXYAve [yObj=0 xObj=0 [findpupil=1]]\n";
		}
		else if( b1 && b2 && b3 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",FNumberObjXYAve(yObj,xObj,findpupil)); s=buf;
		}
		else if( b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",FNumberObjXYAve(yObj,xObj,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",FNumberObjXYAve(0,0,1)); s=buf;
		}
		return s;
	}
	if(s0=="FNumberParaxial" || s0=="fnumberparaxial"){
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="FNumberParaxial [findpupil=1]\n";
		}
		else if( b1 ) {
			findpupil=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",FNumberParaxial(findpupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",FNumberParaxial(1)); s=buf;
		}
		return s;
	}
	if(s0=="FNumberXYAve" || s0=="fnumberxyave"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="FNumberXYAve [yObj=0 xObj=0 [findpupil=1]]\n";
		}
		else if( b1 && b2 && b3 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",FNumberXYAve(yObj,xObj,findpupil)); s=buf;
		}
		else if( b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",FNumberXYAve(yObj,xObj,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",FNumberXYAve(0,0,1)); s=buf;
		}
		return s;
	}
	if(s0=="FocalLengths" || s0=="focallengths"){
		int i1,i2, for_drawing;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="FocalLengths [i1=1 i2=k][for_drawing=0]\n";
		}
		else if(b1 && b2) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			if( b3 ){
				for_drawing=atoi(s3.c_str());
				s=Focallengths(i1,i2,for_drawing);
			}
			else{
				s=Focallengths(i1,i2,0);
			}
		}
		else{
			if(b1){
				for_drawing=atoi(s1.c_str());
				s=Focallengths(1,k,for_drawing);
			}
			else{
				s=Focallengths(1,k,0);
			}
		}
		return s;
	}
	if(s0=="Footprint" || s0=="footprint"){
		double yObj,xObj;
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="Footprint yObj xObj i\n";
		}
		else if(b1 && b2 && b3){
			double ymax,ymin,xmax,xmin,ydia,xdia,ycenter,xcenter;
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			i=atoi(s3.c_str());
			Footprint(ymax,ymin,xmax,xmin,ydia,xdia,ycenter,xcenter,yObj,xObj,i,1);
			sprintf(buf,"no.%d surf footprint yObj=%g xObj=%g\n",i,yObj,xObj); s+=buf;
			sprintf(buf," ymax=%g ymin=%g ydia=%g ycenter=%g\n",ymax,ymin,ydia,ycenter); s+=buf;
			sprintf(buf," xmax=%g xmin=%g xdia=%g xcenter=%g\n",xmax,xmin,xdia,xcenter); s+=buf;
		}
		return s;
	}
	if(s0=="FootprintXdia" || s0=="footprintxdia"){
		double yObj,xObj;
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="FootprintXdia yObj xObj i\n";
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			i=atoi(s3.c_str());
			sprintf(buf,"%.15g\n", FootprintXdia(yObj,xObj,i,1)); s=buf;
		}
		return s;
	}
	if(s0=="FootprintYdia" || s0=="footprintydia"){
		double yObj,xObj;
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="FootprintYdia yObj xObj i\n";
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			i=atoi(s3.c_str());
			sprintf(buf,"%.15g\n", FootprintYdia(yObj,xObj,i,1)); s=buf;
		}
		return s;
	}
	if(s0=="FresnelNumber" || s0=="fresnelnumber"){
		int findpupil;
		double defocus;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="FresnelNumber findpupil defocus\n";
		}
		else if(b1 && b2){
			findpupil=atoi(s1.c_str());
			defocus=atof(s2.c_str());
			sprintf(buf,"%.15g\n", FresnelNumber(findpupil,defocus)); s=buf;
		}
		return s;
	}
	if(s0=="fx"){
		s1=arg(com,1);
		if(s1=="?"){
			s="fx (no argements)\n";
		}
		else {
			sprintf(buf,"%.15g\n",fx()); s=buf;
		}
		return s;
	}
	if(s0=="fy"){
		s1=arg(com,1);
		if(s1=="?"){
			s="fy (no argements)\n";
		}
		else {
			sprintf(buf,"%.15g\n",fy()); s=buf;
		}
		return s;
	}
	if(s0=="g1"){
		s1=arg(com,1);
		if(s1=="?"){
			s="g1 (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",g1()); s=buf;
		}
		return s;
	}
	if(s0=="g1_hat"){
		s1=arg(com,1);
		if(s1=="?"){
			s="g1_hat (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",g1_hat()); s=buf;
		}
		return s;
	}
	if(s0=="gamma"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="gamma [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",gamma(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",gamma(1,k)); s=buf;
		}
		return s;
	}
	if(s0=="GaussBeamTruncatedPower" || s0=="gaussbeamtruncatedpower"){
		double BeamPhiX,BeamPhiY;
		double AperturePhiX,AperturePhiY;
		double ApertureDX,ApertureDY;
		int n;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s+="GaussBeamTruncatedPower";
			s+=" BeamPhiX BeamPhiY AperturePhiX AperturePhiY ApertureDX ApertureDY n\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7){
			double t;
			BeamPhiX=atof(s1.c_str());
			BeamPhiY=atof(s2.c_str());
			AperturePhiX=atof(s3.c_str());
			AperturePhiY=atof(s4.c_str());
			ApertureDX=atof(s5.c_str());
			ApertureDY=atof(s6.c_str());
			n=atoi(s7.c_str());
			t=GaussBeamTruncatedPower(BeamPhiX,BeamPhiY,AperturePhiX,AperturePhiY,ApertureDX,ApertureDY,n);
			sprintf(buf,"T=%g\n",t);
			s=buf;
		}
		return s;
	}
	if(s0=="GDD" || s0=="gdd"){
		int i1,i2,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="GDD i1 i2 j\n";
		}
		else if( b1 && b2 && b3 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			j =atoi(s3.c_str());
			sprintf(buf,"%.15g\n",GDD(i1,i2,j)); s=buf;
		}
		return s;
	}
	if(s0=="gname"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); 
		if(s1=="?"){
			s="gname i [new_value]\n";
		}
		else if(b1 && s2==""){
			i=atoi(s1.c_str());
			sprintf(buf,"%s\n",Get_gname(i).c_str()); s=buf;
		}
		else{
			i=atoi(s1.c_str());
			Set_gname(i,s2);
			sprintf(buf,"%s\n",Get_gname(i).c_str()); s=buf;
		}
		return s;
	}
	if(s0=="g_hat"){
		s1=arg(com,1);
		if(s1=="?"){
			s="g_hat (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",g_hat()); s=buf;
		}
		return s;
	}
	if(s0=="H" || s0=="h"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="H i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",H(i)); s=buf;
		}
		return s;
	}
	if(s0=="Hp" || s0=="hp"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Hp i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Hp(i)); s=buf;
		}
		return s;
	}
	if(s0=="ImageInfinity" || s0=="imageinfinity"){
		s1=arg(com,1);
		if(s1=="?"){
			s="ImageInfinity (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",ImageInfinity()); s=buf;
		}
		return s;
	}
	if(s0=="ImagePlaneTiltX" || s0=="imageplanetiltx"){
		double yObj,xObj,yPupil,xPupil;
		int i,FindPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="ImagePlaneTiltX i [yObj=0 xObj=0 SetRay=\"\" FindPupil=0 yPupil=0 xPupil=0]\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			FindPupil=atoi(s5.c_str());
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			sprintf(buf,"%.15g\n", ImagePlaneTiltX(i,yObj,xObj,SetRay,FindPupil,yPupil,xPupil)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n", ImagePlaneTiltX(i)); s=buf;
		}
		return s;
	}
	if(s0=="IncidentAmplitude" || s0=="incidentamplitude"){
		std::string xyz;
		double yObj,xObj,yPupil,xPupil, pol_phi,pol_ratio;
		int i,j;
		std::string SetRay;
		s1=arg(com,1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		s10=arg(com,10); b10=is_numeric(s10);
		if(s1=="?"){
			s+="IncidentAmplitude xyz i yObj xObj SetRay yPupil xPupil j [pol_phi=90 pol_ratio=0]\n";
		}
		else if(b2 && b3 && b4 && b6 && b7 && b8 && b9 && b10){
			xyz=s1;
			i=atoi(s2.c_str());
			yObj=atof(s4.c_str());
			xObj=atof(s4.c_str());
			SetRay=s5;
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			pol_phi=atof(s9.c_str());
			pol_ratio=atof(s10.c_str());
			sprintf(buf,"%.15g\n",IncidentAmplitude(xyz,i,yObj,xObj,SetRay,yPupil,xPupil,j,pol_phi,pol_ratio));
			s=buf;
		}
		else if(b2 && b3 && b4 && b6 && b7 && b8){
			xyz=s1;
			i=atoi(s2.c_str());
			yObj=atof(s4.c_str());
			xObj=atof(s4.c_str());
			SetRay=s5;
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",IncidentAmplitude(xyz,i,yObj,xObj,SetRay,yPupil,xPupil,j));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentAngle" || s0=="incidentangle"){
		double yObj,xObj,yPupil,xPupil;
		int i,j,InAir;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8);
		if(s1=="?"){
			s+="IncidentAngle i [yObj=yObjctMax xObj=xObjctMax] SetRay yPupil xPupil j InAir\n";
			s+="IncidentAngle i [yObj=yObjctMax xObj=xObjctMax] SetRay yPupil xPupil j (InAir=0)\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			InAir=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",IncidentAngle(i,yObj,xObj,SetRay,yPupil,xPupil,j,InAir));
			s=buf;
		}
		else if(b1 && b3 && b4 && b5){
			i=atoi(s1.c_str());
			yObj=this->yObjectMax;
			xObj=this->xObjectMax;
			SetRay=s2;
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			j=atoi(s5.c_str());
			InAir=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",IncidentAngle(i,yObj,xObj,SetRay,yPupil,xPupil,j,InAir));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentAngleMax" || s0=="incidentanglemax"){
		double yObj,xObj,yPupil,xPupil;
		int j;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="IncidentAngleMax yObj xObj SetRay yPupil xPupil j\n";
		}
		else if(b1 && b2 && b4 && b5 && b6){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3;
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			j=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",IncidentAngleMax(yObj,xObj,SetRay,yPupil,xPupil,j));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentAnglePrinMin" || s0=="incidentangleprinmin"){
		int i,j,FindPupil;
		double yObj1,yObj2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s="IncidentAnglePrinMin i [yObj1=yObjectMax yObj2=-yObjectMax] FindPupil j\n";
		}
		else if(b1 && b2 && b3 && b4 && b5){
			i=atoi(s1.c_str());
			yObj1=atof(s2.c_str());
			yObj2=atof(s3.c_str());
			FindPupil=atoi(s4.c_str());
			j=atoi(s5.c_str());
			sprintf(buf,"%.15g\n",IncidentAnglePrinMin(i,yObj1,yObj2,FindPupil,j));
			s=buf;
		}
		else if(b1 && b2 && b3){
			i=atoi(s1.c_str());
			FindPupil=atoi(s2.c_str());
			j=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",IncidentAnglePrinMin(i,yObjectMax,-yObjectMax,FindPupil,j));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentAngleX" || s0=="incidentanglex"){
		double yObj,xObj,yPupil,xPupil;
		int i,j,InAir;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8);
		if(s1=="?"){
			s+="IncidentAngleX i [yObj=yObjctMax xObj=xObjctMax] SetRay yPupil xPupil j InAir\n";
			s+="IncidentAngleX i [yObj=yObjctMax xObj=xObjctMax] SetRay yPupil xPupil j (InAir=0)\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			InAir=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",IncidentAngleX(i,yObj,xObj,SetRay,yPupil,xPupil,j,InAir));
			s=buf;
		}
		else if(b1 && b3 && b4 && b5){
			i=atoi(s1.c_str());
			yObj=this->yObjectMax;
			xObj=this->xObjectMax;
			SetRay=s2;
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			j=atoi(s5.c_str());
			InAir=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",IncidentAngleX(i,yObj,xObj,SetRay,yPupil,xPupil,j,InAir));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentAngleY" || s0=="incidentangley"){
		double yObj,xObj,yPupil,xPupil;
		int i,j,InAir;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8);
		if(s1=="?"){
			s+="IncidentAngleY i [yObj=yObjctMax xObj=xObjctMax] SetRay yPupil xPupil j InAir\n";
			s+="IncidentAngleY i [yObj=yObjctMax xObj=xObjctMax] SetRay yPupil xPupil j (InAir=0)\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			InAir=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",IncidentAngleY(i,yObj,xObj,SetRay,yPupil,xPupil,j,InAir));
			s=buf;
		}
		else if(b1 && b3 && b4 && b5){
			i=atoi(s1.c_str());
			yObj=this->yObjectMax;
			xObj=this->xObjectMax;
			SetRay=s2;
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			j=atoi(s5.c_str());
			InAir=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",IncidentAngleY(i,yObj,xObj,SetRay,yPupil,xPupil,j,InAir));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentAngleYFanMin" || s0=="incidentangleyfanmin"){
		int i,j,FindPupil;
		double yObj,xObj;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s="IncidentAngleYFanMin i yObj xObj FindPupil j\n";
		}
		else if(b1 && b2 && b4 && b5){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			FindPupil=atoi(s4.c_str());
			j=atoi(s5.c_str());
			sprintf(buf,"%.15g\n",IncidentAngleYFanMin(i,yObj,xObj,FindPupil,j));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentIntensity" || s0=="incidentintensity"){
		double yObj,xObj,yPupil,xPupil, pol_phi,pol_ratio;
		int i,j;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="IncidentAngle i yObj xObj SetRay yPupil xPupil j [pol_phi=90 pol_ratio=0]\n";
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7 && b8 && b9){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			pol_phi=atof(s8.c_str());
			pol_ratio=atof(s9.c_str());
			sprintf(buf,"%.15g\n",IncidentIntensity(i,yObj,xObj,SetRay,yPupil,xPupil,j,pol_phi,pol_ratio));
			s=buf;
		}
		else if(b1 && b2 && b3 && b5 && b6 && b7){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			sprintf(buf,"%.15g\n",IncidentIntensity(i,yObj,xObj,SetRay,yPupil,xPupil,j));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentPhase" || s0=="incidentphase"){
		std::string xyz;
		double yObj,xObj,yPupil,xPupil, pol_phi,pol_ratio;
		int i,j;
		std::string SetRay;
		s1=arg(com,1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		s10=arg(com,10); b10=is_numeric(s10);
		if(s1=="?"){
			s+="IncidentPhase xyz i yObj xObj SetRay yPupil xPupil j [pol_phi=90 pol_ratio=0]\n";
		}
		else if(b2 && b3 && b4 && b6 && b7 && b8 && b9 && b10){
			xyz=s1;
			i=atoi(s2.c_str());
			yObj=atof(s4.c_str());
			xObj=atof(s4.c_str());
			SetRay=s5;
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			pol_phi=atof(s9.c_str());
			pol_ratio=atof(s10.c_str());
			sprintf(buf,"%.15g\n",IncidentPhase(xyz,i,yObj,xObj,SetRay,yPupil,xPupil,j,pol_phi,pol_ratio));
			s=buf;
		}
		else if(b2 && b3 && b4 && b6 && b7 && b8){
			xyz=s1;
			i=atoi(s2.c_str());
			yObj=atof(s4.c_str());
			xObj=atof(s4.c_str());
			SetRay=s5;
			yPupil=atof(s6.c_str());
			xPupil=atof(s7.c_str());
			j=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",IncidentPhase(xyz,i,yObj,xObj,SetRay,yPupil,xPupil,j));
			s=buf;
		}
		return s;
	}
	if(s0=="IncidentPointGlobal" || s0=="incidentpointglobal"){
		int i,FindPupil;
		double yObj,xObj;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s="IncidentPointGlobal i yObj xObj SetRay FindPupil\n";
		}
		else if(b1 && b2 && b3 && b5){
			vector<double> v;
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			FindPupil=atoi(s5.c_str());
			v=IncidentPointGlobal(i,yObj,xObj,SetRay,FindPupil);
			sprintf(buf,"%.15g %.15g %.15g\n",v.x,v.y,v.z); s=buf;
		}
		return s;
	}
	if(s0=="Index" || s0=="index"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Index glassname\n";
		}
		else{
			if(args(com)==1){
				s=cOptics::scmd("Index "+s1+" "+this->color[1],val);
			}
		}
		return s;
	}
	if(s0=="Inflection" || s0=="inflection"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Inflection i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Inflection(i)); s=buf;
		}
		return s;
	}
	if(s0=="InflectionPoint" || s0=="inflectionpoint"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="InflectionPoint i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",InflectionPoint(i)); s=buf;
		}
		return s;
	}
	if(s0=="Koba" || s0=="koba"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Koba i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",koba(i)); s=buf;
		}
		return s;
	}
	if(s0=="LCPar" || s0=="lcpar"){
		int from,to;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="LCPar from to\n";
			s+="  Longitudinal chromatic aberration between no.'from' and no.'to' colors.\n";
		}
		else if( b1 && b2 ) {
			from=atoi(s1.c_str());
			to  =atoi(s2.c_str());
			sprintf(buf,"%.15g\n",LCPar(from,to)); s=buf;
		}
		return s;
	}
	if(s0=="LegendreR0" || s0=="legendrer0"){
		int i;
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="LegendreR0 i [new_value]\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			val=atof(s2.c_str());
			Set_LegendreR0(i,val);
			sprintf(buf,"%.15g\n",Get_LegendreR0(i)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Get_LegendreR0(i)); s=buf;
		}
		return s;
	}
	if(s0=="LegendreToFreeForm" || s0=="legendretofreeform"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="LegendreToFreeForm i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			LegendreToFreeForm(i);
		}
		return s;
	}
	if(s0=="LensDataZemax" || s0=="lensdatazemax"){
		s1=arg(com,1);
		if(s1=="?"){
			s="LensDataZemax (no argument)\n";
		}
		else{
			s=LensDataZemax();
		}
		return s;
	}
	if(s0=="LensCenterPowerWorn" || s0=="lenscenterpowerworn"){
		std::string SCA;
		double So,Co,Ao, WrapAngle,TiltAngle;
		s1=arg(com,1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="LensCenterPowerWorn SCA So Co Ao WrapAngle TiltAngle (SCA = 'S' or 'C' or 'A')\n";
		}
		else if(b2 && b3 && b4 && b5 && b6){
			SCA=s1;
			So=atof(s2.c_str());
			Co=atof(s3.c_str());
			Ao=atof(s4.c_str());
			WrapAngle=atof(s5.c_str());
			TiltAngle=atof(s6.c_str());
			sprintf(buf,"%.15g\n",LensCenterPowerWorn(SCA,So,Co,Ao,WrapAngle,TiltAngle)); s=buf;
			
		}
		return s;
	}
	if(s0=="LSA" || s0=="lsa"){
		double yPupilNormalized;
		int findpupil,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="LSA [yPupilNomalized=1 findpupil=1 [j=1]]\n";
		}
		else if(b1 && b2 && b3){
			yPupilNormalized=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			j=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",LSA(yPupilNormalized,findpupil,j)); s=buf;
		}
		else if(b1 && b2){
			yPupilNormalized=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",LSA(yPupilNormalized,findpupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",LSA()); s=buf;
		}
		return s;
	}
	if(s0=="LSA2nd" || s0=="lsa2nd"){
		s1=arg(com,1);
		if(s1=="?"){
			s="LSA2nd (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",LSA2nd()); s=buf;
		}
		return s;
	}
	if(s0=="LSA3rd" || s0=="lsa3rd"){
		s1=arg(com,1);
		if(s1=="?"){
			s="LSA3rd (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",LSA3rd()); s=buf;
		}
		return s;
	}
	if(s0=="LSA50" || s0=="lsa50"){
		s1=arg(com,1);
		if(s1=="?"){
			s="LSA50 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",LSA50()); s=buf;
		}
		return s;
	}
	if(s0=="LSA70" || s0=="lsa70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="LSA70 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",LSA70()); s=buf;
		}
		return s;
	}
	if(s0=="LSAp" || s0=="lsap"){
		double yObjNormalized;
		int findfield,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="LSAp [yObjNomalized=1 findfield=1 [j=1]]\n";
		}
		else if(b1 && b2 && b3){
			yObjNormalized=atof(s1.c_str());
			findfield=atoi(s2.c_str());
			j=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",LSAp(yObjNormalized,findfield,j)); s=buf;
		}
		else if(b1 && b2){
			yObjNormalized=atof(s1.c_str());
			findfield=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",LSAp(yObjNormalized,findfield)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",LSAp()); s=buf;
		}
		return s;
	}
	if(s0=="LSAp2nd" || s0=="lsap2nd"){
		double yObjNormalized;
		int findfield;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="LSAp2nd [yObjNomalized=1 findfield=1]\n";
		}
		else if(b1 && b2){
			yObjNormalized=atof(s1.c_str());
			findfield=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",LSAp2nd(yObjNormalized,findfield)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",LSAp2nd()); s=buf;
		}
		return s;
	}
	if(s0=="LSAp50" || s0=="lsap50"){
		s1=arg(com,1);
		if(s1=="?"){
			s="LSAp50 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",LSAp50()); s=buf;
		}
		return s;
	}
	if(s0=="LSAp70" || s0=="lsap70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="LSAp70 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",LSAp70()); s=buf;
		}
		return s;
	}
	if(s0=="LSAps" || s0=="lsaps"){
		int findfield,colors;
		double weight;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s+="LSAps findfield colors [weight]\n";
			s+="(findfield=0 is recommended for fast automatic design)\n";
		}
		else if(b1 && b2 && b3){
			findfield=atoi(s1.c_str());
			colors=atoi(s2.c_str());
			weight=atof(s3.c_str());
			s=LSAps(findfield,colors,weight);
		}
		else if(b1 && b2){
			findfield=atoi(s1.c_str());
			colors=atoi(s2.c_str());
			s=LSAps(findfield,colors);
		}
		return s;
	}
	if(s0=="LSAs" || s0=="lsas"){
		int findpupil,colors;
		double weight;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s+="LSAs findpupil colors [weight]\n";
			s+="(findpupil=0 is recommended for fast automatic design)\n";
		}
		else if(b1 && b2 && b3){
			findpupil=atoi(s1.c_str());
			colors=atoi(s2.c_str());
			weight=atof(s3.c_str());
			s=LSAs(findpupil,colors,weight);
		}
		else if(b1 && b2){
			findpupil=atoi(s1.c_str());
			colors=atoi(s2.c_str());
			s=LSAs(findpupil,colors);
		}
		return s;
	}
	if(s0=="M" || s0=="m"){
		int i1,i2;
		double s_; // s は std::string s で使用
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="M [[s=this->s] i1=1 i2=k]\n";
		}
		else if(b1 && b2 && b3){
			s_=atof(s1.c_str());
			i1=atoi(s2.c_str());
			i2=atoi(s3.c_str());
			if(val){
				sprintf(buf,"%.15g\n",M(s_,i1,i2,1)); s=buf;
			}
			else{
				sprintf(buf,"M(%.15g,%d,%d)=%g\n", s_,i1,i2,M(s_,i1,i2,1)); s=buf;
			}
		}
		else if(b1 && b2){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			if(val){
				sprintf(buf,"%.15g\n",M(i1,i2)); s=buf;
			}
			else{
				sprintf(buf,"M(%d,%d)=%g\n", i1,i2,M(i1,i2)); s=buf;
			}
		}
		else{
			if(val){
				sprintf(buf,"%.15g\n",M(1,k)); s=buf;
			}
			else{
				sprintf(buf,"M=%g\n", M(1,k)); s=buf;
			}
		}
		return s;
	}
	if(s0=="MakeDistChart" || s0=="makedistchart"){
		double yObjMax,xObjMax,defocus;
		int FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="MakeDistChart yObjMax xObjMax FindPupil defocus\n";
		}
		else if( b1 && b2 && b3 && b4 ) {
			yObjMax=atof(s1.c_str());
			xObjMax=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			defocus=atof(s4.c_str());
			s=MakeDistChart(yObjMax,xObjMax,FindPupil,defocus)+"\n";
		}
		return s;
	}
	if(s0=="MakeSpot" || s0=="makespot"){
		double yObj,xObj,defocus;
		int FindPupil,footprint,ColorStart,ColorEnd,IsAreaSource,IsLambert,OriginIsGravityCenter,Add;
		s1=arg(com,1);   b1=is_numeric(s1);
		s2=arg(com,2);   b2=is_numeric(s2);
		s3=arg(com,3);   b3=is_numeric(s3);
		s4=arg(com,4);   b4=is_numeric(s4);
		s5=arg(com,5);   b5=is_numeric(s5);
		s6=arg(com,6);   b6=is_numeric(s6);
		s7=arg(com,7);   b7=is_numeric(s7);
		s8=arg(com,8);   b8=is_numeric(s8);
		s9=arg(com,9);   b9=is_numeric(s9);
		s10=arg(com,10); b10=is_numeric(s10);
		s11=arg(com,11); b11=is_numeric(s11);
		if(s1=="?"){
			s+="MakeSpot yObj xObj FindPupil defocus footprint ColorStart ColorEnd\n";
			s+="         IsAreaSource IsLambert OriginIsGravityCenter Add\n";
			s+="MakeSpot yObj xObj FindPupil defocus ColorStart ColorEnd\n";
			s+="         OriginIsGravityCenter (footprint=0 IsAreaSource=0 IsLambert=0 Add=0)\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9 && b10 && b11){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			defocus=atof(s4.c_str());
			footprint=atoi(s5.c_str());
			ColorStart=atoi(s6.c_str());
			ColorEnd=atoi(s7.c_str());
			IsAreaSource=atoi(s8.c_str());
			IsLambert=atoi(s9.c_str());
			OriginIsGravityCenter=atoi(s10.c_str());
			Add=atoi(s11.c_str());
			sprintf(buf,"%d\n",MakeSpot(yObj,xObj,FindPupil,defocus,footprint,ColorStart,ColorEnd,
				                           IsAreaSource,IsLambert,OriginIsGravityCenter,Add));
			s=buf;
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			defocus=atof(s4.c_str());
			ColorStart=atoi(s5.c_str());
			ColorEnd=atoi(s6.c_str());
			OriginIsGravityCenter=atoi(s7.c_str());
			sprintf(buf,"%d\n",MakeSpot(yObj,xObj,FindPupil,defocus,0,ColorStart,ColorEnd,
				                           0,0,OriginIsGravityCenter,0));
			s=buf;
		}
		return s;
	}
	if(s0=="MartinCyl" || s0=="martincyl"){
		double Diopter,th_deg,N;
		s1=arg(com,1);   b1=is_numeric(s1);
		s2=arg(com,2);   b2=is_numeric(s2);
		s3=arg(com,3);   b3=is_numeric(s3);
		if(s1=="?"){
			s+="MartinCyl Diopter th_deg N\n";
		}
		else if(b1 && b2 && b3){
			Diopter=atof(s1.c_str());
			th_deg=atof(s2.c_str());
			N=atof(s3.c_str());
			sprintf(buf,"%.15g\n",MartinCyl(Diopter,th_deg,N));
			s=buf;
		}
		return s;
	}
	if(s0=="MartinEqGeneral" || s0=="martineqgeneral"){
		double Sph0,Cyl0,Axis0_deg,N,rox_deg,roy_deg;
		int order,inverse;
		double Sph,Cyl,Axis_deg;
		s1=arg(com,1);   b1=is_numeric(s1);
		s2=arg(com,2);   b2=is_numeric(s2);
		s3=arg(com,3);   b3=is_numeric(s3);
		s4=arg(com,4);   b4=is_numeric(s4);
		s5=arg(com,5);   b5=is_numeric(s5);
		s6=arg(com,6);   b6=is_numeric(s6);
		s7=arg(com,7);   b7=is_numeric(s7);
		s8=arg(com,8);   b8=is_numeric(s8);
		if(s1=="?"){
			s+="MartinEqGeneral Sph0 Cyl0 Axis0_deg N rox_deg roy_deg order inverse\n";
			s+="return val = Sph Cyl Axis_deg\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8){
			Sph0=atof(s1.c_str());
			Cyl0=atof(s2.c_str());
			Axis0_deg=atof(s3.c_str());
			N=atof(s4.c_str());
			rox_deg=atof(s5.c_str());
			roy_deg=atof(s6.c_str());
			order=atoi(s7.c_str());
			inverse=atoi(s8.c_str());
			MartinEqGeneral(Sph,Cyl,Axis_deg,Sph0,Cyl0,Axis0_deg,N,rox_deg,roy_deg,order,inverse);
			sprintf(buf,"%.15g %.15g %.15g\n", Sph,Cyl,Axis_deg);
			s=buf;
		}
		return s;
	}
	if(s0=="MartinSph" || s0=="martinsph"){
		double Diopter,th_deg,N;
		s1=arg(com,1);   b1=is_numeric(s1);
		s2=arg(com,2);   b2=is_numeric(s2);
		s3=arg(com,3);   b3=is_numeric(s3);
		if(s1=="?"){
			s+="MartinSph Diopter th_deg N\n";
		}
		else if(b1 && b2 && b3){
			Diopter=atof(s1.c_str());
			th_deg=atof(s2.c_str());
			N=atof(s3.c_str());
			sprintf(buf,"%.15g\n",MartinSph(Diopter,th_deg,N));
			s=buf;
		}
		return s;
	}
	if(s0=="MinMTF" || s0=="minmtf"){
		double nu;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="MinMTF nu\n";
		}
		else if( b1 ) {
			nu=atof(s1.c_str());
			sprintf(buf,"%.15g\n",MinMTF(nu)); s=buf;
		}
		return s;
	}
	if(s0=="ModelEyeFundusR" || s0=="modeleyefundusr"){
		double diopter,VD;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ModelEyeFundusR diopter VD\n";
		}
		else if(b1 && b2){
			diopter=atof(s1.c_str());
			VD=atof(s2.c_str());
			sprintf(buf,"%.15g\n",ModelEyeFundusR(diopter,VD)); s=buf;
		}
		return s;
	}
	if(s0=="ModelEyeS1" || s0=="modeleyes1"){
		double diopter,VD;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ModelEyeS1 diopter VD\n";
		}
		else if(b1 && b2){
			diopter=atof(s1.c_str());
			VD=atof(s2.c_str());
			sprintf(buf,"%.15g\n",ModelEyeS1(diopter,VD)); s=buf;
		}
		return s;
	}
	if(s0=="Mpupil" || s0=="mpupil"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="Mpupil [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			if(val){
				sprintf(buf,"%.15g\n",Mpupil(i1,i2)); s=buf;
			}
			else{
				sprintf(buf,"Mpupil(%d,%d)=%g\n", i1,i2,Mpupil(i1,i2)); s=buf;
			}
		}
		else{
			if(val){
				sprintf(buf,"%.15g\n",Mpupil(1,k)); s=buf;
			}
			else{
				sprintf(buf,"Mpupil=%g\n", Mpupil(1,k)); s=buf;
			}
		}
		return s;
	}
	if(s0=="Mpx" || s0=="mpx"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="Mpx [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",Mpx(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Mpx(1,k)); s=buf;
		}
		return s;
	}
	if(s0=="Mpy" || s0=="mpy"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="Mpy [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",Mpy(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Mpy(1,k)); s=buf;
		}
		return s;
	}
	if(s0=="Mstop" || s0=="mstop"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Mstop (no argements)\n";
		}
		else {
			sprintf(buf,"%.15g\n",Mstop()); s=buf;
		}
		return s;
	}
	if(s0=="Mstop1" || s0=="mstop1"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Mstop1 (no argements)\n";
		}
		else {
			sprintf(buf,"%.15g\n",Mstop1()); s=buf;
		}
		return s;
	}
	if(s0=="MTFm" || s0=="mtfm"){
		double nu,yObj,xObj,defocus;
		int Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="MTFm nu yObj xObj [defocus=OptimizedDefocus() Geometrical=0 AberrationFree=0 \n";
			s+="                   ColorStart=1 ColorEnd=cn FindPupil=1] \n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			defocus=atof(s4.c_str());
			Geometrical=atoi(s5.c_str());
			AberrationFree=atoi(s6.c_str());
			ColorStart=atoi(s7.c_str());
			ColorEnd=atoi(s8.c_str());
			FindPupil=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",
				    MTFm(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil));
			s=buf;
		}
		else if(b1 && b2 && b3){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			sprintf(buf,"%.15g\n",MTFm(nu,yObj,xObj));
			s=buf;
		}
		return s;
	}
	if(s0=="MTFs" || s0=="mtfs"){
		double nu,yObj,xObj,defocus;
		int Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="MTFs nu yObj xObj [defocus=OptimizedDefocus() Geometrical=0 AberrationFree=0 \n";
			s+="                   ColorStart=1 ColorEnd=cn FindPupil=1] \n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			defocus=atof(s4.c_str());
			Geometrical=atoi(s5.c_str());
			AberrationFree=atoi(s6.c_str());
			ColorStart=atoi(s7.c_str());
			ColorEnd=atoi(s8.c_str());
			FindPupil=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",
				    MTFs(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil));
			s=buf;
		}
		else if(b1 && b2 && b3){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			sprintf(buf,"%.15g\n",MTFs(nu,yObj,xObj));
			s=buf;
		}
		return s;
	}
	if(s0=="MTFsmave" || s0=="mtfsmave"){
		double nu,yObj,xObj,defocus;
		int Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="MTFsmave nu yObj xObj [defocus=OptimizedDefocus() Geometrical=0 AberrationFree=0 \n";
			s+="                       ColorStart=1 ColorEnd=cn FindPupil=1] \n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			defocus=atof(s4.c_str());
			Geometrical=atoi(s5.c_str());
			AberrationFree=atoi(s6.c_str());
			ColorStart=atoi(s7.c_str());
			ColorEnd=atoi(s8.c_str());
			FindPupil=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",
				    MTFsmave(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil));
			s=buf;
		}
		else if(b1 && b2 && b3){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			sprintf(buf,"%.15g\n",MTFsmave(nu,yObj,xObj));
			s=buf;
		}
		return s;
	}
	if(s0=="MTFx" || s0=="mtfx"){
		double nu,yObj,xObj,defocus;
		int Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="MTFx nu yObj xObj [defocus=OptimizedDefocus() Geometrical=0 AberrationFree=0 \n";
			s+="                   ColorStart=1 ColorEnd=cn FindPupil=1] \n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			defocus=atof(s4.c_str());
			Geometrical=atoi(s5.c_str());
			AberrationFree=atoi(s6.c_str());
			ColorStart=atoi(s7.c_str());
			ColorEnd=atoi(s8.c_str());
			FindPupil=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",
				    MTFx(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil));
			s=buf;
		}
		else if(b1 && b2 && b3){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			sprintf(buf,"%.15g\n",MTFx(nu,yObj,xObj));
			s=buf;
		}
		return s;
	}
	if(s0=="MTFxyave" || s0=="mtfxyave"){
		double nu,yObj,xObj,defocus;
		int Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="MTFxyave nu yObj xObj [defocus=OptimizedDefocus() Geometrical=0 AberrationFree=0 \n";
			s+="                       ColorStart=1 ColorEnd=cn FindPupil=1] \n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			defocus=atof(s4.c_str());
			Geometrical=atoi(s5.c_str());
			AberrationFree=atoi(s6.c_str());
			ColorStart=atoi(s7.c_str());
			ColorEnd=atoi(s8.c_str());
			FindPupil=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",
				    MTFxyave(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil));
			s=buf;
		}
		else if(b1 && b2 && b3){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			sprintf(buf,"%.15g\n",MTFxyave(nu,yObj,xObj));
			s=buf;
		}
		return s;
	}
	if(s0=="MTFy" || s0=="mtfy"){
		double nu,yObj,xObj,defocus;
		int Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s+="MTFy nu yObj xObj [defocus=OptimizedDefocus() Geometrical=0 AberrationFree=0 \n";
			s+="                   ColorStart=1 ColorEnd=cn FindPupil=1] \n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			defocus=atof(s4.c_str());
			Geometrical=atoi(s5.c_str());
			AberrationFree=atoi(s6.c_str());
			ColorStart=atoi(s7.c_str());
			ColorEnd=atoi(s8.c_str());
			FindPupil=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",
				    MTFy(nu,yObj,xObj,defocus,Geometrical,AberrationFree,ColorStart,ColorEnd,FindPupil));
			s=buf;
		}
		else if(b1 && b2 && b3){
			nu=atof(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			sprintf(buf,"%.15g\n",MTFy(nu,yObj,xObj));
			s=buf;
		}
		return s;
	}	
	if(s0=="Mx" || s0=="mx"){
		double yObj,xObj,yPupil,xPupil;
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s= "Mx yObj xObj yPupil xPupil i1 i2\n";
			s+="Mx [i1=1 i2=k]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			i1=atoi(s5.c_str());
			i2=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",Mx(yObj,xObj,yPupil,xPupil,i1,i2,1)); s=buf;
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",Mx(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Mx(1,k)); s=buf;
		}
		return s;
	}
	if(s0=="My" || s0=="my"){
		double yObj,xObj,yPupil,xPupil;
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s= "My yObj xObj yPupil xPupil i1 i2\n";
			s+="My [i1=1 i2=k]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			i1=atoi(s5.c_str());
			i2=atoi(s6.c_str());
			sprintf(buf,"%.15g\n",My(yObj,xObj,yPupil,xPupil,i1,i2,1)); s=buf;
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",My(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",My(1,k)); s=buf;
		}
		return s;
	}
	if(s0=="N" || s0=="n"){
		int i,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="N i [j=1]]\n";
		}
		else if(b1 && b2) {
			i=atoi(s1.c_str());
			j=atoi(s2.c_str());
			if(0<=i && i<=this->k && 1<=j && j<=this->cn){
				sprintf(buf,"%.15g\n",N(i,j)); s=buf;
			}
		}
		else if(b1) {
			i=atoi(s1.c_str());
			if(0<=i && i<=this->k){
				sprintf(buf,"%.15g\n",N(i,1)); s=buf;
			}
		}
		return s;
	}
	if(s0=="NA" || s0=="na"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="NA [yObj=0 xObj=0 [findpupil=1]]\n";
		}
		else if( b1 && b2 && b3 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",NA(yObj,xObj,findpupil)); s=buf;
		}
		else if( b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",NA(yObj,xObj,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",NA(0,0,1)); s=buf;
		}
		return s;
	}
	if(s0=="NAObj" || s0=="naobj"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="NAObj [yObj=0 xObj=0 [findpupil=1]]\n";
		}
		else if( b1 && b2 && b3 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",NAObj(yObj,xObj,findpupil)); s=buf;
		}
		else if( b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",NAObj(yObj,xObj,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",NAObj(0,0,1)); s=buf;
		}
		return s;
	}
	if(s0=="NAObjXYAve" || s0=="naobjxyave"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="NAObjXYAve [yObj=0 xObj=0 [findpupil=1]\n";
		}
		else if( b1 && b2 && b3 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",NAObjXYAve(yObj,xObj,findpupil)); s=buf;
		}
		else if( b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",NAObjXYAve(yObj,xObj,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",NAObjXYAve(0,0,1)); s=buf;
		}
		return s;
	}
	if(s0=="NAParaxial" || s0=="naparaxial"){
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="NAParaxial [findpupil=1]\n";
		}
		else if( b1 ) {
			findpupil=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",NAParaxial(findpupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",NAParaxial(1)); s=buf;
		}
		return s;
	}
	if(s0=="NAXYAve" || s0=="naxyave"){
		double yObj,xObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="NAXYAve [yObj=0 xObj=0 [findpupil=1]]\n";
		}
		else if( b1 && b2 && b3 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",NAXYAve(yObj,xObj,findpupil)); s=buf;
		}
		else if( b1 && b2) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",NAXYAve(yObj,xObj,1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",NAXYAve(0,0,1)); s=buf;
		}
		return s;
	}
	if(s0=="NewtonToR" || s0=="newtontor"){
		double r0,eaphi,FRW_nm,rings;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="NewtonToR r0 eaphi FRW_nm rings\n";
		}
		else if( b1 && b2 && b3 && b4 ){
			r0=atof(s1.c_str());
			eaphi=atof(s2.c_str());
			FRW_nm=atof(s3.c_str());
			rings=atof(s4.c_str());
			sprintf(buf,"r0=%g eaphi=%g FRW=%gnm rings=%g\n", r0,eaphi,FRW_nm,rings);
			s+=buf;
			sprintf(buf,"r=%g %g\n",NewtonToR(r0,eaphi,FRW_nm,rings),NewtonToR(r0,eaphi,FRW_nm,-rings));
			s+=buf;
		}
		return s;
	}
	if(s0=="NNuActual" || s0=="nnuactual"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="NNuActual i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n", NNuActual(i)); s+=buf;
		}
		return s;
	}
	if(s0=="nodal"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="nodal\nnodal i1 i2\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"Δnodal'(%s,%d-%d)=%g\n", color[1].c_str(),i1,i2,nodal1(i1,i2,1)); s+=buf;
			sprintf(buf,"Δnodal (%s,%d-%d)=%g\n", color[1].c_str(),i1,i2,nodal(i1,i2,1) ); s+=buf;
		}
		else{
			sprintf(buf,"Δnodal'(%s)=%g\n", color[1].c_str(),nodal1(1)); s+=buf;
			sprintf(buf,"Δnodal (%s)=%g\n", color[1].c_str(),nodal(1) ); s+=buf;
		}
		return s;
	}
	if(s0=="NotThruSurfaces" || s0=="notthrusurfaces"){
		double yObject,xObject;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="NotThruSurfaces yObject xObject\n";
		}
		else if( b1 && b2 ){
			yObject=atof(s1.c_str());
			xObject=atof(s2.c_str());
			s=NotThruSurfaces(yObject,xObject,0);
		}
		return s;
	}
	if(s0=="OffAxialConicAB" || s0=="offaxialconicab"){
		double a,b,t_deg,LongR,ShortR;
		s1=arg(com,1);  b1=is_numeric(s1); 
		s2=arg(com,2);  b2=is_numeric(s2);
		s3=arg(com,3);  b3=is_numeric(s3);
		if(s1=="?"){
			s="OffAxialConicAB t_deg LongR ShortR\n";
		}
		else if(b1 && b2 && b3){
			t_deg=atof(s1.c_str());
			LongR=atof(s2.c_str());
			ShortR=atof(s3.c_str());
			OffAxialConicAB(a,b,t_deg,LongR,ShortR);
			if(a==0){
				sprintf(buf,"no solution\n");
			}
			else{
				sprintf(buf,"a=%.15g b=%.15g\n",a,b);
			}
			s=buf;
		}
		return s;
	}
	if(s0=="OffAxialMirrorComa" || s0=="offaxialmirrorcoma"){
		int i;
		s1=arg(com,1);  b1=is_numeric(s1); 
		if(s1=="?"){
			s="OffAxialMirrorComa i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",OffAxialMirrorComa(i));
			s=buf;
		}
		return s;
	}
	if(s0=="OPD" || s0=="opd"){
		double yObj,xObj,yPupil,xPupil,yPupil_principal,xPupil_principal,defocus;
		int j0,j, FindPupil;
		std::string SetRay;
		s1 =arg(com,1);  b1= is_numeric(s1); 
		s2 =arg(com,2);  b2= is_numeric(s2);
		s3 =arg(com,3);  b3= is_numeric(s3);
		s4 =arg(com,4);  b4= is_numeric(s4);
		s5 =arg(com,5);  b5= is_numeric(s5); 
		s6 =arg(com,6);  b6= is_numeric(s6);
		s7 =arg(com,7);  b7= is_numeric(s7);
		s8 =arg(com,8);  b8= is_numeric(s8);
		s9 =arg(com,9);  b9= is_numeric(s9);
		s10=arg(com,10); b10=is_numeric(s10);
		if(s1=="?"){
			s="OPD yObj xObj SetRay yPupil xPupil yPupil_principal xPupil_principal defocus j0 j\n";
			s+="OPD yObj xObj FindPupil ypNormalized xpNormalized defocus j0 j\n";
		}
		else if(b1 && b2 && b4 && b5 && b6 && b7 && b8 && b9 && b10){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3.c_str();
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			yPupil_principal=atof(s6.c_str());
			xPupil_principal=atof(s7.c_str());
			defocus=atof(s8.c_str());
			j0=atoi(s9.c_str());
			j=atoi(s10.c_str());
			sprintf(buf,"%.15g\n",OPD(yObj,xObj,SetRay,yPupil,xPupil,yPupil_principal,xPupil_principal,defocus,j0,j));
			s=buf;
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			defocus=atof(s6.c_str());
			j0=atoi(s7.c_str());
			j=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",OPD(yObj,xObj,FindPupil,yPupil,xPupil,defocus,j0,j));
			s=buf;
		}
		return s;
	}
	if(s0=="OPD2" || s0=="opd2"){
		double yObj,xObj,yPupil,xPupil,yPupil_principal,xPupil_principal,defocus,wl_nm0,wl_nm;
		std::string SetRay;
		s1 =arg(com,1);  b1= is_numeric(s1); 
		s2 =arg(com,2);  b2= is_numeric(s2);
		s3 =arg(com,3);
		s4 =arg(com,4);  b4= is_numeric(s4);
		s5 =arg(com,5);  b5= is_numeric(s5); 
		s6 =arg(com,6);  b6= is_numeric(s6);
		s7 =arg(com,7);  b7= is_numeric(s7);
		s8 =arg(com,8);  b8= is_numeric(s8);
		s9 =arg(com,9);  b9= is_numeric(s9);
		s10=arg(com,10); b10=is_numeric(s10);
		if(s1=="?"){
			s="OPD yObj xObj SetRay yPupil xPupil yPupil_principal xPupil_principal defocus wl_nm0 wl_nm\n";
		}
		else if( b1 && b2 && b4 && b5 && b6 && b7 && b8 && b9 && b10 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3.c_str();
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			yPupil_principal=atof(s6.c_str());
			xPupil_principal=atof(s7.c_str());
			defocus=atof(s8.c_str());
			wl_nm0=atof(s9.c_str());
			wl_nm=atof(s10.c_str());
			sprintf(buf,"%.15g\n",
			        OPD2(yObj,xObj,SetRay,yPupil,xPupil,yPupil_principal,xPupil_principal,defocus,wl_nm0,wl_nm));
			s=buf;
		}
		return s;
	}
	if(s0=="OPDPV" || s0=="opdpv"){
		double yObj,xObj,defocus;
		int j,FindPupil,InLambda,AdjustSph,OptimizeDefocusOnAxis;
		s1=arg(com,1); b1=is_numeric(s1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5); 
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="OPDPV [yObj=0 xObj=0 [defocus=0 j=1 FindPupil=1 InLambda=1 AdjustSph=0 OptimizeDefousOnAxis=1]]\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			defocus=atof(s3.c_str());
			j=atoi(s4.c_str());
			FindPupil=atoi(s5.c_str());
			InLambda=atoi(s6.c_str());
			AdjustSph=atoi(s7.c_str());
			OptimizeDefocusOnAxis=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",OPDPV(yObj,xObj,defocus,j,FindPupil,InLambda,AdjustSph,OptimizeDefocusOnAxis,""));
			s=buf;
		}
		else if( b1 && b2 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",OPDPV(yObj,xObj,0,1,1,1,0,1,""));
			s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",OPDPV(0,0,0,1,1,1,0,1,""));
			s=buf;
		}
		return s;
	}
	if(s0=="OPDRMS" || s0=="opdrms"){
		double yObj,xObj,defocus;
		int j,FindPupil,InLambda,AdjustSph,OptimizeDefocusOnAxis;
		s1=arg(com,1); b1=is_numeric(s1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5); 
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="OPDRMS [yObj=0 xObj=0 [defocus=0 j=1 FindPupil=1 InLambda=1 AdjustSph=0 OptimizeDefousOnAxis=1]]\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			defocus=atof(s3.c_str());
			j=atoi(s4.c_str());
			FindPupil=atoi(s5.c_str());
			InLambda=atoi(s6.c_str());
			AdjustSph=atoi(s7.c_str());
			OptimizeDefocusOnAxis=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",OPDRMS(yObj,xObj,defocus,j,FindPupil,InLambda,AdjustSph,OptimizeDefocusOnAxis,""));
			s=buf;
		}
		else if( b1 && b2 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",OPDRMS(yObj,xObj,0,1,1,1,0,1,""));
			s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",OPDRMS(0,0,0,1,1,1,0,1,""));
			s=buf;
		}
		return s;
	}
	if(s0=="OPDRMS0" || s0=="opdrms0"){
		double yObj,xObj;
		int FindPupil;
		s1=arg(com,1); b1=is_numeric(s1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="OPDRMS0 [[yObj=0 xObj=0] FindPupil=1] (AdjustSph=OptimizeDefocus=0)\n";
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",OPDRMS0(yObj,xObj,FindPupil));
			s=buf;
		}
		else if(b1){
			FindPupil=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",OPDRMS0(0,0,FindPupil));
			s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",OPDRMS0(0,0,1));
			s=buf;
		}
		return s;
	}
	if(s0=="OPDs" || s0=="opds"){
		double yObj,xObj,weight;
		int findpupil,colors,n,IgnoreTC;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="OPDs yObj xObj findpupil colors [IgnoreTC=0 [n=0 [weight=1]]]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			n=atoi(s6.c_str());
			weight=atof(s7.c_str());
			s=OPDs(yObj,xObj,findpupil,colors,IgnoreTC,n,weight);
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			n=atoi(s6.c_str());
			s=OPDs(yObj,xObj,findpupil,colors,IgnoreTC,n);
		}
		else if(b1 && b2 && b3 && b4 && b5){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			s=OPDs(yObj,xObj,findpupil,colors,IgnoreTC);
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			s=OPDs(yObj,xObj,findpupil,colors);
		}
		return s;
	}
	if(s0=="OPL" || s0=="opl"){
		double yObj,xObj,yPupil,xPupil,defocus;
		int j,i1,i2;
		std::string SetRay;
		s1 =arg(com,1);  b1= is_numeric(s1);
		s2 =arg(com,2);  b2= is_numeric(s2);
		s3 =arg(com,3);
		s4 =arg(com,4);  b4= is_numeric(s4);
		s5 =arg(com,5);  b5= is_numeric(s5);
		s6 =arg(com,6);  b6= is_numeric(s6);
		s7 =arg(com,7);  b7= is_numeric(s7);
		s8 =arg(com,8);  b8= is_numeric(s8);
		s9 =arg(com,9);  b9= is_numeric(s9);
		if(s1=="?"){
			s="OPL yObj xObj SetRay yPupil xPupil defocus j [i1=0 i2=0]\n";
		}
		else if(b1 && b2 && b4 && b5 && b6 && b7 && b8 && b9){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3.c_str();
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			defocus=atof(s6.c_str());
			j=atoi(s7.c_str());
			i1=atoi(s8.c_str());
			i2=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",OPL(yObj,xObj,SetRay,yPupil,xPupil,defocus,j,i1,i2));
			s=buf;
		}
		else if(b1 && b2 && b4 && b5 && b6 && b7){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3.c_str();
			yPupil=atof(s4.c_str());
			xPupil=atof(s5.c_str());
			defocus=atof(s6.c_str());
			j=atoi(s7.c_str());
			sprintf(buf,"%.15g\n",OPL(yObj,xObj,SetRay,yPupil,xPupil,defocus,j));
			s=buf;
		}
		return s;
	}
	if(s0=="OptimizedDefocus" || s0=="optimizeddefocus"){
		double yObj,xObj;
		int FindPupil, j;
		s1=arg(com,1); b1=is_numeric(s1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s+="OptimizedDefocus [yObj=0 xObj=0 [FindPupil=1 [j=1]]]\n";
			s+="OptimizedDefocus FindPupil (yObj=0 xObj=0 j=1)\n";
		}
		else if( b1 && b2 && b3 && b4 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			j=atoi(s4.c_str());
			sprintf(buf, "%.15g\n", OptimizedDefocus(yObj,xObj,FindPupil,j));
			s=buf;
		}
		else if( b1 && b2 && b3 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			sprintf(buf, "%.15g\n", OptimizedDefocus(yObj,xObj,FindPupil,1));
			s=buf;
		}
		else if( b1 && b2 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf, "%.15g\n", OptimizedDefocus(yObj,xObj,1,1));
			s=buf;
		}
		else if(b1){
			FindPupil=atoi(s1.c_str());
			sprintf(buf, "%.15g\n", OptimizedDefocus(FindPupil));
			s=buf;
		}
		else{
			sprintf(buf, "%.15g\n", OptimizedDefocus(1));
			s=buf;
		}
		return s;
	}
	if(s0=="OptimizeS1fix" || s0=="optimizes1fix"){
		double yObj,xObj;
		int FindPupil, j;
		s1=arg(com,1); b1=is_numeric(s1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s+="OptimizeS1fix [yObj=0 xObj=0 [FindPupil [j=1]]]\n";
			s+="OtiimizeS1fix FindPupil (yObj=0 xObj=0 j=1)\n";
		}
		else if( b1 && b2 && b3 && b4 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			j=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",OptimizeS1fix(yObj,xObj,FindPupil,j)); s=buf;
		}
		else if( b1 && b2 && b3 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			FindPupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",OptimizeS1fix(yObj,xObj,FindPupil,1)); s=buf;
		}
		else if( b1 && b2 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",OptimizeS1fix(yObj,xObj,1,1)); s=buf;
		}
		else if(b1){
			FindPupil=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",OptimizeS1fix(FindPupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",OptimizeS1fix(1)); s=buf;
		}
		return s;
	}
	if(s0=="OSC" || s0=="osc"){
		s1=arg(com,1);
		if(s1=="?"){
			s="OSC (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",OSC()); s=buf;
		}
		return s;
	}
	if(s0=="OSC70" || s0=="osc70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="OSC70 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",OSC70()); s=buf;
		}
		return s;
	}
	if(s0=="OSC50" || s0=="osc50"){
		s1=arg(com,1);
		if(s1=="?"){
			s="OSC50 (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",OSC50()); s=buf;
		}
		return s;
	}
	if(s0=="OSCp" || s0=="oscp"){
		double yObjNormalized;
		int findfield;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="OSCp [yObjNomalized=1 findfield=1]\n";
		}
		else if(b1 && b2){
			yObjNormalized=atof(s1.c_str());
			findfield=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",OSCp(yObjNormalized,findfield)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",OSCp()); s=buf;
		}
		return s;
	}
	if(s0=="ParaxialValues" || s0=="paraxialvalues"){
		s1=arg(com,1);
		if(s1=="?"){
			s="ParaxialValues (no arguments)\n";
		}
		else{
			s=ParaxialValues();
		}
		return s;
	}
	if(s0=="PerturbDxDy" || s0=="perturbdxdy"){
		int i;
		double dr;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="PerturbDxDy i dr\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			dr=atof(s2.c_str());
			sprintf(buf,"%.15g\n",PerturbDxDy(i,dr)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",PerturbDxDy(i)); s=buf;
		}
		return s;
	}
	if(s0=="PerturbXYObjectMaxDxDy" || s0=="perturbxyobjectmax"){
		double dr;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="PerturbXYObjectMax dr\n";
		}
		else if(b1){
			dr=atof(s1.c_str());
			sprintf(buf,"%.15g\n",PerturbXYObjectMax(dr)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",PerturbXYObjectMax()); s=buf;
		}

		return s;
	}
	if(s0=="power"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="power [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",power(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",power(1,k)); s=buf;
		}
		return s;
	}
	if(s0=="PreformKoba" || s0=="preformkoba"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="PreformKoba i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",PreformKoba(i)); s=buf;
		}
		return s;
	}
	if(s0=="PSeriesToDcon" || s0=="pseriestodcon"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="PSeriesToDcon i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			PSeriesToDcon(i);
		}
		return s;
	}
	if(s0=="PSeriesToFreeForm" || s0=="pseriestofreeform"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="PSeriesToFreeForm i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			PSeriesToFreeForm(i);
		}
		return s;
	}
	if(s0=="PupilToPupil" || s0=="pupiltopupil"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="PupilToPupil [i1=1 i2=k]\n";
		}
		else if(b1 && b2) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",PupilToPupil(i1,i2)); s+=buf;
		}
		else{
			sprintf(buf,"%.15g\n",PupilToPupil()); s+=buf;
		}
		return s;
	}
	if(s0=="qBend" || s0=="qbend"){
		int i;
		double q;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="qBend i q\n";
		}
		else if( b1 && b2 ) {
			i=atoi(s1.c_str());
			q=atof(s2.c_str());
			qBend(i,q);
		}
		return s;
	}
	if(s0=="qValue" || s0=="qvalue"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="qValue i\n";
		}
		else if(b1) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",qValue(i)); s=buf;
		}
		return s;
	}
	if(s0=="RayPosX" || s0=="rayposx"){
		int i,j;
		double yObj,xObj,yPupil,xPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="RayPosX i [yObj=yObjectMax xObj=xObjectMax] SetRay yPupil xPupil j\n";
		}
		else if( b1 && b2 && b3 && b5 && b6 && b7 ){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			sprintf(buf,"%.15g\n",RayPosX(i,yObj,xObj,SetRay,yPupil,xPupil,j)); s=buf;
		}
		else if( b1 && b3 && b4 && b5 ){
			i=atoi(s1.c_str());
			SetRay=s2;
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			j=atoi(s5.c_str());
			sprintf(buf,"%.15g\n",RayPosX(i,yObjectMax,xObjectMax,SetRay,yPupil,xPupil,j)); s=buf;
		}
		return s;
	}
	if(s0=="RayPosY" || s0=="rayposy"){
		int i,j;
		double yObj,xObj,yPupil,xPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="RayPosY i [yObj=yObjectMax xObj=xObjectMax] SetRay yPupil xPupil j\n";
		}
		else if( b1 && b2 && b3 && b5 && b6 && b7 ){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			sprintf(buf,"%.15g\n",RayPosY(i,yObj,xObj,SetRay,yPupil,xPupil,j)); s=buf;
		}
		else if( b1 && b3 && b4 && b5 ){
			i=atoi(s1.c_str());
			SetRay=s2;
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			j=atoi(s5.c_str());
			sprintf(buf,"%.15g\n",RayPosY(i,yObjectMax,xObjectMax,SetRay,yPupil,xPupil,j)); s=buf;
		}
		return s;
	}
	if(s0=="RayPosZ" || s0=="rayposz"){
		int i,j;
		double yObj,xObj,yPupil,xPupil;
		std::string SetRay;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="RayPosZ i [yObj=yObjectMax xObj=xObjectMax] SetRay yPupil xPupil j\n";
		}
		else if( b1 && b2 && b3 && b5 && b6 && b7 ){
			i=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			SetRay=s4;
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			j=atoi(s7.c_str());
			sprintf(buf,"%.15g\n",RayPosZ(i,yObj,xObj,SetRay,yPupil,xPupil,j)); s=buf;
		}
		else if( b1 && b3 && b4 && b5 ){
			i=atoi(s1.c_str());
			SetRay=s2;
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			j=atoi(s5.c_str());
			sprintf(buf,"%.15g\n",RayPosZ(i,yObjectMax,xObjectMax,SetRay,yPupil,xPupil,j)); s=buf;
		}
		return s;
	}
	if(s0=="ReduceAsphTerms" || s0=="reduceasphterms"){
		int i,max_order;
		double phi,h_step;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="ReduceAsphTerms i max_order [phi=0 [h_step=1]]\n";
		}
		else if(b1 && b2 && b3 && b4){
			i=atoi(s1.c_str());
			max_order=atoi(s2.c_str());
			phi=atof(s3.c_str());
			h_step=atof(s4.c_str());
			ReduceAsphTerms(i,max_order,phi,h_step);
		}
		else if(b1 && b2 && b3){
			i=atoi(s1.c_str());
			max_order=atoi(s2.c_str());
			phi=atof(s3.c_str());
			ReduceAsphTerms(i,max_order,phi);
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			max_order=atoi(s2.c_str());
			ReduceAsphTerms(i,max_order);
		}
		return s;
	}
	if(s0=="RingsSpecRefSurf" || s0=="ringsspecrefsurf"){
		double th_deg,VirtualAS;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="RingsSpecRefSurf th_deg VirtualAS\n";
			s+="  It is convenient that multiplying VirtualAS by\n";
			s+="  (surface inspection diameter/light diameter)^2 / medium index\n";
		}
		else if( b1 && b2 ){
			th_deg=atof(s1.c_str());
			VirtualAS=atof(s2.c_str());
			s=RingsSpecRefSurf(th_deg,VirtualAS);
		}
		return s;
	}
	if(s0=="RmsPhi" || s0=="rmsphi" || s0=="rmsφ"){
		double yObj,xObj,defocus;
		int ColorStart,ColorEnd,FindPupil,OptimizeDefocus;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s+="RmsPhi yObj xObj [defocus=0 ColorStart=1 ColorEnd=cn FindPupil=1 OptimizeDefocus=1]\n";
			s+="Rmsφ  yObj xObj [defocus=0 ColorStart=1 ColorEnd=cn FindPupil=1 OptimizeDefocus=1]\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6 && b7 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			defocus=atof(s3.c_str());
			ColorStart=atoi(s4.c_str());
			ColorEnd=atoi(s5.c_str());
			FindPupil=atoi(s6.c_str());
			OptimizeDefocus=atoi(s7.c_str());
			sprintf(buf,"%.15g\n",
			        rmsphi(yObj,xObj,defocus,ColorStart,ColorEnd,FindPupil,OptimizeDefocus));
			s=buf;
		}
		else if( b1 && b2 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",rmsphi(yObj,xObj,0,1,this->cn,1,1));
			s=buf;
		}
		return s;
	}
	if(s0=="RmsPhiOff" || s0=="rmsphioff" || s0=="rmsφoff"){
		int OptimizeDefocus;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s+="RmsPhiOff [OptimizeDefocus=1]\n";
			s+="RmsφOff [OptimizeDefocus=1]\n";
		}
		else if( b1 ){
			OptimizeDefocus=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",rmsphiOff(OptimizeDefocus)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",rmsphiOff(1)); s=buf;
		}
		return s;
	}
	if(s0=="RmsPhiOn" || s0=="rmsphion" || s0=="rmsφon"){
		int OptimizeDefocus;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s+="RmsPhiOn [OptimizeDefocus=1]\n";
			s+="RmsφOn [OptimizeDefocus=1]\n";
		}
		else if( b1 ){
			OptimizeDefocus=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",rmsphiOn(OptimizeDefocus)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",rmsphiOn(1)); s=buf;
		}
		return s;
	}
	if(s0=="RotateBlock" || s0=="rotateblock"){
		int i1,i2;
		double Sx,Sy,Sz,Rx,Ry,Rz,th;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s="RotateBlock i1 i2 Sx Sy Sz Rx Ry Rz th\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			Sx=atof(s3.c_str());
			Sy=atof(s4.c_str());
			Sz=atof(s5.c_str());
			Rx=atof(s6.c_str());
			Ry=atof(s7.c_str());
			Rz=atof(s8.c_str());
			th=atof(s9.c_str());			
			RotateBlock(i1,i2,Sx,Sy,Sz,Rx,Ry,Rz,th);
		}
		return s;
	}
	if(s0=="RotateBlockAroundPupil" || s0=="rotateblockaroundpupil"){
		int i1,i2,stop;
		double rox;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="RotateBlockAroundPupil i1 i2 stop rox\n";
		}
		else if( b1 && b2 && b3 && b4 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			stop=atoi(s3.c_str());
			rox=atof(s4.c_str());
			RotateBlockAroundPupil(i1,i2,stop,rox);
		}
		return s;
	}
	if(s0=="RotateBlockX" || s0=="rotateblockx"){
		int i1,i2;
		double th;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="RotateBlockX i1 i2 th\n";
		}
		else if( b1 && b2 && b3 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			th=atof(s3.c_str());			
			RotateBlockX(i1,i2,th);
		}
		return s;
	}
	if(s0=="RotateBlockY" || s0=="rotateblocky"){
		int i1,i2;
		double th;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="RotateBlockY i1 i2 th\n";
		}
		else if( b1 && b2 && b3 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			th=atof(s3.c_str());			
			RotateBlockY(i1,i2,th);
		}
		return s;
	}
	if(s0=="RotateSurface" || s0=="rotatesurface"){
		int i;
		double Sx,Sy,Sz,Rx,Ry,Rz,th;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="RotateSurface i Sx Sy Sz Rx Ry Rz th\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 ) {
			i =atoi(s1.c_str());
			Sx=atof(s2.c_str());
			Sy=atof(s3.c_str());
			Sz=atof(s4.c_str());
			Rx=atof(s5.c_str());
			Ry=atof(s6.c_str());
			Rz=atof(s7.c_str());
			th=atof(s8.c_str());
			RotateSurface(i,Sx,Sy,Sz,Rx,Ry,Rz,th);
		}
		return s;
	}
	if(s0=="RoundA" || s0=="rounda"){
		int n;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="RoundA n\n";
		}
		else if(b1) {
			n=atoi(s1.c_str());
			RoundA(n);
		}
		return s;
	}
	if(s0=="RToNewton" || s0=="rtonewton"){
		double r0,eaphi,FRW_nm,r;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="RToNewton r0 eaphi FRW_nm r\n";
		}
		else if( b1 && b2 && b3 && b4 ){
			r0=atof(s1.c_str());
			eaphi=atof(s2.c_str());
			FRW_nm=atof(s3.c_str());
			r=atof(s4.c_str());
			sprintf(buf,"r0=%g eaphi=%g FRW=%gnm r=%g\n", r0,eaphi,FRW_nm,r);
			s+=buf;
			sprintf(buf,"rings=%g\n", RToNewton(r0,eaphi,FRW_nm,r));
			s+=buf;
		}
		return s;
	}
	if(s0=="s1"){
		int j;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="s1 [j=1]\n";
		}
		else if(b1){
			j=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",this->s1(j)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",this->s1()); s=buf;
		}
		return s;
	}
	if(s0=="s1i"){
		int i;
		double ss;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="s1i [s] i\n";
		}
		else if(b1 && b2){
			ss=atof(s1.c_str());
			i=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",s1i(ss,i,1)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",s1i(i)); s=buf;
		}
		return s;
	}
	if(s0=="S1x" || s0=="s1x"){
		double yObj,xObj,yPupil,xPupil;
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s= "S1x yObj xObj yPupil xPupil i\n";
			s+="S1x [i=k]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			i=atoi(s5.c_str());
			sprintf(buf,"%.15g\n",S1x(yObj,xObj,yPupil,xPupil,i,1)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",S1x(i)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",S1x()); s=buf;
		}
		return s;
	}
	if(s0=="S1y" || s0=="s1y"){
		double yObj,xObj,yPupil,xPupil;
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s= "S1y yObj xObj yPupil xPupil i\n";
			s+="S1y [i=k]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			yPupil=atof(s3.c_str());
			xPupil=atof(s4.c_str());
			i=atoi(s5.c_str());
			sprintf(buf,"%.15g\n",S1y(yObj,xObj,yPupil,xPupil,i,1)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",S1y(i)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",S1y()); s=buf;
		}
		return s;
	}
	if(s0=="SA0" || s0=="sa0"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="SA0 i1 i2 (calculated value)\n";
			s+="SA0 i     (surface property)\n";
			s+="SA0       (calculated value)\n";
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",SA0(i1,i2)); s=buf;
		}
		else if( b1 ){
			i1=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",SA0(i1)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",SA0()); s=buf;
		}
		return s;
	}
	if(s0=="SC" || s0=="sc"){
		s1=arg(com,1);
		if(s1=="?"){
			s="SC (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",SC()); s=buf;
		}
		return s;
	}
	if(s0=="ScheimpflugImagePlaneTiltX" || s0=="scheimpflugimageplanetiltx"){
		s1=arg(com,1);
		if(s1=="?"){
			s="ScheimpflugImagePlaneTiltX (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",ScheimpflugImagePlaneTiltX()); s=buf;
		}
		return s;
	}
	if(s0=="ScheimpflugImagePlaneTiltY" || s0=="scheimpflugimageplanetilty"){
		s1=arg(com,1);
		if(s1=="?"){
			s="ScheimpflugImagePlaneTiltY (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",ScheimpflugImagePlaneTiltY()); s=buf;
		}
		return s;
	}
	if(s0=="SetHyperboloid" || s0=="sethyperboloid"){
		int i;
		double f1,f2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="SetHyperboloid i f1 f2\n";
		}
		else if( b1 && b2 && b3 ){
			i=atoi(s1.c_str());
			f1=atof(s2.c_str());
			f2=atof(s3.c_str());
			SetHyperboloid(i,f1,f2);
		}
		return s;
	}
	if(s0=="SetHyperboloid2" || s0=="sethyperboloid2"){
		int i,i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="SetHyperboloid i i1 i2\n";
		}
		else if( b1 && b2 && b3 ){
			i=atoi(s1.c_str());
			i1=atoi(s2.c_str());
			i2=atoi(s3.c_str());
			SetHyperboloid2(i,i1,i2);
		}
		return s;
	}
	if(s0=="SetM" || s0=="setm"){
		double value;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="SetM value\n";
		}
		else if( b1 ) {
			value=atof(s1.c_str());
			SetM(value);
		}
		return s;
	}
	if(s0=="SetParaboloid" || s0=="setparaboloid"){
		int i;
		double f;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="SetParaboloid i f\n";
		}
		else if( b1 && b2 ){
			i=atoi(s1.c_str());
			f=atof(s2.c_str());
			SetParaboloid(i,f);
		}
		return s;
	}
	if(s0=="SetParaboloid2" || s0=="setparaboloid2"){
		int i,i_f,i_inf;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="SetParaboloid i i_f i_inf\n";
		}
		else if( b1 && b2 && b3 ){
			i=atoi(s1.c_str());
			i_f=atoi(s2.c_str());
			i_inf=atoi(s3.c_str());
			SetParaboloid2(i,i_f,i_inf);
		}
		return s;
	}
	if(s0=="SetSpheroid" || s0=="setspheroid"){
		int i,convex;
		char majoraxis;
		double rshort,ftof;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		if(s1=="?"){
			s+="SetSpheroid i majoraxis convex a ftof\n";
			s+="(ex) SetSpheroid 1 y 0 100 200\n";
		}
		else if( b1 && b3 && b4 && b5 ){
			i=atoi(s1.c_str());
			majoraxis=s2[0];
			convex=atoi(s3.c_str());
			rshort=atof(s4.c_str());
			ftof=atof(s5.c_str());
			sprintf(buf,"%.15g\n",SetSpheroid(i,majoraxis,convex,rshort,ftof));
			s=buf;			
		}
		return s;
	}
	if(s0=="SetSpheroid2" || s0=="setspheroid2"){
		int i,i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="SetSpheroid2 i i1 i2\n";
		}
		else if( b1 && b2 && b3 ){
			i=atoi(s1.c_str());
			i1=atoi(s2.c_str());
			i2=atoi(s3.c_str());
			SetSpheroid2(i,i1,i2);
		}
		return s;
	}
	if(s0=="SetSpheroid3" || s0=="setspheroid3"){
		int i;
		double z1,z2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="SetSpheroid3 i z1 z2\n";
		}
		else if( b1 && b2 && b3 ){
			i=atoi(s1.c_str());
			z1=atof(s2.c_str());
			z2=atof(s3.c_str());
			SetSpheroid3(i,z1,z2);
		}
		return s;
	}
	if(s0=="SetSpline" || s0=="setspline"){
		int i,n;
		double step;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="SetSpline i n step\n";
		}
		else if(b1 && b2 && b3){
			i=atoi(s1.c_str());
			n=atoi(s2.c_str());
			step=atof(s3.c_str());
			SetSpline(i,n,step);
		}
		return s;
	}
	if(s0=="ShapeFactor" || s0=="shapefactor"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="ShapeFactor i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",ShapeFactor(i));
			s=buf;
		}
		return s;
	}
	if(s0=="si"){
		int i;
		double ss;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="si [s] i\n";
		}
		else if(b1 && b2){
			ss=atof(s1.c_str());
			i=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",si(ss,i,1)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",si(i)); s=buf;
		}
		return s;
	}
	if(s0=="SingleRayTrace" || s0=="singleraytrace"){
		double yObj,xObj,yPupil,xPupil,defocus;
		std::string SetRay;
		int findpupil,lastsurf,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s="SingleRayTrace yObj xObj SetRay findpupil yPupil xPupil defocus lastsurf j\n";
		}
		else if( b1 && b2 && b4 && b5 && b6 && b7 && b8 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			SetRay=s3;
			findpupil=atoi(s4.c_str());
			yPupil=atof(s5.c_str());
			xPupil=atof(s6.c_str());
			defocus=atof(s7.c_str());
			lastsurf=atoi(s8.c_str());
			j=atoi(s9.c_str());
			s=SingleRayTrace(yObj,xObj,SetRay,findpupil,yPupil,xPupil,defocus,lastsurf,j);
		}
		return s;
	}
	if(s0=="Sph" || s0=="sph"){
		double yObj,xObj;
		int findpupil,i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Sph [yObj=yObjectMax xObj=xObjectMax [findpupil=0 [i=0]]]\n";
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			i=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",Sph(yObj,xObj,0,findpupil,i)); s=buf;
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",Sph(yObj,xObj,0,findpupil)); s=buf;
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",Sph(yObj,xObj,0)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Sph()); s=buf;
		}
		return s;
	}
	if(s0=="SphEq" || s0=="spheq"){
		double yObj,xObj;
		int findpupil,i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="SphEq [yObj=yObjectMax xObj=xObjectMax [findpupil=0 [i=0]]]\n";
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			i=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",SphEq(yObj,xObj,0,findpupil,i)); s=buf;
		}
		else if(b1 && b2 && b3){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			findpupil=atoi(s3.c_str());
			sprintf(buf,"%.15g\n",SphEq(yObj,xObj,0,findpupil)); s=buf;
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",SphEq(yObj,xObj,0)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",SphEq()); s=buf;
		}
		return s;
	}
	if(s0=="SphEq50" || s0=="spheq50"){
		s1=arg(com,1);
		if(s1=="?"){
			s="SphEq50 (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",SphEq50()); s=buf;
		}
		return s;
	}
	if(s0=="SphEq70" || s0=="spheq70"){
		s1=arg(com,1);
		if(s1=="?"){
			s="SphEq70 (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",SphEq70()); s=buf;
		}
		return s;
	}
	if(s0=="SplineDoubleN" || s0=="splinedoublen"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="SplineDoubleN i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			SplineDoubleN(i);
		}
		return s;
	}
	if(s0=="SplineN" || s0=="splinen"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="SplineN i\n";
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%d\n",SplineN(i)); s=buf;
		}
		return s;
	}
	if(s0=="SplineSetHStep" || s0=="splinesethstep"){
		int i;
		double dh;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="SplineSetHStep i dh\n";
		}
		else if(b1 & b2){
			i=atoi(s1.c_str());
			dh=atof(s2.c_str());
			SplineSetHStep(i,dh);
		}
		return s;
	}
	if(s0=="SPO" || s0=="spo"){
		double yObj,xObj,weight;
		int findpupil,colors,IgnoreTC,yfan;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="SPO yObj xObj findpupil colors [IgnoreTC=0 [weight=1 [yfan=0]]]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			weight=atof(s6.c_str());
			yfan=atoi(s7.c_str());
			s=SPO(yObj,xObj,findpupil,colors,IgnoreTC,weight,yfan);
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			weight=atof(s6.c_str());
			s=SPO(yObj,xObj,findpupil,colors,IgnoreTC,weight);
		}
		else if(b1 && b2 && b3 && b4 && b5){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			s=SPO(yObj,xObj,findpupil,colors,IgnoreTC);
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			s=SPO(yObj,xObj,findpupil,colors);
		}
		return s;
	}
	if(s0=="SPO2" || s0=="spo2"){
		double yObj,xObj,weight;
		int findpupil,colors,n,IgnoreTC;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="SPO2 yObj xObj findpupil colors [IgnoreTC=0 [n=0 [weight=1]]]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			n=atoi(s6.c_str());
			weight=atof(s7.c_str());
			s=SPO2(yObj,xObj,findpupil,colors,IgnoreTC,n,weight);
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			n=atoi(s6.c_str());
			s=SPO2(yObj,xObj,findpupil,colors,IgnoreTC,n);
		}
		else if(b1 && b2 && b3 && b4 && b5){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			s=SPO2(yObj,xObj,findpupil,colors,IgnoreTC);
		}
		else if(b1 && b2 && b3 && b4){
			yObj=atof(s1.c_str());
			xObj=atoi(s2.c_str());
			findpupil=atoi(s3.c_str());
			colors=atoi(s4.c_str());
			IgnoreTC=atoi(s5.c_str());
			s=SPO2(yObj,xObj,findpupil,colors);
		}
		return s;
	}
	if(s0=="SpotXGravityCenter" || s0=="spotxgravitycenter"){
		s1=arg(com,1);
		if(s1=="?"){
			s="SpotXGravityCenter (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",SpotXGravityCenter()); s=buf;
		}
		return s;
	}
	if(s0=="SpotYGravityCenter" || s0=="spotygravitycenter"){
		s1=arg(com,1);
		if(s1=="?"){
			s="SpotYGravityCenter (no argument)\n";
		}
		else{
			sprintf(buf,"%.15g\n",SpotYGravityCenter()); s=buf;
		}
		return s;
	}
	if(s0=="Steepness" || s0=="steepness"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Steepness i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Steepness(i)); s=buf;
		}
		return s;
	}
	if(s0=="StrehlDef" || s0=="strehldef"){
		double yObj,xObj,defocus;
		int j,FindPupil,AdjustSph,OptimizeDefocusOnAxis;
		s1=arg(com,1); b1=is_numeric(s1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5); 
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		if(s1=="?"){
			s="StrehlDef [yObj=xObj=0 [defocus=0 j=1 FindPupil=1 AdjustSph=0 OptimizeDefousOnAxis=1]]\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6 && b7 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			defocus=atof(s3.c_str());
			j=atoi(s4.c_str());
			FindPupil=atoi(s5.c_str());
			AdjustSph=atoi(s6.c_str());
			OptimizeDefocusOnAxis=atoi(s7.c_str());
			sprintf(buf,"%.15g\n",StrehlDef(yObj,xObj,defocus,j,FindPupil,AdjustSph,OptimizeDefocusOnAxis,""));
			s=buf;
		}
		else if( b1 && b2 ){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",StrehlDef(yObj,xObj));
			s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",StrehlDef(0,0));
			s=buf;
		}
		return s;
	}
	if(s0=="SurfaceSag" || s0=="surfacesag"){
		int i;
		double y,x;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="SurfaceSag i y x\n";
		}
		else if( b1 && b2 && b3 ){
			i=atoi(s1.c_str());
			y=atof(s2.c_str());
			x=atof(s3.c_str());
			sprintf(buf, "%.15g\n", surface_sag(i,y,x,0));
			s=buf;
		}
		return s;
	}
	if(s0=="SurfaceSagMax" || s0=="surfacesagmax"){
		int i;
		double hStep;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="SurfaceSagMax i hStep\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			hStep=atof(s2.c_str());
			sprintf(buf, "%.15g\n", SurfaceSagMax(i,hStep));
			s=buf;
		}
		return s;
	}
	if(s0=="SurfaceSagMin" || s0=="surfacesagmin"){
		int i;
		double hStep;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="SurfaceSagMin i hStep\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			hStep=atof(s2.c_str());
			sprintf(buf, "%.15g\n", SurfaceSagMin(i,hStep));
			s=buf;
		}
		return s;
	}
	if(s0=="SurfaceSagTable" || s0=="surfacesagtable"){
		int i;
		double hStep;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="SurfaceSagTable i hStep\n";
		}
		else if( b1 && b2 ){
			i=atoi(s1.c_str());
			hStep=atof(s2.c_str());
			s=SurfaceSagTable(i,hStep);
		}
		return s;
	}
	if(s0=="SurfaceSlope" || s0=="surfaceslope"){
		int i;
		double y,x;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="SurfaceSlope i y x\n";
		}
		else if( b1 && b2 && b3 ){
			i=atoi(s1.c_str());
			y=atof(s2.c_str());
			x=atof(s3.c_str());
			sprintf(buf, "%.15g\n", SurfaceSlope(i,y,x));
			s=buf;
		}
		return s;
	}
	if(s0=="SurfaceSlopeMax" || s0=="surfaceslopemax"){
		int i;
		double hStep;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="SurfaceSlopeMax i hStep\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			hStep=atof(s2.c_str());
			sprintf(buf, "%.15g\n", SurfaceSlopeMax(i,hStep));
			s=buf;
		}
		return s;
	}
	if(s0=="SurfaceSlopeTable" || s0=="surfaceslopetable"){
		int i;
		double hStep;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="SurfaceSlopeTable i hStep\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			hStep=atof(s2.c_str());
			s=SurfaceSlopeTable(i,hStep);
		}
		return s;
	}
	if(s0=="SwapObjPupil" || s0=="swapobjpupil"){
		s1=arg(com,1);
		if(s1=="?"){
			s="SwapObjPupil (no argument)\n";
		}
		else{
			SwapObjPupil();
		}
		return s;
	}
	if(s0=="t1" ){
		int j;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="t1 [j=1]\n";
		}
		if(b1){
			j=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",t1(j)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",t1()); s=buf;
		}
		return s;
	}
	if(s0=="t1i"){
		int i,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="t1i i [j=1]\n";
		}
		else if(b1 && b2) {
			i=atoi(s1.c_str());
			j=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",t1i(i,j));
			s=buf;
		}
		else if(b1) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",t1i(i));
			s=buf;
		}
		return s;
	}
	if(s0=="TangentialCurvature" || s0=="tangentialcurvature"){
		int i;
		double y;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TangentialCurvature i y\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			y=atof(s2.c_str());
			sprintf(buf,"%.15g\n",TangentialCurvature(i,y));
			s=buf;
		}
		return s;
	}
	if(s0=="TangentialCurvatureMax" || s0=="tangentialcurvaturemax"){
		int i;
		double hStep;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TangentialCurvatureMax i hStep\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			hStep=atof(s2.c_str());
			sprintf(buf,"%.15g\n",TangentialCurvatureMax(i,hStep));
			s=buf;
		}
		return s;
	}
	if(s0=="tCalc" || s0=="tcalc"){
		s1=arg(com,1);
		if(s1=="?"){
			s="tCalc (no argements)\n";
		}
		else{
			tCalc();
		}
		return s;
	}
	if(s0=="TelecentricityObj" || s0=="telecentricityobj"){
		s1=arg(com,1);
		if(s1=="?"){
			s="TelecentricityObj (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",TelecentricityObj()); s=buf;
		}
		return s;
	}
	if(s0=="TelephotoF1" || s0=="telephotof1"){
		double L,e,f;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="TelephotoF1 L e f\n";
		}
		else if(b1 && b2 && b3){
			L=atof(s1.c_str());
			e=atof(s2.c_str());
			f=atof(s3.c_str());
			sprintf(buf,"%.15g\n", TelephotoF1(L,e,f)); s=buf;
		}
		return s;
	}
	if(s0=="TelephotoF2" || s0=="telephotof2"){
		double L,e,f;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		if(s1=="?"){
			s="TelephotoF2 L e f\n";
		}
		else if(b1 && b2 && b3){
			L=atof(s1.c_str());
			e=atof(s2.c_str());
			f=atof(s3.c_str());
			sprintf(buf,"%.15g\n", TelephotoF2(L,e,f)); s=buf;
		}
		return s;
	}
	if(s0=="ThinT" || s0=="thint"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ThinT i1 i2\n";
		}
		else if(b1 && b2){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n", ThinT(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n", ThinT()); s=buf;
		}
		return s;
	}
	if(s0=="ti"){
		int i,j;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ti i [j=1]\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			j=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",ti(i,j));
			s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",ti(i));
			s=buf;
		}
		return s;
	}
	if(s0=="tLack" || s0=="tlack"){
		s1=arg(com,1);
		if(s1=="?"){
			s="tLack (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",tLack()); s=buf;
		}
		return s;
	}
	if(s0=="ToAplanaticSurf" || s0=="toaplanaticsurf"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="ToAplanaticSurf i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",ToAplanaticSurf(i));
			s=buf;
		}
		return s;
	}
	if(s0=="ToBlock" || s0=="toblock"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ToBlock i1 i2\n";
		}
		else if(b1 && b2) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			ToBlock(i1,i2);
		}
		return s;
	}
	if(s0=="ToConcentricSurf" || s0=="toconcentricsurf"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="ToConcentricSurf i\n";
		}
		else if( b1 ) {
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",ToConcentricSurf(i));
			s=buf;
		}
		return s;
	}
	if(s0=="Toff" || s0=="toff"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Toff (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",Toff()); s=buf;
		}
		return s;
	}
	if(s0=="ToIdealLens" || s0=="toideallens"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ToIdealLens i1 i2\n";
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			ToIdealLens(i1,i2);
		}
		return s;
	}
	if(s0=="ToLMTestLens" || s0=="tolmtestlens"){
		double sph,cyl;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ToLMTestLens sph cyl\n";
		}
		else if( b1 && b2 ){
			sph=atof(s1.c_str());
			cyl=atof(s2.c_str());
			ToLMTestLens(sph,cyl);
		}
		return s;
	}
	if(s0=="Ton" || s0=="ton"){
		s1=arg(com,1);
		if(s1=="?"){
			s="Ton (no arguments)\n";
		}
		else{
			sprintf(buf,"%.15g\n",Ton()); s=buf;
		}
		return s;
	}
	if(s0=="ToroidRa" || s0=="toroidra"){
		double x,y,rx,ry,kp;
		int IsXToroid;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="ToricRa x y rx ry kp IsXToroid\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6) {
			x=atof(s1.c_str());
			y=atof(s2.c_str());
			rx=atof(s3.c_str());
			ry=atof(s4.c_str());
			kp=atof(s5.c_str());
			IsXToroid=atoi(s6.c_str());
			sprintf(buf,"%.10g\n",ToroidRa(x,y,rx,ry,kp,IsXToroid)); s=buf;
		}
		return s;
	}
	if(s0=="ToroidZ" || s0=="toroidz"){
		double x,y,rx,ry,kp;
		int IsXToroid;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5);
		s6=arg(com,6); b6=is_numeric(s6);
		if(s1=="?"){
			s="ToroidZ x y rx ry kp IsXToroid\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6) {
			x=atof(s1.c_str());
			y=atof(s2.c_str());
			rx=atof(s3.c_str());
			ry=atof(s4.c_str());
			kp=atof(s5.c_str());
			IsXToroid=atoi(s6.c_str());
			sprintf(buf,"%.10g\n",ToroidZ(x,y,rx,ry,kp,IsXToroid)); s=buf;
		}
		return s;
	}
	if(s0=="ToScheimpflugImagePlane" || s0=="toscheimpflugimageplane"){
		s1=arg(com,1);
		if(s1=="?"){
			s="ToScheimpflugImagePlane (no arguments)\n";
		}
		else{
			ToScheimpflugImagePlane();
		}
		return s;
	}
	if(s0=="TotalInAirThickness" || s0=="totalinairthickness"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TotalInAirThickness [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",TotalInAirThickness(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",TotalInAirThickness()); s=buf;
		}
		return s;
	}
	if(s0=="TotalOpticalThickness" || s0=="totalopticalthickness"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TotalOpticalThickness [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",TotalOpticalThickness(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",TotalOpticalThickness()); s=buf;
		}
		return s;
	}
	if(s0=="TotalThickness" || s0=="totalthickness"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TotalThickness [i1=1 i2=k]\n";
		}
		else if( b1 && b2 ) {
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			sprintf(buf,"%.15g\n",TotalThickness(i1,i2)); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",TotalThickness()); s=buf;
		}
		return s;
	}
	if(s0=="ToThinLens" || s0=="tothinlens"){
		int i1,i2;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ToThinLens i1 i2\n";
		}
		else if( b1 && b2 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			ToThinLens(i1,i2);
		}
		return s;
	}
	if(s0=="TransformACoefficients" || s0=="transformacoefficients"){
		int i;
		double newNormH;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="TransformACoefficients i newNormH\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			newNormH=atof(s2.c_str());
			TransformACoefficients(i,newNormH);
		}
		return s;
	}
	if(s0=="Transmittance" || s0=="transmittance"){
		double yObj,xObj;
		int IsAreaSource,IsLambert;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Transmittance yObj xObj IsAreaSource IsLambert\n";
		}
		else if( b1 && b2 && b3 && b4 ) {
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			IsAreaSource=atoi(s3.c_str());
			IsLambert=atoi(s4.c_str());
			sprintf(buf,"%.15g\n",Transmittance(yObj,xObj,IsAreaSource,IsLambert)); s=buf;
		}
		return s;
	}
	if(s0=="vdiopter1i"){
		double ss;
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="vdiopter1i [s] i\n";
		}
		else if(b1 && b2){
			ss=atof(s1.c_str());
			i=atoi(s2.c_str());
			sprintf(buf,"%.15g\n", vdiopter1i(ss,i,1)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",vdiopter1i(i,1)); s=buf;
		}
		return s;
	}
	if(s0=="vertex"){
		int i;
		double si;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="vertex i [si]\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			si=atof(s2.c_str());
			sprintf(buf,"%.15g\n",vertex(i,si)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",vertex(i)); s=buf;
		}
		return s;
	}
	if(s0=="vertex1"){
		int i;
		double si;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="vertex1 i [si]\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			si=atof(s2.c_str());
			sprintf(buf,"%.15g\n",vertex1(i,si)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",vertex1(i)); s=buf;
		}
		return s;
	}
	if(s0=="VertexGlobal" || s0=="vertexglobal"){
		int i;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="VertexGlobal i\n";
		}
		else if(b1){
			vector<double> v;
			i=atoi(s1.c_str());
			v=VertexGlobal(i);
			sprintf(buf,"%.15g %.15g %.15g\n",v.x,v.y,v.z); s=buf;
		}
		return s;
	}
	if(s0=="Wl" || s0=="wl"){
		int j;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="Wl j\n";
		}
		else if(b1){
			j=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Wl(j)); s=buf;
		}
		return s;
	}
	if(s0=="xObjectMaxAng" || s0=="xobjectmaxang"){
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="xObjectMaxAng [new_value]\n";
		}
		else if(b1){
			val=atof(s1.c_str());
			Set_xObjectMaxAng(val);
			sprintf(buf,"%.15g\n",Get_xObjectMaxAng()); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Get_xObjectMaxAng()); s=buf;
		}
		return s;
	}
	if(s0=="xScan" || s0=="xscan"){
		std::string RotateAxisXY;
		double th_deg;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="xScan [RotateAxisXY='y'] th_deg\n";
			s+="xScan [RotateAxisXY='y']\n";
			s+="    (RotateAxisXY='x' or 'y')\n";
		}
		else if(b1){
			th_deg=atof(s1.c_str());
			xScan(th_deg);
		}
		else if(b2){
			RotateAxisXY=s1;
			th_deg=atof(s2.c_str());
			xScan(RotateAxisXY,th_deg);
		}
		else{
			RotateAxisXY=s1;
			if(RotateAxisXY==""){
				sprintf(buf,"%.15g\n",xScan()); s=buf;
			}
			else{
				sprintf(buf,"%.15g\n",xScan(RotateAxisXY)); s=buf;
			}
		}
		return s;
	}
	if(s0=="xUsedRange" || s0=="xusedrange"){
		double xObjMax,xObjMin;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="xUsedRange xObjMax xObjMin\n";
		}
		else if( b1 & b2 ){
			xObjMax=atof(s1.c_str());
			xObjMin=atof(s2.c_str());
			s=xUsedRange(xObjMax,xObjMin);
		}
		return s;
	}
	if(s0=="yObjectMaxAng" || s0=="yobjectmaxang"){
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="yObjectMaxAng [new_value]\n";
		}
		else if(b1){
			val=atof(s1.c_str());
			Set_yObjectMaxAng(val);
			sprintf(buf,"%.15g\n",Get_yObjectMaxAng()); s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",Get_yObjectMaxAng()); s=buf;
		}
		return s;
	}
	if(s0=="yScan" || s0=="yscan"){
		std::string RotateAxisXY;
		double th_deg;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s+="yScan [RotateAxisXY='x'] th_deg\n";
			s+="yScan [RotateAxisXY='x']\n";
			s+="    (RotateAxisXY='x' or 'y')\n";
		}
		else if(b1){
			th_deg=atof(s1.c_str());
			yScan(th_deg);
		}
		else if(b2){
			RotateAxisXY=s1;
			th_deg=atof(s2.c_str());
			yScan(RotateAxisXY,th_deg);
		}
		else{
			RotateAxisXY=s1;
			if(RotateAxisXY==""){
				sprintf(buf,"%.15g\n",yScan()); s=buf;
			}
			else{
				sprintf(buf,"%.15g\n",yScan(RotateAxisXY)); s=buf;
			}
		}
		return s;
	}
	if(s0=="yUsedRange" || s0=="yusedrange"){
		double yObjMax,yObjMin;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="yUsedRange yObjMax yObjMin\n";
		}
		else if( b1 & b2 ){
			yObjMax=atof(s1.c_str());
			yObjMin=atof(s2.c_str());
			s=yUsedRange(yObjMax,yObjMin);
		}
		return s;
	}
	if(s0=="yVignetting" || s0=="yvignetting"){
		double yObj;
		int findpupil;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="yVignetting [yObj=yObjectmax findpupil=1]\n";
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			findpupil=atoi(s2.c_str());
			sprintf(buf,"%.15g",yVignetting(yObj,findpupil)); s=buf;
		}
		else{
			sprintf(buf,"%.15g",yVignetting()); s=buf;
		}
		return s;
	}
	if(s0=="ZernikeC" || s0=="zernikec"){
		double yObj,xObj,defocus;
		int term_no,j,FindPupil,InLambda,AdjustSph,OptimizeDefocus;
		s1=arg(com,1); b1=is_numeric(s1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5); 
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		s9=arg(com,9); b9=is_numeric(s9);
		if(s1=="?"){
			s="ZernikeC term_no [yObj=0 xObj=0 [defocus=0 j=1 FindPupil=1 InLambda=1 AdjustSph=0 OptimizeDefous=1]]\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9 ){
			term_no=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			defocus=atof(s4.c_str());
			j=atoi(s5.c_str());
			FindPupil=atoi(s6.c_str());
			InLambda=atoi(s7.c_str());
			AdjustSph=atoi(s8.c_str());
			OptimizeDefocus=atoi(s9.c_str());
			sprintf(buf,"%.15g\n",ZernikeC(term_no,yObj,xObj,defocus,j,FindPupil,InLambda,AdjustSph,OptimizeDefocus));
			s=buf;
		}
		else if( b1 && b2 && b3 ){
			term_no=atoi(s1.c_str());
			yObj=atof(s2.c_str());
			xObj=atof(s3.c_str());
			sprintf(buf,"%.15g\n",ZernikeC(term_no,yObj,xObj,0,1,1,1,0,1));
			s=buf;
		}
		else if( b1 ){
			term_no=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",ZernikeC(term_no,0,0,0,1,1,1,0,1));
			s=buf;
		}
		return s;
	}
	if(s0=="ZernikeCHigh" || s0=="zernikechigh"){
		double yObj,xObj,defocus;
		int j,FindPupil,InLambda,AdjustSph,OptimizeDefocus;
		s1=arg(com,1); b1=is_numeric(s1); 
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		s5=arg(com,5); b5=is_numeric(s5); 
		s6=arg(com,6); b6=is_numeric(s6);
		s7=arg(com,7); b7=is_numeric(s7);
		s8=arg(com,8); b8=is_numeric(s8);
		if(s1=="?"){
			s="ZernikeC [yObj=0 xObj=0 [defocus=0 j=1 FindPupil=0 InLambda=0 AdjustSph=0 OptimizeDefous=0]]\n";
		}
		else if(b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			defocus=atof(s3.c_str());
			j=atoi(s4.c_str());
			FindPupil=atoi(s5.c_str());
			InLambda=atoi(s6.c_str());
			AdjustSph=atoi(s7.c_str());
			OptimizeDefocus=atoi(s8.c_str());
			sprintf(buf,"%.15g\n",ZernikeCHigh(yObj,xObj,defocus,j,FindPupil,InLambda,AdjustSph,OptimizeDefocus));
			s=buf;
		}
		else if(b1 && b2){
			yObj=atof(s1.c_str());
			xObj=atof(s2.c_str());
			sprintf(buf,"%.15g\n",ZernikeCHigh(yObj,xObj,0,1,0,0,0,0));
			s=buf;
		}
		else{
			sprintf(buf,"%.15g\n",ZernikeCHigh(0,0,0,1,0,0,0,0));
			s=buf;
		}
		return s;
	}
	if(s0=="ZernikeMaxOrder" || s0=="zernikemaxorder"){
		int i,val;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ZernikeMaxOrder i [new_value]\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			val=atoi(s2.c_str());
			Set_ZernikeMaxOrder(i,val);
			sprintf(buf,"%d\n",Get_ZernikeMaxOrder(i)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%d\n",Get_ZernikeMaxOrder(i)); s=buf;
		}
		return s;
	}
	if(s0=="ZernikeR0" || s0=="zerniker0"){
		int i;
		double val;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ZernikeR0 i [new_value]\n";
		}
		else if(b1 && b2){
			i=atoi(s1.c_str());
			val=atof(s2.c_str());
			Set_ZernikeR0(i,val);
			sprintf(buf,"%.15g\n",Get_ZernikeR0(i)); s=buf;
		}
		else if(b1){
			i=atoi(s1.c_str());
			sprintf(buf,"%.15g\n",Get_ZernikeR0(i)); s=buf;
		}
		return s;
	}
	if(s0=="Zoom" || s0=="zoom"){
		int i1,i2,i3;
		double MorF;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		s3=arg(com,3); b3=is_numeric(s3);
		s4=arg(com,4); b4=is_numeric(s4);
		if(s1=="?"){
			s="Zoom i1 i2 i3 m(f)\n";
		}
		else if( b1 && b2 && b3 && b4 ){
			i1=atoi(s1.c_str());
			i2=atoi(s2.c_str());
			i3=atoi(s3.c_str());
			MorF=atof(s4.c_str());
			if( !Zoom(i1,i2,i3,MorF) ) s="Can't zoom, no solution etc..\n";
		}
		return s;
	}
	if(s0=="ZoomCamTable" || s0=="zoomcamtable"){
		double fv,mo,th1,m_th1,th2,m_th2,th_start,th_end,th_step;
		int variator_leading;
		s1 =arg(com,1);  b1 =is_numeric(s1);
		s2 =arg(com,2);  b2 =is_numeric(s2);
		s3 =arg(com,3);  b3 =is_numeric(s3);
		s4 =arg(com,4);  b4 =is_numeric(s4);
		s5 =arg(com,5);  b5 =is_numeric(s5);
		s6 =arg(com,6);  b6 =is_numeric(s6);
		s7 =arg(com,7);  b7 =is_numeric(s7);
		s8 =arg(com,8);  b8 =is_numeric(s8);
		s9 =arg(com,9);  b9 =is_numeric(s9);
		s10=arg(com,10); b10=is_numeric(s10);
		if(s1=="?"){
			s+="ZoomCamTable fv mo th1 m_th1 th2 m_th2 th_start th_end th_step variator_leading\n";
			s+="    fv       : variator focallength\n";
			s+="    mo       : magnification(or focallength) without variator\n";
			s+="    th1      : a position of com\n";
			s+="    m_th1    : total magnification(or focallength) at th=th1\n";
			s+="    th2      : another position of cam\n";
			s+="    m_th2    : total magnification(or focallength) at th=th2\n";
			s+="    th_start : start of th\n";
			s+="    th_end   : end of th\n";
			s+="    th_step  : step of th\n";
			s+="    variator_leading : Is variator leading compensator?\n";
			s+="    <note> x=y=0 when th=0\n";
		}
		else if( b1 && b2 && b3 && b4 && b5 && b6 && b7 && b8 && b9 && b10 ){
			fv =atof(s1.c_str());
			mo =atof(s2.c_str());
			th1=atof(s3.c_str());
			m_th1=atof(s4.c_str());
			th2  =atof(s5.c_str());
			m_th2=atof(s6.c_str());
			th_start=atof(s7.c_str());
			th_end  =atof(s8.c_str());
			th_step =atof(s9.c_str());
			variator_leading=atoi(s10.c_str());
			s=ZoomCamTable(fv,mo,th1,m_th1,th2,m_th2,th_start,th_end,th_step,variator_leading);
		}
		return s;
	}
	if(s0=="zScan" || s0=="zscan"){
		double dz;
		s1=arg(com,1); b1=is_numeric(s1);
		if(s1=="?"){
			s="zScan dz\n";
		}
		else if(b1){
			dz=atof(s1.c_str());
			zScan(dz);
		}
		else{
			sprintf(buf,"%.15g\n",zScan()); s=buf;
		}
		return s;
	}
	if(s0=="ZValue" || s0=="zvalue"){
		int i,i1;
		s1=arg(com,1); b1=is_numeric(s1);
		s2=arg(com,2); b2=is_numeric(s2);
		if(s1=="?"){
			s="ZValue i\nZValue i1 i2\n";
		}
		else if( b1 ){
			i=atoi(s1.c_str());
			i1=i+1;
			if( b2 ){
				i1=atoi(s2.c_str());
			}
			sprintf(buf,"%d-%d ZValue=%g\n", i,i1,ZValue(i,i1));
			s=buf;
		}
		return s;
	}
	
	return s="";
}

double* cLens1::property(std::string& s,int &args){
	// arg(s,0)のプロパティのデータのアドレスを返す．
	// さらに，arg(s,0)とその引数をsから削除する（optimize()で変化量を表す引数を抜き出すのに必要）．
	// arg(s,0)がプロパティを示すものでなければ0を返す．
	//
	// argsにより，引数の数を代入する．sに対応するプロパティがないときは-1を代入する．

	std::string s0,s1,s2,s3;
	bool b1,b2,b3;
	double *x=0;
	int *ix=0;

	args=-1;
	s0=arg(s,0);
		
	{   // 0変数プロパティ /////////////////////////////
		double *p=0;
		
		if     (s0=="EPD"  || s0=="epd" ) p=&EPD;
		else if(s0=="EPDx" || s0=="epdx") p=&EPDx;
		else if(s0=="EPx" || s0=="epx") p=&EPx;
		else if(s0=="EPy" || s0=="epy") p=&EPy;
		else if(s0=="s") p=&(this->s);
		else if(s0=="t") p=&t;
		else if(s0=="rImage" || s0=="rimage") p=&rImage();
		else if(s0=="rObj" || s0=="robj") p=&rObj();
		else if(s0=="s1fix") p=&s1fix;
		else if(s0=="xObjectMax" || s0=="xobjectmax") p=&xObjectMax;
		else if(s0=="yObjectMax" || s0=="yobjectmax") p=&yObjectMax;
		else if(s0=="var") p=&var;
		else if(s0=="var1") p=&var1;
		else if(s0=="var2") p=&var2;
		else if(s0=="var3") p=&var3;
		else if(s0=="var4") p=&var4;
		else if(s0=="var5") p=&var5;

		if(p!=0){
			x=p;
			s=remove_arg(s,0,0);
			args=0;
			return x;
		}
	}

	{   // 1変数プロパティ //////////////////////////////
		double& (cLens1::*fp)(int)=0;

		if     (s0=="a10")              fp=&cLens1::a10;  // 古いコンパイラ(VC6など)だと fp=a10; でも通る
		else if(s0=="a11")              fp=&cLens1::a11;
		else if(s0=="a12")              fp=&cLens1::a12;
		else if(s0=="a13")              fp=&cLens1::a13;
		else if(s0=="a14")              fp=&cLens1::a14;
		else if(s0=="a15")              fp=&cLens1::a15;
		else if(s0=="a16")              fp=&cLens1::a16;
		else if(s0=="a18")              fp=&cLens1::a18;
		else if(s0=="a20")              fp=&cLens1::a20;
		else if(s0=="a1")               fp=&cLens1::a1;
		else if(s0=="a2")               fp=&cLens1::a2;
		else if(s0=="a3")               fp=&cLens1::a3;
		else if(s0=="a4")               fp=&cLens1::a4;
		else if(s0=="a5")               fp=&cLens1::a5;
		else if(s0=="a6")               fp=&cLens1::a6;
		else if(s0=="a7")               fp=&cLens1::a7;
		else if(s0=="a8")               fp=&cLens1::a8;
		else if(s0=="a9")               fp=&cLens1::a9;
		else if(s0=="As0"   || s0=="as0"  ) fp=&cLens1::As0;
		else if(s0=="As45"  || s0=="as45" ) fp=&cLens1::As45;
		else if(s0=="AsTol" || s0=="astol") fp=&cLens1::AsTol;
		else if(s0=="c")                fp=&cLens1::c;
		else if(s0=="caCOA" || s0=="cacoa") fp=&cLens1::caCOA;
		else if(s0=="cbCOA" || s0=="cbcoa") fp=&cLens1::cbCOA;
		else if(s0=="CHMx" || s0=="chmx") fp=&cLens1::CHMx;
		else if(s0=="CHMy" || s0=="chmy") fp=&cLens1::CHMy;
		else if(s0=="CM0" || s0=="cm0") fp=&cLens1::CM0;
		else if(s0=="kp")               fp=&cLens1::kp;
		else if(s0=="d" )               fp=&cLens1::d;
		else if(s0=="DconA4"  || s0=="dcona4"  ) fp=&cLens1::DconA4;
		else if(s0=="DconA6"  || s0=="dcona6"  ) fp=&cLens1::DconA6;
		else if(s0=="DconA8"  || s0=="dcona8"  ) fp=&cLens1::DconA8;
		else if(s0=="DconA10" || s0=="dcona10" ) fp=&cLens1::DconA10;
		else if(s0=="DconA12" || s0=="dcona12" ) fp=&cLens1::DconA12;
		else if(s0=="DconA14" || s0=="dcona14" ) fp=&cLens1::DconA14;
		else if(s0=="DconA16" || s0=="dcona16" ) fp=&cLens1::DconA16;
		else if(s0=="DconA18" || s0=="dcona18" ) fp=&cLens1::DconA18;
		else if(s0=="DconA20" || s0=="dcona20" ) fp=&cLens1::DconA20;
		else if(s0=="DconRn"  || s0=="dconrn"  ) fp=&cLens1::DconRn;
		else if(s0=="dN"      || s0=="dn" )     fp=&cLens1::dN;
		else if(s0=="dx" )              fp=&cLens1::dx;
		else if(s0=="dx1" )             fp=&cLens1::dx1;
		else if(s0=="dy" )              fp=&cLens1::dy;
		else if(s0=="dy1" )             fp=&cLens1::dy1;
		else if(s0=="dz" )              fp=&cLens1::dz;
		else if(s0=="dz1" )             fp=&cLens1::dz1;
		else if(s0=="EAdx" || s0=="eadx") fp=&cLens1::EAdx;
		else if(s0=="EAdy" || s0=="eady") fp=&cLens1::EAdy;
		else if(s0=="EAx" || s0=="eax") fp=&cLens1::EAx;
		else if(s0=="EAy" || s0=="eay") fp=&cLens1::EAy;
		else if(s0=="fideal")           fp=&cLens1::fideal;
		else if(s0=="aCOA" || s0=="acoa") fp=&cLens1::aCOA;
		else if(s0=="bCOA" || s0=="bcoa") fp=&cLens1::bCOA;
		else if(s0=="tCOA" || s0=="tcoa") fp=&cLens1::tCOA;
		else if(s0=="gpitch")           fp=&cLens1::gpitch;
		else if(s0=="grx")              fp=&cLens1::grx;
		else if(s0=="gry")              fp=&cLens1::gry;
		else if(s0=="grz")              fp=&cLens1::grz;
		else if(s0=="kpx")              fp=&cLens1::kpx;
		else if(s0=="kpy")              fp=&cLens1::kpy;
		else if(s0=="Nd" || s0=="nd")   fp=&cLens1::Nd;
		else if(s0=="Newton"    || s0=="newton"    ) fp=&cLens1::Newton;
		else if(s0=="NewtonTol" || s0=="newtontol" ) fp=&cLens1::NewtonTol;
		else if(s0=="NormH" || s0=="normh" ) fp=&cLens1::NormH;
		else if(s0=="Nud" || s0=="nud") fp=&cLens1::Nud;
		else if(s0=="pideal")           fp=&cLens1::pideal;
		else if(s0=="r")                fp=&cLens1::r;
		else if(s0=="rox")              fp=&cLens1::rox;
		else if(s0=="rox1")             fp=&cLens1::rox1;
		else if(s0=="roy")              fp=&cLens1::roy;
		else if(s0=="roy1")             fp=&cLens1::roy1;
		else if(s0=="roz")              fp=&cLens1::roz;
		else if(s0=="roz1")             fp=&cLens1::roz1;
		else if(s0=="rx")               fp=&cLens1::rx;
		else if(s0=="ry")               fp=&cLens1::ry;
		else if(s0=="SA0" || s0=="sa0") fp=&cLens1::SA0;

		if(fp!=0){
			args=1;
			s1=arg(s,1); b1=is_numeric(s1);
			if(b1){
				// fpは double& (cLens::*)(int) 型であり，
				// メンバ関数をfpで呼ぶときはインスタンス(ここでは*this)の指定を省略できない．
				// 演算子 ->* を使う．
				x=&(this->*fp)(atoi(s1.c_str()));
				s=remove_arg(s,0,1);
				return x;
			}
		}
	}

	{   // 2変数プロパティ //////////////////////////////
		double& (cLens1::*fp)(int,int)=0;

		if     (s0=="ZC"  || s0=="zc" ) fp=&cLens1::ZC;
		else if(s0=="SplineH" || s0=="splineh" ) fp=&cLens1::SplineH;
		else if(s0=="SplineZ" || s0=="splinez" ) fp=&cLens1::SplineZ;

		if(fp!=0){
			args=2;
			s1=arg(s,1); b1=is_numeric(s1);
			s2=arg(s,2); b2=is_numeric(s2);
			if(b1 && b2){
				// fpは double& (cLens::*)(int,int) 型であり，
				// メンバ関数をfpで呼ぶときはインスタンス(ここでは*this)の指定を省略できない．
				// 演算子 ->* を使う．
				x=&(this->*fp)(atoi(s1.c_str()),atoi(s2.c_str()));
				s=remove_arg(s,0,2);
				return x;
			}
		}
	}

	{   // 3変数プロパティ //////////////////////////////
		double& (cLens1::*fp)(int,int,int)=0;

		if     (s0=="b")  fp=&cLens1::b;
		else if(s0=="LeC" || s0=="lec") fp=&cLens1::LeC;

		if(fp!=0){
			args=3;
			s1=arg(s,1); b1=is_numeric(s1);
			s2=arg(s,2); b2=is_numeric(s2);
			s3=arg(s,3); b3=is_numeric(s3);
			if(b1 && b2 && b3){
				// fpは double& (cLens::*)(int,int,int) 型であり，
				// メンバ関数をfpで呼ぶときはインスタンス(ここでは*this)の指定を省略できない．
				// 演算子 ->* を使う．
				x=&(this->*fp)(atoi(s1.c_str()),atoi(s2.c_str()),atoi(s3.c_str()));
				s=remove_arg(s,0,3);
				return x;
			}
		}
	}

	return 0;
}

double* cLens1::property(std::string& s){
	int dummy;
	return property(s,dummy);
}