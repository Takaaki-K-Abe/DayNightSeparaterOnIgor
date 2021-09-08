#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

////////////////////////////////////////////////////////////////////////////////
/// @file 			DayNightDivider.ipf
/// @breif 			昼夜を分けるパネルマクロ
/// @author         
/// @date           
/// Version:       $
/// Revision:      $
/// @note           ファイルに備考などを明記する場合はここへ書き込む
/// @attention      ファイルに注意書きなどを明記する場合はここへ書き込む
///
////////////////////////////////////////////////////////////////////////////////

// === Interface ===
// パネル作成
	//PauseUpdate; Silent 1
Macro DayNightDivider()
	NewPanel/N=SelectLocationType/W=(150,50,350,160)
	Variable/G	LocType=0
	CheckBox	check0,	pos={52,25},	size={78,20},	title="Fixed",	value=0,	mode=1,	proc=LocTypeSel
	CheckBox	check1,	pos={52,45},	size={78,20},	title="Moving",	value=0,	mode=1,	proc=LocTypeSel
	Button 	btn1,	pos={52, 65},	title="Select",	proc=StartSunAlt, size={90,20}
EndMacro

// Macro: DayNightDivider()のlocation type を選択する
Function LocTypeSel(name,value)

	String name
	Variable value
	NVAR LocType= root:LocType
	
	strswitch (name)
		case "check0":
			LocType= 0
		break
		case "check1":
			LocType= 1
		break
	endswitch	
	CheckBox check0, value = LocType==0
	CheckBox check1, value = LocType==1
End

////////////////////////////////////////////////////////////////////////////////
/// @brief     DayNightDivider()の関数
///            ポジションが時系列で変わるかどうかで変える
/// @fn        StartSunAlt
/// @param     name
/// @return         
///
////////////////////////////////////////////////////////////////////////////////
Function StartSunAlt(name)
	string name
	NVAR LocType

	If (LocType==0) ///< ロケーションが変わらない時
		Execute "SunAlterFix()"
	elseif (LocType==1) ///< ロケーションが変わる時
		Execute "SunAlterVar()"
	endif
	DoWindow/K SelectLocationType
End

// 計算するロケーションが変わらない時
Proc SunAlterfix(Mold, Year, GMTDiff, BaseLat, BaseLon)
	string Mold
	variable Year, GMTDiff, BaseLat, BaseLon
	Prompt Mold, " Wave for mold (second-scale): ", popup, WaveList("*",";","")
	Prompt Year, " What year?: "
	Prompt GMTDiff, "Time difference from GMT (hour):"
	Prompt BaseLat, " Latitude (°): "
	Prompt BaseLon, " Longitude (°): "

	SunAltitudefix($Mold, Year, GMTDiff, BaseLat, BaseLon)
End

// 計算するロケーションが変わる時
Proc SunAlterVar(Mold, Year, GMTDiff, BaseLat, BaseLon)
    string Mold, BaseLat, BaseLon
    variable Year, GMTDiff
    Prompt Mold, " Wave for mold (second-scale): ", popup, WaveList("*",";","")
    Prompt Year, " What year?: "
    Prompt GMTDiff, "Time difference from GMT (hour):"
    Prompt NameOfLatitudeWave, " Latitude (degree): ", popup, WaveList("*",";","")
    Prompt NameOfLongitudeWave, " Longitude (degree): ", popup, WaveList("*",";","")

    SunAltitudeVar($Mold, Year, GMTDiff, $BaseLat, $BaseLon)
End

////////////////////////////////////////////////////////////////////////////////
/// @brief     ある緯度経度の地点の日の出日の入りの時間を計算する
/// @fn        SunAltitudeFix()
// -- waves ---
/// @param[in] Mold: テンプレートに用いるデータセット
// -- variables ---
/// @param[in] Year: 年
/// @param[in] GMTDiff: GMTからの時差
/// @param[in] BaseLat: 緯度
/// @param[in] BaseLon: 経度
/// @return         
///
////////////////////////////////////////////////////////////////////////////////
Function SunAltitudeFix(Mold, Year, GMTDiff, BaseLat, BaseLon)

	wave     Mold
	variable Year, GMTDiff, BaseLat, BaseLon

	Duplicate/O  Mold Times, ThetaO, SunDec, SunRR, SunEq, SunOri, SunH, SunAlt
	Variable Dates

	Dates=floor(leftx(Mold)/24/60/60-365.25*(Year-1904))
	Times[0]=(leftx(Mold)/24/60/60-floor(leftx(Mold)/24/60/60-1))*24

	Variable i = 0
	for(i = 0; i<numpnts(Mold); i += 1)
		ThetaO[i]=2*pi*(Dates+Times[i]/24-1)/365

		SunDec[i]=0.006918-0.399912*cos(ThetaO[i])+0.070257*sin(ThetaO[i])-0.006758*cos(2*ThetaO[i])+0.000907*sin(2*ThetaO[i])-0.002697*cos(3*ThetaO[i])+0.001480*sin(3*ThetaO[i])
		SunRR[i]=1/(1.000110+0.034221*cos(ThetaO[i])+0.001280*sin(ThetaO[i])+0.000719*cos(2*ThetaO[i])+0.000077*sin(2*ThetaO[i]))^0.5
		SunEq[i]=0.000075+0.001868*cos(ThetaO[i])-0.032077*sin(ThetaO[i])-0.014615*cos(2*ThetaO[i])-0.040849*sin(2*ThetaO[i])

		SunH[i]=(Times[i]-12)*pi/12+(BaseLon-GMTDiff*15)*pi/180+SunEq[i]

		SunAlt[i]=asin(sin(BaseLat*pi/180)*sin(SunDec[i])+cos(BaseLat*pi/180)*cos(SunDec[i])*cos(SunH[i]))
		SunOri[i]=atan(cos(BaseLat*pi/180)*cos(SunDec[i])*sin(SunH[i])/(sin(BaseLat*pi/180)*sin(SunAlt[i])-sin(SunDec[i])))

		// 
		if ( SunAlt[i] - SunAlt[i-1] > 0 && SunOri[i] < 0)
			SunOri[i] += pi
		elseif (SunAlt[i] - SunAlt [i-1] < 0 && SunOri[i] > 0)
			SunOri[i] -= pi
		endif

		// 
		if (SunOri[i]-SunOri[i-1]>pi/2&&SunOri[i]-SunOri[i-1]<3*pi/2)
			SunOri[i]-=pi
		endif
		
		// 
		if (SunOri[i]-SunOri[i-1]<-pi/2&&SunOri[i]-SunOri[i-1]>-3*pi/2)
			SunOri[i]+=pi
		endif 
		
		Times[i+1]=Times[i]+1/60/60
	EndFor
	
	Duplicate/O SunAlt DayNight
		SunAlt *= 180 / pi
		SunOri *= 180 / pi
		DayNight = SunAlt>0 ? 1 : 0
	// Display SunAlt
	// ModifyGraph zColor(SunAlt)={DayNight,*,*,RedWhiteBlue,1}

	// 余分なデータの削除
	Killwaves Times, ThetaO, SunDec, SunRR, SunEq, SunH
 
End

// equation of time

////////////////////////////////////////////
Function SunAltitudeVar(Mold, Year, GMTDiff, BaseLat, BaseLon)

	wave     Mold, BaseLat, BaseLon
	variable Year, GMTDiff

	Duplicate/O  Mold Times, ThetaO, SunDec, SunRR, SunEq, SunOri, SunH, SunAlt
	Make/O/D/N=(numpnts(Mold)) SunDeclination // 赤緯
	Make/O/D/N=(numpnts(Mold)) EquationOfTime // 均時差
	Make/O/D/N=(numpnts(Mold)) HourAngle // 時角
	Make/O/D/N=(numpnts(Mold)) SunAltitude // 太陽高度
	
	Variable Dates
	Dates = floor( leftx(Mold)/24/60/60 - 365.25*( Year - 1904 ) )

	variable startTime = ( leftx(Mold)/24/60/60 - floor( leftx(Mold)/24/60/60 - 1 ) ) * 24
	Times = startTime + x

	variable i
	for(i = 1; i < numpnts(Mold); i += 1)

		// Times[i+1] = Times[i] + 1 / 60 / 60

		ThetaO[i] = 2 * pi * ( Dates + Times[i]/24 -1) / 365

		SunDec[i] = 0.006918-0.399912*cos(ThetaO[i]) + 0.070257*sin(ThetaO[i])
			SunDec[i] += -0.006758*cos(2*ThetaO[i]) + 0.000907*sin(2*ThetaO[i])
			SunDec[i] += -0.002697*cos(3*ThetaO[i]) + 0.001480*sin(3*ThetaO[i])
		
//		SunRR[i] = 1/(1.000110+0.034221*cos(ThetaO[i]) + 0.001280*sin(ThetaO[i])
//			SunRR[i] += 0.000719*cos(2*ThetaO[i]) + 0.000077*sin(2*ThetaO[i]))^0.5

		SunEq[i] = 0.000075+0.001868*cos(ThetaO[i]) - 0.032077*sin(ThetaO[i])
			SunEq[i] += -0.014615*cos(2*ThetaO[i]) - 0.040849*sin(2*ThetaO[i])

		SunH[i] = (Times[i]-12)*pi/12 + (BaseLon[i] - GMTDiff*15)*pi/180 + SunEq[i]

		SunAlt[i] = asin( sin(BaseLat[i] * pi/180) * sin(SunDec[i]) + cos(BaseLat[i]*pi/180) * cos(SunDec[i]) * cos(SunH[i]) )
		SunOri[i] = atan( cos(BaseLat[i] * pi/180) * cos(SunDec[i]) * sin(SunH[i]) / ( sin(BaseLat[i] * pi / 180)*sin(SunAlt[i])-sin(SunDec[i]) ) )

		if(SunAlt[i]-SunAlt[i-1]>0 && SunOri[i]<0)
			SunOri[i] += pi
		elseif (SunAlt[i]-SunAlt[i-1]<0 && SunOri[i]>0)
			SunOri[i] -= pi
		endif

		if (SunOri[i]-SunOri[i-1]>pi/2&&SunOri[i]-SunOri[i-1]<3*pi/2)
			SunOri[i] -= pi
		endif

		if (SunOri[i]-SunOri[i-1]< -pi/2 && SunOri[i]-SunOri[i-1] > -3*pi/2)
			SunOri[i] += pi
		endif

	endfor

 	Duplicate/O SunAlt DayNight
		SunAlt *= 180 / pi
		SunOri *= 180 / pi
		DayNight = SunAlt>0 ? 1 : 0

	DoWindow SunAltitude
	if(V_flag)
		DoWindow/K SunAltitude
	endif

	Display SunAlt
	DoWindow/C SunAltitude

	ModifyGraph zColor(SunAlt)={DayNight,*,*,RedWhiteBlue,1}

	Killwaves Times, ThetaO, SunDec, SunRR, SunEq, SunH

End
///////////////////////////////////////


Function CalculateSolarAltitude(ReferenceWave, Year, TimeDiffFromUTC, Latitude, Longitude)

	wave ReferenceWave
	variable Latitude, Longitude
	variable Year, TimeDiffFromUTC
	
	string saveDF = GetDataFolder(1)
	
	if(!DataFolderExists("root:DayNightSeparater"))
		NewDataFolder root:DayNightSeparater
	endif
	
	SetDataFolder root:DayNightSeparater
	

	Make/O/D/N=(numpnts(ReferenceWave)) TimeValue //
	Make/O/D/N=(numpnts(ReferenceWave)) ThetaO // 赤緯
	Make/O/D/N=(numpnts(ReferenceWave)) SunDeclination // 赤緯
	Make/O/D/N=(numpnts(ReferenceWave)) EquationOfTime // 均時差
	Make/O/D/N=(numpnts(ReferenceWave)) HourAngle // 時角

	Duplicate/O ReferenceWave SolarAltitude
	Duplicate/O ReferenceWave SolarOrientation

	// Make/O/D/N=(numpnts(ReferenceWave)) SolarAltitude // 太陽高度
	// Make/O/D/N=(numpnts(ReferenceWave)) SolarOrientation // 太陽高度
	
	Variable JulianDate
	JulianDate = floor( leftx(ReferenceWave)/24/60/60 - 365.25*( Year - 1904 ) )

	variable startTime = ( leftx(ReferenceWave)/24/60/60 - floor( leftx(ReferenceWave)/24/60/60 - 1 ) ) * 24
	TimeValue = startTime + x/60/60

	ThetaO = 2 * pi * ( JulianDate + TimeValue/24 -1) / 365

	SunDeclination = 0.006918-0.399912*cos(ThetaO) + 0.070257*sin(ThetaO) - 0.006758*cos(2*ThetaO)
		SunDeclination += 0.000907*sin(2*ThetaO) -0.002697*cos(3*ThetaO) + 0.001480*sin(3*ThetaO)
	
	EquationOfTime = 0.000075+0.001868*cos(ThetaO) - 0.032077*sin(ThetaO)
		EquationOfTime += -0.014615*cos(2*ThetaO) - 0.040849*sin(2*ThetaO)

	HourAngle = (TimeValue-12)*pi/12 + (Longitude - TimeDiffFromUTC*15)*pi/180 + EquationOfTime

	SolarAltitude = asin( sin(Latitude * pi/180) * sin(SunDeclination) + cos(Latitude*pi/180) * cos(SunDeclination) * cos(HourAngle) )
	SolarOrientation = atan( cos(Latitude * pi/180) * cos(SunDeclination) * sin(HourAngle) / ( sin(Latitude * pi / 180)*sin(SolarAltitude)-sin(SunDeclination) ) )
	
	duplicate SolarOrientaition SolOlnotIF

	variable i
	for(i = 1; i<numpnts(ReferenceWave); i += 1)

		if(SolarAltitude[i]-SolarAltitude[i-1]>0 && SolarOrientation[i]<0)
			SolarOrientation[i] += pi
		elseif (SolarAltitude[i]-SolarAltitude[i-1]<0 && SolarOrientation[i]>0)
			SolarOrientation[i] -= pi
		endif

		if (SolarOrientation[i]-SolarOrientation[i-1] > pi/2 && SolarOrientation[i]-SolarOrientation[i-1] < 3*pi/2)
			SolarOrientation[i] -= pi
		endif
		if (SolarOrientation[i]-SolarOrientation[i-1]< -pi/2 && SolarOrientation[i]-SolarOrientation[i-1] > -3*pi/2)
			SolarOrientation[i] += pi
		endif

	endfor

	SolarAltitude *= 180 / pi
	SolarOrientation *= 180 / pi
 	
	Duplicate/O SolarAltitude DayNight
		DayNight = SolarAltitude>0 ? 1 : 0
		note DayNight, "Type:Mask;"

	DoWindow SolAltView
	if(V_flag)
		DoWindow/K SolAltView
	endif
	Display SolarAltitude
	DoWindow/C SolAltView

	ModifyGraph zColor(SolarAltitude)={DayNight,*,*,RedWhiteBlue,1}

	Killwaves TimeValue, ThetaO, SunDeclination, EquationOfTime, HourAngle
	
	Duplicate/O DayNight $(saveDF + "DayNight")
	
	SetDataFolder saveDF

End
///////////////////////////////////////



