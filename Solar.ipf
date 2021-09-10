#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

////////////////////////////////////////////////////////////////////////////////
/// @file 			Solar.ipf
/// @breif 			Compute solar altitude orientation, daytime, twilight, and night
/// @author         Takaaki K Abe
/// @date           
/// Version:       	1.0
/// Revision:      	0
/// @note           ファイルに備考などを明記する場合はここへ書き込む
/// @attention      ファイルに注意書きなどを明記する場合はここへ書き込む
///
////////////////////////////////////////////////////////////////////////////////

Menu "Astronomy"
	Submenu "Solar"
		"Solar Altitude at one point",/Q, Panel_CalculateSolarAltitude()
		"Solar Altitude at multi points",/Q, Panel_CalculateSolarAltitude()
	End
End


Function Panel_CalculateSolarAltitude() : Panel

	SetDefaultParametersForCalclateSolarAltitude()
	
	variable pos = 0
	
	variable Width, Height, Vertical, Horizontal	
	Horizontal = 703; Vertical = 56
	Width = 330; Height = 365
	
	DoWindow SolarAltitudeMenu
	if(V_flag == 1)
		DoWindow/K SolarAltitudeMenu
	endif
	
	PauseUpdate; Silent 1		// building window...
	NewPanel/K=1/W=(Horizontal,Vertical,Horizontal + Width,Vertical + Height )/K=1
	DoWindow/C SolarAltitudeMenu
	
	// title 
	pos += 9
	TitleBox title_ParamSet, pos={100,pos},size={116,14},title="\\f01Calculate Solar Altitude"
	TitleBox title_ParamSet, fSize=14,frame=0

	// Profiles
	pos += 25
	GroupBox group_parameters, frame = 0, title = "\Z11Parameters for Calculation"
	GroupBox group_parameters, pos = {13, pos}, size = {310, 130}

	SetVariable Latitude, title = "Latitude (degree)", pos = {25, pos + 21}, size = {137, 30}
	SetVariable Latitude, limits = {-90,90,0.1}, value = root:Astronomic:Solar:Parameters:Latitude
	SetVariable Longitude, title = "Longitude (degree)", pos = {173, pos + 21}, size = {137, 30}
	SetVariable Longitude, limits = {-180,180,0.1}, value = root:Astronomic:Solar:Parameters:Longitude
	SetVariable TimeDiffFromUTC, title = "Time Difference from UTC time (h)", pos = {25, pos + 39}, size = {200, 30}
	SetVariable TimeDiffFromUTC, limits = {-12,12,1}, value = root:Astronomic:Solar:Parameters:TimeDiffFromUTC
	SetVariable IntervalOfTime, title = "Time-intervel for solar altitude (s)", pos = {25, pos + 57}, size = {200, 30}
	SetVariable IntervalOfTime, limits = {0,600,1}, value = root:Astronomic:Solar:Parameters:IntervalOfTime

	PopupMenu popup_ReferenceWaveSelect, title="Select reference wave", proc=ReferenceWaveSelect
	PopupMenu popup_ReferenceWaveSelect, pos={27, pos + 75}, size={98,20}
	PopupMenu popup_ReferenceWaveSelect, value= WaveList("*", ";", ""), popvalue="   -   "
	
	Button button_StartCalc, title="Calculate Solar Altitude", proc = Button_Launch_SolarAltCalculation
	Button button_StartCalc, pos = {25, pos + 100}, size={180, 20}

	// Button button_CloseParamWindow, title = "Close", proc = Button_CloseParamWindow
	// Button button_CloseParamWindow, pos = {105, pos}, size = {90,20}
	
	variable PanelHeight = pos + 25
	
end

////////////////////////////////////////////////////////////////////////////////
/// @brief		S
/// @param		
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
Function Button_Launch_SolarAltCalculation(ctrlName)
	string ctrlName
	
	SVAR RefWaveName = root:Astronomic:Solar:Parameters:ReferenceWaveName
	wave ReferenceWave = $RefWaveName

	NVAR Latitude = root:Astronomic:Solar:Parameters:Latitude
	NVAR Longitude = root:Astronomic:Solar:Parameters:Longitude
	NVAR TimeDiffFromUTC = root:Astronomic:Solar:Parameters:TimeDiffFromUTC
	NVAR intervalOfTime = root:Astronomic:Solar:Parameters:intervalOfTime

	Calculate_SolarAltitude_From_Coordinate(ReferenceWave, TimeDiffFromUTC, Latitude, Longitude, intervalOfTime=intervalOfTime)

End Function

////////////////////////////////////////////////////////////////////////////////
/// @brief		name
/// @param		
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
Function 	ReferenceWaveSelect(ctrlName,popNum,popStr) : PopupMenuControl

	String ctrlName
	Variable popNum
	String popStr

	SVAR ReferenceWaveName = root:Astronomic:Solar:Parameters:ReferenceWaveName
	ReferenceWaveName = popStr

End Function


STRUCTURE SolarAltitude

	wave ReferenceWave
	variable TimeDiffFromUTC

EndSTRUCTURE

Function SetDefaultParametersForCalclateSolarAltitude()

	string saveDF = GetDataFolder(1)
	
	if(!DataFolderExists("root:Astronomic"))
		NewDataFolder root:Astronomic
	endif

	if(!DataFolderExists("root:Astronomic:Solar"))
		NewDataFolder root:Astronomic:Solar
	endif
	
	if(!DataFolderExists("root:Astronomic:Solar:Parameters"))
		NewDataFolder root:Astronomic:Solar:Parameters
	endif

	SetDataFolder root:Astronomic:Solar:Parameters

	SVAR/Z ReferenceWaveName
	if(!SVAR_Exists(ReferenceWaveName))
		string/G ReferenceWaveName
	endif

	//! in the case for calculating one point
	NVAR/Z Latitude
	if(!NVAR_Exists(Latitude))
		variable/G Latitude = 35
	endif
	
	NVAR/Z Longitude
	if(!NVAR_Exists(Longitude))
		variable/G Longitude = 135
	endif

	//! the time difference in standard time from UTC time
	NVAR/Z TimeDiffFromUTC
	if(!NVAR_Exists(TimeDiffFromUTC))
		variable/G TimeDiffFromUTC = 9
	endif

	//! the time interval for calculating solar altitude (seconds)
	NVAR/Z IntervalOfTime
	if(!NVAR_Exists(IntervalOfTime))
		variable/G IntervalOfTime = 60
	endif

	SetDataFolder saveDF

End Function

////////////////////////////////////////////////////////////////////////////////
/// @brief		
/// @param		
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
Function StartCalculationOfSolarAltitude(RefWaveName)

	string RefWaveName

	wave ReferenceWave = $RefWaveName
	NVAR Latitude = root:Astronomic:Solar:Parameters:Latitude
	NVAR Longitude = root:Astronomic:Solar:Parameters:Longitude
	NVAR TimeDiffFromUTC = root:Astronomic:Solar:Parameters:TimeDiffFromUTC
	NVAR intervalOfTime = root:Astronomic:Solar:Parameters:intervalOfTime

	Calculate_SolarAltitude_From_Coordinate(ReferenceWave, TimeDiffFromUTC, Latitude, Longitude, intervalOfTime=intervalOfTime)

End Function

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
Function Calculate_SolarAltitude_From_Coordinate(ReferenceWave, TimeDiffFromUTC, Latitude, Longitude, [intervalOfTime])

	wave ReferenceWave
	variable Latitude, Longitude
	variable TimeDiffFromUTC
	variable intervalOfTime
	
	if(ParamIsDefault(intervalOfTime))
		intervalOfTime = 1
	endif
	
	string saveDF = GetDataFolder(1)
	
	// move to working directory
	if(!DataFolderExists("root:Astronomic:Solar"))
		NewDataFolder root:Astronomic:Solar
	endif
	
	SetDataFolder root:Astronomic:Solar
	
	variable year, month, day
	variable startDateTime = leftx(ReferenceWave)
	string dateString = secs2Date(startDateTime, -2, "/")
	sscanf dateString, "%d/%d/%d", year, month, day

	variable NumberOfWavePoints = ceil(numpnts(ReferenceWave)*deltax(ReferenceWave)/intervalOfTime)

	// preperation for calculation
	Make/O/D/N=(NumberOfWavePoints) HourValue // 通し時間
	Make/O/D/N=(NumberOfWavePoints) ThetaO // 
	Make/O/D/N=(NumberOfWavePoints) SolarDeclination // 赤緯
	Make/O/D/N=(NumberOfWavePoints) EquationOfTime // 均時差
	Make/O/D/N=(NumberOfWavePoints) HourAngle // 時角

	Make/O/D/N=(NumberOfWavePoints) SolarAltitude
		SetScale/P x startDateTime, intervalOfTime, "dat", SolarAltitude
	Make/O/D/N=(NumberOfWavePoints) SolarOrientation
		SetScale/P x startDateTime, intervalOfTime, "dat", SolarOrientation

	// Make/O/D/N=(numpnts(ReferenceWave)) SolarAltitude // 太陽高度
	// Make/O/D/N=(numpnts(ReferenceWave)) SolarOrientation // 太陽高度
	
	// Compute the parameters for solar altitude
	Variable DateValue
	DateValue = floor( leftx(ReferenceWave)/24/60/60 - 365.25*( Year - 1904 ) )

	variable startTime = ( leftx(ReferenceWave)/24/60/60 - floor( leftx(ReferenceWave)/24/60/60 - 1 ) ) * 24
	HourValue = startTime + x/60/60*intervalOfTime

	ThetaO = 2 * pi * ( DateValue + HourValue/24 -1) / 365

	SolarDeclination =  0.006918-0.399912*cos(ThetaO) + 0.070257*sin(ThetaO) - 0.006758*cos(2*ThetaO)
	SolarDeclination += 0.000907*sin(2*ThetaO) -0.002697*cos(3*ThetaO) + 0.001480*sin(3*ThetaO)
	
	EquationOfTime =  0.000075+0.001868*cos(ThetaO) - 0.032077*sin(ThetaO)
	EquationOfTime += -0.014615*cos(2*ThetaO) - 0.040849*sin(2*ThetaO)

	HourAngle = (HourValue-12)*pi/12 + ( Longitude - TimeDiffFromUTC * 15) * pi/180
	HourAngle += EquationOfTime
	// calculate solar altitude
	SolarAltitude =  sin(Latitude * pi/180) * sin(SolarDeclination) 
	SolarAltitude += cos(Latitude*pi/180) * cos(SolarDeclination) * cos(HourAngle)
	SolarAltitude =  asin( SolarAltitude )
	
	// calculate solar orientation
	SolarOrientation =  cos(Latitude * pi/180) * cos(SolarDeclination) * sin(HourAngle)
	SolarOrientation /= (sin(Latitude * pi / 180) * sin(SolarAltitude) - sin(SolarDeclination))
	SolarOrientation =  atan( SolarOrientation )
	Duplicate/D/O SolarAltitude SolAlt_Diff
	Differentiate/Meth=1 SolAlt_Diff

	SolarOrientation += SolAlt_Diff > 0 && SolarOrientation < 0 ? pi : 0
	SolarOrientation -= SolAlt_Diff < 0 && SolarOrientation > 0 ? pi : 0

	variable i
	for(i = 1; i<NumberOfWavePoints; i += 1)

		if (SolarOrientation[i]-SolarOrientation[i-1] > pi/2 && SolarOrientation[i]-SolarOrientation[i-1] < 3*pi/2)
			SolarOrientation[i] -= pi
		endif
		if (SolarOrientation[i]-SolarOrientation[i-1]< -pi/2 && SolarOrientation[i]-SolarOrientation[i-1] > -3*pi/2)
			SolarOrientation[i] += pi
		endif

	endfor

	SolarAltitude *= 180 / pi
	SolarOrientation *= 180 / pi

	// Separate daytime, twilight nighttime
	variable Solar_Threshold = -1/60*50
	Duplicate/O SolarAltitude Daytime
		Daytime = SolarAltitude > Solar_Threshold ? 1 : 0
		note Daytime, "Type:Mask; \n 1: Daytime \n 0: Nighttime"

	Duplicate/O SolarAltitude Twilight
		Twilight = SolarAltitude <= Solar_Threshold && SolarAltitude > -18 ? 1 : 0
		note Twilight, "Type:Mask; \n 1: Twilight \n 0: Nighttime"
	
	Duplicate/O SolarAltitude Nighttime
		Nighttime = SolarAltitude <= -18 ? 1 : 0
		note Nighttime, "Type:Mask; \n 1: Night"

	DisplaySolarAltitude(SolarAltitude, Daytime, Twilight)
		
	Duplicate/O Daytime $(saveDF + "Daytime")
	Duplicate/O Twilight $(saveDF + "Twilight")
	Duplicate/O Nighttime $(saveDF + "Nighttime")
	Duplicate/O SolarAltitude $(saveDF + "SolarAltitude")
	
	Killwaves ThetaO, SolarDeclination, EquationOfTime, HourAngle
	Killwaves SolAlt_Diff, HourValue

	SetDataFolder saveDF

End
///////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// @brief          
/// @param[in]      
/// @param[out]     
/// @return         
/// @author         
/// @date           
/// @version        
/// @note           
/// @attention      
/// @par            
///                 
///
////////////////////////////////////////////////////////////////////////////////
Function DisplaySolarAltitude(SolarAltitude, Daytime, Twilight)

	wave SolarAltitude, Daytime, Twilight
	// string PathTwilight = GetWavesDataFolder(Twilight, 2)
	// string PathDaytime = GetWavesDataFolder(Daytime, 2)
	// string PathSolarAltitude = GetWavesDataFolder(SolarAltitude, 2)

	DoWindow SolAltView
	if(V_flag)
		DoWindow/K SolAltView
	endif
	Display/K=1/R Twilight
	AppendtoGraph SolarAltitude
	DoWindow/C SolAltView

	// ModifyGraph zColor(SolarAltitude)={DayNight,*,*,RedWhiteBlue,1}
	ModifyGraph zColor(SolarAltitude)={Daytime,*,*,BlueGreenOrange,0}
	ModifyGraph manTick(left)={0,30,0,0},manMinor(left)={2,0}
	ModifyGraph mode(Twilight)=7, rgb(Twilight)=(65535,54607,32768)
	ModifyGraph hbFill(Twilight)=4, lsize(Twilight)=0
	ModifyGraph tick(right)=3,noLabel(right)=2,axRGB(right)=(65535,65535,65535)
	Label left "Solor altitude (°)"
	Label bottom "Date and Time"
	SetAxis left -90,90

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
Function Calculate_SolarAltitude_From_MovingPoint(TimeDiffFromUTC, Latitude, Longitude, [intervalOfTime])

	wave Latitude, Longitude
	variable TimeDiffFromUTC
	variable intervalOfTime
	
	wave ReferenceWave = Latitude

	if(numpnts(Latitude) != numpnts(Longitude))
		print "The number of points in latitude and longitude waves is not equal."
		Abort
	endif
	
	if(ParamIsDefault(intervalOfTime))
		intervalOfTime = 1
	endif
	
	string saveDF = GetDataFolder(1)
	
	// move to working directory
	if(!DataFolderExists("root:Astronomic:Solar"))
		NewDataFolder root:Astronomic:Solar
	endif
	
	SetDataFolder root:Astronomic:Solar
	
	// read start year
	variable year, month, day
	variable startDateTime = leftx(ReferenceWave)
	string dateString = secs2Date(startDateTime, -2, "/")
	sscanf dateString, "%d/%d/%d", year, month, day

	variable NumberOfWavePoints = ceil(numpnts(ReferenceWave)*deltax(ReferenceWave)/intervalOfTime)

	// preperation for calculation
	Make/O/D/N=(NumberOfWavePoints) HourValue // 通し時間
	Make/O/D/N=(NumberOfWavePoints) ThetaO // 
	Make/O/D/N=(NumberOfWavePoints) SolarDeclination // 赤緯
	Make/O/D/N=(NumberOfWavePoints) EquationOfTime // 均時差
	Make/O/D/N=(NumberOfWavePoints) HourAngle // 時角
	
	Make/O/D/N=(NumberOfWavePoints) LatitudeForCalc // 時角
		// SetScale/P x startDateTime, intervalOfTime, "dat", LatitudeForCalc
		LatitudeForCalc = Latitude[ p/(deltax(Latitude)/intervalOfTime) ]
	Make/O/D/N=(NumberOfWavePoints) LongitudeForCalc // 時角
		// SetScale/P x startDateTime, intervalOfTime, "dat", LongitudeForCalc
		LongitudeForCalc = Longitude[ p/(deltax(Longitude)/intervalOfTime) ]

	Make/O/D/N=(NumberOfWavePoints) SolarAltitude
		SetScale/P x startDateTime, intervalOfTime, "dat", SolarAltitude
	Make/O/D/N=(NumberOfWavePoints) SolarOrientation
		SetScale/P x startDateTime, intervalOfTime, "dat", SolarOrientation
	
	// Compute the parameters for solar altitude
	Variable DateValue
	DateValue = floor( leftx(ReferenceWave)/24/60/60 - 365.25*( Year - 1904 ) )

	variable startTime = ( leftx(ReferenceWave)/24/60/60 - floor( leftx(ReferenceWave)/24/60/60 - 1 ) ) * 24
	HourValue = startTime + x/60/60*intervalOfTime

	ThetaO = 2 * pi * ( DateValue + HourValue/24 -1) / 365

	SolarDeclination =  0.006918-0.399912*cos(ThetaO) + 0.070257*sin(ThetaO) - 0.006758*cos(2*ThetaO)
	SolarDeclination += 0.000907*sin(2*ThetaO) -0.002697*cos(3*ThetaO) + 0.001480*sin(3*ThetaO)
	
	EquationOfTime = 0.000075+0.001868*cos(ThetaO) - 0.032077*sin(ThetaO)
	EquationOfTime += -0.014615*cos(2*ThetaO) - 0.040849*sin(2*ThetaO)

	HourAngle = (HourValue-12)*pi/12 + (LongitudeForCalc[p] - TimeDiffFromUTC*15)*pi/180 + EquationOfTime

	// calculate solar altitude
	SolarAltitude = asin( sin(LatitudeForCalc * pi/180) * sin(SolarDeclination) + cos(LatitudeForCalc*pi/180) * cos(SolarDeclination) * cos(HourAngle) )

	// calculate solar orientation
	SolarOrientation = atan( cos(LatitudeForCalc * pi/180) * cos(SolarDeclination) * sin(HourAngle) / ( sin(LatitudeForCalc * pi / 180)*sin(SolarAltitude)-sin(SolarDeclination) ) )
	
	Duplicate/D/O SolarAltitude SolAlt_Diff
	Differentiate/Meth=1 SolAlt_Diff

	SolarOrientation += SolAlt_Diff > 0 && SolarOrientation < 0 ? pi : 0
	SolarOrientation -= SolAlt_Diff < 0 && SolarOrientation > 0 ? pi : 0

	variable i
	for(i = 1; i<NumberOfWavePoints; i += 1)

		if (SolarOrientation[i]-SolarOrientation[i-1] > pi/2 && SolarOrientation[i]-SolarOrientation[i-1] < 3*pi/2)
			SolarOrientation[i] -= pi
		endif
		if (SolarOrientation[i]-SolarOrientation[i-1]< -pi/2 && SolarOrientation[i]-SolarOrientation[i-1] > -3*pi/2)
			SolarOrientation[i] += pi
		endif

	endfor

	SolarAltitude *= 180 / pi
	SolarOrientation *= 180 / pi

	variable Solar_Threshold = -1/60*50
	Duplicate/O SolarAltitude Daytime
		Daytime = SolarAltitude > Solar_Threshold ? 1 : 0
		note Daytime, "Type:Mask; \n 1: Daytime \n 0: Nighttime"
	
	Duplicate/O SolarAltitude Twilight
		Twilight = SolarAltitude <= Solar_Threshold && SolarAltitude > -18 ? 1 : 0
		note Twilight, "Type:Mask; \n 1: Twilight \n 0: Nighttime"
	
	Duplicate/O SolarAltitude Nighttime
		Nighttime = SolarAltitude <= -18 ? 1 : 0
		note Nighttime, "Type:Mask; \n 1: Night"

	DisplaySolarAltitude(SolarAltitude, Daytime, Twilight)
	
	Duplicate/O Daytime $(saveDF + "Daytime")
	Duplicate/O SolarAltitude $(saveDF + "SolarAltitude")

	Killwaves ThetaO, SolarDeclination, EquationOfTime, HourAngle
	Killwaves SolAlt_Diff, HourValue, LatitudeForCalc, LongitudeForCalc

	SetDataFolder saveDF

End
///////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// @brief		ReturnTwilights_and_Nights
/// @param		
/// @param		type (variable) 1:civil, 2:nautical, 3:astronomical twilights, 4: midnights
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
Function ReturnTwilights_and_Nights(SolarAltitude, type)

	wave SolarAltitude
	variable type
	string NameOfOutput
	variable minSolAlt, maxSolAlt

	variable Solar_Threshold = 1/60*50
	switch(type)
		case 0:
			NameOfOutput = "Twilight"
			maxSolAlt = Solar_Threshold
			minSolAlt = -18
			break
		case 1:
			NameOfOutput = "Twilight_civil"
			maxSolAlt = Solar_Threshold
			minSolAlt = -6
			break
		case 2:
			NameOfOutput = "Twilight_nautical"
			maxSolAlt = -6
			minSolAlt = -12
			break
		case 3:
			NameOfOutput = "Twilights_astro"
			maxSolAlt = -12
			minSolAlt = -18
			break
		case 4:
			NameOfOutput = "MidNight"
			maxSolAlt = -18
			minSolAlt = -120
			break
	endswitch

	Duplicate/O SolarAltitude $NameOfOutput
	wave outputWave = $NameOfOutput
	outputWave = SolarAltitude <= maxSolAlt && SolarAltitude > minSolAlt ? 1 : 0

End Function
