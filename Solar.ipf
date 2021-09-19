#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

////////////////////////////////////////////////////////////////////////////////
/// @file		Solar.ipf
/// @breif		Compute solar altitude orientation, daytime, twilight, and night
/// @author		Takaaki K Abe
/// @date           
/// Version:	1.0
/// Revision:	0
/// @note		parameters are derived from (http://www.es.ris.ac.jp/~nakagawa/met_cal/solar.html)
/// 						
/// @attention      
///
////////////////////////////////////////////////////////////////////////////////

Menu "Astronomy"
	Submenu "Solar"
		"Solar Altitude at one point",/Q, Panel_CalculateSolarAltitude_One()
		"Solar Altitude at moving points",/Q, Panel_CalculateSolarAltitude_Moving()
	End
End


//	================================================================================================
//										Solar Altitude at one point
//	================================================================================================


////////////////////////////////////////////////////////////////////////////////
/// @brief		Create Panel for Calculation of Solar Altitude
/// @param		
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
Function Panel_CalculateSolarAltitude_One() : Panel

	SetDefaultParametersForCalclateSolarAltitude()
		
	variable Width = 250
	variable Height = 200
	variable Vertical = 56
	variable Horizontal = 703

	variable vertical_pos = 0
	variable vertical_spacing = 21
	string PathOfVariables = "root:Astronomic:Solar:Parameters:"
	
	DoWindow SolarAltitudeMenu
	if(V_flag == 1)
		DoWindow/K SolarAltitudeMenu
	endif
	
	PauseUpdate; Silent 1		// building window...
	NewPanel/K=1/W=(Horizontal,Vertical,Horizontal + Width,Vertical + Height )/K=1
	DoWindow/C SolarAltitudeMenu
	
	// title 
	vertical_pos += 9
	TitleBox title_ParamSet, pos={30, vertical_pos},size={116,14},title="\\f01Calculate Solar Altitude"
	TitleBox title_ParamSet, fSize=14,frame=0

	// Profiles
	vertical_pos += 25
	GroupBox group_parameters, frame = 0, title = "\Z11Parameters"
	GroupBox group_parameters, pos = {13, vertical_pos}, size = {220, 5*vertical_Spacing}

	vertical_pos += vertical_spacing
	SetVariable Latitude, title = "Latitude (degree)", pos = {25, vertical_pos}, size = {137, 30}
	SetVariable Latitude, limits = {-90,90,0.1}, value = $( pathOfVariables+"Latitude" )
	
	vertical_pos += vertical_spacing
	SetVariable Longitude, title = "Longitude (degree)", pos = {25, vertical_pos}, size = {137, 30}
	SetVariable Longitude, limits = {-180,180,0.1}, value = $( pathOfVariables+"Longitude" )
	
	vertical_pos += vertical_spacing
	SetVariable TimeDiffFromUTC, title = "Time Difference from UTC time (h)", pos = {25, vertical_pos}, size = {200, 30}
	SetVariable TimeDiffFromUTC, limits = {-12,12,1}, value = $( pathOfVariables+"TimeDiffFromUTC" )
	
	vertical_pos += vertical_spacing
	SetVariable IntervalOfTime, title = "Time-intervel for solar altitude (s)", pos = {25, vertical_pos}, size = {200, 30}
	SetVariable IntervalOfTime, limits = {0,600,1}, value = $( pathOfVariables+"IntervalOfTime" )

	vertical_pos += 1.5*vertical_spacing 
	PopupMenu popup_ReferenceWaveSelect, title="Select reference wave", proc=ReferenceWaveSelect
	PopupMenu popup_ReferenceWaveSelect, pos={27, vertical_pos}, size={98,20}
	PopupMenu popup_ReferenceWaveSelect, value= WaveList("*", ";", ""), popvalue="   -   "
	
	vertical_pos += 1.25*vertical_spacing 
	Button button_StartCalc, title="Calculate Solar Altitude", proc = Button_Launch_SolarAltCalculation
	Button button_StartCalc, pos = {25, vertical_pos}, size={180, 20}

	// Button button_CloseParamWindow, title = "Close", proc = Button_CloseParamWindow
	// Button button_CloseParamWindow, pos = {105, pos}, size = {90,20}
	
	variable PanelHeight = vertical_pos + 25
	
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


////////////////////////////////////////////////////////////////////////////////
/// @brief		name
/// @param		
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
static Function SetDefaultParametersForCalclateSolarAltitude()

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

	//! in the case for calculating moving point
	SVAR/Z LatitudeWaveName
	if(!SVAR_Exists(LatitudeWaveName))
		string/G LatitudeWaveName
	endif

	SVAR/Z LongitudeWaveName
	if(!SVAR_Exists(LongitudeWaveName))
		string/G LongitudeWaveName
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
		
	Duplicate/O Daytime $(saveDF + "Daytime")
	Duplicate/O Twilight $(saveDF + "Twilight")
	Duplicate/O Nighttime $(saveDF + "Nighttime")
	Duplicate/O SolarAltitude $(saveDF + "SolarAltitude")
	
	DisplaySolarAltitude($(saveDF+"SolarAltitude"), $(saveDF+"Daytime"), $(saveDF+"Twilight"))
	
	Killwaves ThetaO, SolarDeclination, EquationOfTime, HourAngle
	Killwaves SolAlt_Diff, HourValue
	Killwaves/Z Daytime, Nighttime, Twilight

	SetDataFolder saveDF

End
///////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
/// @brief          
/// @param[in]      
/// @param[out]     
/// @return                      
///
////////////////////////////////////////////////////////////////////////////////
Function DisplaySolarAltitude(SolarAltitude, Daytime, Twilight)

	wave SolarAltitude, Daytime, Twilight

	DoWindow SolAltView
	if(V_flag)
		DoWindow/F SolAltView
		return 0
		// DoWindow/K SolAltView
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

//	================================================================================================
//										Solar Altitude at moving points
//	================================================================================================


Function Panel_CalculateSolarAltitude_Moving() : Panel

	SetDefaultParametersForCalclateSolarAltitude()

	variable Width = 250
	variable Height = 200
	variable Vertical = 56
	variable Horizontal = 703

	variable vertical_pos = 0
	variable vertical_spacing = 19
	string pathOfVariables = "root:Astronomic:Solar:Parameters:"

	DoWindow SolarAltitudeMenu
	if(V_flag == 1)
		DoWindow/K SolarAltitudeMenu
	endif
		
	PauseUpdate; Silent 1		// building window...
	NewPanel/K=1/W=(Horizontal,Vertical,Horizontal + Width,Vertical + Height )
	DoWindow/C SolarAltitudeMenu
	
	// title 
	vertical_pos += 9
	TitleBox title_ParamSet, pos={30,vertical_pos},size={116,14},title="\\f01Calculate Solar Altitude"
	TitleBox title_ParamSet, fSize=14,frame=0

	// Profiles
	vertical_pos += 25
	GroupBox group_Startparameters, frame = 0, title = "\Z11Waves and Parameters"
	GroupBox group_Startparameters, pos = {13, vertical_pos}, size = {220, 5*vertical_spacing}

	vertical_pos += vertical_spacing
	PopupMenu popup_LatitudeWaveSelect, title="Latitude wave", proc=LatitudeWaveSelect
	PopupMenu popup_LatitudeWaveSelect, pos={27, vertical_pos}, size={98,20}
	PopupMenu popup_LatitudeWaveSelect, value= WaveList("*", ";", ""), popvalue="   -   "
	vertical_pos += vertical_spacing
	PopupMenu popup_LongitudeWaveSelect, title="Longitude wave", proc=LongitudeWaveSelect
	PopupMenu popup_LongitudeWaveSelect, pos={27, vertical_pos}, size={98,20}
	PopupMenu popup_LongitudeWaveSelect, value= WaveList("*", ";", ""), popvalue="   -   "

	vertical_pos += vertical_spacing
	SetVariable TimeDiffFromUTC, title = "Time Difference from UTC time (h)", pos = {25, vertical_pos}, size = {200, 30}
	SetVariable TimeDiffFromUTC, limits = {-12,12,1}, value = $(pathOfVariables+"TimeDiffFromUTC")
	vertical_pos += vertical_spacing
	SetVariable IntervalOfTime, title = "Time-intervel for solar altitude (s)", pos = {25, vertical_pos}, size = {200, 30}
	SetVariable IntervalOfTime, limits = {0,600,1}, value = $(pathOfVariables+"IntervalOfTime")

	vertical_pos += vertical_spacing + 10
	Button button_StartCalc, title="Calculate Solar Altitude", proc = Button_Launch_SolarAltCalculation_Moving
	Button button_StartCalc, pos = {25, vertical_pos}, size={180, 20}

	variable PanelHeight = vertical_pos + 25
	
end

////////////////////////////////////////////////////////////////////////////////
/// @brief		PopupMenuControl of popup_LatitudeWaveSelect in Panel_CalculateSolarAltitude_Moving()
////////////////////////////////////////////////////////////////////////////////
Function 	LatitudeWaveSelect(ctrlName,popNum,popStr) : PopupMenuControl

	String ctrlName
	Variable popNum
	String popStr

	SVAR LatitudeWaveName = root:Astronomic:Solar:Parameters:LatitudeWaveName
	LatitudeWaveName = popStr

End Function

////////////////////////////////////////////////////////////////////////////////
/// @brief		PopupMenuControl of popup_LongitudeWaveSelect in Panel_CalculateSolarAltitude_Moving()
////////////////////////////////////////////////////////////////////////////////
Function 	LongitudeWaveSelect(ctrlName,popNum,popStr) : PopupMenuControl

	String ctrlName
	Variable popNum
	String popStr

	SVAR LongitudeWaveName = root:Astronomic:Solar:Parameters:LongitudeWaveName
	LongitudeWaveName = popStr

End Function

////////////////////////////////////////////////////////////////////////////////
/// @brief		Proc of button_StartCalc in Panel_CalculateSolarAltitude_Moving()
////////////////////////////////////////////////////////////////////////////////
Function Button_Launch_SolarAltCalculation_Moving(ctrlName)
	string ctrlName
	
	string pathofVariables = "root:Astronomic:Solar:Parameters:"


	SVAR Latitude = $( pathOfVariables + "LatitudeWaveName")
	wave LatitudeWave =$Latitude
	SVAR Longitude = $( pathOfVariables + "LongitudeWaveName")
	wave LongitudeWave =$Longitude
	NVAR TimeDiffFromUTC = $( pathOfVariables + "TimeDiffFromUTC")
	NVAR intervalOfTime = $( pathOfVariables + "intervalOfTime")

	if(numpnts(LatitudeWave)!=numpnts(LongitudeWave))
		print "The data length of Latitude is not equal to that of Longitude!"
		print "Please confirm the data length"
		Abort
	endif

	Calculate_SolarAltitude_From_MovingPoint(TimeDiffFromUTC, LatitudeWave, LongitudeWave, intervalOfTime=intervalOfTime)

End Function

////////////////////////////////////////////////////////////////////////////////
/// @brief     ある緯度経度の地点の日の出日の入りの時間を計算する
// -- waves ---
/// @param[in] Latitude (wave) time series of Latitude
/// @param[in] Longitude (wave) time series of Longitude
// -- variables ---
/// @param[in] TimeDiffFromUTC (variable) Time difference between standard time and UTC
/// @param[in] intervalOfTime (variable) time interval for solar altitude wave (second)
/// @return SolarAltitude
/// @return Daytime
/// @return Twilight
/// @return Nighttime
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
	
	Duplicate/O Daytime $(saveDF + "Daytime")
	Duplicate/O Twilight $(saveDF + "Twilight")
	Duplicate/O Nighttime $(saveDF + "Nighttime")
	Duplicate/O SolarAltitude $(saveDF + "SolarAltitude")

	Killwaves ThetaO, SolarDeclination, EquationOfTime, HourAngle
	Killwaves SolAlt_Diff, HourValue, LatitudeForCalc, LongitudeForCalc
	Killwaves/Z Daytime, Nighttime, Twilight
	
	DisplaySolarAltitude($(saveDF+"SolarAltitude"), $(saveDF+"Daytime"), $(saveDF+"Twilight"))

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

