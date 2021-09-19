#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Menu "Astronomy"
	Submenu "miscs"
		"Create track between two points",/Q, LanchCreateTimeSeriesPath()
	End
End

Function LanchCreateTimeSeriesPath()

	SetDefaultParameters()
	Panel_ComputeTrack()

End

Function Panel_ComputeTrack() : Panel

	variable pos = 0
	variable spacing = 19
	
	variable Width = 250
	variable Height = 345

	variable Vertical = 56
	variable Horizontal = 100

	DoWindow ComputeTrack
	if(V_flag == 1)
		DoWindow/K ComputeTrack
	endif
		
	PauseUpdate; Silent 1		// building window...
	NewPanel/K=1/W=(Horizontal,Vertical,Horizontal + Width,Vertical + Height )
	DoWindow/C ComputeTrack
	
	// title 
	pos += 9
	TitleBox title_ParamSet, pos={50,pos},size={116,14},title="\\f01Compute line track"
	TitleBox title_ParamSet, fSize=14,frame=0

	// Profiles
	pos += 25
	GroupBox group_Startparameters, frame = 0, title = "\Z11Coordinates at start point"
	GroupBox group_Startparameters, pos = {13, pos}, size = {220, 5*spacing}
	pos += spacing
		SetVariable startLatitude, title = "Latitude (degree)", pos = {25, pos}, size = {182, 30}
		SetVariable startLatitude, limits = {-90,90,0.1}, value = root:Astronomic:miscs:track:startLatitude
	pos += spacing
		SetVariable startLongitude, title = "Longitude (degree)", pos = {25, pos}, size = {182, 30}
		SetVariable startLongitude, limits = {-180,180,0.1}, value = root:Astronomic:miscs:track:startLongitude
	pos += spacing
		SetVariable startDateString, title = "Date (yyyy/mm/dd)"
		SetVariable startDateString, pos = {25, pos}, size = {182, 30}
		SetVariable startDateString, value = root:Astronomic:miscs:track:startDateString
	pos += spacing
		SetVariable startTimeString, pos = {25, pos}, size = {182, 30}
		SetVariable startTimeString, title = "Time (hh:mm:ss)"
		SetVariable startTimeString, value = root:Astronomic:miscs:track:startTimeString

	pos += 25
	GroupBox group_Endparameters, frame = 0, title = "\Z11Coordinates at end point"
	GroupBox group_Endparameters, pos = {13, pos}, size = {220, 5*spacing}

	pos += spacing
		SetVariable endLatitude, title = "Latitude (degree)", pos = {25, pos}, size = {182, 30}
		SetVariable endLatitude, limits = {-90,90,0.1}, value = root:Astronomic:miscs:track:endLatitude
	pos += spacing
		SetVariable endLongitude, title = "Longitude (degree)", pos = {25, pos}, size = {182, 30}
		SetVariable endLongitude, limits = {-180,180,0.1}, value = root:Astronomic:miscs:track:endLongitude

	pos += spacing
		SetVariable endDateString, title = "Date (yyyy/mm/dd)"
		SetVariable endDateString, pos = {25, pos}, size = {182, 30}
		SetVariable endDateString, value = root:Astronomic:miscs:track:endDateString
	pos += spacing
		SetVariable endTimeString, pos = {25, pos}, size = {182, 30}
		SetVariable endTimeString, title = "Time (hh:mm:ss)"
		SetVariable endTimeString, value = root:Astronomic:miscs:track:endTimeString
	
	pos += 30
		GroupBox group_WaveCalc, frame = 0, title = "\Z11Create Lat and Lon waves"
		GroupBox group_WaveCalc, pos = {13, pos}, size = {220, 4.5*spacing}
	
	pos += spacing
		SetVariable intervalOfTime, title = "Time-intervel for solar altitude (s)"
		SetVariable intervalOfTime, limits = {0,3600,1}, value = root:Astronomic:miscs:track:intervalOfTime
		SetVariable intervalOfTime, pos = {25, pos}, size = {200, 30}
	pos += spacing
		SetVariable waveNote, pos = {25, pos}, size = {182, 30}
		SetVariable waveNote, title = "Comment"
		SetVariable waveNote, value = $"root:Astronomic:miscs:track:waveNote"

	pos += spacing
		Button button_StartCalc, title="Compute", proc = Button_Launch_ComputePath
		Button button_StartCalc, pos = {25, pos}, size={100, 20}

	// PopupMenu popup_ReferenceWaveSelect, title="Select reference wave", proc=ReferenceWaveSelect
	// PopupMenu popup_ReferenceWaveSelect, pos={27, pos + 75}, size={98,20}
	// PopupMenu popup_ReferenceWaveSelect, value= WaveList("*", ";", ""), popvalue="   -   "
	

	// // Button button_CloseParamWindow, title = "Close", proc = Button_CloseParamWindow
	// // Button button_CloseParamWindow, pos = {105, pos}, size = {90,20}
	
	print pos
	
	variable PanelHeight = pos + 25
	
end

////////////////////////////////////////////////////////////////////////////////
/// @brief		S
/// @param		
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
Function Button_Launch_ComputePath(ctrlName)
	string ctrlName

	CreateTimeSeriesPath()

End Function


////////////////////////////////////////////////////////////////////////////////
/// @brief		
/// @param		
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
Function CreateTimeSeriesPath()

	string saveDF = GetDataFolder(1)

	string pathToParameter = "root:Astronomic:miscs:track:"
	
	SetDataFolder pathToParameter

	SVAR startDateString = $(pathToParameter + "startDateString")
	SVAR startTimeString = $(pathToParameter + "startTimeString")
	SVAR endDateString = $(pathToParameter + "endDateString")
	SVAR endTimeString = $(pathToParameter + "endTimeString")
	SVAR waveNote = $(pathToParameter + "waveNote")
	
	NVAR intervalOfTime = $(pathToParameter + "intervalOfTime")

	NVAR startLatitude = $(pathToParameter + "startLatitude")
	NVAR startLongitude = $(pathToParameter + "startLongitude")
	NVAR endLatitude = $(pathToParameter + "endLatitude")
	NVAR endLongitude = $(pathToParameter + "endLongitude")

	variable startDateTime = ReadDateTimeString(startDateString, startTimeString)
	variable endDateTime = ReadDateTimeString(endDateString, endTimeString)

	if(startDateTime > endDateTime)
		print "The end date and time is earlier than start date and time!"
		Abort
	endif
	
	variable numberOfPoints = ceil((endDateTime - startDateTime)/intervalOfTime) + 1

	variable LatDelta = (endLatitude - startLatitude)/numberOfPoints
	variable LonDelta = (endLongitude - startLongitude)/numberOfPoints

	make/O/N=(numberOfPoints) Latitude
	SetScale/P x startDateTime,intervalOfTime,"dat", Latitude
	Latitude = startLatitude + p*LatDelta
	note Latitude waveNote

	make/O/N=(numberOfPoints) Longitude
	SetScale/P x startDateTime,intervalOfTime,"dat", Longitude
	Longitude = startLongitude + p*LonDelta
	note Longitude waveNote

	Duplicate/O Latitude $(saveDF + NameOfWave(Latitude))
	Duplicate/O Longitude $(saveDF + NameOfWave(Longitude))

	killwaves Latitude, Longitude

	SetDataFolder saveDF

End

Function ReadDateTimeString(dateString, timeString)

	string dateString
	string timeString
	
	variable year, month, day
		sscanf dateString, "%d/%d/%d", year, month, day
	variable dateValue = date2secs(year, month, day)

	variable hour, minute, second
		sscanf timeString, "%d:%d:%d", hour, minute, second
	variable timevalue = hour*60*60 + minute*60 + second

	variable dateTimeValue = dateValue + timevalue
	return dateTimeValue
End

////////////////////////////////////////////////////////////////////////////////
/// @brief		
/// @param		
/// @return	
///
////////////////////////////////////////////////////////////////////////////////
static Function SetDefaultParameters()

	string saveDF = GetDataFolder(1)

	string pathToBase = "root:Astronomic"
	string pathToMisc = pathToBase + ":miscs"
	string pathToTrack = pathToMisc + ":track"
	if(!DataFolderExists(pathToBase))
		NewDataFolder $pathToBase
	endif
	if(!DataFolderExists(pathToMisc))
		NewDataFolder $pathToMisc
	endif
	if(!DataFolderExists(pathToTrack))
		NewDataFolder $pathToTrack
	endif	

	SetDataFolder $pathToTrack

	SVAR/Z startDateString
	if(!SVAR_Exists(startDateString))
		string/G startDateString = "2020/1/1"
	endif
	SVAR/Z startTimeString
	if(!SVAR_Exists(startTimeString))
		string/G startTimeString = "00:00:00"
	endif
	SVAR/Z endDateString
	if(!SVAR_Exists(endDateString))
		string/G endDateString = "2020/2/1"
	endif
	SVAR/Z endTimeString
	if(!SVAR_Exists(endTimeString))
		string/G endTimeString = "00:00:00"
	endif
	SVAR/Z waveNote
	if(!SVAR_Exists(waveNote))
		string/G waveNote = ""
	endif


	NVAR/Z intervalOfTime
	if(!NVAR_Exists(intervalOfTime))
		variable/G intervalOfTime = 35
	endif
	
	NVAR/Z startLatitude
	if(!NVAR_Exists(startLatitude))
		variable/G startLatitude = 35
	endif
	
	NVAR/Z startLongitude
	if(!NVAR_Exists(startLongitude))
		variable/G startLongitude = 135
	endif

	NVAR/Z endLatitude
	if(!NVAR_Exists(endLatitude))
		variable/G endLatitude = 40
	endif
	
	NVAR/Z endLongitude
	if(!NVAR_Exists(endLongitude))
		variable/G endLongitude = 140
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