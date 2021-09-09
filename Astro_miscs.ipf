#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


Menu "Astronomical calculation"
	Submenu "miscs"
		"Create track between two points"
	End
End

Function LanchCreateTimeSeriesPath()

	if(!DataFolderExists("root:Astronomic"))
		NewDataFolder "root:Astronomic"
	endif

	if(!DataFolderExists("root:Astronomic:miscs"))
		NewDataFolder "root:Astronomic:miscs"
	endif

	if(!DataFolderExists("root:Astronomic:miscs:track"))
		NewDataFolder "root:Astronomic:miscs:track"
	endif

	CreateTimeSeriesPath()
End


Function CreateTimeSeriesPath()

	string pathToParameter = "root:Astronomic:miscs:track:"

	SVAR startDateString = $(pathToParameter + "startDateTimeString")
	SVAR startTimeString = $(pathToParameter + "startTimeString")
	SVAR endDateString = $(pathToParameter + "endDateTimeString")
	SVAR endTimeString = $(pathToParameter + "endTimeString")
	
	NVAR intervalOfTime = $(pathToParameter + "intervalOfTime")

	NVAR startLatitude = $(pathToParameter + "startLatitude")
	NVAR startLongitude = $(pathToParameter + "startLongitude")
	NVAR endLatitude = $(pathToParameter + "endLatitude")
	NVAR endLongitude = $(pathToParameter + "endLongitude")

	variable startDateTime = ReadDateTimeString(startDateString, startTimeString)
	variable endDateTime = ReadDateTimeString(endDateString, endTimeString)

	variable numberOfPoints = (endDateTime - startDateTime)/intervalOfTime
	
	make/O/N=numpnts(numberOfPoints) Track

End

Function ReadDateTimeString(dateString, timeString)

	string dateString
	string timeString
	
	variable year, month, day
		sscanf dateString, "%d/%d/%d", year, month, day
	variable dateValue = date2secs(year, month, day)

	variable hour, minute, second
		sscanf timString, "%d:%d:%d", hour, minute, second
	variable timevalue = hour*60*60 + minute*60 + second

	variable dateTimeValue = dateValue + timevalue
	return dateTimeValue
End
