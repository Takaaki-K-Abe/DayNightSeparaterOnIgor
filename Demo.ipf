#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

Macro Runtest()

variable year = 1904
variable TimeDiffFromUTC = 9
variable Latitude = 35
variable Longitude = 135
variable NumDays = 60

make/O/N=(24*60*60*NumDays) eq1
SetScale/P x 2678400,1,"dat", eq1

CalculateSolarAltitude(root:Eq1, Year, TimeDiffFromUTC, Latitude, Longitude)

End