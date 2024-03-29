# LENGTH
m2km = 1e-3 
km2m = 1e3

au2km = 149597870700*m2km
km2au = 1/au2km

# TIME
min2sec = 60
sec2min = 1/min2sec

hr2min = 60
min2hr = 1/hr2min

day2hr = 24
hr2day = 1/day2hr

sec2hr = sec2min*min2hr
hr2sec = hr2min*min2sec

sec2day = sec2hr*hr2day
day2sec = day2hr*hr2sec

min2day = min2hr*hr2day
day2min = day2hr*hr2min

yr2day = 365.25 # Julian year
day2yr = 1/yr2day

yr2hr = yr2day*day2hr
hr2yr = hr2day*day2yr

yr2min = yr2hr*hr2min
min2yr = min2hr*hr2yr

yr2sec = yr2min*min2sec
sec2yr = sec2min*min2yr
