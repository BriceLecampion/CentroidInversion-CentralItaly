Point centroid inversion of a Earthquake
======
## EPFL doctoral course Introduction to Earthquake source physics (2022)
### Brice Lecampion - June 16, 2022 (6 hours of effective work)

In this repo, you will find examples of scripts for the processing, forward modeling (using axitra) and grid-search for the inversion of a Mw5.1 EQ that occured in Central Italy in 2017

I just took few broadband stations (6) from the INGV network

EQ details 
---------
 A magnitude Mw 5.1 earthquake occured in region: 3 km NW Capitignano (AQ), on
 18-01-2017 09:25:40 (UTC) 5 years ago
    18-01-2017 10:25:40 (UTC +01:00) Italian time
and geographic coordinates (lat, lon) 42.5450, 13.2770 at 10 km depth.

The earthquake was located by: Bollettino Sismico Italiano INGV. 

Earthquake with magnitude of Mw 5.1 on date 18-01-2017 and time 10:25:40 (Italy) in region 3 km NW Capitignano (AQ)

Wilber 3: Select Stations
2017-01-18 Mww5.3 Central Italy
Latitude 	Longitude 	Date 	Depth 	Magnitude 	Description 	Related Pages
42.6591° N 	13.2086° E 	2017-01-18 09:25:41 UTC 	10.0 km 	Mww5.3 	Central Italy


Folders
-------
+ axitra - contain the source of axitra (with the makefile corrected to work with gcc)  - NOT commited/Pushed to this /
git repo - git clone axitra, make the changes & compile yourself (with axitra folder in that location)
+ data  - contain miniseed data + network inventory xml 

Scripts
-------
- DataProcessing-prelim.py  : first prelim data-processing - helpful 
- DataProcessing.py   : data processing for the 6 chosen stations (detrend, filtering, padding, resampling)
- ForwardSimulation-prelim.py : simulation for only one station ensuring things works as needed
- ForwardSimulation-AllTraces.py : 1 simulation for all stations & all traces defining some useful functions etc.
- CoarseGridSeach.py : a coarse grid search 

What is still missing
----------
- geo-plot of all stations and hypocenter 
- ... refine the domain space for the search
- ... compare with ingv results..., add "better" stations
- ... report etc.
