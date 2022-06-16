#
# This file is part of CentroidInversion-CentralItaly.
#
# Created by Brice Lecampion on 15.06.22.
# Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2022.  All rights reserved.
# See the LICENSE.TXT file for more details. 
#
#

# this is a preliminary simulation for all stations in the selected list in order to iron out all the proper formatting of
# the forward problem via axitra and checks that base hypocenter location reported from INGV is not that bad as first guess

from obspy import read, read_inventory,Trace,Stream
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
import sys
axitra_path='./axitra/MOMENT_DISP_F90_OPENMP/src/'
sys.path.extend([axitra_path])
from axitra import *

#
# EQ metadata taken from terromoti.ingv catalog
eq_hc= {"time":UTCDateTime("2017-01-18 09:25:41Z"), # slightly different origin time than from IRIS
    "Latitude" : 42.6591,"Longitude" : 13.2086,"Depth" : 6000,"Mw" : 5.3,
        "Mo":6.39e23,"Dip":58,"Strike":331,"Rake":-91}

tb=UTCDateTime("2017-01-18T09:25:40.00Z") # this is the origin time from IRIS

# loading up the processed data (see script DataProcessing.py )
st=read('data/ProcessedData.miniseed')
inv=read_inventory('./data/ingv.xml')

# take only ONE single station for now
# select between
station_list=['MOMA', 'PTQR', 'SNTG', 'SRES', 'T0110', 'TRTR']
# note MOMA N& Z zeros.
station_coor=np.zeros([len(station_list),4])
for i in range(len(station_list)):
    station_name=station_list[i]
    st2=st.select(station=station_name)
    st2.plot()
    inv2=inv.select(station=station_name)
    station=inv2.networks[0].stations[0]
    station_coor[i]=np.array([int(i),station.latitude,station.longitude,0.*station.elevation])

# change station coordinate into a planar coordinate system
origin_y = np.min(station_coor, axis = 0)[1]
origin_x = np.min(station_coor, axis = 0)[2]
stations_xy = station_coor #xy coordinates of stations (unit in m) #stations to stations_xy is deep copy!
stations_xy.astype(int)
for i in range(station_coor.shape[0]):
    stations_xy[i][1] = (station_coor[i][1] - origin_y) * 100000
    stations_xy[i][2] = (station_coor[i][2] - origin_x) * 100000

############################################################
# Regional of Aquila & surroundings
############################################################
# thickness (or top), Vp, Vs, rho, Qp, Qs
vel_model = np.array([[0., 3750., 2140., 2270., 350., 175.],
                  [1500., 4940., 2820., 2485., 640., 320.],
                  [4500., 6010., 3430., 2706., 800., 400.],
                  [7500., 5550., 3150., 2609., 840., 420.],
                  [14500., 5880., 3360., 2677., 860., 430.],
                  [29500., 7110., 4010., 3010., 900., 450.],
                  [35500., 7100., 3990., 3012., 1000., 500.],
                  [43500., 7900., 4400., 3276., 1500., 750.]])

############################################################
# earthquake hypocenter
# in the planar system
hyp_y = (eq_hc["Latitude"] - origin_y) * 100000
hyp_x = (eq_hc["Longitude"] - origin_x) * 100000
hyp_z = eq_hc["Depth"] # depth
sources=np.array([[1,hyp_x, hyp_y,hyp_z]] )

M0=7.25e17 #eq_hc["Mo"] # Seismic moment
strike = eq_hc["Strike"]
dip = eq_hc["Dip"]
rake =  eq_hc["Rake"]
t_rise=0.2

hist = np.array([[1,M0,strike,dip,rake,0.,0.,t_rise]])

dist = np.linalg.norm(sources[0,1:4]-station_coor[0,1:4])
print('distance is ',dist, 'meters')
rho = np.mean(vel_model[:,3])  #density
beta = np.mean(vel_model[:,2])#Vs

LowFreq_asympt = M0* 1./(4.*np.pi*rho*beta**3) / dist
print('Low frequency asymptote ',LowFreq_asympt)


##### FORWARD SIMULATION
# Fill in the instance of Axitra Class
ap = Axitra(vel_model, station_coor, sources, fmax=2.5, duration=90, latlon=False,axpath=axitra_path)
ap = moment.green(ap)# Compute green's function
t, vx, vy, vz = moment.conv(ap,hist,source_type=4,t0=0.05,unit=2) # convolve with DC moment tensor source fction
delta_t= t[2]-t[1]

# SOME USEFUL FUNCTIONS
# An utility function to process the axitra results and format them into a ObsPy Trace object
from obspy import Trace
def createSyntheticTraceFromAxitra(vs,delta_s,fmin=0.03,fmax=0.1,re_sample=1) :
    """
    :param vs: 1D np.array containing the signal sample at reg. inter
    :param delta_s: the DT
    :param fmin: min frequency of the band pass filter
    :param fmax: max frequency of the band pass filter
    :param re_sample: new sampling rate of the final trace.
    :return:
        an obsPy Trace object properly filtered and resampled
    """
    synt_tr = Trace(vs)
    synt_tr.stats.delta=delta_s # to set the proper sampling rate of the time trace
    synt_tr.filter('bandpass',freqmin=fmin,freqmax=fmax,corners=4) # low-pass to the same frequency band than the data
    synt_tr.resample(re_sample)
    return synt_tr

import scipy
# function for plotting prediction and data at a given station ....
def DatavsPrediction(station_name,plot=True):
    """
    :param station_name: string with station name

    :return: rms
    """
    i=station_list.index(station_name)
    tr_x=createSyntheticTraceFromAxitra(vx[i,:],delta_t)
    tr_y=createSyntheticTraceFromAxitra(vy[i,:],delta_t)
    tr_z=createSyntheticTraceFromAxitra(vz[i,:],delta_t)
    st2=st.select(station=station_name)
    rms_x = np.sum((st2[0].data-tr_x.data)**2.0) / ( np.sum((st2[0].data)**2.0) + np.sum((tr_x.data)**2.0) )
    rms_y = np.sum((st2[1].data-tr_y.data)**2.0) / ( np.sum((st2[1].data)**2.0) + np.sum((tr_y.data)**2.0) )
    rms_z = np.sum((st2[2].data-tr_z.data)**2.0) / ( np.sum((st2[2].data)**2.0) + np.sum((tr_z.data)**2.0) )
    rms=np.array([rms_x,rms_y,rms_z])
    rms_t=0.
    n_tr=0
    isN=np.isnan(rms)
    for c in range(3):
        if isN[c]==False:
            rms_t+=rms[c]
            n_tr=n_tr+1
    if n_tr>0:
        rms_t = rms_t /n_tr

    if plot==True :
        print("rms station "+station_name,rms_t)
        plt.plot(tr_x.data,'.-r')
        plt.plot(st2[0].data,'.-k')
        plt.xlabel("Samples")
        plt.ylabel('x-velocity (m/s)')
        plt.title('x-velocity (m/s) - station :'+ station_name)
        plt.legend(['pred','data'])
        plt.show()
        plt.plot(tr_y.data,'.-r')
        plt.plot(st2[1].data,'.-k')
        plt.xlabel("Samples")
        plt.ylabel('y-velocity (m/s)')
        plt.title('y-velocity (m/s) - station :'+ station_name)
        plt.legend(['pred','data'])
        plt.show()
        plt.plot(tr_z.data, '.-r')
        plt.plot(st2[2].data, '.-k')
        plt.xlabel("Samples")
        plt.ylabel('z-velocity (m/s)')
        plt.legend(['pred', 'data'])
        plt.title('z-velocity (m/s) - station :' + station_name)
        plt.show()

    return rms_t

def rmsAllStations(plot=False):
    rms=0.
    for j in range(len(station_list)):
        rms +=DatavsPrediction(station_list[j],plot=plot)
    rms=rms/len(station_list)
    return rms


# Check all stations
rmsAllStations(plot=True)

# now that we see that the predicted and measured traces make sense
# we can move onto the grid search...
