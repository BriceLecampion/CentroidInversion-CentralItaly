#
# This file is part of CentroidInversion-CentralItaly.
#
# Created by Brice Lecampion on 15.06.22.
# Copyright (c) ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, Geo-Energy Laboratory, 2016-2022.  All rights reserved.
# See the LICENSE.TXT file for more details. 
#


#
from obspy import read, read_inventory
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt

tb=UTCDateTime("2017-01-18T09:25:40.Z")

te=UTCDateTime("2017-01-18T09:27:40Z")
#t1=np.arange(0,102.4,0.2)
st=read('data/2017-01-18-mww53-central-italy.miniseed')
st.resample(5.0)

stations=[]
network =[]
for i in range(len(st)):
    stations.append(st[i].stats["station"])
    network.append(st[i].stats["network"])
# network up to 158 is 'IV'
iv_stations=stations[0:158][0::3]
iv_stations = list( dict.fromkeys(iv_stations) )
#
from obspy.clients.fdsn import Client

# Put the correct network name
netw = 'IV'
client_ingv = Client('INGV')

ingv_net = client_ingv.get_stations(network=netw,
                            station="ARRO,ASSB,ATCC,CADA,CESI,CESX,CSP1,EL6,FAGN,FEMA,FIAM,FIU1,GAG1,GIGS,MDAR,MNTP,MOMA,MTL1,OFFI,PIO1,PTQR,SEF1,SNTG,SRES,T0110,T1217,T1220,T1247,TRTR",
                            level='response')

ingv_net.write('ingv.xml', 'STATIONXML')

inv=read_inventory('./data/ingv.xml')

inv_new = inv.select(channel='HH*')
inv_new=inv_new.remove(station='A*')
inv_new=inv_new.remove(station='CE*')
inv_new=inv_new.remove(station='F*')
inv_new=inv_new.remove(station='G*')
inv_new=inv_new.remove(station='OF*')
inv_new=inv_new.remove(station='T1*')

st_new=st.select(inventory=inv_new)

# 'SRES'  'SNTG' 'MOMA' 'T0110' 'TRTR' 'PTQR'

#looks at one trace

for i in range(len(st_new)):
    tr=st_new[i]
    tr.detrend('linear')
    pre_filt=[0.001, 0.005, 45, 50]
    tr.remove_response(inventory=inv,pre_filt=pre_filt,output="VEL",water_level=60)
    tr.filter('bandpass',freqmin=0.03,freqmax=0.1,corners=4)
    tr.trim(tb,te)
    tr.plot()

# SNTG

