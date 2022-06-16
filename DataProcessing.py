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

tb=UTCDateTime("2017-01-18T09:25:40.00Z")

te=UTCDateTime("2017-01-18T09:27:09.90Z")
#t1=np.arange(0,102.4,0.2)
st=read('data/6stats-2017-01-18-mww53-central-italy.miniseed')
st.resample(5.0,strict_length=False)

inv=read_inventory('./data/ingv.xml')

# inv_new = inv.select(channel='HH*')
# inv_new=inv_new.remove(station='A*')
# inv_new=inv_new.remove(station='CE*')
# inv_new=inv_new.remove(station='F*')
# inv_new=inv_new.remove(station='G*')
# inv_new=inv_new.remove(station='OF*')
# inv_new=inv_new.remove(station='T1*')

# 'SRES'  'SNTG' 'MOMA' 'T0110' 'TRTR' 'PTQR'

for i in range(len(st)):
    tr=st[i]
    tr.detrend('linear')
    pre_filt=[0.001, 0.005, 45, 50]
    tr.remove_response(inventory=inv,pre_filt=pre_filt,output="VEL",water_level=60)
    tr.filter('bandpass',freqmin=0.03,freqmax=0.1,corners=4)
    tr.trim(tb,te,pad=True)

st.rotate(method="->ZNE", inventory=inv)

st.resample(1.0,strict_length=False)

st.write("./data/ProcessedData.miniseed",format="MSEED")

