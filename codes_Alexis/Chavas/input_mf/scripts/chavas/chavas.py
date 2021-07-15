# Script principal de Chavas
# Appel les autres fonctions
# Morgane Guiziou 2019

import numpy as np
import sys
import copy
from shapely.geometry import LineString
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import CLE15
import matplotlib.pyplot as plt
import glob
from glob import glob
from datetime import datetime
import time
import pandas as pd
import h5py
from palette import *
palette = '/net/diva/home/mcms/guizioum/high_wind_speed.pal'
cwnd = getColorMap( rgbFile = palette )

# Debut du decompte du temps
start_time = time.time()

########################## TEST VALUE ##########################

#Az=np.arange(370)
#Az=Az[0:370:10] 		# angle 0 to 360 
#Rfit=200+40*np.cos(Az*np.pi/180) 		# R34
#Vfit=17		 # wind speed 34 kts (17m/s)
#Vmax=60 		# max wind speed in best track
#Vfm=5 
#Dfm=70
#Lat=12
#Rlim=800

###############################################################

#j='/net/diva/home/mcms/guizioum/chavas/SMOS_201802011243_SH07_CEBILE_FIX.mat'
#j='/net/diva/home/mcms/guizioum/chavas/SMOS_201808181505_EP14_LANE_FIX.mat'
j='/net/diva/home/mcms/guizioum/chavas/SMOS_201808190329_EP14_LANE_FIX.mat'

f=h5py.File(j,'r')
Az=np.arange(370)
Az=Az[0:361:1]		# angle 0 to 360
file_name= j.split('/')[-1]
date=file_name.split('_')[1]
date1=datetime.strptime(date[0:8], '%Y%m%d')
date1=date1.strftime('%Y/%m/%d')
name_cycl = file_name.split('_')[3]

					# VARIABLES #

date_file=f['date'][:]
date_orig=f['track_orig']['date'][:]
Xdate=np.abs(date_file-date_orig)
IDXdate=np.where(Xdate==Xdate.min())
Lat=f['clat'][:]
Vmax=f['max_ws_kph'][:]/3.6              		 # Wind speed m/s
destlat=f['track_orig']['clat'][:,36]
sourcelat=f['track_orig']['clat'][:,35]
destlon=f['track_orig']['clon'][:,36]
sourcelon=f['track_orig']['clon'][:,35]
destT=f['track_orig']['date'][:,36]*24*3600		 #time second
sourceT=f['track_orig']['date'][:,35]*24*3600		 #time second


					# CALCUL VFIT #
					
dlon=6371000*(np.pi/180)*(destlon-sourcelon)
dlat=6371000*(np.pi/180)*(destlat-sourcelat)
d=np.sqrt(dlon**2+dlat**2)

dt=destT-sourceT
Vfm=d/dt

#Dfm=258 # Cebile
#Dfm=28 # Lane 18/08
Dfm=109 # Lane 19/08

Vfit=17
Rlim=800


#df=glob('/net/diva/home/mcms/guizioum/chavas/data_r34_CEBILE_201802011243.txt')
#df=glob('/net/diva/home/mcms/guizioum/chavas/data_r34_LANE_201808181505.txt')
df=glob('/net/diva/home/mcms/guizioum/chavas/data_r34_LANE_201808190329.txt')

df=str(df)
df=df.replace("['","")
df=df.replace("']","")
df=pd.read_csv(df)
Rfit = df.iloc[:,0]

RRF, WWF, xf,yf=CLE15.glochv(Rfit,Vfit,Az,Vmax,Vfm,Dfm,Lat,Rlim)

fig = figure(figsize = (12,12))
im=plt.pcolormesh(xf,yf,WWF,cmap=cwnd,vmin=0,vmax=80)
plt.xlim([-200,200])
plt.ylim([-200,200])
plt.title('Chavas modelisaton of wind speed field '+name_cycl+' '+date1 , fontsize=20)
plt.xlabel('Longitude',fontsize=15)
plt.ylabel('Latitude',fontsize=15)
cb=plt.colorbar(im,orientation='horizontal')
cb.set_label(label='Wind speed (m/s)',fontsize=15)
savefig('/net/merzhin1/vol/safo/morgane/data/figures/CHAVAS/chavas_modelisation_wind_speed_field_'+name_cycl+'_'+date+'.png')

fig1 = figure(figsize = (12,12))
plt.plot(RRF,WWF)
plt.title('Wind speed profil '+name_cycl+' '+date1 , fontsize=20)
plt.xlabel('Radius (km)',fontsize=15)
plt.ylabel('Wind speed (m/s)',fontsize=15)
savefig('/net/merzhin1/vol/safo/morgane/data/figures/CHAVAS/wind_speed_profil_'+name_cycl+'_'+date+'.png')

plt.close('all')

np.savetxt('WWF_'+name_cycl+'_'+date+'.txt',WWF)
np.savetxt('RRF_'+name_cycl+'_'+date+'.txt',RRF)                                                
np.savetxt('xf_'+name_cycl+'_'+date+'.txt',xf)                                                  
np.savetxt('yf_'+name_cycl+'_'+date+'.txt',yf)  

##############################################


#RRF=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/data/RRF_'+name_cycl+'_'+date+'.txt')
#WWF=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/data/WWF_'+name_cycl+'_'+date+'.txt')
#xf=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/data/xf_'+name_cycl+'_'+date+'.txt')
#yf=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/data/yf_'+name_cycl+'_'+date+'.txt')




# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))

