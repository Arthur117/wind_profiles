# Headers

import h5py
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from datetime import datetime
import glob
from glob import glob
import time
#import fonc_data
from palette import *
palette = '/net/diva/home/mcms/guizioum/high_wind_speed.pal'
cwnd = getColorMap( rgbFile = palette )
import csv

# Debut du decompte du temps
start_time = time.time()

clat = 12.2552 #LANE_20180818
clon = -138.2162 #LANE_20180818

#clat = 12.7092 #LANE_20180819
#clon = -141.0608 #LANE_20180819

#clat = -15.7321 #CEBILE_20180201
#clon = 76.0757 #CEBILE_20180201

######################################################################################################
#################################			CHAVAS			##########################################
######################################################################################################



#chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201802011243_SH07_CEBILE_FIX.mat'
#chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201808181505_EP14_LANE_FIX.mat'
#chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201808190329_EP14_LANE_FIX.mat'


f=h5py.File(chv,'r')
file_name= chv.split('/')[-1]
date=file_name.split('_')[1]
date1=datetime.strptime(date[0:8], '%Y%m%d')
date1=date1.strftime('%Y/%m/%d')
name_cycl = file_name.split('_')[3]

RRF=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/data/RRF_'+name_cycl+'_'+date+'.txt')
WWF=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/data/WWF_'+name_cycl+'_'+date+'.txt')
xf=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/data/xf_'+name_cycl+'_'+date+'.txt')
yf=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/data/yf_'+name_cycl+'_'+date+'.txt')

WWF_mean=[]
for mchv in range(len(RRF)):
	WWF_mean=np.append(WWF_mean,WWF[(RRF==mchv)].mean())

######################################################################################################
#################################			SMOS 		##############################################
######################################################################################################



#smos = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/CEBILE_20180201T124400_SAR_s1a_SMOS.nc'
#smos = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/LANE_20180818T150600_SAR_s1a_SMOS.nc'
#smos = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/LANE_20180819T032900_SAR_s1a_SMOS.nc'


file_smos = nc.Dataset(smos,'r')
delta = np.squeeze(file_smos.variables['delta_time'][:])
lon = np.squeeze(file_smos.variables['Lon'][:])
lat = np.squeeze(file_smos.variables['Lat'][:])
vent_sar = np.squeeze(file_smos.variables['WindSpeed_sar_mean'][:])
vent_sar = np.ma.masked_where((vent_sar>= 999),vent_sar)
vent_smos = np.squeeze(file_smos.variables['WindSpeed_smos'][:]) 
vent_smos_a = np.ma.masked_where((vent_smos >= 999)*(vent_smos.mask == False)*(vent_sar.mask == False)+(vent_sar.mask == True) , vent_smos)

dim=vent_smos.shape
xsmos=np.zeros(dim)
ysmos=np.zeros(dim)

for k in range(len(xsmos)): 
	xsmos[k,:]=np.arange(dim[1]) 

for l in range(len(ysmos[1])): 
	ysmos[:,l]=np.arange(dim[0])

X=np.abs(lon-clon)
xsmosc=np.ma.where(X==X.min())
Y= np.abs(lat-clat)
ysmosc=np.ma.where(Y==Y.min())

lon_smos,lat_smos=np.meshgrid(lon,lat)
dsmos=np.sqrt((lon_smos-lon[xsmosc])**2+(lat_smos-lat[ysmosc])**2)
dsmos_km=dsmos*111 

Az=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/angle_'+name_cycl+'_'+date+'_smos.txt',delimiter=',')
Z=Az.size

x_smos=(dsmos_km)*np.cos(Az.T*np.pi/180);  
y_smos=(dsmos_km)*np.sin(Az.T*np.pi/180);

######################################################################################################
#################################			SAR			##############################################
######################################################################################################


#sar = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/CEBILE_20180201T124400_SAR_s1a_SMOS.nc'
#sar = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/LANE_20180818T150600_SAR_s1a_SMOS.nc'
#sar = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/LANE_20180819T032900_SAR_s1a_SMOS.nc'

file_sar = nc.Dataset(sar,'r')
delta = np.squeeze(file_sar.variables['delta_time'][:])
lon = np.squeeze(file_sar.variables['Lon'][:])
lat = np.squeeze(file_sar.variables['Lat'][:])
vent_sar = np.squeeze(file_sar.variables['WindSpeed_sar_mean'][:])
vent_sar = np.ma.masked_where((vent_sar>= 999),vent_sar)


dimsar=vent_sar.shape
xsar=np.zeros(dimsar)
ysar=np.zeros(dimsar)

for i in range(len(xsar)): 
	xsar[i,:]=np.arange(dimsar[1]) 

for l in range(len(ysar[1])): 
	ysar[:,l]=np.arange(dimsar[0])

Xsar=np.abs(lon-clon)
xsarc=np.ma.where(Xsar==Xsar.min())
Ysar= np.abs(lat-clat)
ysarc=np.ma.where(Ysar==Ysar.min())

lon_sar,lat_sar=np.meshgrid(lon,lat)
dsar=np.sqrt((lon_sar-lon[xsarc])**2+(lat_sar-lat[ysarc])**2)
dsar_km=dsar*111 




x1=np.linspace(-100,100,1000)
y1=np.linspace(-100,100,1000)
X,Y=np.meshgrid(x1,y1)

F=X**2+Y**2-100**2
F2=X**2+Y**2-50**2
F3=X**2+Y**2-25**2

[x,y]=np.where(WWF==WWF.max())

fig1=plt.figure(figsize=(12,12))
im=plt.pcolormesh(xf,yf,WWF,cmap=cwnd,vmin=0,vmax=80)
plt.contour(X,Y,F,[0])
plt.contour(X,Y,F2,[0])
plt.contour(X,Y,F3,[0])
plt.scatter(xf[x,y],yf[x,y],marker='+',color='k',label='position Vmax')
plt.xlim([-200,200])
plt.ylim([-200,200])
plt.title('Chavas modelisaton of wind speed field '+name_cycl+' '+date1 , fontsize=20)
plt.xlabel('Distance from center',fontsize=15)
plt.ylabel('Distance from center',fontsize=15)
plt.legend()
cb=plt.colorbar(im,orientation='horizontal')
cb.set_label(label='Wind speed (m/s)',fontsize=15)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/CHAVAS/chavas_modelisation_wind_speed_field_'+name_cycl+'_'+date+'.png')

fig2=plt.figure(figsize=(12,12))
im=plt.pcolormesh(lon,lat,vent_sar,cmap=cwnd,vmin=0,vmax=80)
plt.scatter(clon,clat,marker='+',color='k',label='TC center')
plt.scatter(lon[a],lat[b],marker='*',color='k',label='position Vmax')
plt.title('SAR wind speed field '+name_cycl+' '+date1 , fontsize=20)
plt.xlabel('Longitude',fontsize=15)
plt.ylabel('Latitude',fontsize=15)
plt.xlim([72,79])
plt.ylim([-20,-13])
cb=plt.colorbar(im,orientation='horizontal')
cb.set_label(label='Wind speed (m/s)',fontsize=15)

