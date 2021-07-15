# Script permettant de tracer le profil radial du vent à basse résolution



# Headers

import h5py
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from datetime import datetime



######################################################################################################
#################################			CHAVAS			##########################################
######################################################################################################



chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201802011243_SH07_CEBILE_FIX.mat'
#chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201808181505_EP14_LANE_FIX.mat'
#chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201808190329_EP14_LANE_FIX.mat'


f=h5py.File(chv,'r')
file_name= chv.split('/')[-1]
date=file_name.split('_')[1]
date1=datetime.strptime(date[0:8], '%Y%m%d')
date1=date1.strftime('%Y/%m/%d')
name_cycl = file_name.split('_')[3]
clat=f['clat'][:]
clon=f['clon'][:]


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



smos = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/CEBILE_20180201T124400_SAR_s1a_SMOS.nc'
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



######################################################################################################
#################################			SAR			##############################################
######################################################################################################


sar = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/CEBILE_20180201T124400_SAR_s1a_SMOS.nc'
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








######################################################################################################
#################################			SMAP			##########################################
#################################      (seulement LANE) 	##########################################
######################################################################################################


smap = '/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/LANE_20180818T151100_SAR_s1a_SMAP.nc'


file_smap = nc.Dataset(smap,'r')
delta = np.squeeze(file_smap.variables['delta_time'][:])
lon_smap = np.squeeze(file_smap.variables['Lon'][:])
lat_smap = np.squeeze(file_smap.variables['Lat'][:])
vent_sar = np.squeeze(file_smap.variables['WindSpeed_sar_mean'][:])
vent_sar = np.ma.masked_where((vent_sar>= 999),vent_sar)
vent_smap = np.squeeze(file_smap.variables['WindSpeed_smap'][:]) 
vent_smap_a = np.ma.masked_where((vent_smap >= 999)*(vent_smap.mask == False)*(vent_sar.mask == False)+(vent_sar.mask == True) , vent_smap)

dimsmap=vent_smap.shape
xsmap=np.zeros(dimsmap)
ysmap=np.zeros(dimsmap)

for k in range(len(xsmap)): 
	xsmap[k,:]=np.arange(dimsmap[1]) 

for l in range(len(ysmap[1])): 
	ysmap[:,l]=np.arange(dimsmap[0])

Xsmap=np.abs(lon_smap-clon)
xsmapc=np.ma.where(Xsmap==Xsmap.min())
Ysmap= np.abs(lat_smap-clat)
ysmapc=np.ma.where(Ysmap==Ysmap.min())

lon_smap1,lat_smap1=np.meshgrid(lon_smap,lat_smap)
dsmap=np.sqrt((lon_smap1-lon_smap[xsmapc])**2+(lat_smap1-lat_smap[ysmapc])**2)
dsmap_km=dsmap*111



vent_sar_mean=[]
dsar_mean=[]
vent_smap_mean=[]
dsmap_mean=[]
vent_smos_mean=[]
dsmos_mean=[]
vent_chv_mean=[]
dchv_mean=[]

p=25

r=np.arange(1700)
rmin=r[0:1500:p]
rmax=r[p:1600:p]

for i in range(len(rmin)):
	vent_smos_mean=np.append(vent_smos_mean,vent_smos[(dsmos_km>=rmin[i])*(dsmos_km<=rmax[i])].mean())
	dsmos_mean=np.append(dsmos_mean,dsmos_km[(dsmos_km>=rmin[i])*(dsmos_km<=rmax[i])].mean())
	vent_sar_mean=np.append(vent_sar_mean,vent_sar[(dsar_km>=rmin[i])*(dsar_km<=rmax[i])].mean())
	dsar_mean=np.append(dsar_mean,dsar_km[(dsar_km>=rmin[i])*(dsar_km<=rmax[i])].mean())
	vent_sar_mean[vent_sar_mean==0]=np.nan
	vent_smap_mean=np.append(vent_smap_mean,vent_smap[(dsmap_km>=rmin[i])*(dsmap_km<=rmax[i])].mean())
	dsmap_mean=np.append(dsmap_mean,dsmap_km[(dsmap_km>=rmin[i])*(dsmap_km<=rmax[i])].mean())
	vent_smap_mean[vent_smap_mean==0]=np.nan
	vent_chv_mean=np.append(vent_chv_mean,WWF[(RRF>=rmin[i])*(RRF<=rmax[i])].mean())
	dchv_mean=np.append(dchv_mean,RRF[(RRF>=rmin[i])*(RRF<=rmax[i])].mean())

fig=plt.figure(figsize=(12,12))
plt.plot(dsmos_mean,vent_smos_mean,'r',label='SMOS')
#plt.plot(dsmap_mean,vent_smap_mean,'b',label='SMAP')
plt.plot(dchv_mean,vent_chv_mean,'g',label='CHAVAS')
plt.plot(dsar_mean,vent_sar_mean,'orange',label='SENTINEL-1')
plt.title('Wind speed profil global '+name_cycl+' '+date+' pas = '+str(p)+'km',fontsize=20)
plt.xlabel('Radius (km)',fontsize=16)
plt.ylabel('Wind speed',fontsize=16)
plt.xlim([0,800])
plt.ylim([0,60])
plt.grid()
plt.legend()
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/CHAVAS/wind_profil_'+name_cycl+'_'+date+'_(p='+str(p)+')_1.png')


#plt.close('all')





