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
import fonc_data
from palette import *
palette = '/net/diva/home/mcms/guizioum/high_wind_speed.pal'
cwnd = getColorMap( rgbFile = palette )


# Debut du decompte du temps
start_time = time.time()


######################################################################################################
#################################			CHAVAS		##############################################
######################################################################################################



#chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201802011243_SH07_CEBILE_FIX.mat'
#chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201808181505_EP14_LANE_FIX.mat'
chv='/net/diva/home/mcms/guizioum/chavas/SMOS_201808190329_EP14_LANE_FIX.mat'


f=h5py.File(chv,'r')
file_name= chv.split('/')[-1]
date=file_name.split('_')[1]
date1=datetime.strptime(date[0:8], '%Y%m%d')
date1=date1.strftime('%Y/%m/%d')
name_cycl = file_name.split('_')[3]
clat=f['clat'][:]
clon=f['clon'][:]

RRF=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/RRF_'+name_cycl+'_'+date+'.txt')
WWF=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/WWF_'+name_cycl+'_'+date+'.txt')
xf=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/xf_'+name_cycl+'_'+date+'.txt')
yf=np.loadtxt('/net/diva/home/mcms/guizioum/chavas/yf_'+name_cycl+'_'+date+'.txt')

WWF_mean=[]
for mchv in range(len(RRF)):
	WWF_mean=np.append(WWF_mean,WWF[(RRF==mchv)].mean())
	
	






######################################################################################################
#################################			SAR			##############################################
######################################################################################################


#j= '/net/merzhin1/vol/safo/morgane/data/CEBILE/sentinel1/s1a-ew-owi-cm-20180201t132219-20180201t132527-000003-022E2F.nc'
#j = '/net/merzhin1/vol/safo/morgane/data/LANE/sentinel1/s1a-ew-owi-cm-20180818t150300-20180818t150504-000003-0288C3.nc'
j = '/net/merzhin1/vol/safo/morgane/data/LANE/sentinel1/s1a-ew-owi-cm-20180819t032959-20180819t033303-000003-0288F7.nc'



if clon <=0:	
	clon = clon+360 
else:
	clon=clon

file_vent=nc.Dataset(j,'r')
grp = file_vent.groups['owiInversionTables_UV']
lon_sar=np.squeeze(file_vent.variables['owiLon'][:])
lat_sar=np.squeeze(file_vent.variables['owiLat'][:])
ws_sar=np.squeeze(grp.variables['owiWindSpeed_Tab_dualpol_2steps'][:])
mask=np.squeeze(file_vent.variables['owiMask'][:])
angle=np.squeeze(file_vent.variables['owiIncidenceAngle'][:])
name_cycl = j.split('/')[7]
sat = j.split('/')[-2]
name_sar, sat_sar, x, y, z, w, xlim, ylim, day_sar,day_sar_1, date_sar, hour_sar_debut, hour_sar_fin = fonc_data.data_sar(lon_sar,lat_sar,mask,ws_sar,angle,j)


dimsar=z.shape
xsar=np.zeros(dimsar)
ysar=np.zeros(dimsar)

for i in range(len(xsar)): 
	xsar[i,:]=np.arange(dimsar[1]) 

for l in range(len(ysar[1])): 
	ysar[:,l]=np.arange(dimsar[0])

Xsar=np.abs(lon_sar-clon)
xsarc=np.ma.where(Xsar==Xsar.min())
Ysar= np.abs(lat_sar-clat)
ysarc=np.ma.where(Ysar==Ysar.min())


#lon,lat=np.meshgrid(lon_sar,lat_sar)
dsar=np.sqrt((lon_sar-lon_sar[xsarc[0][0],xsarc[1][0]])**2+(lat_sar-lat_sar[ysarc])**2)
dsar_km=dsar*111 

vent_sar_mean=[]
dsar_mean=[]
vent_chavas_mean=[]
dchv_mean=[]

p=1

r=np.arange(1700)
rmin=r[0:1500:p]
rmax=r[p:1600:p]



for i in range(len(rmin)):
	vent_sar_mean=np.append(vent_sar_mean,z[(dsar_km>=rmin[i])*(dsar_km<=rmax[i])].mean())
	dsar_mean=np.append(dsar_mean,dsar_km[(dsar_km>=rmin[i])*(dsar_km<=rmax[i])].mean())
	vent_sar_mean[vent_sar_mean==0]=np.nan

	vent_chavas_mean=np.append(vent_chavas_mean,WWF[(RRF>=rmin[i])*(RRF<=rmax[i])].mean())
	dchv_mean=np.append(dchv_mean,RRF[(RRF>=rmin[i])*(RRF<=rmax[i])].mean())



fig=plt.figure(figsize = (10,10))
plt.plot(dchv_mean,vent_chavas_mean,'g',label='CHAVAS')
plt.plot(dsar_mean,vent_sar_mean,'orange',label='SENTINEL-1')
plt.title('Wind speed profil global '+name_cycl+' '+date1+' resolution = '+str(p)+'km',fontsize=20)
plt.xlabel('Radius (km)',fontsize=16)
plt.ylabel('Wind speed',fontsize=16)
plt.grid()
plt.xlim([0,800])
plt.ylim([0,60])
plt.legend()
plt.savefig('/net/merzhin1/vol/safo/morgane/data/figures/CHAVAS/wind_profil_'+name_cycl+'_'+date+'_(p='+str(p)+')_HD_1.png')

#plt.close('all')





