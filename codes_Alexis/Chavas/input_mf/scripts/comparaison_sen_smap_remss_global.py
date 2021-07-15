# Script permettant la comparaison des données vents entre Sentinel-1 et SMAP (REMSS)
# 29/03/2019 M.G

# Headers

import h5py
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime as dt
from datetime import datetime
from datetime import timedelta
from os.path import basename, splitext
from statistics import stdev
import fonc_data , fonc_remap, fonc_create_netcdf
import glob
from glob import glob
import pandas as pd
import os
import os.path
import sys
from palette import *
import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
from mpl_toolkits.basemap import Basemap
import time
from pyresample import load_area, save_quicklook, SwathDefinition
from pyresample.geometry import GridDefinition
from pyresample.kd_tree import resample_nearest

# Debut du decompte du temps
start_time = time.time()


					# Lecture fichier csv

df=pd.read_csv('/net/diva/home/mcms/guizioum/multiResReport.csv') # Sentinel-1


# Déclaration/Définition "variables"

R=df['Resolution']
L2M=df['L2M file']
name=df['Comments']
L2M_m=[]
cycl=[]


# Catégorie

cats = {}

for j,f in enumerate(L2M):
    r=R[j]
    if r == 3 :
        file_name=f.split('/')[-1]
        cycl_name=name[j].split(' ')[0]
        if 'cat-' in name[j]:
            cat = int(name[j].split('cat-')[1].split(' ')[0])
        else:
            cat = 0
        cats[file_name] = cat

				# Données provenant de Sentinel-1 (SAR)
	 

path_sar=glob('/net/merzhin1/vol/safo/morgane/data/*/sentinel1/*.nc')

for i,j in enumerate(path_sar):
	print(i,j)
	cat = cats[j.split('/')[-1]] 
	file_sar=nc.Dataset(j,'r')
	name_cycl = j.split('/')[7]
	name_sar = j.split('/')[-1]
	sat_sar = name_sar.split('-')[0]
	sat = j.split('/')[-2]
	grp = file_sar.groups['owiInversionTables_UV']
	lon_sar=np.squeeze(file_sar.variables['owiLon'][:])
	lat_sar=np.squeeze(file_sar.variables['owiLat'][:])
	ws_sar=np.squeeze(grp.variables['owiWindSpeed_Tab_dualpol_2steps'][:])
	mask=np.squeeze(file_sar.variables['owiMask'][:])
	angle=np.squeeze(file_sar.variables['owiIncidenceAngle'][:])

	name_sar, sat_sar, x, y, z, w, xlim, ylim, day_sar,day_sar_1, date_sar, hour_sar_debut, hour_sar_fin = fonc_data.data_sar(lon_sar,lat_sar,mask,ws_sar,angle,j)


				# Données provenant de SMAP (REMSS)
	

	src = '/net/merzhin1/vol/safo/morgane/data/'+name_cycl+'/smap_remss'
	path_smap_remss = glob('%s/RSS_smap_wind_daily_%s_v01.0.nc' % (src,day_sar_1))
	path_smap_remss = str(path_smap_remss)
	path_smap_remss = path_smap_remss.replace("['","")
	path_smap_remss = path_smap_remss.replace("']","")
	print(path_smap_remss)

	file_smap_remss=nc.Dataset(path_smap_remss,'r')
	file_smap_remss.set_auto_maskandscale(False)
	file_smap_remss.set_auto_scale(True)

	name_smap = path_smap_remss.split('/')[-1]
	day_smap=name_smap.split('_')[6]
	month_smap=name_smap.split('_')[5]
	year_smap=name_smap.split('_')[4]

	lon_smap=np.squeeze(file_smap_remss.variables['lon'][:])
	lat_smap=np.squeeze(file_smap_remss.variables['lat'][:])
	ws_smap=np.ma.squeeze(file_smap_remss.variables['wind'][:])
	time_smap = np.ma.squeeze(file_smap_remss.variables['minute'][:])
	sat_rad = name_smap.split('_')[1]

	try:
		vent_smap, dt, a,b,date1, a1, b1= fonc_data.data_smap_remss(lon_smap,lat_smap,ws_smap, time_smap, x, y, year_smap, month_smap, day_smap, date_sar)

							# Remapping

		lon_smap[lon_smap>180] = lon_smap[lon_smap>180]-360
		x[x>180] = x[x>180]-360
		lon, lat = meshgrid(lon_smap,lat_smap)
		swath_sar = SwathDefinition(x.data,y.data)
		swath_smap = SwathDefinition(lon,lat)

		i = np.arange(ws_smap.shape[0]) # SMAP shape
		j = np.arange(ws_smap.shape[1]) # SMAP shape
		j_smap,i_smap = np.meshgrid(j,i) # SMAP shape

		result_i = resample_nearest(swath_smap, i_smap, swath_sar, radius_of_influence=21500, fill_value=None) #SAR shape
		result_j = resample_nearest(swath_smap, j_smap, swath_sar, radius_of_influence=21500, fill_value=None) #SAR shape

		if result_i.mask.all() == False:
			ws_smapss = vent_smap[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			ws_smapss = np.ma.masked_where(ws_smapss==-99.99, ws_smapss)
			i_smapss = i_smap[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			j_smapss = j_smap[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			lon_smapss = lon[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			lat_smapss = lat[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]

			dim = i_smapss.shape
			vent_sarmoy = np.ma.zeros(dim)-999
			vent_sarmed = np.ma.zeros(dim)-999
			vent_sarstd = np.ma.zeros(dim)-999
			n_sarvalid = np.ma.zeros(dim)-999
			angle_sarmoy = np.ma.zeros(dim)-999

			for n in np.arange(i_smapss.shape[0]):
				for m in np.arange(j_smapss.shape[1]):
					s = ((result_i == i_smapss[n,m]) & (result_j == j_smapss[n,m])) # mask booléen
					vent_sarmoy[n,m] = z[s].mean()
					vent_sarmed[n,m] = np.median(z[s])
					vent_sarstd[n,m] = z[s].std()
					n_sarvalid[n,m] = len(s)-z[s].mask.sum()
					angle_sarmoy[n,m] = w[s].mean()


			angle_sarmoy = np.ma.masked_where(angle_sarmoy == -999 , angle_sarmoy)
			vent_sarmoy = np.ma.masked_where(vent_sarmoy == -999 , vent_sarmoy)
			vent_sarmed = np.ma.masked_where(vent_sarmed == -999 , vent_sarmed)
			vent_sarstd = np.ma.masked_where(vent_sarstd == -999 , vent_sarstd)
			n_sarvalid = np.ma.masked_where((n_sarvalid <= 0)*(ws_smapss==-99.99)*(vent_sarmoy == -999), n_sarvalid)

					# Graphiques comparant les vitesses des vents selon l'origine

			date_smap_moyenne = date1.strftime('%Y/%m/%d %H:%M:%S ')
			date_smap = date1.strftime('%Y%m%dT%H%M%S')
			hour_smap_moyenne = date_smap_moyenne.split(' ')[1]
			smap_date = date1.strftime('%Y/%m/%d')

			palette = '/net/diva/home/mcms/guizioum/high_wind_speed.pal'
			cwnd = getColorMap( rgbFile = palette )

			m = Basemap(llcrnrlon=x.min(),llcrnrlat=y.min(),urcrnrlon=x.max(),urcrnrlat=y.max(), 			resolution='i', projection='cyl', lat_0 =y.mean() , lon_0 =x.mean() )

			x1,y1 = m(x,y)
			a2,b2 = m(lon_smapss,lat_smapss)

			fig = figure(figsize = (16,8))
			fig.suptitle(name_cycl+' '+ smap_date+' SAR '+hour_sar_debut+' SMAP (REMSS) '+hour_smap_moyenne+ ' cat: '+str(cat),fontsize = 18)

			ax = fig.add_axes((0.05,0.1,0.68,0.02))
			ax1 = fig.add_axes((0.05,0.22,0.18,0.65))
			ax2 = fig.add_axes((0.30,0.22,0.18,0.65))
			ax3 = fig.add_axes((0.55,0.22,0.18,0.65))
			ax4 = fig.add_axes((0.78,0.22,0.18,0.65))
			ax5 = fig.add_axes((0.78,0.1,0.18,0.02))

			im = m.pcolormesh(x1.filled(x1.mean()),y1.filled(y1.mean()),z,cmap=cwnd,vmin = 0 , vmax = 80,ax = ax1)
			m.drawcoastlines(ax = ax1)
			m.fillcontinents(color='grey',ax = ax1)
			ax1.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
			ax1.set_ylabel('Latitude')
			ax1.set_title('Sentinel-1 '+sat_sar)
			m.drawparallels(np.arange(int(y.min()),int(y.max()+1),1),labels=[0,1,0,0],ax = ax1)
			m.drawmeridians(np.arange(int(x.min()),int(x.max()+1),1),labels=[0,0,0,1],ax = ax1)

			im = m.pcolormesh(a2,b2,ws_smapss,cmap=cwnd, vmin = 0 , vmax = 80, ax = ax2)
			m.drawcoastlines(ax = ax2)
			m.fillcontinents(color='grey',ax = ax2)
			ax2.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
			ax2.set_title('SMAP')
			m.drawparallels(np.arange(int(y.min()),int(y.max()+1),1),labels=[0,1,0,0],ax = ax2)
			m.drawmeridians(np.arange(int(x.min()),int(x.max()+1),1),labels=[0,0,0,1],ax = ax2)

			im = m.pcolormesh(a2,b2,vent_sarmoy,cmap=cwnd,vmin = 0 , vmax = 80, ax = ax3)
			m.drawcoastlines(ax = ax3)
			m.fillcontinents(color='grey',ax = ax3)
			plt.colorbar(im, cax = ax, orientation='horizontal', label='Wind speed (m/s)')
			ax3.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
			ax3.set_title('Sentinel-1 -> SMAP')
			m.drawparallels(np.arange(int(y.min()),int(y.max()+1),1),labels=[0,1,0,0],ax = ax3)
			m.drawmeridians(np.arange(int(x.min()),int(x.max()+1),1),labels=[0,0,0,1],ax = ax3)

			im = m.pcolormesh(a2,b2,n_sarvalid,cmap='jet',vmin = 0 , vmax = 1000,ax = ax4)
			m.drawcoastlines(ax = ax4)
			m.fillcontinents(color='grey',ax = ax4)
			plt.colorbar(im, cax = ax5, orientation='horizontal', label='Number of points')
			ax4.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
			ax4.set_title('Number of points by pixel SMAP')
			m.drawparallels(np.arange(int(y.min()),int(y.max()+1),1),labels=[0,1,0,0],ax = ax4)
			m.drawmeridians(np.arange(int(x.min()),int(x.max()+1),1),labels=[0,0,0,1],ax = ax4)

			savefig('/net/merzhin1/vol/safo/morgane/data/figures/SAR_SMAP(REMSS)/SAR_SMAP_(REMSS)_'+name_cycl+'_'+day_sar+'_compare_pyre.png')
			plt.close('all')

					# Create netcdf

			day_sar1 = datetime.strptime(name_sar.split('-')[4][0:8],'%Y%m%d')
			day_sar1 = datetime.strftime(day_sar1, '%Y%m%d')
			day_sar1=str(day_sar1)
			file_net = nc.Dataset('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/'+name_cycl+'_'+date_smap+'_SAR_'+sat_sar+'_SMAP.nc', 'w', format='NETCDF4')

			fonc_create_netcdf.create_netcdf_smap(file_net,name_cycl, name_sar, name_smap, hour_sar_debut, hour_sar_fin, date_smap_moyenne, lon_smapss, lat_smapss, vent_sarmoy, vent_sarmed, vent_sarstd, ws_smapss, n_sarvalid, dt, angle_sarmoy,sat,cat,sat_sar,sat_rad  )

			file_net.close()
		
	except:
		print('Pas de colocalisation '+name_cycl+' '+day_sar)

	file_sar.close()
	file_smap_remss.close()
# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))

