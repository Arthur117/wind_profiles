# Script permettant la comparaison des données vents entre Sentinel-1 et SMOS via netcdf intermediaire
# 01/04/2019 M.G

# Headers

import os
import conda
conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib
import h5py
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
palette = '/net/diva/home/mcms/guizioum/high_wind_speed.pal'
from palette import *
cwnd = getColorMap( rgbFile = palette )
from datetime import datetime
import glob
from glob import glob
import scipy
from scipy import stats
import time

# Debut du decompte du temps
start_time = time.time()

path = glob('/net/merzhin1/vol/safo/morgane/data/*/netcdf_intermediaire/*_SMAP_JPL.nc')

for i,j in enumerate(path[0:1]):
#	j='/net/merzhin1/vol/safo/morgane/data/IRMA/netcdf_intermediaire/IRMA_20170907_SAR_s1a_SMAP.nc'
	print(i,j)
	file_stat = nc.Dataset(j,'r')
	file_name =  j.split('/')[-1]
	name_cycl = j.split('/')[7]
	date = datetime.strptime(file_name.split('_')[1],'%Y%m%d')
	date1=date.strftime('%Y%m%d')
	date = date.strftime('%Y/%m/%d')
	delta = np.squeeze(file_stat.variables['delta_time'][:])
	lon = np.squeeze(file_stat.variables['Lon'][:])
	lat = np.squeeze(file_stat.variables['Lat'][:])
	nb_points = np.squeeze(file_stat.variables['number_points'][:])
	vent_sar = np.squeeze(file_stat.variables['WindSpeed_sar_mean'][:])
	angle = np.squeeze(file_stat.variables['sar_mean_incidence_angle'][:])
	hour_smap = file_stat.hour_smap
	hour_sar = file_stat.hour_sar
											# SMAP #

	vent_smap = np.squeeze(file_stat.variables['WindSpeed_smap_jpl'][:]) # SMAP
	vent_smap_a = np.ma.masked_where((vent_smap >= 999)*(vent_smap.mask == False)*(vent_sar.mask == False)+(vent_sar.mask == True)  , vent_smap)

										# MASK #

	nb_points = np.ma.masked_where(nb_points <= 600,nb_points)
	vent_sar.mask = nb_points.mask.copy()+vent_smap_a.mask.copy()
	nb_points.mask = nb_points.mask.copy()+vent_smap_a.mask.copy()
	angle.mask = nb_points.mask.copy()+vent_smap_a.mask.copy()
	vent_smap.mask = nb_points.mask.copy()+vent_smap_a.mask.copy()
	dim = vent_smap.shape
	dt = np.ma.zeros(dim)


#	if delta < 20000:
	for i in range(len(dt)):
		dt[i] = delta
	try:
		vent_sar_1 = np.ma.append(vent_sar_1,vent_sar)
		angle_1 = np.ma.append(angle_1,angle)
		vent_smap_1 = np.ma.append(vent_smap_1,vent_smap)
		dt_1 = np.ma.append(dt_1,dt)
		nb_1 = np.ma.append(nb_1,nb_points)
	except:
		vent_sar_1 = vent_sar
		angle_1 = angle
		vent_smap_1 = vent_smap
		dt_1 = dt
		nb_1 = nb_points
										# STATISTIQUES

	biais_sar = (vent_sar.mean()-vent_sar)/vent_sar
	biais_smap = (vent_smap.mean()-vent_smap)/vent_smap
	
	diff = vent_sar - vent_smap

	lim_max = abs(diff.max())
	lim_min = abs(diff.min())
	if lim_max < lim_min:
		lim = lim_min
	else:
		lim = lim_max


										# FIGURE 

	try:

		m = Basemap(llcrnrlon=lon.min(),llcrnrlat=lat.min(),urcrnrlon=lon.max(),urcrnrlat=lat.max(), 			resolution='i', projection='tmerc', lat_0 =lat.mean() , lon_0 =lon.mean() )

		x1,y1 = m(lon,lat)


		fig = figure(figsize = (15,7))
		fig.suptitle(name_cycl+' SMAP(JPL) '+hour_smap+'SAR '+hour_sar,fontsize = 18)

		ax = fig.add_axes((0.1,0.1,0.47,0.02)) 
		ax1 = fig.add_axes((0.1,0.22,0.22,0.65))
		ax2 = fig.add_axes((0.35,0.22,0.22,0.65))
		ax3 = fig.add_axes((0.66,0.22,0.24,0.65))
		ax4 = fig.add_axes((0.66,0.1,0.24,0.02))

		im = m.pcolormesh(x1,y1,vent_sar,cmap=cwnd,vmin = 0 , vmax = 80, ax = ax1)
		m.drawcoastlines(ax = ax1)
		ax1.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
		ax1.set_ylabel('Latitude', fontsize = 12)
		ax1.set_title('Sentinel-1')		
		m.fillcontinents(color='grey',ax = ax1)
		m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax1)
		m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax1)

		im = m.pcolormesh(x1,y1,vent_smap,cmap=cwnd, vmin = 0 , vmax = 80,ax = ax2)
		m.drawcoastlines(ax = ax2)
		m.fillcontinents(color='grey',ax = ax2)
		ax2.set_xlabel('Longitude', labelpad = 20, fontsize = 12)
		ax2.set_title('SMAP(JPL)')
		plt.colorbar(im, cax = ax, orientation='horizontal', label='Wind speed (m/s)',extend='both')
		m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax2)
		m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax2)

		im = m.pcolormesh(x1,y1,nb_points,cmap='jet', ax = ax3)
		m.drawcoastlines(ax = ax3)
		m.fillcontinents(color='grey',ax = ax3)
		ax3.set_xlabel('Longitude', labelpad = 20, fontsize = 12)
		ax3.set_title('Number of points by pixel SMAP(JPL)')
		plt.colorbar(im, cax = ax4, orientation='horizontal', label='Number of points',extend='both')
		m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax3)
		m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax3)
#		savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_RESUME/STATS_SAR_SMAP(JPL)/'+name_cycl+'_'+date1+'_SAR_SMAP(JPL)_compare(@SMAP(JPL)_res).png')


		fig_1 = figure(figsize = (10,7))
		fig_1.suptitle(name_cycl+' '+date, fontsize = 18)
		plot([0,100],[0,100],c='k')
		im = scatter(vent_smap,vent_sar,alpha = 0.3)
		plt.grid()
		xlim(0,70)
		ylim(0,70)
		xlabel('Vent SMAP(JPL) (m/s)', fontsize = 12)
		ylabel('Vent Sentinel-1 (m/s)', fontsize = 12)
#		savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_RESUME/STATS_SAR_SMAP(JPL)/'+name_cycl+'_'+date1+'_SAR_vs_SMAP(JPL).png')


		fig_2 = figure(figsize = (10,7))
		fig_2.suptitle(name_cycl+' '+date+' '+'Difference between Sentinel-1 and SMAP(JPL) wind speed', fontsize = 18)
		im = plt.hist(diff.compressed(),bins=50)
		xlim(-lim-1,lim+1)
		xlabel('Wind speed difference (SAR-SMAP(JPL)) (m/s)',fontsize=12)
		ylabel('Number of points',fontsize=12)
		plt.grid()
#		savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_RESUME/STATS_SAR_SMAP(JPL)/'+name_cycl+'_'+date1+'_SAR_SMAP(JPL)_WindSpeed_difference(hist).png')


		fig_3 = figure(figsize = (10,7))
		fig_3.suptitle(name_cycl+' '+date, fontsize = 18)
		im = scatter(angle,diff)
		plt.grid()
		plt.xlabel('SAR mean angle incidence', fontsize = 12)
		plt.ylabel('Difference between SAR & SMAP(JPL) (m/s)', fontsize = 12)
#		savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_RESUME/STATS_SAR_SMAP(JPL)/'+name_cycl+'_'+date1+'_SAR-SMAP(JPL)_vs_angle.png')


#		plt.close('all')
	except:
		pass


# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))






