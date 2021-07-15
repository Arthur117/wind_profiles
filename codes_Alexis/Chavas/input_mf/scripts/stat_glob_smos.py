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
import time

# Debut du decompte du temps
start_time = time.time()


path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/LESTER*20160831*_SMOS.nc')

for i,j in enumerate(path):
	print(i,j)
	file_stat = nc.Dataset(j,'r')
	file_name =  j.split('/')[-1]
	name_cycl = file_name.split('_')[0]
	date = datetime.strptime(file_name.split('_')[1],'%Y%m%dT%H%M%S')
	date1 = date.strftime('%Y%m%d')
	date = date.strftime('%Y/%m/%d')
	delta = np.squeeze(file_stat.variables['delta_time'][:])
	lon = np.squeeze(file_stat.variables['Lon'][:])
	lat = np.squeeze(file_stat.variables['Lat'][:])
	nb_points = np.squeeze(file_stat.variables['number_points'][:])
	vent_sar = np.squeeze(file_stat.variables['WindSpeed_sar_mean'][:])
	angle = np.squeeze(file_stat.variables['sar_mean_incidence_angle'][:])
	hour_smos = file_stat.hour_smos
	hour_sar = file_stat.hour_sar
			# SMOS #

	vent_smos = np.squeeze(file_stat.variables['WindSpeed_smos'][:]) 
	vent_smos_a = np.ma.masked_where((vent_smos >= 999)*(vent_smos.mask == False)*(vent_sar.mask == False) , vent_smos)

			# MASK #
	
	nb_points = np.ma.masked_where(nb_points <= 300,nb_points)
	vent_sar.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	nb_points.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	angle.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	vent_smos.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	dim = vent_smos.shape
	dt = np.ma.zeros(dim)

#	if delta<10000:

	for i in range(len(dt)):
		dt[i] = delta

			# STATISTIQUES

	biais_sar = (vent_sar.mean()-vent_sar)/vent_sar
	biais_smos = (vent_smos.mean()-vent_smos)/vent_smos
	diff = vent_sar - vent_smos

	lim_max = abs(diff.max())
	lim_min = abs(diff.min())
	if lim_max < lim_min:
		lim = lim_min
	else:
		lim = lim_max

											# FIGURE 

	m = Basemap(llcrnrlon=lon.min(),llcrnrlat=lat.min(),urcrnrlon=lon.max(),urcrnrlat=lat.max(), 			resolution='i', projection='tmerc', lat_0 =lat.mean() , lon_0 =lon.mean() )
	lon, lat = np.meshgrid(lon,lat)
	x1,y1 = m(lon,lat)

	fig = figure(figsize = (15,7))
	fig.suptitle(name_cycl+' SMOS '+hour_smos+' SAR '+hour_sar,fontsize = 18)

	ax = fig.add_axes((0.1,0.1,0.47,0.02)) 
	ax1 = fig.add_axes((0.1,0.22,0.22,0.65))
	ax2 = fig.add_axes((0.35,0.22,0.22,0.65))
	ax3 = fig.add_axes((0.66,0.22,0.24,0.65))
	ax4 = fig.add_axes((0.66,0.1,0.24,0.02))

	im = m.pcolormesh(x1,y1,vent_sar,cmap=cwnd,vmin = 0 , vmax = 80, ax = ax1)
	m.drawcoastlines(ax = ax1)
	ax1.set_xlabel('Longitude',labelpad = 20,fontsize=14)
	ax1.set_ylabel('Latitude',fontsize=14)
	ax1.set_title('Sentinel-1')		
	m.fillcontinents(color='grey',ax = ax1)
	m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax1)
	m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax1)

	im = m.pcolormesh(x1,y1,vent_smos,cmap=cwnd, vmin = 0 , vmax = 80,ax = ax2)
	m.drawcoastlines(ax = ax2)
	m.fillcontinents(color='grey',ax = ax2)
	ax2.set_xlabel('Longitude', labelpad = 20,fontsize=14)
	ax2.set_title('SMOS',fontsize=14)
	plt.colorbar(im, cax = ax, orientation='horizontal', label='Wind speed (m/s)')
	m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax2)
	m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax2)

	im = m.pcolormesh(x1,y1,nb_points,cmap='jet', ax = ax3)
	m.drawcoastlines(ax = ax3)
	m.fillcontinents(color='grey',ax = ax3)
	ax3.set_xlabel('Longitude', labelpad = 20,fontsize=14)
	ax3.set_title('Number of points by pixel SMOS',fontsize=14)
	plt.colorbar(im, cax = ax4, orientation='horizontal', label='Number of points')
	m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax3)
	m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax3)
	savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMOS_RESUME/STATS_SAR_SMOS/'+name_cycl+'_'+date1+'_SAR_SMOS_compare(@SMOS_res).png')

	fig_1 = figure(figsize = (15,9))
	fig_1.suptitle(name_cycl+' SMOS '+hour_smos+' SAR '+hour_sar,fontsize = 20)
	plot([0,100],[0,100],c='k')
	im = scatter(vent_smos,vent_sar,alpha = 0.3)
	plt.grid()
	plt.tick_params(axis = 'both', labelsize = 16)
	xlim(0,70)
	ylim(0,70)
	xlabel('Wind SMOS (m/s)',fontsize=16)
	ylabel('Wind Sentinel-1 (m/s)',fontsize=16)
	savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMOS_RESUME/STATS_SAR_SMOS/'+name_cycl+'_'+date1+'_SAR_vs_SMOS.png')

	fig_2 = figure(figsize = (15,9))
	fig_2.suptitle(name_cycl+' SMOS '+hour_smos+' SAR '+hour_sar,fontsize = 20)
	im = plt.hist(diff.compressed(),bins = 50)
	xlabel('Wind speed difference (SAR-SMOS) (m/s)',fontsize=16)
	ylabel('Number of points',fontsize=16)
	plt.tick_params(axis = 'both', labelsize = 16)
	xlim(-30,30,5)
	plt.grid()
	savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMOS_RESUME/STATS_SAR_SMOS/'+name_cycl+'_'+date1+'_SAR_SMOS_WindSpeed_difference(hist).png')

	fig_3 = figure(figsize = (15,9))
	fig_3.suptitle(name_cycl+' SMOS '+hour_smos+' SAR '+hour_sar,fontsize = 20)
	im = scatter(angle,diff)
	plt.tick_params(axis = 'both', labelsize = 16)
	plt.grid()
	plt.xlabel('SAR mean angle incidence',fontsize=16)
	plt.ylabel('Difference between SAR & SMOS (m/s)',fontsize=16)
	savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMOS_RESUME/STATS_SAR_SMOS/'+name_cycl+'_'+date1+'_SAR-SMOS_vs_angle.png')

	#plt.close('all')


# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))


