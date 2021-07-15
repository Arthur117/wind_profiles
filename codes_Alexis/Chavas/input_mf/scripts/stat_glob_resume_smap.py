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

# Début du décompte du temps
start_time = time.time()

path = glob('/net/merzhin1/vol/safo/morgane/data/*/netcdf_intermediaire/*_SMAP.nc')

for i,j in enumerate(path):

	print(i,j)
	file_stat = nc.Dataset(j,'r')
	file_name =  j.split('/')[-1]
	name_cycl = j.split('/')[7]
	date = datetime.strptime(file_name.split('_')[1],'%Y%m%d')
	date1=date.strftime('%Y%m%dT')
	date = date.strftime('%Y/%m/%d ')
	delta = np.squeeze(file_stat.variables['delta_time'][:])
	lon = np.squeeze(file_stat.variables['Lon'][:])
	lat = np.squeeze(file_stat.variables['Lat'][:])
	nb_points = np.squeeze(file_stat.variables['number_points'][:])
	vent_sar = np.squeeze(file_stat.variables['WindSpeed_sar_mean'][:])
	angle = np.squeeze(file_stat.variables['sar_mean_incidence_angle'][:])
	vent_diffu = np.squeeze(file_stat.variables['WindSpeed_smap'][:]) # SMAP
	hour_smap = file_stat.hour_smap
	hour_sar = file_stat.hour_sar

	
	vent_diffu_a = np.ma.masked_where((vent_diffu >= 999)*(vent_diffu.mask == False)*(vent_sar.mask == False)+(vent_sar.mask == True) , vent_diffu)

	nb_points = np.ma.masked_where(nb_points <= 600,nb_points)
	nb_points.mask = nb_points.mask.copy()+vent_diffu_a.mask.copy()
	vent_sar.mask = nb_points.mask.copy()+vent_diffu_a.mask.copy()
	angle.mask = nb_points.mask.copy()+vent_diffu_a.mask.copy()
	vent_diffu.mask = nb_points.mask.copy()+vent_diffu_a.mask.copy()

	diff = vent_sar - vent_diffu
	biais_sar = (vent_sar.mean()-vent_sar)/vent_sar
	biais_diff = (vent_diffu.mean()-vent_diffu)/vent_diffu

	lim_max = abs(diff.max())
	lim_min = abs(diff.min())
	if lim_max < lim_min:
		lim = lim_min
	else:
		lim = lim_max

	R = np.ma.corrcoef(vent_diffu,vent_sar)
	diff_moy = diff.mean()
	diff_med = np.ma.median(diff)
	nb_sum = vent_diffu.shape-vent_diffu.mask.sum()
	nb_moy = nb_points.mean()

	try:

		m = Basemap(llcrnrlon=lon.min(),llcrnrlat=lat.min(),urcrnrlon=lon.max(),urcrnrlat=lat.max(), 			resolution='i', projection='tmerc', lat_0 =lat.mean() , lon_0 =lon.mean() )

		lon, lat = np.meshgrid(lon,lat)
		x1,y1 = m(lon,lat)


		fig = figure(figsize = (20,10))
		fig.suptitle(name_cycl+' SMAP '+hour_smap+'SAR '+hour_sar,fontsize = 18)

		ax = fig.add_axes((0.27,0.52,0.01,0.35)) 
		ax1 = fig.add_axes((0.1,0.52,0.15,0.35))
		ax2 = fig.add_axes((0.30,0.52,0.15,0.35))
		ax3 = fig.add_axes((0.48,0.52,0.15,0.35))
		ax9 = fig.add_axes((0.65,0.52,0.01,0.35))
		ax4 = fig.add_axes((0.71,0.52,0.15,0.35))
		ax8 = fig.add_axes((0.87,0.52,0.01,0.35))
		ax5 = fig.add_axes((0.1,0.1,0.25,0.35))
		ax6 = fig.add_axes((0.38,0.1,0.25,0.35))
		ax7 = fig.add_axes((0.66,0.1,0.25,0.35))

		im = m.pcolormesh(x1,y1,vent_sar,cmap=cwnd,vmin = 0 , vmax = 80, ax = ax1)
		m.drawcoastlines(ax = ax1)
		ax1.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
		ax1.set_ylabel('Latitude', fontsize = 12)
		ax1.set_title('Sentinel-1', fontsize = 12)		
		m.fillcontinents(color='grey',ax = ax1)
		m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax1)
		m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax1)

		im = m.pcolormesh(x1,y1,vent_diffu,cmap=cwnd, vmin = 0 , vmax = 80,ax = ax2)
		m.drawcoastlines(ax = ax2)
		m.fillcontinents(color='grey',ax = ax2)
		ax2.set_xlabel('Longitude', labelpad = 20, fontsize = 12)
		ax2.set_title('SMAP(REMSS)', fontsize = 12)
		plt.colorbar(im, cax = ax, orientation='vertical', label='Wind speed (m/s)')
		m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax2)
		m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax2)

		im = m.pcolormesh(x1,y1,diff,cmap='seismic',ax = ax3, vmax=-np.min(diff))
		m.drawcoastlines(ax = ax3)
		m.fillcontinents(color='grey',ax = ax3)
		ax3.set_xlabel('Longitude', labelpad = 20, fontsize = 12)
		ax3.set_title('Difference (SAR-SMAP(REMSS)) (m/s)', fontsize = 12)
		plt.colorbar(im, cax = ax9, orientation='vertical', label='SAR-SMAP(REMSS) (m/s)')
		m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax3)
		m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax3)

		im = m.pcolormesh(x1,y1,nb_points,cmap='jet', ax = ax4)
		m.drawcoastlines(ax = ax4)
		m.fillcontinents(color='grey',ax = ax4)
		ax4.set_xlabel('Longitude', labelpad = 20, fontsize = 12)
		ax4.set_title('Number of points by pixel SMAP(REMSS)', fontsize = 12)
		plt.colorbar(im, cax = ax8, orientation='vertical', label='Number of points')
		m.drawparallels(np.arange(int(lat.min()),int(lat.max()+1),1),labels=[0,1,0,0],ax = ax4)
		m.drawmeridians(np.arange(int(lon.min()),int(lon.max()+1),1),labels=[0,0,0,1],ax = ax4)

		ax5.scatter(vent_diffu,vent_sar,alpha = 0.3)
		ax5.plot([0,100],[0,100],c='k')
		ax5.grid()
		ax5.set_xlim(0,70)
		ax5.set_ylim(0,70)
		ax5.set_xlabel('Wind speed SMAP(REMSS) (m/s)', fontsize = 12)
		ax5.set_ylabel('Wind speed SAR (m/s)', fontsize = 12)
		ax5.set_title('Wind speed SMAP(REMSS) vs wind speed SAR ', fontsize = 12)

		ax6.hist(diff.compressed(),bins=50)
		ax6.grid()
		ax6.set_xlim(-lim-1,lim+1)
		ax6.set_xlabel('Difference between wind speed (m/s)', fontsize = 12)	
		ax6.set_ylabel('Number of points', fontsize = 12)

		ax7.scatter(angle,diff,alpha = 0.5)
		ax7.grid()
		ax7.set_xlabel('SAR mean angle incidence (degrees)', fontsize = 12)
		ax7.set_ylabel('Difference between wind speed SAR & SMAP(REMSS) (m/s)', fontsize = 12)
		ax7.set_title('Wind speed difference vs angle incidence', fontsize = 12)

		savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(REMSS)_RESUME/'+name_cycl+'_'+date1+'_SAR_SMAP(REMSS)_stat_pyre.png')

		plt.close('all')
	except:
		pass


# Affichage du temps d execution
print("Temps d'exécution : %s secondes ---" % (time.time() - start_time))

