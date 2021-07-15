# Headers


import h5py
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
palette = '/net/diva/home/mcms/guizioum/high_wind_speed.pal'
from palette import *
cwnd = getColorMap( rgbFile = palette )
from datetime import datetime
import glob
from glob import glob
import scipy
from scipy import stats
import time
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER



# Début du décompte du temps
start_time = time.time()


path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*_SMAP(JPL).nc')

for i,j in enumerate(path):
	print(i,j)
	file_stat = nc.Dataset(j,'r')
	file_name =  j.split('/')[-1]
	name_cycl = file_name.split('_')[0]
	date = datetime.strptime(file_name.split('_')[1],'%Y%m%dT%H%M%S')
	date1=date.strftime('%Y%m%dT%H%M%S')
	date = date.strftime('%Y/%m/%d %H:%M:%S')
	delta = np.squeeze(file_stat.variables['delta_time'][:])
	lon = np.squeeze(file_stat.variables['Lon'][:])
	lat = np.squeeze(file_stat.variables['Lat'][:])
	nb_points = np.squeeze(file_stat.variables['number_points'][:])
	vent_sar = np.squeeze(file_stat.variables['WindSpeed_sar_mean'][:])
	angle = np.squeeze(file_stat.variables['sar_mean_incidence_angle'][:])
	vent_diffu = np.squeeze(file_stat.variables['WindSpeed_smap_jpl'][:]) # SMAP JPL
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

		fig = figure(figsize = (20,10))
		ax = fig.add_axes((0.27,0.52,0.01,0.35)) 
		ax1 = fig.add_axes((0.1,0.52,0.15,0.35),projection=ccrs.PlateCarree())
		ax2 = fig.add_axes((0.30,0.52,0.15,0.35),projection=ccrs.PlateCarree())
		ax3 = fig.add_axes((0.48,0.52,0.15,0.35),projection=ccrs.PlateCarree())
		ax9 = fig.add_axes((0.65,0.52,0.01,0.35))
		ax4 = fig.add_axes((0.71,0.52,0.15,0.35),projection=ccrs.PlateCarree())
		ax8 = fig.add_axes((0.87,0.52,0.01,0.35))
		ax5 = fig.add_axes((0.1,0.1,0.25,0.35))
		ax6 = fig.add_axes((0.38,0.1,0.25,0.35))
		ax7 = fig.add_axes((0.66,0.1,0.25,0.35))

		fig.suptitle(name_cycl+' SMAP (JPL) '+hour_smap+'SAR '+hour_sar,fontsize = 18)
		im = ax1.pcolor(lon,lat,vent_sar,cmap=cwnd,vmin = 0 , vmax = 80)
		try:
			ax1.set_extent([lon.min()-1,lon.max()+1,lat.min()-1,lat.max()+1])
		except:
			pass
		ax1.coastlines('50m')
		plt.colorbar(im, cax = ax, orientation='vertical', label='Wind speed (m/s)')
		ax1.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
		ax1.set_title('Sentinel-1 -> SMAP', fontsize = 12)
		ax1.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor = 'black' , facecolor='grey'))
		ax1.gridlines(linewidth = 0.5 , linestyle='--')
		gl = ax1.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
		gl.xlabels_top = False
		gl.ylabels_left = False
		gl.ylabels_right=True
		gl.xlines = True
		gl.xformatter = LONGITUDE_FORMATTER
		gl.yformatter = LATITUDE_FORMATTER
		
	
		im = ax2.pcolor(lon,lat,vent_diffu,cmap=cwnd, vmin = 0 , vmax = 80)
		try:
			ax2.set_extent([lon.min()-1,lon.max()+1,lat.min()-1,lat.max()+1])
		except:
			pass
		ax2.coastlines('50m')
		ax2.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
		ax2.set_title('SMAP', fontsize = 12)
		ax2.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor = 'black' , facecolor='grey'))
		ax2.gridlines(linewidth = 0.5 , linestyle='--')
		gl = ax2.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
		ax2.text(0.5, 5, 'Longitude', va='bottom', ha='center', rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes)
		gl.xlabels_top = False
		gl.ylabels_left = False
		gl.ylabels_right=True
		gl.xlines = True
		gl.xformatter = LONGITUDE_FORMATTER
		gl.yformatter = LATITUDE_FORMATTER

	
		im = ax3.pcolor(lon,lat,diff,cmap='seismic', vmax=-np.min(diff))
		try:
			ax3.set_extent([lon.min()-1,lon.max()+1,lat.min()-1,lat.max()+1])
		except:
			pass
		ax3.coastlines('50m')
		plt.colorbar(im, cax = ax9, orientation='vertical', label='Wind speed difference (m/s)')
		ax3.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
		ax3.set_title('Difference between (SAR-SMAP(JPL))', fontsize = 12)
		ax3.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor = 'black' , facecolor='grey'))
		ax3.gridlines(linewidth = 0.5 , linestyle='--')
		gl = ax3.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
		gl.xlabels_top = False
		gl.ylabels_left = False
		gl.ylabels_right=True
		gl.xlines = True
		gl.xformatter = LONGITUDE_FORMATTER
		gl.yformatter = LATITUDE_FORMATTER

		im = ax4.pcolor(lon,lat,nb_points,cmap='jet' , vmax = 1000)
		try:
			ax4.set_extent([lon.min()-1,lon.max()+1,lat.min()-1,lat.max()+1])
		except:
			pass
		ax4.coastlines('50m')
		plt.colorbar(im, cax = ax8, orientation='vertical', label='Number of points')
		ax4.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
		ax4.set_title('Number of points by pixel SMAP', fontsize = 12)
		ax4.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor = 'black' , facecolor='grey'))
		ax4.gridlines(linewidth = 0.5 , linestyle='--')
		gl = ax4.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
		gl.xlabels_top = False
		gl.ylabels_left = False
		gl.ylabels_right=True
		gl.xlines = True
		gl.xformatter = LONGITUDE_FORMATTER
		gl.yformatter = LATITUDE_FORMATTER


		ax5.scatter(vent_diffu,vent_sar,alpha = 0.3)
		ax5.plot([0,100],[0,100],c='k')
		ax5.grid()
		ax5.set_xlim(0,70)
		ax5.set_ylim(0,70)
		ax5.set_xlabel('Wind speed SMAP (JPL) (m/s)', fontsize = 12)
		ax5.set_ylabel('Wind speed SAR (m/s)', fontsize = 12)
		ax5.set_title('Wind speed SMAP (JPL) vs wind speed SAR ', fontsize = 12)

		ax6.hist(diff.compressed(),bins=50)
		ax6.grid()
		ax6.set_xlim(-lim-1,lim+1)
		ax6.set_xlabel('Difference between wind speed', fontsize = 12)	
		ax6.set_ylabel('Number of points', fontsize = 12)

		ax7.scatter(angle,diff,alpha = 0.5)
		ax7.grid()
		ax7.set_xlabel('SAR mean angle incidence', fontsize = 12)
		ax7.set_ylabel('Difference between wind speed SAR & SMAP (JPL) (m/s)', fontsize = 12)
		ax7.set_title('Wind speed difference vs angle incidence', fontsize = 12)

		savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_RESUME/'+name_cycl+'_'+date1+'_SAR_SMAP(JPL)_stat.png')

		plt.close('all')
	except:
		pass


# Affichage du temps d execution
print("Temps d'exécution : %s secondes ---" % (time.time() - start_time))

