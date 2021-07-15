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
import statistics

# Debut du decompte du temps
start_time = time.time()

path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*SMAP(JPL).nc')

for i,j in enumerate(path):
#j='/net/merzhin1/vol/safo/morgane/data/IRMA/netcdf_intermediaire/IRMA_20170907_SAR_s1a_SMAP.nc'
	print(i,j)
	file_stat = nc.Dataset(j,'r')
	file_name =  j.split('/')[-1]
	name_cycl = j.split('/')[7]
	date = datetime.strptime(file_name.split('_')[1],'%Y%m%dT%H%M%S')
	date1=date.strftime('%Y%m%d')
	date = date.strftime('%Y/%m/%d')
	delta = np.squeeze(file_stat.variables['delta_time'][:])
	lon = np.squeeze(file_stat.variables['Lon'][:])
	lat = np.squeeze(file_stat.variables['Lat'][:])
	nb_points = np.squeeze(file_stat.variables['number_points'][:])
	vent_sar = np.squeeze(file_stat.variables['WindSpeed_sar_mean'][:])
	angle = np.squeeze(file_stat.variables['sar_mean_incidence_angle'][:])

	
											# SMAP #

	vent_smap = np.squeeze(file_stat.variables['WindSpeed_smap_jpl'][:]) # SMAP
	vent_smap_a = np.ma.masked_where((vent_smap >= 999)*(vent_smap.mask == False)*(vent_sar.mask == False)+(vent_sar.mask == True) , vent_smap)

										# MASK #

	nb_points = np.ma.masked_where(nb_points <= 600,nb_points)
	vent_sar.mask = nb_points.mask.copy()+vent_smap_a.mask.copy()
	nb_points.mask = nb_points.mask.copy()+vent_smap_a.mask.copy()
	angle.mask = nb_points.mask.copy()+vent_smap_a.mask.copy()
	vent_smap.mask = nb_points.mask.copy()+vent_smap_a.mask.copy()
	dim = vent_smap.shape
	dt = np.ma.zeros(dim)


	if delta < 10000:
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


		angle_1 = np.ma.masked_where(angle_1>99 , angle_1)
		biais_sar = (vent_sar.mean()-vent_sar)/vent_sar
		biais_smap = (vent_smap.mean()-vent_smap)/vent_smap

		diff = vent_sar - vent_smap
		diff_smap = vent_sar_1 - vent_smap_1
		std_diff=diff.std()
		

		lim_max = abs(diff.max())
		lim_min = abs(diff.min())
		if lim_max < lim_min:
			lim = lim_min
		else:
			lim = lim_max

		R_1 = np.ma.corrcoef(vent_sar_1,vent_smap_1)
		diff_moy = diff_smap.mean()
		diff_med = np.ma.median(diff_smap)
		nb_sum = vent_smap_1.shape-vent_smap_1.mask.sum()
		nb_moy = nb_1.mean()

		nbins = 10
		x = angle_1.compressed()
		y = diff_smap.compressed()

		n, _ = np.histogram(x, bins=nbins)
		sy, _ = np.histogram(x, bins=nbins, weights=y)
		sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
		mean = sy / n
		std = np.sqrt(sy2/n - mean*mean)



fig_1 = figure(figsize = (10,7))
fig_1.suptitle('Wind speed function of mean sar angle incidence (dt<2h45min)', fontsize = 18)
im = scatter(x,y,marker='.',alpha=0.2) 
plt.xlabel('SAR mean angle incidence',fontsize=16)
plt.ylabel('Wind speed difference (SAR-SMAP(JPL)) (m/s)',fontsize=16)
plt.tick_params(axis = 'both', labelsize = 16)
plt.grid()
plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, fmt='r-',capthick = 2)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_GLOB/Wind_speed_SAR-SMAP(JPL)=f(SAR_mean_angle_incidence)_dt<10000_pyre.png')


fig_2 = figure(figsize = (10,7))
fig_2.suptitle('Wind speed (dt<2h45min)', fontsize = 18)
plt.plot([0,100],[0,100],c='k')
im = scatter(vent_smap_1,vent_sar_1,marker='.',alpha = 0.2, color = c)
plt.grid()
plt.text(1,65,'r² = '+str(round(R_1[0,1],3)),fontsize=16)
plt.text(1,61,'diff_moy = '+str(round(diff_moy,3))+' m/s',fontsize=16)
plt.text(1,57,'diff_med = '+str(round(diff_med,3))+' m/s',fontsize=16)
plt.text(1,53,'nb_tot = '+str(nb_sum),fontsize=16)
plt.text(1,49,'nb_moy = '+str(round(nb_moy,0)),fontsize=16)
#plt.text(1,44,'std = '+str(round(std_diff,2)),fontsize=16)
plt.xlabel('SMAP(JPL) Wind speed (m/s)',fontsize=16)
plt.ylabel('SAR Wind speed (m/s)',fontsize=16)
plt.tick_params(axis = 'both', labelsize = 16)
xlim(0,70)
ylim(0,70)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_GLOB/Wind_speed(SAR&SMAP(JPL))_dt<10000_pyre.png')

fig_3 = figure(figsize = (10,7))
fig_3.suptitle('Difference between Sentinel-1 and SMAP(JPL) wind speed (dt<2h45min)', fontsize = 18)
im = plt.hist(diff_smap.compressed(),bins=50)
plt.grid()
plt.xlabel('Wind speed difference',fontsize=16)
plt.xlim(-20,20)
plt.tick_params(axis = 'both', labelsize = 16)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_GLOB/Hist(SAR-SMAP(JPL))_dt<10000_pyre.png')


fig_4 = figure(figsize = (10,7))
fig_4.suptitle('Wind speed difference (SAR-SMAP(JPL)) vs delta T (dt<2h45min)', fontsize = 18)
im = scatter(dt_1,diff_smap,marker='.' )
plt.xlabel('Difference between time acquisition ( s )',fontsize=16)
plt.ylabel('Wind speed difference (SAR-SMAP(JPL)) (m/s)',fontsize=16)
plt.ylim(-20,20)
plt.grid()
plt.tick_params(axis = 'both', labelsize = 16)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMAP(JPL)_GLOB/Diff(SAR-SMAP_R)_vs_dt<10000_pyre.png')

#plt.close('all')



file_stat.close()

# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))
