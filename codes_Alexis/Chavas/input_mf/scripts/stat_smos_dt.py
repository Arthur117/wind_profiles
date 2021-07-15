# Script permettant la comparaison des données vents entre Sentinel-1 et SMOS via netcdf intermediaire
# 01/04/2019 M.G

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

# Debut du decompte du temps
start_time = time.time()


path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*_SMOS.nc')

for i,j in enumerate(path):
	print(i,j)
	file_stat = nc.Dataset(j,'r')
	delta = np.squeeze(file_stat.variables['delta_time'][:])
	lon = np.squeeze(file_stat.variables['Lon'][:])
	lat = np.squeeze(file_stat.variables['Lat'][:])
	nb_points = np.squeeze(file_stat.variables['number_points'][:])
	vent_sar = np.squeeze(file_stat.variables['WindSpeed_sar_mean'][:])
	vent_sar = np.ma.masked_where((vent_sar>= 999),vent_sar)
	angle = np.squeeze(file_stat.variables['sar_mean_incidence_angle'][:])
	vent_smos = np.squeeze(file_stat.variables['WindSpeed_smos'][:]) 
	vent_smos_a = np.ma.masked_where((vent_smos >= 999)*(vent_smos.mask == False)*(vent_sar.mask == False)+(vent_sar.mask == True) , vent_smos)
	
	




			# MASK #

	nb_points = np.ma.masked_where(nb_points <= 400,nb_points)
	vent_sar.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	nb_points.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	angle.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	vent_smos.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()

	dim = vent_smos.shape
	dt = np.ma.zeros(dim)

	if delta<10000:
		diff = vent_sar - vent_smos

		for i in range(len(dt)):
			dt[i] = delta

		try:
			vent_sar_1 = np.append(vent_sar_1,vent_sar)
			angle_1 = np.ma.append(angle_1,angle)
			vent_smos_1 = np.append(vent_smos_1,vent_smos)
			dt_1 = np.ma.append(dt_1,dt)
			nb_1 = np.ma.append(nb_1,nb_points)
		except:
			vent_sar_1 = vent_sar
			angle_1 = angle
			vent_smos_1 = vent_smos
			dt_1 = dt
			nb_1 = nb_points	


					# STATISTIQUES #
		vent_sar_1=np.ma.masked_where((vent_sar_1 >=999),vent_sar_1)
		vent_smos_1=np.ma.masked_where((vent_smos_1 >=999),vent_smos_1)
		diff_smos = vent_sar_1 - vent_smos_1
		
		if vent_sar_1<30:
			diff_smos_inf=vent_sar_1[vent_sar_1<30]-vent_smos_1[vent_sar_1<30]


		R_1 = np.ma.corrcoef(vent_sar_1,vent_smos_1)
		diff_moy = diff_smos.mean()
		diff_med = np.ma.median(diff_smos)
		nb_sum = vent_smos_1.shape-vent_smos_1.mask.sum()
		nb_moy = nb_1.mean()
		R = R_1[0,1]
	
		lim_max = abs(diff.max())
		lim_min = abs(diff.min())
		if lim_max < lim_min:
			lim = lim_min
		else:
			lim = lim_max

		nbins = 10
		x = nan_to_num(angle_1)
		y = nan_to_num(diff_smos)
		n, _ = np.histogram(x, bins=nbins)
		sy, _ = np.histogram(x, bins=nbins, weights=y)
		sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
		mean = sy / n
		std = np.sqrt(sy2/n - mean*mean)


fig_1 = figure(figsize = (10,7))
fig_1.suptitle('Wind speed vs mean sar angle incidence dt<2h45', fontsize = 18)
im = scatter(x,y,marker='.',alpha=0.2,s=50)
plt.xlabel('SAR mean angle incidence (degrees)',fontsize=16)
plt.ylabel('Wind speed difference (SAR-SMOS) (m/s)',fontsize=16)
plt.grid()
xlim(0,60)
plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, fmt='r-',capthick = 2)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMOS_GLOB/Wind_speed_SAR-SMOS=f(SAR_mean_angle_incidence)_dt<10000.png')

fig_2 = figure(figsize = (10,7))
fig_2.suptitle('Wind speed dt<2h45', fontsize = 18)
plt.plot([0,100],[0,100],c='k')
im = scatter(vent_smos_1,vent_sar_1,marker='.',alpha=0.2,s=50 )
plt.grid()
plt.text(1,60,'r² = '+str(round(R,3)),fontsize=16)
plt.text(1,57,'diff_moy = '+str(round(diff_moy,3))+' (m/s)',fontsize=16)
plt.text(1,54,'diff_med = '+str(round(diff_med,3))+' (m/s)',fontsize=16)
plt.text(1,51,'nb_tot = '+str(nb_sum),fontsize=16)
plt.text(1,48,'nb_moy = '+str(round(nb_moy,0)),fontsize=16)
xlim(0,70)
ylim(0,70)
plt.xlabel('SMOS Wind speed (m/s)',fontsize=16)
plt.ylabel('SAR Wind speed (m/s)',fontsize=16)
plt.tick_params(axis = 'both', labelsize = 16)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMOS_GLOB/Wind_speed(SAR&SMOS)_dt<10000.png')

fig_3 = figure(figsize = (10,7))
fig_3.suptitle('Difference between Sentinel-1 and SMOS wind speed dt<2h45', fontsize = 18)
im = plt.hist(diff_smos.compressed(),bins=50)
plt.grid()
plt.xlabel('Wind speed difference (m/s)',fontsize=16)
plt.ylabel('Number of points',fontsize=16)
plt.tick_params(axis = 'both', labelsize = 16)
plt.xlim(-40,40,5)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMOS_GLOB/Hist(SAR-SMOS)_dt<10000.png')

fig_4 = figure(figsize = (10,7))
fig_4.suptitle('Wind speed difference (SAR-SMOS) vs delta T dt<2h45', fontsize = 18)
im = scatter(dt_1,diff_smos,marker='.' )
plt.grid()
plt.tick_params(axis = 'both', labelsize = 16)
plt.xlabel('Delta t between time acquisition (s)',fontsize=16)
plt.ylabel('Wind speed difference (SAR-SMOS) (m/s)',fontsize=16)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_SAR_SMOS_GLOB/Diff(SAR-SMOS)_vs_dt<10000.png')
file_stat.close()
# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))


