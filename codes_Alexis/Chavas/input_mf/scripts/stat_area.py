# Script permettant la comparaison des données vents entre Sentinel-1 et SMOS via netcdf intermediaire par zone
# 21/05/2019 M.G

# Headers


import h5py
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from datetime import datetime
import glob
from glob import glob
import scipy
from scipy import stats
import time
from scipy.stats import binned_statistic
import fonc_poly
import pandas as pd

# Debut du decompte du temps
start_time = time.time()


df=pd.read_csv('/net/merzhin1/vol/safo/morgane/data/cycl_area.txt')
name = df.iloc[:,0]
date = df.iloc[:,1]
hour = df.iloc[:,2]
rad = df.iloc[:,4]
area = df.iloc[:,7]
sen = df.iloc[:,3]

a1 = np.where(area == '01')
a2 = np.where(area == '02')
a3 = np.where(area == '03')
a4 = np.where(area == '04')
a5 = np.where(area == '05')
a6 = np.where(area == '06')
a7 = np.where(area == '07')
a11= np.where((area == '11')|(area == '11b'))
a12 = np.where((area == '12')|(area == '12b'))
a0 = np.where((area != '01')&(area != '02')&(area != '03')&(area != '04')&(area != '05')&(area != '06')&(area != '07')&(area != '11')&(area != '11b')&(area != '12')&(area != '12b'))

name1 = [''.join(map(str, l)) for l in name[a1[0]]]
for i,j in enumerate(np.unique(name1)):
	path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/'+j.split(' ')[-1]+'_*_SMOS.nc')
	for k,l in enumerate(path):
		print(k,l)
		area_fix = 'area_6'
		rad_fix = 'SMAP'
		sen_fix = 'SAR (Sentinel 1)'
		file_stat = nc.Dataset(l,'r')
		delta = np.squeeze(file_stat.variables['delta_time'][:])
		lon = np.squeeze(file_stat.variables['Lon'][:])
		lat = np.squeeze(file_stat.variables['Lat'][:])
		nb_points = np.squeeze(file_stat.variables['number_points'][:])
		vent_sar = np.squeeze(file_stat.variables['WindSpeed_sar_mean'][:])
		angle = np.squeeze(file_stat.variables['sar_mean_incidence_angle'][:])
						
		vent_smos = np.squeeze(file_stat.variables['WindSpeed_smos'][:])
		vent_smos_a = np.ma.masked_where((vent_smos >= 999)*(vent_smos.mask == False)*(vent_sar.mask == False)+(vent_sar.mask == True)  , vent_smos)

											# MASK #
		nb_points = np.ma.masked_where(nb_points <= 400,nb_points)
		vent_sar.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
		print(vent_sar.min())
		nb_points.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
		angle.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
		vent_smos.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
		dim = vent_smos.shape
		dt = np.ma.zeros(dim)

		if delta < 10000:
			for i in range(len(dt)):
				dt[i] = delta

			try:
				vent_sar_1 = np.ma.append(vent_sar_1,vent_sar)
				angle_1 = np.ma.append(angle_1,angle)
				vent_smos_1 = np.ma.append(vent_smos_1,vent_smos)
				dt_1 = np.ma.append(dt_1,dt)
				nb_1 = np.ma.append(nb_1,nb_points)
			except:
				vent_sar_1 = vent_sar
				angle_1 = angle
				vent_smos_1 = vent_smos
				dt_1 = dt
				nb_1 = nb_points
												# STATISTIQUES

			diff = vent_sar - vent_smos
			diff_smos = vent_sar_1 - vent_smos_1
			R= np.corrcoef(vent_sar,vent_smos)
			R_1 = np.corrcoef(vent_sar_1.compressed(),vent_smos_1.compressed())
			diff_moy = diff_smos.mean()
			std_diff = diff_smos.std()
			diff_med = np.ma.median(diff_smos)
			nb_sum = vent_smos_1.shape-vent_smos_1.mask.sum()
			nb_moy = nb_1.mean()

			lim_max = abs(diff.max())
			lim_min = abs(diff.min())
			if lim_max < lim_min:
				lim = lim_min
			else:
				lim = lim_max

			try:
				nbins = 10
				x = angle_1.compressed()
				y = diff_smos.compressed()
				n, _ = np.histogram(x, bins=nbins)
				sy, _ = np.histogram(x, bins=nbins, weights=y)
				sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
				mean = sy / n
				std = np.sqrt(sy2/n - mean*mean)
			except:
				pass
try:
	fig_1 = figure(figsize = (10,7))
	fig_1.suptitle('Wind speed function of mean '+sen_fix+' angle incidence '+area_fix, fontsize = 18)
	im = scatter(x,y,marker='.',alpha=0.2) 
	plt.xlabel(sen_fix+' mean angle incidence',fontsize=12)
	plt.ylabel('Wind speed difference ('+sen_fix+'-'+rad_fix+') (m/s)',fontsize=16)
	plt.grid()
	plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, fmt='r-',capthick = 2)
	plt.tick_params(axis = 'both', labelsize = 16)
#	savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/Wind_speed_(SAR-'+rad_fix+')_VS_SAR_mean_angle_incidence)_'+area_fix+'.png')
except:
	pass

fig_2 = figure(figsize = (10,7))
fig_2.suptitle('Wind speed '+area_fix, fontsize = 18)
plt.plot([0,100],[0,100],c='k')
im = scatter(vent_sar_1,vent_smos_1,marker='.',alpha = 0.2,s=25)
plt.grid()
try:
	plt.text(1,65,'r² = '+str(round(R_1[0,1],3)),fontsize=16)
	plt.text(1,61,'diff_moy = '+str(round(diff_moy,3)),fontsize=16)
	plt.text(1,57,'diff_med = '+str(round(diff_med,3)),fontsize=16)
	plt.text(1,53,'nb_tot = '+str(nb_sum),fontsize=16)
	plt.text(1,49,'nb_moy = '+str(round(nb_moy,0)),fontsize=16)
except:
	pass
xlim(0,70)
ylim(0,70)
plt.xlabel(''+rad_fix+' Wind speed (m/s)',fontsize=16)
plt.ylabel(sen_fix+' Wind speed (m/s)',fontsize=16)
plt.tick_params(axis = 'both', labelsize = 16)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/Wind_speed(SAR&'+rad_fix+')_'+area_fix+'.png')

fig_3 = figure(figsize = (10,7))
fig_3.suptitle('Difference between '+sen_fix+' and '+rad_fix+' wind speed '+area_fix, fontsize = 18)
im = plt.hist(diff_smos.compressed(),bins=50)
plt.grid()
plt.xlabel('Wind speed difference (m/s)',fontsize=16)
plt.tick_params(axis = 'both', labelsize = 16)
plt.xlim(-20,20)
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/Hist(SAR-'+rad_fix+')_'+area_fix+'.png')

fig_4 = figure(figsize = (10,7))
fig_4.suptitle('Wind speed difference ('+sen_fix+'-'+rad_fix+') vs delta T '+area_fix, fontsize = 18)
im = scatter(dt_1,diff_smos,marker='.' )
plt.xlabel('Difference between time acquisition ( s )',fontsize=16)
plt.ylabel('Wind speed difference ('+sen_fix+'-'+rad_fix+') (m/s)',fontsize=16)
plt.tick_params(axis = 'both', labelsize = 16)
plt.grid()
#savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/Diff(SAR-'+rad_fix+')_vs_dt_'+area_fix+'.png')
file_stat.close()
#plt.close('all')
print(R_1,diff_moy,std_diff)
# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))
