# Script permettant la comparaison des données vents entre Sentinel-1 et SMOS/SMAP via netcdf intermediaire par zone
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
a11 = np.where((area == '11')|(area == '11b'))
a12 = np.where((area == '12')|(area == '12b'))
a0 = np.where((area != '01')&(area != '02')&(area != '03')&(area != '04')&(area != '05')&(area != '06')&(area != '07')&(area != '11')&(area != '11b')&(area != '12')&(area != '12b'))

sum1=0
sum2=0
sum3=0
c='all'
a='all'

name1 = [''.join(map(str, l)) for l in name[a6[0]]] # a6 à modifier pour chaque zone citée plus haut
name1 = np.unique(name1)
#for i,j in enumerate(name1):  # à decommenter pour faire par zone
	#print(i,j)
#	path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/'+j.split(' ')[-1]+'_*_SMOS.nc')     #  à decommenter pour faire par zone
path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*_SMOS.nc')                     # à commenter pour faire par zone
for k,l in enumerate(path):
	#print(k,l)
	area_fix = 'area_1'
	rad_fix = 'SMOS'
	sen_fix = 'SAR (Sentinel 1)'
	file_stat = nc.Dataset(l,'r')
	delta = np.squeeze(file_stat.variables['delta_time'][:])
	lon = np.squeeze(file_stat.variables['Lon'][:])
	lat = np.squeeze(file_stat.variables['Lat'][:])
	nb_points = np.squeeze(file_stat.variables['number_points'][:])
	vent_sar = np.squeeze(file_stat.variables['WindSpeed_sar_mean'][:])
	angle = np.squeeze(file_stat.variables['sar_mean_incidence_angle'][:])
	cat=file_stat.category
	vent_smos = np.squeeze(file_stat.variables['WindSpeed_smos'][:])


	vent_smos_a = np.ma.masked_where((vent_smos >= 999)*(vent_smos.mask == False)*(vent_sar.mask == False)+(vent_sar.mask == True)  , vent_smos)

										# MASK #
	nb_points = np.ma.masked_where(nb_points <= 400,nb_points)
	vent_sar.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	nb_points.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	angle.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	vent_smos.mask = nb_points.mask.copy()+vent_smos_a.mask.copy()
	dim = vent_smos.shape
	dt = np.ma.zeros(dim)

	if delta < 10000:
		for i in range(len(dt)):
			dt[i] = delta

		#if cat==c:
		sum1=sum1+1	
		try:
			vent_sar_1 = np.ma.append(vent_sar_1,vent_sar)
			print(vent_sar_1.shape)
			angle_1 = np.ma.append(angle_1,angle)
			vent_smos_1 = np.ma.append(vent_smos_1,vent_smos)
			print(vent_smos_1.shape)
			dt_1 = np.ma.append(dt_1,dt)
			nb_1 = np.ma.append(nb_1,nb_points)
		except:
			vent_sar_1 = vent_sar
			print(vent_sar_1.shape)
			angle_1 = angle
			vent_smos_1 = vent_smos
			print(vent_smos_1.shape)
			dt_1 = dt
			nb_1 = nb_points
												# STATISTIQUES

try:
	vent_sar_1=np.ma.masked_invalid(vent_sar_1)
	diff_smos = vent_sar_1 - vent_smos_1
	print('smos ',diff_smos.shape)
	diff_moy = diff_smos.mean()
	std_diff = diff_smos.std()
	diff_med = np.ma.median(diff_smos)

	nbins = 10
	x = vent_sar_1.compressed()
	y = diff_smos.compressed()
	n, _ = np.histogram(x, bins=nbins)
	sy, _ = np.histogram(x, bins=nbins, weights=y)
	sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)
except:
	pass


####################################################################################################
#												SMAP(RSS)										   #
####################################################################################################


#for i,j in enumerate(name1):
#	path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/'+j.split(' ')[-1]+'_*_SMAP.nc')            #  à decommenter pour faire par zone
path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*_SMAP.nc')        # à commenter pour faire par zone
for k,l in enumerate(path):
	#print(k,l)
	area_fix_smap = 'area_1'
	rad_fix_smap = 'SMAP(RSS)'
	sen_fix_smap = 'SAR (Sentinel 1)'
	file_stat_smap = nc.Dataset(l,'r')
	delta_smap = np.squeeze(file_stat_smap.variables['delta_time'][:])
	lon_smap = np.squeeze(file_stat_smap.variables['Lon'][:])
	lat_smap = np.squeeze(file_stat_smap.variables['Lat'][:])
	nb_points_smap = np.squeeze(file_stat_smap.variables['number_points'][:])
	vent_sar_smap = np.squeeze(file_stat_smap.variables['WindSpeed_sar_mean'][:])
	angle_smap = np.squeeze(file_stat_smap.variables['sar_mean_incidence_angle'][:])
	cat_smap=file_stat_smap.category
	vent_smap = np.squeeze(file_stat_smap.variables['WindSpeed_smap'][:])


	vent_smap_a = np.ma.masked_where((vent_smap >= 999)*(vent_smap.mask == False)*(vent_sar_smap.mask == False)+(vent_sar_smap.mask == True)  , vent_smap)

										# MASK #
	nb_points_smap = np.ma.masked_where(nb_points_smap <= 400,nb_points_smap)
	vent_sar_smap.mask = nb_points_smap.mask.copy()+vent_smap_a.mask.copy()
	nb_points_smap.mask = nb_points_smap.mask.copy()+vent_smap_a.mask.copy()
	angle_smap.mask = nb_points_smap.mask.copy()+vent_smap_a.mask.copy()
	vent_smap.mask = nb_points_smap.mask.copy()+vent_smap_a.mask.copy()
	dim_smap = vent_smap.shape
	dt_smap = np.ma.zeros(dim_smap)

	if delta_smap < 10000:
		for i in range(len(dt_smap)):
			dt_smap[i] = delta_smap

		#if cat_smap==c:
		sum2=sum2+1	
		try:
			vent_sar_1_smap = np.ma.append(vent_sar_1_smap,vent_sar_smap)
			angle_1_smap = np.ma.append(angle_1_smap,angle_smap)
			vent_smap_1= np.ma.append(vent_smap_1,vent_smap)
			dt_1_smap = np.ma.append(dt_1_smap,dt_smap)
			nb_1_smap = np.ma.append(nb_1_smap,nb_points_smap)
		except:
			vent_sar_1_smap = vent_sar_smap
			angle_1_smap = angle_smap
			vent_smap_1 = vent_smap
			dt_1_smap = dt_smap
			nb_1_smap = nb_points_smap
													# STATISTIQUES
try:
	vent_sar_1_smap=np.ma.masked_invalid(vent_sar_1_smap)
	diff_smap = vent_sar_smap - vent_smap
	diff_smap_smap = vent_sar_1_smap - vent_smap_1
	print('smap ',diff_smap_smap.shape)
	diff_moy_smap = diff_smap_smap.mean()
	std_diff_smap = diff_smap_smap.std()
	diff_med_smap = np.ma.median(diff_smap_smap)

	nbins_smap = 10
	x_smap = vent_sar_1_smap.compressed()
	y_smap = diff_smap_smap.compressed()
	n_smap, __smap = np.histogram(x_smap, bins=nbins_smap)
	sy_smap, __smap = np.histogram(x_smap, bins=nbins_smap, weights=y_smap)
	sy2_smap, __smap = np.histogram(x_smap, bins=nbins_smap, weights=y_smap*y_smap)
	mean_smap = sy_smap / n_smap
	std_smap = np.sqrt(sy2_smap/n_smap - mean_smap*mean_smap)
except:
	pass


####################################################################################################
#												SMAP(JPL)										   #
####################################################################################################


#for i,j in enumerate(name1): 
#	path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/'+j.split(' ')[-1]+'_*_SMAP(JPL).nc')    #  à decommenter pour faire par zone
path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*_SMAP(JPL).nc')          # à commenter pour faire par zone
for k,l in enumerate(path):
	#print(k,l)
	area_fix_smapjpl = 'area_1'
	rad_fix_smapjpl = 'SMAP(JPL)'
	sen_fix_smapjpl = 'SAR (Sentinel 1)'
	file_stat_smapjpl = nc.Dataset(l,'r')
	delta_smapjpl = np.squeeze(file_stat_smapjpl.variables['delta_time'][:])
	lon_smapjpl = np.squeeze(file_stat_smapjpl.variables['Lon'][:])
	lat_smapjpl = np.squeeze(file_stat_smapjpl.variables['Lat'][:])
	nb_points_smapjpl = np.squeeze(file_stat_smapjpl.variables['number_points'][:])
	vent_sar_smapjpl = np.squeeze(file_stat_smapjpl.variables['WindSpeed_sar_mean'][:])
	angle_smapjpl = np.squeeze(file_stat_smapjpl.variables['sar_mean_incidence_angle'][:])
	cat_smapjpl=file_stat_smapjpl.category
	vent_smapjpl = np.squeeze(file_stat_smapjpl.variables['WindSpeed_smap_jpl'][:])


	vent_smapjpl_a = np.ma.masked_where((vent_smapjpl >= 999)*(vent_smapjpl.mask == False)*(vent_sar_smapjpl.mask == False)+(vent_sar_smapjpl.mask == True)  , vent_smapjpl)

										# MASK #
	nb_points_smapjpl = np.ma.masked_where(nb_points_smapjpl <= 400,nb_points_smapjpl)
	vent_sar_smapjpl.mask = nb_points_smapjpl.mask.copy()+vent_smapjpl_a.mask.copy()
	nb_points_smapjpl.mask = nb_points_smapjpl.mask.copy()+vent_smapjpl_a.mask.copy()
	angle_smapjpl.mask = nb_points_smapjpl.mask.copy()+vent_smapjpl_a.mask.copy()
	vent_smapjpl.mask = nb_points_smapjpl.mask.copy()+vent_smapjpl_a.mask.copy()
	dim_smapjpl = vent_smapjpl.shape
	dt_smapjpl = np.ma.zeros(dim_smapjpl)

	if delta_smapjpl < 10000:
		for i in range(len(dt_smapjpl)):
			dt_smapjpl[i] = delta_smapjpl

#			if cat_smapjpl==c:
		sum3=sum3+1	
		try:
			vent_sar_1_smapjpl = np.ma.append(vent_sar_1_smapjpl,vent_sar_smapjpl)
			angle_1_smapjpl = np.ma.append(angle_1_smapjpl,angle_smapjpl)
			vent_smapjpl_1= np.ma.append(vent_smapjpl_1,vent_smapjpl)
			dt_1_smapjpl = np.ma.append(dt_1_smapjpl,dt_smapjpl)
			nb_1_smapjpl = np.ma.append(nb_1_smapjpl,nb_points_smapjpl)
		except:
			vent_sar_1_smapjpl = vent_sar_smapjpl
			angle_1_smapjpl = angle_smapjpl
			vent_smapjpl_1 = vent_smapjpl
			dt_1_smapjpl = dt_smapjpl
			nb_1_smapjpl = nb_points_smapjpl
											
	# STATISTIQUES
try:
	vent_sar_1_smapjpl=np.ma.masked_invalid(vent_sar_1_smapjpl)
	diff_smapjpl = vent_sar_smapjpl - vent_smapjpl
	diff_smapjpl_smapjpl = vent_sar_1_smapjpl - vent_smapjpl_1
	print('smapjpl ',diff_smapjpl_smapjpl.shape)
	diff_moy_smapjpl = diff_smapjpl_smapjpl.mean()
	std_diff_smapjpl = diff_smapjpl_smapjpl.std()
	diff_med_smapjpl = np.ma.median(diff_smapjpl_smapjpl)


	nbins_smapjpl = 10
	x_smapjpl = vent_sar_1_smapjpl.compressed()
	y_smapjpl = diff_smapjpl_smapjpl.compressed()
	n_smapjpl, __smapjpl = np.histogram(x_smapjpl, bins=nbins_smapjpl)
	sy_smapjpl, __smapjpl = np.histogram(x_smapjpl, bins=nbins_smapjpl, weights=y_smapjpl)
	sy2_smapjpl, __smapjpl = np.histogram(x_smapjpl, bins=nbins_smapjpl, weights=y_smapjpl*y_smapjpl)
	mean_smapjpl = sy_smapjpl / n_smapjpl
	std_smapjpl = np.sqrt(sy2_smapjpl/n_smapjpl - mean_smapjpl*mean_smapjpl)
except:
	pass
				
try:					
	print(diff_moy,std_diff,diff_med,size(diff_smos))
except:
	print('0 0 0 0')	
try:			
	print(diff_moy_smap,std_diff_smap,diff_med_smap,size(diff_smap_smap))
except:
	print('0 0 0 0')	
try:
	print(diff_moy_smapjpl,std_diff_smapjpl,diff_med_smapjpl,size(diff_smapjpl_smapjpl))
except:
	print('0 0 0 0')

fig = figure(figsize(16,12))
ax1 = fig.add_axes((0.07,0.5,0.38,0.35))
ax2 = fig.add_axes((0.55,0.5,0.38,0.35))
ax3 = fig.add_axes((0.07,0.05,0.38,0.35))

fig.suptitle('Wind speed difference function of '+sen_fix_smapjpl+' wind speed  (area '+a+' category '+c+')', fontsize = 18)
try:
	im1 = ax1.scatter(x,y,c='green',marker='.',alpha=0.2)
	ax1.set_title('SMOS', fontsize = 16) 
	ax1.set_ylim([-40,40])
	ax1.set_xlim([0,80])
	ax1.set_xlabel(sen_fix_smapjpl+' Wind speed',fontsize=12)
	ax1.set_ylabel('WS difference ('+sen_fix+'-'+rad_fix+') (m/s)',fontsize=14)
	ax1.grid()
	ax1.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, fmt='k-',capthick = 2,alpha=0.6)
	ax1.tick_params(axis = 'both', labelsize = 14)
except:
	pass

try:
	im2 = ax2.scatter(x_smap,y_smap,c='green',marker='.',alpha=0.2) 
	ax2.set_title('SMAP(RSS)', fontsize = 16) 
	ax2.set_ylim([-40,40])
	ax2.set_xlim([0,80])
	ax2.set_xlabel(sen_fix_smapjpl+' Wind speed',fontsize=12)
	ax2.set_ylabel('WS difference ('+sen_fix_smap+'-'+rad_fix_smap+') (m/s)',fontsize=14)
	ax2.grid()
	ax2.errorbar((__smap[1:] + __smap[:-1])/2, mean_smap, yerr=std_smap, fmt='k-',capthick = 2,alpha=0.6)
	ax2.tick_params(axis = 'both', labelsize = 14)
except:
	pass

try:
	im3 = ax3.scatter(x_smapjpl,y_smapjpl,c='green',marker='.',alpha=0.2) 
	ax3.set_title('SMAP(JPL)', fontsize = 16)
	ax3.set_ylim([-40,40])
	ax3.set_xlim([0,80])
	ax3.set_xlabel(sen_fix_smapjpl+' Wind speed',fontsize=12)
	ax3.set_ylabel('WS difference ('+sen_fix_smapjpl+'-'+rad_fix_smapjpl+') (m/s)',fontsize=14)
	ax3.grid()
	ax3.errorbar((__smapjpl[1:] + __smapjpl[:-1])/2, mean_smapjpl, yerr=std_smapjpl, fmt='k-',capthick = 2,alpha=0.6)
	ax3.tick_params(axis = 'both', labelsize = 14)
except:
	pass
savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/STAT_WS/Wind_speed_difference_VS_WS_SAR_area'+a+'_cat'+c+'.png')

try:
	figure(figsize(12,12))
	plt.scatter(x,y,c='green',marker='.',alpha=0.2)
	plt.ylim([-40,40])
	plt.xlim([0,80])
	plt.xlabel(sen_fix+' Wind speed',fontsize=14)
	plt.ylabel('WS difference ('+sen_fix+'-'+rad_fix+') (m/s)',fontsize=14)
	plt.grid()
	plt.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, fmt='k-',capthick = 2,alpha=0.6)
	plt.tick_params(axis = 'both', labelsize = 14)
	plt.title('Wind speed difference function of '+sen_fix+' wind speed  (area '+a+' category '+c+')', fontsize = 18)
	savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/STAT_WS/Wind_speed_(SAR-'+rad_fix+')_VS_WS_SAR_area'+a+'_cat'+c+'.png')
except:
	pass

try:
	figure(figsize(12,12))
	plt.scatter(x_smap,y_smap,c='green',marker='.',alpha=0.2)
	plt.ylim([-40,40])
	plt.xlim([0,80])
	plt.xlabel(sen_fix_smap+' Wind speed',fontsize=14)
	plt.ylabel('WS difference ('+sen_fix_smap+'-'+rad_fix_smap+') (m/s)',fontsize=14)
	plt.grid()
	plt.errorbar((__smap[1:] + __smap[:-1])/2, mean_smap, yerr=std_smap, fmt='k-',capthick = 2,alpha=0.6)
	plt.tick_params(axis = 'both', labelsize = 14)
	plt.title('Wind speed difference function of '+sen_fix_smapjpl+' wind speed (area '+a+' category '+c+') ', fontsize = 18)
	savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/STAT_WS/Wind_speed_(SAR-'+rad_fix_smap+')_VS_WS_SAR_area'+a+'_cat'+c+'.png')
except:
	pass

try:
	figure(figsize(12,12))
	plt.scatter(x_smapjpl,y_smapjpl,c='green',marker='.',alpha=0.2)
	plt.ylim([-40,40])
	plt.xlim([0,80])
	plt.xlabel(sen_fix_smapjpl+' Wind speed',fontsize=14)
	plt.ylabel('WS difference ('+sen_fix_smapjpl+'-'+rad_fix_smapjpl+') (m/s)',fontsize=14)
	plt.grid()
	plt.errorbar((__smapjpl[1:] + __smapjpl[:-1])/2, mean_smapjpl, yerr=std_smapjpl, fmt='k-',capthick = 2,alpha=0.6)
	plt.tick_params(axis = 'both', labelsize = 14)
	plt.title('Wind speed difference function of '+sen_fix_smapjpl+' wind speed (area '+a+' category '+c+')', fontsize = 18)
	savefig('/net/merzhin1/vol/safo/morgane/data/figures/STATS_ZONE/STAT_WS/Wind_speed_(SAR-'+rad_fix_smapjpl+')_VS_WS_SAR_area'+a+'_cat'+c+'.png')
except:
	pass
plt.close('all')
