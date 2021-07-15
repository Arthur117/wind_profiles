# Script permettant la projection des données vents entre Sentinel-1 et SMAP JPL (podaac) 
# 01/04/2019 M.G

# Headers

import netCDF4 as nc
import numpy as np
from pyresample import load_area, save_quicklook, SwathDefinition
from pyresample.kd_tree import resample_nearest
import h5py
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime
from datetime import datetime
import pandas as pd
from palette import *
import matplotlib.pyplot as plt
import glob
from glob import glob
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import time

# Debut du decompte du temps
start_time = time.time()



					# Lecture fichier csv #

df=pd.read_csv('/net/diva/home/mcms/guizioum/multiResReport.csv') # Sentinel-1


				# Déclaration/Définition "variables" #

R=df['Resolution']
L2M=df['L2M file']
name=df['Comments']
L2M_m=[]
cycl=[]


						# Catégorie #

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

			# Données provenant de  SAR #

path_sar=glob('/net/merzhin1/vol/safo/morgane/data/*/sentinel1/*.nc')

for k,l in enumerate(path_sar):
	print(k,l)
	cat = cats[l.split('/')[-1]] 
	file_sar = nc.Dataset(l,'r')
	name_sar = l.split('/')[-1]
	sat_sar = name_sar.split('-')[0]
	sat = l.split('/')[-2]
	name_cycl = l.split('/')[7]

	grp = file_sar.groups['owiInversionTables_UV']
	lon_sar=np.squeeze(file_sar.variables['owiLon'][:])
	lat_sar=np.squeeze(file_sar.variables['owiLat'][:])
	ws_sar=np.squeeze(grp.variables['owiWindSpeed_Tab_dualpol_2steps'][:])
	mask=np.squeeze(file_sar.variables['owiMask'][:])
	angle=np.squeeze(file_sar.variables['owiIncidenceAngle'][:])

	x = np.ma.array(lon_sar)
	y = np.ma.array(lat_sar)
	z = np.ma.masked_where(mask!=0 , ws_sar)
	w = np.ma.masked_where(mask!=0 , angle)
	xlim = (x.min(), x.max())
	ylim = (y.min(), y.max())
	day_sar = datetime.strptime(name_sar.split('-')[4][0:8],'%Y%m%d')

	day_sar_1 = day_sar.strftime('%Y_%m_%d')    
	day_sar = day_sar.strftime('%Y%m%d')
	date_sar = datetime.strptime(name_sar.split('-')[4][0:8]+name_sar.split('-')[4][9:15],'%Y%m%d%H%M%S')
	hour_sar1 = datetime.strptime(name_sar.split('-')[4][9:15],'%H%M%S')
	hour_sar_debut = hour_sar1.strftime('%H:%M:%S')
	hour_sar2 = datetime.strptime(name_sar.split('-')[5][9:15],'%H%M%S')
	hour_sar_fin = hour_sar2.strftime('%H:%M:%S')

				#  Données provenant de SMAP (JPL) #

	src = '/net/merzhin1/vol/safo/morgane/data/'+name_cycl+'/smap_podaac'
	path_smap = glob('%s/SMAP_L2B_SSS_*_%sT*_*_V4.2.h5' % (src,day_sar))
	for g in path_smap:
		print(g)
		name_smap = g.split('/')[-1]
		date_smap = datetime.strptime(name_smap.split('_')[4],'%Y%m%dT%H%M%S')
		hour_smap = datetime.strptime(name_smap.split('_')[4][9:15],'%H%M%S')
		sat_rad = name_smap.split('_')[0]
		date_smap1 = datetime.strftime(date_smap,'%Y%m%dT%H%M%S')
		file_smap = h5py.File(g, 'r')

		lon_smap = file_smap['/lon'][:]
		lon_smap = np.array(lon_smap)
		lon_smap = np.ma.masked_where(lon_smap == -9999, lon_smap)

		lat_smap = file_smap['/lat'][:]
		lat_smap = np.ma.masked_where(lat_smap == -9999, lat_smap)

		ws_smap = file_smap['/smap_spd'][:]
		ws_smap = np.ma.masked_where(ws_smap == -9999, ws_smap)

		swath_sar = SwathDefinition(x,y)
		swath_smap = SwathDefinition(lon_smap,lat_smap)

		i = np.arange(ws_smap.shape[0]) # SMAP shape
		j = np.arange(ws_smap.shape[1]) # SMAP shape
		j_smap,i_smap = np.meshgrid(j,i) # SMAP shape

						# Resample / Remapping #

		result_i = resample_nearest(swath_smap, i_smap, swath_sar, radius_of_influence=21500, fill_value=None) #SAR shape
		result_j = resample_nearest(swath_smap, j_smap, swath_sar, radius_of_influence=21500, fill_value=None) #SAR shape
		result = resample_nearest(swath_smap, ws_smap,swath_sar, radius_of_influence=21500, fill_value=None)

		if result_i.mask.all() == False:
			ws_smapss = ws_smap[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			i_smapss = i_smap[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			j_smapss = j_smap[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			lon_smapss = lon_smap[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			lat_smapss = lat_smap[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]

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
			n_sarvalid = np.ma.masked_where(n_sarvalid <= 0 , n_sarvalid)

			dt = abs((hour_sar1 - hour_smap).total_seconds())


						# Graphiques comparant les vitesses des vents selon l'origine #

			a1 = lon_smapss.max()
			b1 = lon_smapss.min()
			a2 = x.max()
			b2 = x.min()
			if a1 > a2:
				lon_max = a1
			else:
				lon_max = a2
			if b1 > b2:
				lon_min = b1
			else:
				lon_min = b2

			a3 = lat_smapss.max()
			b3 = lat_smapss.min()
			a4 = y.max()
			b4 = y.min()
			if a3 > a4:
				lat_max = a3
			else:
				lat_max = a4
			if b3 > b4:
				lat_min = b3
			else:
				lat_min = b4

			palette = '/net/diva/home/mcms/guizioum/high_wind_speed.pal'
			cwnd = getColorMap( rgbFile = palette )


			fig = figure(figsize = (16,8))
			fig.suptitle(name_cycl+' SAR '+hour_sar_debut+' SMAP (JPL) '+ str(date_smap)+' cat: '+str(cat),fontsize = 18)

			ax = fig.add_axes((0.05,0.1,0.68,0.02))
			ax1 = fig.add_axes((0.05,0.22,0.18,0.65),projection=ccrs.PlateCarree())
			ax2 = fig.add_axes((0.30,0.22,0.18,0.65),projection=ccrs.PlateCarree())
			ax3 = fig.add_axes((0.55,0.22,0.18,0.65),projection=ccrs.PlateCarree())
			ax4 = fig.add_axes((0.78,0.22,0.18,0.65),projection=ccrs.PlateCarree())
			ax5 = fig.add_axes((0.78,0.1,0.18,0.02))

			im = ax1.pcolor(x.filled(x.mean()),y.filled(y.mean()),z,cmap=cwnd,vmin = 0 , vmax = 80)
			try:
				ax1.set_extent([lon_min-1,lon_max+1,lat_min-1,lat_max+1])
			except:
				pass
			ax1.coastlines('50m')
			ax1.set_title('Sentinel-1 '+sat_sar)
			ax1.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor = 'black' , facecolor='grey'))
			ax1.gridlines(linewidth = 0.5 , linestyle='--')
			gl = ax1.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
			gl.xlabels_top = False
			gl.ylabels_left = False
			gl.ylabels_right=True
			gl.xlines = True
			gl.xformatter = LONGITUDE_FORMATTER
			gl.yformatter = LATITUDE_FORMATTER
			ax1.text(-0.02, 23, 'Latitude', va='bottom', ha='center', rotation='vertical', rotation_mode='anchor', transform=ax.transAxes,fontsize = 12)
			ax1.text(0.13, 5, 'Longitude', va='bottom', ha='center', rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes,fontsize = 12)

			im = ax2.pcolor(lon_smapss,lat_smapss,ws_smapss,cmap=cwnd, vmin = 0 , vmax = 80)
			try:
				ax2.set_extent([lon_min-1,lon_max+1,lat_min-1,lat_max+1])
			except:
				pass
			ax2.coastlines('50m')
			ax2.set_title('SMAP')
			ax2.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor = 'black' , facecolor='grey'))
			ax2.gridlines(linewidth = 0.5 , linestyle='--')
			gl = ax2.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
			ax2.text(0.5, 5, 'Longitude', va='bottom', ha='center', rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes,fontsize = 12)
			gl.xlabels_top = False
			gl.ylabels_left = False
			gl.ylabels_right=True
			gl.xlines = True
			gl.xformatter = LONGITUDE_FORMATTER
			gl.yformatter = LATITUDE_FORMATTER

			im = ax3.pcolor(lon_smapss,lat_smapss,vent_sarmoy,cmap=cwnd,vmin = 0 , vmax = 80)
			try:
				ax3.set_extent([lon_min-1,lon_max+1,lat_min-1,lat_max+1])
			except:
				pass
			ax3.coastlines('50m')
			plt.colorbar(im, cax = ax, orientation='horizontal', label='Wind speed (m/s)')
			ax3.set_xlabel('Longitude',labelpad = 20)
			ax3.set_title('Sentinel-1 '+sat_sar+' @ SMAP resolution')
			ax3.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor = 'black' , facecolor='grey'))
			ax3.gridlines(linewidth = 0.5 , linestyle='--')
			gl = ax3.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
			gl.xlabels_top = False
			gl.ylabels_left = False
			gl.ylabels_right=True
			gl.xlines = True
			gl.xformatter = LONGITUDE_FORMATTER
			gl.yformatter = LATITUDE_FORMATTER
			ax3.text(0.87, 5, 'Longitude', va='bottom', ha='center', rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes,fontsize = 12)

			im = ax4.pcolor(lon_smapss,lat_smapss,n_sarvalid,cmap='jet',vmin = 0 , vmax = 1000)
			try:
				ax4.set_extent([lon_min-1,lon_max+1,lat_min-1,lat_max+1])
			except:
				pass
			ax4.coastlines('50m')
			plt.colorbar(im, cax = ax5, orientation='horizontal', label='Number of points')
			ax4.set_xlabel('Longitude',labelpad = 20)
			ax4.set_title('Number of points by pixel SMAP')
			ax4.add_feature(cfeature.NaturalEarthFeature('physical','land','50m',edgecolor = 'black' , facecolor='grey'))
			ax4.gridlines(linewidth = 0.5 , linestyle='--')
			gl = ax4.gridlines(crs=ccrs.PlateCarree(), linewidth=0.5, color='black', alpha=0.5, linestyle='--', draw_labels=True)
			gl.xlabels_top = False
			gl.ylabels_left = False
			gl.ylabels_right=True
			gl.xlines = True
			gl.xformatter = LONGITUDE_FORMATTER
			gl.yformatter = LATITUDE_FORMATTER
			ax4.text(1.22, 5, 'Longitude', va='bottom', ha='center', rotation='horizontal', rotation_mode='anchor', transform=ax.transAxes,fontsize = 12)

			fig.savefig('/net/merzhin1/vol/safo/morgane/data/figures/SAR_SMAP_(JPL)/SAR_SMAP_JPL_compare_'+name_cycl+'_'+date_smap1+'.png')
			plt.close('all')

								# Création netcdf #


			file_net = nc.Dataset('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/'+name_cycl+'_'+date_smap1+'_SAR_'+sat_sar+'_SMAP(JPL).nc', 'w', format='NETCDF4')
			file_net.description = 'Data processed from Sentinel-1 to SMAP (JPL) resolution for tropical hurricane '+name_cycl
			file_net.source_sar = name_sar
			file_net.source_smap = name_smap
			file_net.cyclone = name_cycl
			file_net.hour_sar = 'Begin : '+hour_sar_debut+' End : '+hour_sar_fin
			file_net.hour_smap = str(date_smap)
			file_net.category = str(cat)
			file_net.SAR = sat+' '+sat_sar
			file_net.SMAP = sat_rad

			file_net.createDimension('lon', len(lon_smapss[0,:]))
			file_net.createDimension('lat', len(lat_smapss[:,0]))
			file_net.createDimension('date', None)


			longitude = file_net.createVariable('Lon', 'f4', ('lat','lon'))
			longitude.units = "degrees_east"
			longitude[:] = lon_smapss

			latitude = file_net.createVariable('Lat', 'f4', ('lat','lon'))  
			latitude.units = "degrees_north"
			latitude[:] = lat_smapss

			ventsarmoy = file_net.createVariable('WindSpeed_sar_mean', 'f4', ('lat','lon'))
			ventsarmoy.units = "m/s"
			ventsarmoy.long_name = " Mean value of wind speed "
			ventsarmoy[:]=vent_sarmoy

			mediane = file_net.createVariable('WindSpeed_sar_median', 'f4', ('lat','lon') )
			mediane.units = "m/s"
			mediane.long_name = " Median value of wind speed "
			mediane[:]=vent_sarmed

			ecart_type = file_net.createVariable('WindSpeed_sar_std', 'f4', ('lat','lon'))
			ecart_type.units = "m/s"
			ecart_type.long_name = " Standard deviation "
			ecart_type[:] = vent_sarstd

			ventsmap = file_net.createVariable('WindSpeed_smap_jpl', 'f4', ('lat','lon'))
			ventsmap.units = " m/s "
			ventsmap.long_name = " Wind speed by smap acquisition (jpl)"
			ventsmap[:] = ws_smapss

			nsarvalid = file_net.createVariable('number_points', 'f4',  ('lat','lon'))
			nsarvalid.units = " no dimension "
			nsarvalid.long_name = " Number of points "
			nsarvalid[:] = n_sarvalid

			date = file_net.createVariable('delta_time','f4','date')
			date.units = "s"
			date.long_name = "Time difference between sar and smap (jpl) acquisition"
			date[:] = dt

			anglesarmoy = file_net.createVariable('sar_mean_incidence_angle', 'f4',  ('lat','lon'))
			anglesarmoy.units = "degrees"
			anglesarmoy.long_name = " Incidence angle "
			anglesarmoy[:]=angle_sarmoy
			file_net.close()

		else:
			print('Pas de colocalisation')
	file_smap.close()
	file_sar.close()
# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))

