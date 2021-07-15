# Script permettant la comparaison des données vents entre Sentinel-1 et SMOS
# 15/03/2019 M.G

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
	#j = '/net/merzhin1/vol/safo/morgane/data/CEBILE/sentinel1/s1b-ew-owi-cm-20180207t131831-20180207t131931-000003-01120C.nc'
	cat = cats[j.split('/')[-1]] 
	file_vent=nc.Dataset(j,'r')
	grp = file_vent.groups['owiInversionTables_UV']
	lon_sar=np.squeeze(file_vent.variables['owiLon'][:])
	lat_sar=np.squeeze(file_vent.variables['owiLat'][:])
	ws_sar=np.squeeze(grp.variables['owiWindSpeed_Tab_dualpol_2steps'][:])
	mask=np.squeeze(file_vent.variables['owiMask'][:])
	angle=np.squeeze(file_vent.variables['owiIncidenceAngle'][:])
	name_cycl = j.split('/')[7]
	sat = j.split('/')[-2]
	name_sar, sat_sar, x, y, z, w, xlim, ylim, day_sar,day_sar_1, date_sar, hour_sar_debut, hour_sar_fin = fonc_data.data_sar(lon_sar,lat_sar,mask,ws_sar,angle,j)
	print(ws_sar.min(), ws_sar.max())


				# Données provenant de SMOS



	src='/net/merzhin1/vol/safo/morgane/data/'+name_cycl+'/smos'
	path_smos = glob('%s/%s.01d.mat' % (src,day_sar))
	path_smos = str(path_smos)
	path_smos = path_smos.replace("['","")
	path_smos = path_smos.replace("']","")
	print(path_smos)
	sat_rad = src.split('/')[-1]
	fid = h5py.File(path_smos, 'r')
	name_smos = path_smos.split('/')[-1]
	smos = {}
	smos['date_smos'] = datetime.strptime(name_smos.split('.')[0],'%Y%m%d')
	smos['lon'] = fid['morning']['lon'][::]
	smos['lat'] = fid['morning']['lat'][::]
	smos['date_mo']=fid['morning']['grid_latlon_date'][::] 
	smos['spd_mo'] = np.array(fid['morning']['grid_latlon_rws10'][::])
	smos['date_ev']=fid['evening']['grid_latlon_date'][::] 
	smos['spd_ev'] = np.array(fid['evening']['grid_latlon_rws10'][::]) 
	smos['dist_mo'] = np.array(fid['morning']['grid_latlon_xtrack_dist_km'][::])
	smos['dist_ev'] = np.array(fid['evening']['grid_latlon_xtrack_dist_km'][::])
	try:
		date1, a1, b1, vent_smos, dt, dist = fonc_data.data_smos(x,y,smos,date_sar)

					# Remapping

		smos['lon'][smos['lon']>180] = smos['lon'][smos['lon']>180]-360
		x[x>180] = x[x>180]-360

		swath_sar = SwathDefinition(x,y)
		swath_smos = SwathDefinition(smos['lon'],smos['lat'])

		i = np.arange(smos['lon'].shape[0]) # SMAP shape
		j = np.arange(smos['lon'].shape[1]) # SMAP shape
		j_smos,i_smos = np.meshgrid(j,i) # SMAP shape

		result_i = resample_nearest(swath_smos, i_smos, swath_sar, radius_of_influence=30000, fill_value=None) #SAR shape
		result_j = resample_nearest(swath_smos, j_smos, swath_sar, radius_of_influence=30000, fill_value=None) #SAR shape

		if result_i.mask.all() == False:
			ws_smosss = vent_smos[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			print(ws_smosss.min(),ws_smosss.max())
	except:
		pass
			i_smosss = i_smos[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			j_smosss = j_smos[result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			lon_smosss = smos['lon'][result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			lat_smosss = smos['lat'][result_i.min():result_i.max()+1,result_j.min():result_j.max()+1]
			dim = i_smosss.shape
			vent_sarmoy = np.ma.zeros(dim)-999
			vent_sarmed = np.ma.zeros(dim)-999
			vent_sarstd = np.ma.zeros(dim)-999
			n_sarvalid = np.ma.zeros(dim)-999
			angle_sarmoy = np.ma.zeros(dim)-999
			for n in np.arange(i_smosss.shape[0]):
				for m in np.arange(j_smosss.shape[1]):
					s = ((result_i == i_smosss[n,m]) & (result_j == j_smosss[n,m])) # mask booléen
					vent_sarmoy[n,m] = z[s].mean()
					vent_sarmed[n,m] = np.median(z[s])
					vent_sarstd[n,m] = z[s].std()
					n_sarvalid[n,m] = len(s)-z[s].mask.sum()
					angle_sarmoy[n,m] = w[s].mean()

			angle_sarmoy = np.ma.masked_where(angle_sarmoy == -999 , angle_sarmoy)
			vent_sarmoy = np.ma.masked_where(vent_sarmoy == -999 , vent_sarmoy)
			vent_sarmed = np.ma.masked_where(vent_sarmed == -999 , vent_sarmed)
			vent_sarstd = np.ma.masked_where(vent_sarstd == -999 , vent_sarstd)
			n_sarvalid = np.ma.masked_where((n_sarvalid <= 0)*(vent_sarmoy == -999), n_sarvalid)
			

			# Graphiques comparant les vitesses des vents selon l'origine

			d=date1.strftime('%Y%m%dT%H%M%S')

			date_smos_moyenne = date1.strftime('%Y/%m/%d %H:%M:%S %p')
			hour_smos_moyenne = date_smos_moyenne.split(' ')[1]
			smos_date = date1.strftime('%Y/%m/%d')

			palette = '/net/diva/home/mcms/guizioum/high_wind_speed.pal'
			cwnd = getColorMap( rgbFile = palette )

			m = Basemap(llcrnrlon=x.min(),llcrnrlat=y.min(),urcrnrlon=x.max(),urcrnrlat=y.max(), 			resolution='i', projection='tmerc', lat_0 =y.mean() , lon_0 =x.mean() )

			x1,y1 = m(x,y)
			a2,b2 = m(lon_smosss,lat_smosss)

			fig = figure(figsize = (16,8))
			fig.suptitle(name_cycl+' '+ smos_date+' SAR '+sat_sar+' '+hour_sar_debut+' SMOS '+hour_smos_moyenne+'  cat: '+str(cat) ,fontsize = 18)

			ax = fig.add_axes((0.05,0.1,0.68,0.02))
			ax1 = fig.add_axes((0.05,0.22,0.18,0.65))
			ax2 = fig.add_axes((0.30,0.22,0.18,0.65))
			ax3 = fig.add_axes((0.55,0.22,0.18,0.65))
			ax4 = fig.add_axes((0.78,0.22,0.18,0.65))
			ax5 = fig.add_axes((0.78,0.1,0.18,0.02))

			im = m.pcolormesh(x1.filled(x1.mean()),y1.filled(y1.mean()),z,cmap=cwnd,vmin = 0 , vmax = 80, ax = ax1)
			m.drawcoastlines(ax = ax1)
			m.fillcontinents(color='grey',ax = ax1)
			ax1.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
			ax1.set_ylabel('Latitude', fontsize = 12)
			ax1.set_title('Sentinel-1 '+sat_sar, fontsize = 12)
			m.drawparallels(np.arange(int(y.min()),int(y.max()+1),1),labels=[0,1,0,0],ax = ax1)
			m.drawmeridians(np.arange(int(x.min()),int(x.max()+1),1),labels=[0,0,0,1],ax = ax1)


			im = m.pcolormesh(a2,b2,ws_smosss,cmap=cwnd, vmin = 0 , vmax = 80,ax = ax2)
			m.drawcoastlines(ax = ax2)
			m.fillcontinents(color='grey',ax = ax2)
			ax2.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
			ax2.set_title('SMOS', fontsize = 12)
			m.drawparallels(np.arange(int(y.min()),int(y.max()+1),1),labels=[0,1,0,0],ax = ax2)
			m.drawmeridians(np.arange(int(x.min()),int(x.max()+1),1),labels=[0,0,0,1],ax = ax2)

			im = m.pcolormesh(a2,b2,vent_sarmoy,cmap=cwnd,vmin = 0 , vmax = 80,ax = ax3)
			m.drawcoastlines(ax = ax3)
			m.fillcontinents(color='grey',ax = ax3)
			plt.colorbar(im, cax = ax, orientation='horizontal', label='Wind speed (m/s)')
			ax3.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
			ax3.set_title('Sentinel-1 -> SMOS', fontsize = 12)
			m.drawparallels(np.arange(int(y.min()),int(y.max()+1),1),labels=[0,1,0,0],ax = ax3)
			m.drawmeridians(np.arange(int(x.min()),int(x.max()+1),1),labels=[0,0,0,1],ax = ax3)

			im = m.pcolormesh(a2,b2,n_sarvalid,cmap='jet',vmin = 0 , vmax = 1000,ax = ax4)
			m.drawcoastlines(ax = ax4)
			m.fillcontinents(color='grey',ax = ax4)
			plt.colorbar(im, cax = ax5, orientation='horizontal', label='Number of points')
			ax4.set_xlabel('Longitude',labelpad = 20, fontsize = 12)
			ax4.set_title('Number of points by pixel SMOS', fontsize = 12)
			m.drawparallels(np.arange(int(y.min()),int(y.max()+1),1),labels=[0,1,0,0],ax = ax4)
			m.drawmeridians(np.arange(int(x.min()),int(x.max()+1),1),labels=[0,0,0,1],ax = ax4)

			savefig('/net/merzhin1/vol/safo/morgane/data/figures/SAR_SMOS/SAR_SMOS_'+name_cycl+'_'+d+'_compare_pyre.png')
	
			plt.close('all')


	# Création fichier netCDF

		file_net = nc.Dataset('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/'+name_cycl+'_'+d+'_SAR_'+sat_sar+'_SMOS.nc', 'w', format='NETCDF4')
	
		fonc_create_netcdf.create_netcdf_smos(file_net,name_cycl, name_sar, name_smos, hour_sar_debut, hour_sar_fin, date_smos_moyenne, lon_smosss, lat_smosss, vent_sarmoy, vent_sarmed, vent_sarstd, ws_smosss, n_sarvalid, dt, angle_sarmoy,sat,cat,sat_sar,sat_rad )
		
		file_net.close()

	except UnboundLocalError:
		print('Pas de Colocalisation SMOS '+name_cycl+' '+day_sar)

fid.close()

file_vent.close()

# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))


######################################################################################################
####################################### fin programme ################################################
######################################################################################################

# Permet de voir le nom des variables de smos

# def printname(name): 
#     print (name)
# fid.visit(printname)  


# Variables potentiellement utiles

# smos['error_mo'] = fid['morning']['grid_latlon_var_rws10'][::]
# smos['xtd_mo'] = fid['morning']['grid_latlon_xtrack_dist_km'][::]
# smos['error_ev'] = fid['evening']['grid_latlon_var_rws10'][::]
# smos['xtd_ev'] = fid['evening']['grid_latlon_xtrack_dist_km'][::]
