# Script delimitant zone cyclonique
# 16/05/2019 M.G

# Headers

import matplotlib.pyplot as plt
from shapely.geometry.polygon import LinearRing, Polygon
from shapely.geometry import Point
import time
import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.ticker as mticker
import glob
from glob import glob
from datetime import datetime
import time
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import fonc_poly
import pandas as pd

# Debut du decompte du temps
start_time = time.time()

# Delimitation zones

(area_1,area_2,area_3,area_4,area_5,area_6,area_7,area_11,area_11b,area_12,area_12b)=fonc_poly.area()

figure(figsize(20,12))
fonc_poly.poly_fig(area_1,area_2,area_3,area_4,area_5,area_6,area_7,area_11,area_11b,area_12,area_12b)

path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*.nc')
lat_moy = []
lon_moy = []

df=pd.read_csv('/net/merzhin1/vol/safo/morgane/data/cycl_area.txt')
e = np.where(df.iloc[:,3] == 'SMOS') 



for i,j in enumerate(path):
	print(i,j)
	file_stat = nc.Dataset(j,'r')
	file_name =  j.split('/')[-1]
	name_cycl = file_name.split('_')[0]
	sat_sar = file_name.split('_')[2]
	sat_rad = (file_name.split('_')[-1]).split('.')[0]
	date = file_name.split('_')[1][0:8]
	lat = np.squeeze(file_stat.variables['Lat'][:])
	lon = np.squeeze(file_stat.variables['Lon'][:])

	lat_moy.append(lat.mean())
	lon_moy.append(lon.mean())

	if sat_rad == 'SMOS':
		plt.plot([lon.mean()], [lat.mean()],marker='o',color='b',alpha = 0.5,label="SMOS " if i == 317 else "")
	if sat_rad == 'SMAP' : 
		plt.plot([lon.mean()], [lat.mean()],marker='o',color='r',alpha = 0.5)
	if sat_rad == 'SMAP(JPL)' :
		plt.plot([lon.mean()], [lat.mean()],marker='o',color='r',alpha = 0.5,label="SMAP " if i == 169 else "")
plt.text(-100,50,'1',fontsize='16')
plt.text(-135,50,'2',fontsize='16')
plt.text(-175,50,'3',fontsize='16')
plt.text(105,50,'4',fontsize='16')
plt.text(55,50,'5',fontsize='16')
plt.text(25,-10,'6',fontsize='16')
plt.text(95,-5,'7,8,9,10',fontsize='16')
plt.text(-130,-5,'11',fontsize='16')
plt.text(165,-5,'11',fontsize='16')
plt.text(-130,-30,'12',fontsize='16')
plt.text(165,-30,'12',fontsize='16')








plt.title('Tropical cyclone localisation with WMO Regional Specialized Meteorological Centres and their areas of responsibility ', fontsize = 18)
plt.legend()
savefig('/net/merzhin1/vol/safo/morgane/data/figures/map_areas_localisation_tropical_cyclone.png')


# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))
