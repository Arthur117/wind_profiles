# Script delimitant zone cyclonique et les plaçant dans fichier texte
# 21/05/2019 M.G

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


# Debut du decompte du temps
start_time = time.time()

# Delimitation zones

(area_1,area_2,area_3,area_4,area_5,area_6,area_7,area_11,area_11b,area_12,area_12b)=fonc_poly.area()

path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*.nc')
lat_moy = []
lon_moy = []

with open('/net/merzhin1/vol/safo/morgane/data/cycl_area.txt','a') as f: 
	f.write('%10s,%12s, %12s,%6s,%8s,%10s,%10s,%s\n' % ('Name','Date','Hour','SAR','RAD','Lon mean','Lat mean','Area'))
	for i,j in enumerate(path):
		print(i,j)
		file_stat = nc.Dataset(j,'r')
		file_name =  j.split('/')[-1]
		name_cycl = file_name.split('_')[0]
		sat_sar = file_name.split('_')[3]
		sat_rad = (file_name.split('_')[-1]).split('.')[0]
		date = file_name.split('_')[1][0:8]
		hour = file_name.split('_')[1][9:15]
		lat = np.squeeze(file_stat.variables['Lat'][:])
		lon = np.squeeze(file_stat.variables['Lon'][:])

		point = Point(lon.mean(),lat.mean())

		if area_1.contains(point):
			a = '01'
		elif area_2.contains(point):
			a = '02'
		elif area_3.contains(point):
			a = '03'
		elif area_4.contains(point):
			a = '04'
		elif area_5.contains(point):
			a = '05'
		elif area_6.contains(point):
			a = '06'
		elif area_7.contains(point):
			a = '07'
		elif area_11.contains(point):
			a = '11'
		elif area_11b.contains(point):
			a = '11b'
		elif area_12.contains(point):
			a = '12'
		elif area_12b.contains(point):
			a = '12b'
		else:
			a = ' '

		f.write('%10s,%12s,%12s,%6s,%8s,%10.2f,%10.2f,%s\n' % (name_cycl,date,hour,sat_sar,sat_rad,lon.mean(),lat.mean(),a))

# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))
