# Script permettant la localisation des cyclones sur une carte
# 15/05/2019 M.G.

# Headers

import numpy as np
import netCDF4 as nc
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import glob
from glob import glob
import time
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Debut du decompte du temps
start_time = time.time()

path = glob('/net/merzhin1/vol/safo/morgane/data/netcdf_intermediaire/*.nc')

figure(figsize = (20,12))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.stock_img()
ax.coastlines()
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

for i,j in enumerate(path):
	print(i,j)
	file_stat = nc.Dataset(j,'r')
	name_cycl = j.split('/')[7]
	lat = np.squeeze(file_stat.variables['Lat'][:])	
	lon = np.squeeze(file_stat.variables['Lon'][:])

	lat_moy = lat.mean()
	lon_moy = lon.mean()

	plt.plot([lon_moy], [lat_moy],marker='o',color='b')
	plt.title('Tropical cyclone localisation', fontsize = 18)
	savefig('/net/merzhin1/vol/safo/morgane/data/figures/map_localisation_tropical_cyclone.png')





# Affichage du temps d execution
print("Temps d execution : %s secondes ---" % (time.time() - start_time))
