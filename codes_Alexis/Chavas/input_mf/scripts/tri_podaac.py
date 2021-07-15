# Script permettant le tri des fichiers netcdf par cyclone pour fichier SMAP Podaac
# 08/04/2019 M.G

# Header

import pandas as pd
import os
import os.path
import shutil
import sys
import glob
from glob import glob
from datetime import datetime
from datetime import timedelta
import netCDF4 as nc
from netCDF4 import Dataset
import numpy as np
import h5py


					# Lecture fichier csv

df=pd.read_csv('/net/diva/home/mcms/guizioum/multiResReport.csv')



path_sar=glob('/net/merzhin1/vol/safo/morgane/data/IRMA/sentinel1/*.nc')
src='/net/merzhin1/vol/safo/morgane/temp_data/SMAP/*/*'

for i,j in enumerate(path_sar):
	name_sar = j.split('/')[-1]
	cycl_name = j.split('/')[7]
	day_sar = datetime.strptime(name_sar.split('-')[4][0:8],'%Y%m%d')
	day_sar = datetime.strftime(day_sar,'%Y%m%d')
	time_sar = datetime.strptime(name_sar.split('-')[4],'%Y%m%dt%H%M%S')
	deltaT = timedelta(hours=6)
	deltaT1 = timedelta(minutes=45)
	dst='/net/merzhin1/vol/safo/morgane/data/'+cycl_name+'/smap_podaac/'
	smap_file=glob('%s/SMAP_L2B_SSS_*_%sT*_R16010_V4.2.h5' % (src,day_sar))	
	# print(smap_file)
	for l, m in enumerate(smap_file):
		name_smap = m.split('/')[-1]
		#print(name_smap)
		time_smap = datetime.strptime(name_smap.split('_')[4],'%Y%m%dT%H%M%S')
		#print(time_smap)
		if (time_smap+deltaT1 >= time_sar  - deltaT) & (time_smap+deltaT1 <= time_sar + deltaT):
			print(time_smap,time_sar,name_smap)
			if not os.path.isdir(dst):
				os.makedirs(dst)
			shutil.copy(m,dst)








####################
nc_attrs = .ncattrs()  
smap = {}
file_smap = h5py.File(m, 'r')
list_elmts = [key for key in file_smap['/'].keys()]
for key in list_elmts:
	print(key)

