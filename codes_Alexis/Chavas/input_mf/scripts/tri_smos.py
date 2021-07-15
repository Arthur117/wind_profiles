# Script permettant le tri des fichiers netcdf par cyclone pour fichier SMOS
# 14/03/2019 M.G

# Header

import pandas as pd
import os
import os.path
import shutil
import sys
import glob
from glob import glob


# Lecture fichier csv

# df=pd.read_csv('/net/diva/home/mcms/guizioum/multiResReport.csv') # Sentinel-1
df=pd.read_csv('/net/diva/home/mcms/guizioum/multiResReport_rs2.csv') # RS2

# Déclaration/Définition "variables"

R=df['Resolution']
L2M=df['L2M file']
name=df['Comments']
smos=[]
L2M_m=[]
cycl=[]
sar=[]



# Déplacement des fichiers nécessaires dans dossier correspondant

src='/net/merzhin1/vol/safo/morgane/temp_data/smos_output_test_eaf.new_gmf.merged/'

dirs=os.listdir(src)
for i,f in enumerate(dirs):
	smos_name=f.split('.')[0]
	smos.append(smos_name)

for j,k in enumerate(L2M):
    r=R[j]
    if r == 3 :
       file_name=k.split('/')[-1]
       L2M_m.append(file_name)
       sar_date=k.split('-')[4][0:8]
       sar.append(k.split('-')[4][0:8])
       cycl_name=name[j].split(' ')[0]
       cycl.append(cycl_name)
       dst='/net/merzhin1/vol/safo/morgane/data/'+cycl_name+'/smos/'
       smos_file=glob('%s/%s.01d.mat' % (src,sar_date))
       print(smos_file)
       if not os.path.isdir(dst):
           os.makedirs(dst)
       assert len(smos_file) == 1
       shutil.copy(smos_file,dst)






