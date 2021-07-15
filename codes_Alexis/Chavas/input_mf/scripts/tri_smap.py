# Script permettant le tri des fichiers netcdf par cyclone pour fichier SMAP
# 28/03/2019 M.G

# Header

import pandas as pd
import os
import os.path
import shutil
import sys
import glob
from glob import glob
import datetime
from datetime import datetime

# Lecture fichier csv

df=pd.read_csv('/net/diva/home/mcms/guizioum/multiResReport.csv')


# Déclaration/Définition "variables"

R=df['Resolution']
L2M=df['L2M file']
name=df['Comments']
smos=[]
L2M_m=[]
cycl=[]
sar=[]



# Déplacement des fichiers nécessaires dans dossier correspondant

# src='/net/merzhin1/vol/safo/morgane/temp_data/smap_remss/*/' # remss

src='/net/merzhin1/vol/safo/morgane/temp_data/SMAP/*/*/'   # podaac

for j,k in enumerate(L2M):
    r=R[j]
    if r == 3 :
       file_name=k.split('/')[-1]
       L2M_m.append(file_name)
       sar_date=k.split('-')[5][0:8]
       sar_hour=k.split('-')[5][9:15]
       sar_hour=datetime.strptime(k.split('-')[5][9:15],'%H%M%S')
#       sar_date=datetime.strptime(sar_date,'%Y%m%d') 
#       sar_date2 = datetime.strftime(sar_date,'%Y_%m_%d')
       sar.append(sar_date)
       cycl_name=name[j].split(' ')[0]
       cycl.append(cycl_name)
#       dst='/net/merzhin1/vol/safo/morgane/data/'+cycl_name+'/smap_remss/'
       dst='/net/merzhin1/vol/safo/morgane/data/'+cycl_name+'/smap_podaac/'
       smap_file=glob('%s/SMAP_L2B_SSS_*_%sT*_R16010_V4.2.h5' % (src,sar_date))
       print(smap_file)
       try:
           if not os.path.isdir(dst):
               os.makedirs(dst)
           assert len(smap_file) == 1
           shutil.copy(smap_file[0],dst)
       except:
           pass


