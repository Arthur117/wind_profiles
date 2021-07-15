# Script permettant le tri des fichiers netcdf par cyclone pour fichier SAR
# 14/03/2019 M.G

# Header

import pandas as pd
import os
import os.path
import shutil
import sys


# Lecture fichier csv

# df=pd.read_csv('/net/diva/home/mcms/guizioum/multiResReport.csv') # Sentinel-1
df=pd.read_csv('/net/diva/home/mcms/guizioum/multiResReport_rs2.csv') # Radarsat-2

# Déclaration/Définition "variables"

R=df['Resolution']
L2M=df['L2M file']
name=df['Comments']
L2M_m=[]
cycl=[]


#Création dossier et déplacement des fichiers nécessaires dans dossier correspondant

    
for j,f in enumerate(L2M):
    r=R[j]
    if r == 3 :
       file_name=f.split('/')[-1]
       L2M_m.append(file_name)
       cycl_name=name[j].split(' ')[0]
       cat = name[j].split(' ')[-1]
       print(cat)
       cycl.append(cycl_name)
#       dst='/net/merzhin1/vol/safo/morgane/data/'+cycl_name+'/sentinel1/'
       dst='/net/merzhin1/vol/safo/morgane/data/'+cycl_name+'/radarsat2/'
       if not os.path.isdir(dst):
           os.makedirs(dst)
#       shutil.copy("/net/merzhin1/vol/safo/morgane/temp_data/sentinel1/"+file_name,dst)
       shutil.copy("/net/merzhin1/vol/safo/morgane/temp_data/radarsat2/"+file_name,dst)
  





