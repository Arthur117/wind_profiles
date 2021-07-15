# Script permettant le remapping des données SAR à la résolution SMOS
# 27/03/2019 M.G.

import numpy as np

def remap(a1,b1,z,w,mask,x,y):
	res = 0.25 # resolution produit
	dim = a1.shape
	angle_sarmoy = np.zeros(dim)-999
	vent_sarmoy = np.zeros(dim)-999
	vent_sarmed = np.zeros(dim)-999
	vent_sarstd = np.zeros(dim)-999
	n_sarvalid = np.zeros(dim)-999
	for i in np.arange(dim[0]):
	    for j in np.arange(dim[1]):
	         k = np.ma.where((x>a1[i,j]-(res/2)) * (x<a1[i,j]+(res/2)) * (y>b1[i,j]-(res/2)) * (y<b1[i,j]+(res/2)))
	         n_sarvalid[i,j] = len(k[0])-z[k[0],k[1]].mask.sum()
	         if n_sarvalid[i,j] != 0:
	             angle_sarmoy[i,j] = w[k[0],k[1]].mean()
	             vent_sarmoy[i,j]=z[k[0],k[1]].mean()
	             vent_sarmed[i,j] = np.ma.median(z[k[0],k[1]])
	             vent_sarstd[i,j] = z[k[0],k[1]].std()

	angle_sarmoy = np.ma.masked_where(angle_sarmoy == -999 , angle_sarmoy)
	vent_sarmoy = np.ma.masked_where(vent_sarmoy == -999 , vent_sarmoy)
	vent_sarmed = np.ma.masked_where(vent_sarmed == -999 , vent_sarmed)
	vent_sarstd = np.ma.masked_where(vent_sarstd == -999 , vent_sarstd)
	n_sarvalid = np.ma.masked_where(n_sarvalid <= 0 , n_sarvalid)
	return angle_sarmoy, vent_sarmoy, vent_sarmed, vent_sarstd, n_sarvalid
