from glob import glob
import numpy as np

import pdb
import pickle

import matplotlib.pyplot as plt

from ER11E04_nondim_rfitinput import ER11E04_nondim_rfitinput
from adjust_chavas_on_sar_data import adjust_chavas_on_sar_data

##=======================================
## Storm parameters from SAR
##=======================================

Vmax = 57                      #[ms-1] {50} maximum azimuthal-mean wind speed
Rfit = 100*1000                #[m] {300*1000} a wind radius
Vfit = 17                      #[m] {12} wind speed at rfit
Lat  = 12.3

_file = 'polar_analysis_synthesis_s1b-ew-owi-cm-20180818t083657-20180818t084001-000003-016B1D_gd.pkl'
_file = 'polar_analysis_synthesis_s1a-ew-owi-cm-20180818t150300-20180818t150604-000003-0288C3_gd.pkl'
_path = '/home/amouche/science/hurricanes/CyclObs/'
res = pickle.load(open(_path+_file,'rb'))


files = glob(_path+'polar*pkl')
for _file in files:

    

    ret = adjust_chavas_on_sar_data(_file)

