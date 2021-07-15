import sys

import numpy as np

sys.path.append('/home/amouche/PythonTool/palette')
from palette import *


import CLE15_2020_06_23


Rfit = np.array([350,350,350])
Vfit = 17
Az   = np.array([0,90,180])
Vmax = 25 
Vfm  = 5 
Dfm  = 0
Lat  = 20.
Rlim = 800 # outer core limit



## Storm parameters
Vmax = 39                      #[ms-1] {50} maximum azimuthal-mean wind speed
Rfit = 93*1000                #[m] {300*1000} a wind radius
Vfit = 17                      #[m] {12} wind speed at rfit
Lat  =12.3

RRF, WWF, xf,yf = CLE15.chv(Rfit,Vfit,Vmax,Lat)
