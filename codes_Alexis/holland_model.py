import numpy as np

def bs_function(dp,Lat,Vtrans,dpt=0):
    xa = 0.6*(1-dp/215)
    return -4.4*10**(-5)*dp**2. + 0.01*dp + 0.03*dpt - 0.014*Lat + 0.15*Vtrans**xa + 1.0



def holland_model(radius, Rmax, dp, Lat, Router, Vouter, rho=None, Vtrans=0, dpt=0):

    bs = bs_function(dp,Lat,Vtrans,dpt)
    
    if rho is None:
        rho = 1.2041 # Constant. Need for SST to refine


    Rnorm_Router = (Rmax/Router)
    Ratio = (100*bs*dp*Rnorm_Router**bs)/(rho*np.exp(Rnorm_Router))
    Xouter = np.log(Vouter)/np.log(Ratio)

    # X fonction
    x = radius*0+0.5
    x[radius>Rmax]= 0.5 + (radius[radius>Rmax]-Rmax)*(Xouter-0.5)/(Router-Rmax)

    # wind Fonction
    rnorm = (Rmax/radius)**bs
    u10_Holland = ((100*bs*dp*rnorm)/(rho*np.exp(rnorm)))**x 

    return u10_Holland,bs
