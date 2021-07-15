import copy

import numpy as np
from scipy.interpolate import interp1d

def ER11_radprof_raw(Vmax,r_in,rmax_or_r0,fcor,CkCd,rr_ER11,eyealpha=1):
    fcor = np.abs(fcor)
    if rmax_or_r0 == 'rmax':
        rmax = r_in
    else:
        print('rmax_or_r0 must be set to "rmax"')
        
    ## CALCULATE Emanuel and Rotunno (2011) theoretical profile
    V_ER11 = (1./rr_ER11)*(Vmax*rmax + .5*fcor*rmax**2)*((2*(rr_ER11/rmax)**2)/(2-CkCd+CkCd*(rr_ER11/rmax)**2))**(1/(2-CkCd)) - .5*fcor*rr_ER11
    ## make V=0 at r=0
    V_ER11[rr_ER11==0] = 0
   
    if rmax_or_r0 == 'rmax':
        i_rmax = np.argwhere(V_ER11==np.max(V_ER11))[0,0]
        f = interp1d( V_ER11[i_rmax+1:],rr_ER11[i_rmax+1:],fill_value='extrapolate')
        r0_profile = f(0.)
        r_out = r0_profile.tolist() #use value from profile itself
    else:
        print('rmax_or_r0 must be set to"rmax"')
    return V_ER11,r_out

def ER11_radprof(Vmax,r_in,rmax_or_r0,fcor,CkCd,rr_ER11):
    dr = rr_ER11[1]-rr_ER11[0]
    ##  Call ER11_radprof_raw
    V_ER11,r_out = ER11_radprof_raw(Vmax,r_in,rmax_or_r0,fcor,CkCd,rr_ER11)
    
    if rmax_or_r0 == 'rmax':
        drin_temp = r_in-rr_ER11[np.argwhere(V_ER11==np.max(V_ER11))[0,0]]
    elif rmax_or_r0 =='r0':
        f = interp1d( V_ER11[2:],rr_ER11[2:])
        drin_temp = r_in-f(0).tolist() 
    ## Calculate error in Vmax
    dVmax_temp = Vmax - np.max(V_ER11)

    ## Check is errors are too large and adjust accordingly
    r_in_save = copy.copy(r_in)
    Vmax_save = copy.copy(Vmax)

    n_iter = 0
    #print('while1')
    while((np.abs(drin_temp)>dr/2) or (np.abs(dVmax_temp/Vmax_save)>=10**-2)): #if error is sufficiently large NOTE: FIRST ARGUMENT MUST BE ">" NOT ">=" or else rmax values at exactly dr/2 intervals (e.g. 10.5 for dr=1 km) will not converge

	#drin_temp/1000

        n_iter = n_iter + 1
        if n_iter>20:
	    #sprintf('ER11 CALCULATION DID NOT CONVERGE TO INPUT (RMAX,VMAX) = (#3.1f km,#3.1f m/s) Ck/Cd = #2.2f!',r_in_save/1000,Vmax_save,CkCd)
            V_ER11 = np.float('NaN')*np.zeros(rr_ER11.size)
            r_out = np.float('NaN')
            break

        ##Adjust estimate of r_in according to error
        r_in = r_in + drin_temp

        ##Vmax second
        #print('while2')
        while(np.abs(dVmax_temp/Vmax)>=10**-2): #if error is sufficiently large

            #          dVmax_temp

            ##Adjust estimate of Vmax according to error
            Vmax = Vmax + dVmax_temp

            [V_ER11,r_out] = ER11_radprof_raw(Vmax,r_in,rmax_or_r0,fcor,CkCd,rr_ER11)
            Vmax_prof = np.max(V_ER11)
            dVmax_temp = Vmax_save-Vmax_prof

        [V_ER11,r_out] = ER11_radprof_raw(Vmax,r_in,rmax_or_r0,fcor,CkCd,rr_ER11)
        Vmax_prof = np.max(V_ER11)
        dVmax_temp = Vmax_save-Vmax_prof
        if rmax_or_r0=='rmax':
            drin_temp = r_in_save-rr_ER11[np.argwhere(V_ER11==Vmax_prof)[0,0]]
        elif rmax_or_r0=='r0':
            f = interp1d( V_ER11[2:],rr_ER11[2:])
            pdb.set_trace()
	    #drin_temp = r_in_save-f(0)

    return V_ER11,r_out
