import copy
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from shapely.geometry import Point, LineString
import os.path

#================================= PROFILE FUNCTIONS =====================================

#Chavas
def ER11E04_nondim_rfitinput(Vmax, rfit, Vfit, fcor, Cdvary, C_d, w_cool, CkCdvary, CkCd, eye_adj, alpha_eye): 

    ## Initialization
    fcor = np.abs(fcor)
    if CkCdvary == 1:
        CkCd_coefquad = 5.5041e-04
        CkCd_coeflin  = -0.0259
        CkCd_coefcnst = 0.7627
        CkCd = CkCd_coefquad * Vmax ** 2 + CkCd_coeflin * Vmax + CkCd_coefcnst
    CkCd = np.min((1.9,CkCd))
    Mfit = rfit * Vfit + .5 * fcor * rfit ** 2
    soln_converged = 0
    
    while soln_converged == 0:
        rmaxrfit_min     = .01
        rmaxrfit_max     = 1.
        rmaxrfit_new     = (rmaxrfit_max + rmaxrfit_min) / 2 #first guess -- in the middle
        rmaxrfit         = rmaxrfit_new                      #initialize
        drmaxrfit        = rmaxrfit_max - rmaxrfit           #initialize
        drmaxrfit_thresh = .0001
        jterN            = 0

        while np.abs(drmaxrfit) >= drmaxrfit_thresh: #keep looping til changes in estimate are very small
            jterN = jterN + 1
            
            rmax  = rmaxrfit_new * rfit
            ## Step 1: Calculate ER11 M/Mm vs. r/rm
            #[~, ~, rrfracrm_ER11, MMfracMm_ER11] = ER11_radprof_nondim(Vmax, rmax, fcor, CkCdvary, CkCd)
            drfracrm = .01
            if(rmax > 100 * 1000):
                drfracrm = drfracrm / 10. # extra precision for large storms
                
            rfracrm_min    = 0.   #[-] start at r=0
            rfracrm_max    = 50.  #[-] extend out to many rmaxs
            rrfracrm_ER11  = np.arange(rfracrm_min, rfracrm_max + drfracrm, drfracrm) #[] r/r0 vector
            rr_ER11        = rrfracrm_ER11 * rmax
            rmax_or_r0     = 'rmax'
            VV_ER11, dummy = ER11_radprof(Vmax, rmax, rmax_or_r0, fcor, CkCd, rr_ER11)
            #fig=plt.figure()
            #plt.plot(rr_ER11, VV_ER11);
            #plt.grid()
            #plt.savefig('ER11.png')

            if not np.isnan(np.max(VV_ER11)):    # ER11_radprof converged
                Mm            = rmax * Vmax + .5 * fcor * rmax ** 2
                MMfracMm_ER11 = (rr_ER11 * VV_ER11 + .5 * fcor * rr_ER11 ** 2) / Mm
                rmaxr0_min    = .01
                rmaxr0_max    = .75
                rmaxr0_new    = (rmaxr0_max + rmaxr0_min) / 2 #first guess -- in the middle
                rmaxr0        = rmaxr0_new                    #initialize
                drmaxr0       = rmaxr0_max - rmaxr0           #initialize
                drmaxr0_thresh= .000001
                
                iterN = 0
                while np.abs(drmaxr0) >= drmaxr0_thresh:
                    iterN = iterN + 1
                    ##Calculate E04 M/M0 vs r/r0
                    r0 = rmax / rmaxr0_new   #[m]
                    Nr = 100000

                    
                    rrfracr0_E04, MMfracM0_E04 = E04_outerwind_r0input_nondim_MM0(r0, fcor, Cdvary, C_d, w_cool, Nr)
                                   
                    
                    ##Convert ER11 to M/M0 vs. r/r0 space
                    rrfracr0_ER11 = rrfracrm_ER11 * (rmaxr0_new)
                    M0_E04        = .5 * fcor * r0 ** 2
                    MMfracM0_ER11 = MMfracMm_ER11 * (Mm / M0_E04)
                    
                    #fig=plt.figure()
                    #plt.plot(rrfracr0_E04,MMfracM0_E04)
                    #plt.plot(rrfracr0_ER11,MMfracM0_ER11)
                    #plt.grid();plt.ylim([0,1]);plt.xlim([0,1])
                    #plt.savefig('moment_'+str(iterN)+'.png')
                    
                    l1 = LineString(zip(rrfracr0_E04 , MMfracM0_E04))
                    l2 = LineString(zip(rrfracr0_ER11 , MMfracM0_ER11))
                    intersection = l1.intersection(l2)
                    #if intersection.wkt == 'GEOMETRYCOLLECTION EMPTY':   ##no intersections r0 too large --> rmaxr0 too small   
                    if intersection.wkt == 'LINESTRING EMPTY':   ##no intersections r0 too large --> rmaxr0 too small
                        drmaxr0 = np.abs(drmaxr0) / 2
                    else:			##at least one intersection -- r0 too small --> rmaxr0 too large
                        if intersection.wkt.split(' ')[0] == 'POINT':
                            X0,Y0 = intersection.coords[0]
                        elif intersection.wkt.split(' ')[0] == 'MULTIPOINT':
                            #X0,Y0 = intersection[0].coords[0]
                            X0 = [];Y0 =[]
                            for ikt in np.arange(len(intersection)):
                                X0.append(intersection[ikt].coords[0][0])
                                Y0.append(intersection[ikt].coords[0][1])
                            X0 = np.array(X0); Y0 = np.array(Y0)

                        #at least one intersection -- rmaxr0 too large
                        drmaxr0  = -np.abs(drmaxr0) / 2
                        rmerger0 = np.mean(X0)
                        MmergeM0 = np.mean(Y0)
                        
	            ###update value of rmaxr0
                    rmaxr0     = rmaxr0_new    #this is the final one
                    rmaxr0_new = rmaxr0_new + drmaxr0

                    #print(iterN,rmaxr0,rmaxr0_new,drmaxr0)
                    #print(X0,Y0)
                    #print(intersection.wkt, rmerger0, MmergeM0)

                ## Calculate some things
                M0   = .5 * fcor * r0 ** 2
                Mm   = .5 * fcor * rmax ** 2 + rmax * Vmax
                MmM0 = Mm / M0

                ## Define merged solution
                ii_ER11       = np.argwhere((rrfracr0_ER11 < rmerger0) & (MMfracM0_ER11 < MmergeM0))[:, 0]
                ii_E04        = np.argwhere((rrfracr0_E04 >= rmerger0) & (MMfracM0_E04 >= MmergeM0))[:, 0]
                MMfracM0_temp = np.hstack((MMfracM0_ER11[ii_ER11], MMfracM0_E04[ii_E04]))
                rrfracr0_temp = np.hstack((rrfracr0_ER11[ii_ER11], rrfracr0_E04[ii_E04]))
                del ii_ER11
                del ii_E04

                ## Check to see how close solution is to input value of (rfitr0, MfitM0)
                rfitr0 = rfit / r0
                MfitM0 = Mfit / M0

                ## 2020-06-23 fixed bug returning NaN (and wrong solution) if rfit > current r0
                if rfitr0 <= 1:
                    f           = interp1d(rrfracr0_temp, MMfracM0_temp)
                    #print MMfracM0_temp
                    MfitM0_temp = f(rfitr0)
                    MfitM0_err  = MfitM0 - MfitM0_temp
                else:    #true rfit exceeds current r0, so doesnt exist in profile!
                    MfitM0_err  = 1000000;  #simply need smaller r0 -- set MfitM0_err to any positive number
                #print(rfitr0, MfitM0, MfitM0_err)
                if MfitM0_err > 0:    #need smaller rmax (r0)
                    drmaxrfit = np.abs(drmaxrfit) / 2.
                else:    #need larger rmax (r0)
                    drmaxrfit = -np.abs(drmaxrfit) / 2.
                #pdb.set_trace()

                
            else:    #ER11_radprof did not converge -- convergence fails for low CkCd and high Ro = Vm/(f*rm)
                ##Must reduce rmax (and thus reduce Ro)
                drmaxrfit = -abs(drmaxrfit) / 2

            ##update value of rmaxrfit
            rmaxrfit     = rmaxrfit_new    #this is the final one
            rmaxrfit_new = rmaxrfit + drmaxrfit
             
        ##Check if solution converged
        if not np.isnan(np.max(VV_ER11)):
            soln_converged = 1
        else:
            soln_converged = 0
            CkCd           = CkCd + .1
            #print('Adjusting CkCd to find convergence')

    ## Finally: Interpolate to a grid
    # drfracr0 = .0001
    # rfracr0_min = 0    #[-] r=0
    # rfracr0_max = 1   #[-] r=r0
    # rrfracr0 = rfracr0_min:drfracr0:rfracr0_max #[] r/r0 vector
    # MMfracM0 = interp1(rrfracr0_temp,MMfracM0_temp,rrfracr0,'pchip',NaN)
    drfracrm    = .01                                           #calculating VV at radii relative to rmax ensures no smoothing near rmax!
    rfracrm_min = 0.                                            #[-] r=0
    rfracrm_max = r0 / rmax                                     #[-] r=r0
    rrfracrm    = np.arange(rfracrm_min, rfracrm_max, drfracrm) #[] r/r0 vector
    f           = interp1d(rrfracr0_temp * (r0 / rmax), MMfracM0_temp * (M0 / Mm))
    MMfracMm    = f(rrfracrm)
    rrfracr0    = rrfracrm * rmax / r0                          #save this as output
    MMfracM0    = MMfracMm * Mm / M0

    ## Calculate dimensional wind speed and radii
    # VV = (M0/r0)*((MMfracM0./rrfracr0)-rrfracr0)  #[ms-1]
    # rr = rrfracr0*r0   #[m]
    # rmerge = rmerger0*r0
    # Vmerge = (M0/r0)*((MmergeM0./rmerger0)-rmerger0)  #[ms-1]
    VV = (Mm / rmax) * (MMfracMm / rrfracrm) - .5 * fcor * rmax * rrfracrm #[ms-1]
    rr = rrfracrm * rmax #[m]

    ## Make sure V=0 at r=0
    VV[rr==0] = 0

    rmerge = rmerger0 * r0
    Vmerge = (M0 / r0) * ((MmergeM0 / rmerger0) - rmerger0)  #[ms-1]
    ## Adjust profile in eye, if desired
    #if(eye_adj==1)
    #    r_eye_outer = rmax
    #    V_eye_outer = Vmax
    #    [VV] = radprof_eyeadj(rr,VV,alpha_eye,r_eye_outer,V_eye_outer)
    #
    #    sprintf('EYE ADJUSTMENT: eye alpha = #3.2f',alpha_eye)
    return rr, VV, rmax, r0, rmerge, Vmerge
    
def ER11_radprof(Vmax, r_in, rmax_or_r0, fcor, CkCd, rr_ER11):
    '''Given a guess of Vmax and rmax as input (r_in), returns the associated INNER radial profile of V following ER11 theory (V_ER11) as well as the radius of 0 wind speed (r_out).
    Then uses this estimated profile V_ER11 and update both the estimates of rmax (using drin_temp the gap between the previous and current estimate);
                                                                         and Vmax (using dVmax_temp the gap between the previous and current estimate).
    Then he re-computes the theoretical profile V_ER11 (again using ER11_radprof_raw() until the rmax and Vmax values converge (i.e is they don't change that muchat each update).'''
    dr = rr_ER11[1] - rr_ER11[0]
    ## Call ER11_radprof_raw
    V_ER11, r_out = ER11_radprof_raw(Vmax, r_in, rmax_or_r0, fcor, CkCd, rr_ER11)
    
    if rmax_or_r0 == 'rmax':
        drin_temp = r_in - rr_ER11[np.argwhere(V_ER11 == np.max(V_ER11))[0,0]]
    elif rmax_or_r0 == 'r0':
        f         = interp1d( V_ER11[2:], rr_ER11[2:])
        drin_temp = r_in - f(0).tolist() 
    ## Calculate error in Vmax
    dVmax_temp = Vmax - np.max(V_ER11)

    ## Check if errors are too large and adjust accordingly
    r_in_save = copy.copy(r_in)
    Vmax_save = copy.copy(Vmax)

    n_iter = 0
    while((np.abs(drin_temp) > (dr / 2)) or (np.abs(dVmax_temp / Vmax_save) >= 10 ** -2)): #if error is sufficiently large NOTE: FIRST ARGUMENT MUST BE ">" NOT ">=" or else rmax values at exactly dr/2 intervals (e.g. 10.5 for dr=1 km) will not converge

	# drin_temp/1000

        n_iter = n_iter + 1
        if n_iter>20:
	    # sprintf('ER11 CALCULATION DID NOT CONVERGE TO INPUT (RMAX,VMAX) = (#3.1f km,#3.1f m/s) Ck/Cd = #2.2f!',r_in_save/1000,Vmax_save,CkCd)
            V_ER11 = np.float('NaN') * np.zeros(rr_ER11.size)
            r_out  = np.float('NaN')
            break

        ## Adjust estimate of r_in according to error
        r_in = r_in + drin_temp

        ##Vmax second
        while(np.abs(dVmax_temp / Vmax) >= 10 ** -2): # if error is sufficiently large

            #          dVmax_temp

            ## Adjust estimate of Vmax according to error
            Vmax            = Vmax + dVmax_temp

            [V_ER11, r_out] = ER11_radprof_raw(Vmax, r_in, rmax_or_r0, fcor, CkCd, rr_ER11)
            Vmax_prof       = np.max(V_ER11)
            dVmax_temp      = Vmax_save - Vmax_prof

        [V_ER11, r_out] = ER11_radprof_raw(Vmax, r_in, rmax_or_r0, fcor, CkCd, rr_ER11)
        Vmax_prof       = np.max(V_ER11)
        dVmax_temp      = Vmax_save - Vmax_prof
        if rmax_or_r0 == 'rmax':
            drin_temp = r_in_save - rr_ER11[np.argwhere(V_ER11 == Vmax_prof)[0,0]]
        elif rmax_or_r0=='r0':
            f = interp1d( V_ER11[2:], rr_ER11[2:])
	    #drin_temp = r_in_save-f(0)

    return V_ER11,r_out

def ER11_radprof_raw(Vmax, r_in, rmax_or_r0, fcor, CkCd, rr_ER11, eyealpha=1):
    '''Given a guess of rmax as input (r_in), returns the associated INNER radial profile of V following ER11 theory (V_ER11).
    Also returns the estimated radius of 0 wind speed (r_out), estimated using the found profile.'''
    fcor = np.abs(fcor)
    if rmax_or_r0 == 'rmax':
        rmax = r_in
    else:
        print('rmax_or_r0 must be set to "rmax"')
        
    ## CALCULATE Emanuel and Rotunno (2011) theoretical profile
    V_ER11 = (1. / rr_ER11) * (Vmax * rmax + .5 * fcor * rmax ** 2) * ((2 * (rr_ER11 / rmax) ** 2) / (2 - CkCd + CkCd * (rr_ER11 / rmax) ** 2)) ** (1 / (2 - CkCd)) - .5 * fcor * rr_ER11
    ## make V=0 at r=0
    V_ER11[rr_ER11 == 0] = 0
   
    if rmax_or_r0 == 'rmax':
        i_rmax     = np.argwhere(V_ER11 == np.max(V_ER11))[0,0]
        f          = interp1d(V_ER11[i_rmax + 1:], rr_ER11[i_rmax + 1:], fill_value='extrapolate')
        r0_profile = f(0.)
        r_out      = r0_profile.tolist() # use value from profile itself
    else:
        print('rmax_or_r0 must be set to"rmax"')
    return V_ER11, r_out

def E04_outerwind_r0input_nondim_MM0(r0, fcor, Cdvary, C_d, w_cool, Nr):
    
    ## Initialization
    fcor = abs(fcor)
    M0   = .5 * fcor * r0 ** 2 #[m2/s] M at outer radius

    drfracr0 = .001
    if ((r0 > 2500 * 1000) | (r0 < 200 * 1000)):
        drfracr0 = drfracr0/10 
        #extra precision for very large storm to avoid funny bumps near r0 (though rest of solution is stable!)
        #or for tiny storm that requires E04 extend to very small radii to match with ER11

    if Nr > 1 / drfracr0:
        Nr = 1 / drfracr0 #grid radii must be > 0

    rfracr0_max  = 1                                                        #[-] start at r0, move radially inwards
    rfracr0_min  = rfracr0_max - (Nr - 1) * drfracr0                        #[-] inner-most node
    rrfracr0     = np.arange(rfracr0_min, rfracr0_max + drfracr0, drfracr0) #[] r/r0 vector
    MMfracM0     = np.float('NaN') * np.zeros(rrfracr0.size)                #[] M/M0 vector initialized to 1 (M/M0 = 1 at r/r0=1)
    MMfracM0[-1] = 1

    ## First step inwards from r0: d(M/M0)/d(r/r0) = 0 by definition
    rfracr0_temp = rrfracr0[-2] #one step inwards from r0
    #dMfracM0_drfracr0_temp = 0  #[] d(M/M0)/d(r/r0) = 0 at r/r0 = 1
    MfracM0_temp = MMfracM0[-1]
    MMfracM0[-2] = MfracM0_temp
    ##################################################################
    ##Variable C_d: code from Cd_Donelan04.m (function call is slow) ######
    ##Piecewise linear fit parameters estimated from Donelan2004_fit.m
    C_d_lowV     = 6.2e-4
    V_thresh1    = 6    #m/s transition from constant to linear increasing
    V_thresh2    = 35.4 #m/s transition from linear increasing to constant
    C_d_highV    = 2.35e-3
    linear_slope = (C_d_highV - C_d_lowV) / (V_thresh2 - V_thresh1)
    ##################################################################

    ## Integrate inwards from r0 to obtain profile of M/M0 vs. r/r0
    for ii in range(0, np.int(Nr) - 2, 1): #first two nodes already done above
        
        ## Calculate C_d varying with V, if desired
        if Cdvary==1:
            
            ##Calculate V at this r/r0 (for variable C_d only)
            V_temp = (M0 / r0) * ((MfracM0_temp / rfracr0_temp) - rfracr0_temp)
            
            ##Calculate C_d
            if V_temp <= V_thresh1:
                C_d = C_d_lowV
            elif V_temp > V_thresh2:
                C_d = C_d_highV
            else:
                C_d = C_d_lowV + linear_slope * (V_temp - V_thresh1)

        ## Calculate model parameter, gamma
        gam = C_d * fcor * r0 / w_cool   #[] non-dimensional model parameter

        ## Update dMfracM0_drfracr0 at next step inwards
        dMfracM0_drfracr0_temp = gam * ((MfracM0_temp - rfracr0_temp ** 2) ** 2) / (1 - rfracr0_temp ** 2)

        ## Integrate M/M0 radially inwards
        MfracM0_temp = MfracM0_temp - dMfracM0_drfracr0_temp * drfracr0

        ## Update r/r0 to follow M/M0
        rfracr0_temp = rfracr0_temp - drfracr0 #[] move one step inwards

        ## Save updated values
        MMfracM0[MMfracM0.shape[0] - 1 - ii - 2] = MfracM0_temp

    return rrfracr0,MMfracM0

def initialize_chavas(spdm, Lat, print_values):
    '''Initialize the values of Vmax, Rfit, Vfit, fcor, Cdvary, Cd, w_cool, CkCdvary, CkCd, eye_adj, alpha_eye for the Chavas profile.'''
    ### STORM PARAMETERS (Vmax, Vfit, Rfit, and Lat)
    Vmax = np.max(spdm[:200])                      #[ms-1] {50} maximum azimuthal-mean wind speed
    Vfit = 17                                      #[ms-1] {12} wind speed at Rfit
    Rmax = np.argmax(spdm[:200])                   #To compute Rfit
    diff = np.abs(spdm - Vfit)                     #To compute Rfit
    Rfit = (np.argmin(diff[Rmax:]) + Rmax) * 1000  #[m] {300*1000} a wind radius
    Lat  = Lat                                     #[°]         latitude of storm center

    ### DEFAULT PARAMETERS
    fcor = coriolis(Lat)         #[s-1]  {5e-5} Coriolis parameter at storm center

    ### ENVIRONMENTAL PARAMETERS
    # Outer region
    Cdvary    = 1                  #[-]    {1} 0 : Outer region Cd = constant (defined on next line) 1 : Outer region Cd = f(V) (empirical 	Donelan et al. 2004)
    Cd        = 1.5e-3             #[-]    {1.5e-3} ignored if Cdvary = 1 surface momentum exchange (i.e. drag) coefficient
    w_cool    = 2/1000             #[ms-1] {2/1000 Chavas et al 2015} radiative-subsidence rate in the rain-free tropics above the boundary layer top
    # Inner region
    CkCdvary  = 1                  #[-]    {1} 0 : Inner region Ck/Cd = constant (defined on next line) 1 : Inner region Ck/Cd = f(Vmax) (empirical Chavas et al. 2015)
    CkCd      = 1                  #[-]    {1} ignored if CkCdvary = 1 ratio of surface exchange coefficients of enthalpy and momentum capped at 1.9 (things get weird >=2)
    # Eye adjustment
    eye_adj   = 0                  #[-]    {1} 0 = use ER11 profile in eye 1 = empirical adjustment
    alpha_eye = .15                #[-]    {.15 empirical Chavas et al 2015} V/Vm in eye is reduced by factor (r/rm)^alpha_eye ignored if eye_adj=0
    if print_values:
        print(
            'CHAVAS - Initial values',
            '\n Vmax_ini  =', "{:.2f}".format(Vmax),
            '\n Rfit_ini  =', "{:.2f}".format(Rfit),
            '\n Vfit_ini  =', "{:.2f}".format(Vfit),
            '\n fcor      =', "{:.2f}".format(fcor),
            '\n Cdvary    =', "{:.2f}".format(Cdvary),
            '\n Cd        =', "{:.2f}".format(Cd),
            '\n w_cool    =', "{:.2f}".format(w_cool),
            '\n CkCdvary  =', "{:.2f}".format(CkCdvary),
            '\n CkCd      =', "{:.2f}".format(CkCd),
            '\n eye_adj   =', "{:.2f}".format(eye_adj),
            '\n alpha_eye =', "{:.2f}".format(alpha_eye)
        )
    return Vmax, Rfit, Vfit, fcor, Cdvary, Cd, w_cool, CkCdvary, CkCd, eye_adj, alpha_eye

def fit_chavas(Vmax, Rfit, Vfit, fcor, Cdvary, Cd, w_cool, CkCdvary, CkCd, eye_adj, alpha_eye, print_values):
    '''Returns the VV and rr associated to Chavas (different than other profiles), as well as optimal Rmax, r0, Rmerge and Vmerge.'''
    rr, VV, rmax, r0, rmerge, Vmerge = ER11E04_nondim_rfitinput(Vmax, Rfit, Vfit, fcor, Cdvary, Cd, w_cool, CkCdvary, CkCd, eye_adj, alpha_eye)
    if print_values:
        print(
            'CHAVAS - Fit values',
            '\n Rmax_fit   =', "{:.2f}".format(Rmax),
            '\n R0_fit     =', "{:.2f}".format(r0),
            '\n Rmerge_fit =', "{:.2f}".format(rmerge),
            '\n Vmerge_fit =', "{:.2f}".format(Vmerge)

        )
    return rr, VV, rmax, r0, rmerge, Vmerge

# Rankine
def rankine_profile(r, x, alpha, Vmin, Rmax):
    '''We assume V = alpha * r inside and V * r ^ x = alpha outside'''
    V   = r * 0.
    Vin = Vmin + alpha * r
    Vou = Vmin + (alpha * Rmax) * (Rmax ** x) / (r ** x)
    # Vou = (Vmin + alpha * Rmax) * (Rmax ** x) / (r ** x)
    V[r <=Rmax] = Vin[r <=Rmax]
    V[r > Rmax] = Vou[r > Rmax]
    return V

def initialize_rankine(spdm, x, alpha, print_values):
    '''Initialize the values of x, alpha, Vmin and Rmax for the Rankine profile.'''
    Vmin   = spdm[0] # 14
    Rmax   = np.argmax(spdm[:200]) # 40
    if Rmax < 5: # sometimes the TC is not well-centered and the center of the grid is on the eyewall. In this case we set Rmax = 5 because in the optimizing process we constrain Rmax > 5.
        Rmax = 5.01
    if print_values:
        print(
            'RANKINE - Initial values',
            '\n x_ini     =', "{:.2f}".format(x),
            '\n alpha_ini =', "{:.2f}".format(alpha),
            '\n Vmin_ini  =', "{:.2f}".format(Vmin),
            '\n Rmax_ini  =', "{:.2f}".format(Rmax)        
        )
    return x, alpha, Vmin, Rmax

def fit_rankine(r, spdm, x, alpha, Vmin, Rmax, print_values):
    '''Fit the rankine profile given initial values of x, alpha, Vmin and Rmax.
    Returns the optimal parameters found with curve_fit()'''
    popt, pcov                           = curve_fit(rankine_profile, r, spdm, p0=[x, alpha, Vmin, Rmax], bounds=((0, 0, 0, 5), (1, 50, 50, 500)))
    x_fit, alpha_fit, Vmin_fit, Rmax_fit = popt[0], popt[1], popt[2], popt[3]
    if print_values:
        print(
            'RANKINE - Fit values',
            '\n x_fit     =', "{:.2f}".format(x_fit),
            '\n alpha_fit =', "{:.2f}".format(alpha_fit),
            '\n Vmin_fit  =', "{:.2f}".format(Vmin_fit),
            '\n Rmax_fit  =', "{:.2f}".format(Rmax_fit)        
        )
    return x_fit, alpha_fit, Vmin_fit, Rmax_fit

# Holland
def holland_profile(r, Lat, pn, pc, Vmin, Rmax, Vmax):
    '''We assume that rho is constant and equals 1.15 kg.m-3'''
    fcor  = coriolis(Lat)
    rho   = 1.15
    B     = (Vmax ** 2) * rho * np.e / (pn - pc)
    A     = Rmax ** B
    
    V   = r * 0.
    V   = Vmin + np.sqrt(A * B * (pn - pc) * np.exp((-1) * A / (r ** B)) / (rho * r ** B) + (r ** 2 * fcor ** 2) / 4) - (r * fcor / 2)
    return V

def initialize_holland(spdm, Lat, pn, pc, print_values):
    '''The user chooses the values for Lat, pn (in mbar) and pc (mbar).'''
    Lat    = Lat                   #[°] latitude of storm center
    pn     = pn * 100.             # to convert mbar pressures in Pa pressures
    pc     = pc * 100.             # to convert mbar pressures in Pa pressures
    Vmin   = spdm[0]               # 14
    Rmax   = np.argmax(spdm[:200]) # 40
    if Rmax < 5: # sometimes the TC is not well-centered and the center of the grid is on the eyewall. In this case we set Rmax = 5 because in the optimizing process we constrain Rmax > 5.
        Rmax = 5.01
    Vmax   = np.max(spdm[:200])    # 48
    if print_values:
        print(
            'HOLLAND - Initial values',
            '\n Lat       =', "{:.2f}".format(Lat),
            '\n pn_ini    =', "{:.2f}".format(pn),
            '\n pc_ini    =', "{:.2f}".format(pc),
            '\n Vmin_ini  =', "{:.2f}".format(Vmin),
            '\n Rmax_ini  =', "{:.2f}".format(Rmax),
            '\n Vmax_ini  =', "{:.2f}".format(Vmax)
        )
    return Lat, pn, pc, Vmin, Rmax, Vmax

def fit_holland(r, spdm, Lat, pn, pc, Vmin, Rmax, Vmax, print_values):
    '''Fit the Holland profile given initial values of Lat, pn, pc, Vmin, Rmax and Vmax.
    Returns the optimal parameters found with curve_fit()'''
    popt, pcov = curve_fit(lambda r, pn, pc, Vmin, Rmax, Vmax: holland_profile(r, Lat, pn, pc, Vmin, Rmax, Vmax), r, spdm, p0=[pn, pc, Vmin, Rmax, Vmax], bounds=((1000 * 100, 850 * 100, 0, 5, 0), (1100 * 100, 1000 * 100, 50, 500, 200))) # Lat is fixed
    pn, pc, Vmin, Rmax, Vmax = popt[0], popt[1], popt[2], popt[3], popt[4]
    if print_values:
        print(
            'HOLLAND - Fit values',
            '\n Lat       =', "{:.2f}".format(Lat),
            '\n pn_fit    =', "{:.2f}".format(pn),
            '\n pc_fit    =', "{:.2f}".format(pc),
            '\n Vmin_fit  =', "{:.2f}".format(Vmin),
            '\n Rmax_fit  =', "{:.2f}".format(Rmax),
            '\n Vmax_fit  =', "{:.2f}".format(Vmax)
        )
    return Lat, pn, pc, Vmin, Rmax, Vmax

# Willoughby
def willoughby_profile_no_smooth(r, n, X1, Vmin, Rmax, Vmax):
    '''No polynomial ramp smoothing here.
    We assume V(0) = Vmin != 0 to fit SAR data'''
    V    = r * 0.
    Vinf = (Vmax - Vmin) * ((r / Rmax) ** n) + Vmin
    Vsup = (Vmax - Vmin) * np.exp((-1.) * ((r - Rmax) / X1)) + Vmin
    V[r <=Rmax] = Vinf[r <=Rmax]
    V[r > Rmax] = Vsup[r > Rmax]
    return V

def initialize_willoughby(spdm, n, print_values):
    '''Initialize the values of n, X1, Vmin, Rmax and Vmax for the Willoughby profile.
    By default X1 = 3 * Rmax'''
    n      = n
    Vmin   = spdm[0] # 14
    Rmax   = np.argmax(spdm[:200]) # 40
    X1     = 3 * Rmax
    if Rmax < 5: # sometimes the TC is not well-centered and the center of the grid is on the eyewall. In this case we set Rmax = 5 because in the optimizing process we constrain Rmax > 5.
        Rmax = 5.01
    Vmax   = np.max(spdm[:200])    # 48
    if print_values:
        print(
            'WILLOUGHBY - Initial values',
            '\n n_ini    =', "{:.2f}".format(n),
            '\n X1_ini   =', "{:.2f}".format(X1),
            '\n Vmin_ini =', "{:.2f}".format(Vmin),
            '\n Rmax_ini =', "{:.2f}".format(Rmax),
            '\n Vmax_ini =', "{:.2f}".format(Vmax),
        )
    return n, X1, Vmin, Rmax, Vmax

def fit_willoughby_no_smooth(r, spdm, n, X1, Vmin, Rmax, Vmax, print_values):
    '''Fit the Willoughby profile given initial values of n, X1, Vmin, Rmax, Vmax.
    Returns the optimal parameters found with curve_fit()'''
    popt, pcov              = curve_fit(willoughby_profile_no_smooth, r, spdm, p0=[n, X1, Vmin, Rmax, Vmax], bounds=((0, 0, 0, 5, 0), (50, 5000, 50, 500, 200)))
    n, X1, Vmin, Rmax, Vmax = popt[0], popt[1], popt[2], popt[3], popt[4]
    if print_values:
        print(
            'WILLOUGHBY - Fit values',
            '\n n_fit    =', "{:.2f}".format(n),
            '\n X1_fit   =', "{:.2f}".format(X1),
            '\n Vmin_fit =', "{:.2f}".format(Vmin),
            '\n Rmax_fit =', "{:.2f}".format(Rmax),
            '\n Vmax_fit =', "{:.2f}".format(Vmax)
        )
    return n, X1, Vmin, Rmax, Vmax

#================================= OTHER FUNCTIONS =====================================
def printWS_ggd(path):
    '''Given the path of a _gd file, print the wind speed of the image'''
    ds = xr.open_dataset(path)
    # No need for np.array()?
    ws = np.array(ds['wind_speed'])
    xx = np.array(ds['x'])
    yy = np.array(ds['y'])
    plt.pcolormesh(xx, yy, ws[0, :, :])
    
def printWS(path):
    '''Given the path of a _rotated file, print the wind speed of the image'''
    ds = xr.open_dataset(path)
    plt.pcolormesh(ds['x'], ds['y'], ds['wind_speed'])
    
def coriolis(lat):
    '''Latitude must be in degrees.'''
    Omega = 7.2921e-5                            # Earth rotation vector
    f     = 2 * Omega * np.sin(lat * np.pi / 180) # Coriolis parameter at 20° latitude and assuming it's constant 
    return f

def compute_mean_wind_spd(ds, r_window_len):
    '''Returns azimuthal-mean total (not azimuthal) wind speed'''
    # Define (r, theta) grid
    r     = np.arange(r_window_len)
    th    = np.arange(361)
    r, th = np.meshgrid(r, th)
    ds_r  = np.array(ds['r_polar'])

    ds_th = np.mod(np.array(ds['theta']) * 180. / np.pi, 360) # convert theta from radians to degrees
    ds_ws = np.array(ds['wind_speed'])
    # Possible to call griddata() without using meshgrid() before? 
    spd   = griddata((ds_r.flatten(), ds_th.flatten()), ds_ws.flatten(), (r, th), method='nearest')
    spdm   = np.nanmean(spd, axis=0)  
    return spdm

def compute_mean_tangential_wind_spd(ds, r_window_len):
    '''Returns azimuthal-mean tangential wind speed'''
    # Define (r, theta) grid
    r      = np.arange(r_window_len)
    th     = np.arange(361)
    r, th  = np.meshgrid(r, th)
    ds_r   = np.array(ds['r_polar'])

    ds_th  = np.mod(np.array(ds['theta']) * 180. / np.pi, 360) # convert theta from radians to degrees
    ds_tws = np.array(ds['wind_speed'])
    ds_aws = np.abs(np.array(ds['tangential_wind'])) # normed azimuthal wind
    ds_ws  = np.multiply(ds_tws, ds_aws)             # azimuthal wind
    # Possible to call griddata() without using meshgrid() before? 
    spd    = griddata((ds_r.flatten(), ds_th.flatten()), ds_ws.flatten(), (r, th), method='nearest')
    spdm_ch= np.nanmean(spd, axis=0)  
    return spdm_ch

def initialize_radius(spdm):
    '''Given the spdm, returns the largest radius (and asociated spdm) on which the profile can be fitted. 
    Indeed, sometimes the spdm isn't defined from r = 0 to r = 500, in this case the largest domain is taken instead.'''
    if np.isnan(np.min(spdm)):
        # last_valid_index = (~np.isnan(spdm)).cumsum(0).argmax(0)
        last_valid_index = np.min(np.where(np.isnan(spdm))[0])# returns the index of the last valid value before the first nan
        r                = np.arange(last_valid_index) + .0001 # to avoid having both r = 0 and n < 0 during fitting process
    else:
        last_valid_index = 501
        r                = np.arange(last_valid_index) + .0001 # to avoid having both r = 0 and n < 0 during fitting process
    spdm = spdm[:last_valid_index]
    return r, spdm

def plot_curves(i, file, r, spdm, FIT):
    # Compute fitted profiles
    V_rankine              = rankine_profile(r, *FIT['Rankine'])
    V_holland              = holland_profile(r, *FIT['Holland'])
    V_willoughby_no_smooth = willoughby_profile_no_smooth(r, *FIT['Willoughby'])
    V_chavas               = FIT['Chavas'][1]         # Different from other profiles: here V_chavas has already been computed before, and is stored in FIT['Chavas'][1]
    r_chavas               = FIT['Chavas'][0] / 1000. # Convert from m to km
    
    Rmax          = np.argmax(spdm[:200]) # center on Rmax
    # Compute indexes to print in the right window (for Chavas profile only)
    upper_bound   = np.int(np.floor(len(r_chavas) * len(spdm) / (FIT['Chavas'][3] / 1000.))) # to not plot the entire Chavas profile (because it goes until r0 which is larger than 500 km in the general case)
    index25       = np.int(np.floor(len(r_chavas) * 25 / (FIT['Chavas'][3] / 1000.)))
    upper_bound50 = np.int(np.floor(len(r_chavas) * 50 / (FIT['Chavas'][3] / 1000.)))
    Rmax_chavas   = np.int(np.floor(len(r_chavas) * Rmax / (FIT['Chavas'][3] / 1000.)))
    
    # Title
    plt.figure(figsize=(18, 8))
    plt.suptitle('N°' + '%d'%i + ": " + os.path.basename(file), fontsize=14)
    
    # Large scale
    plt.subplot(1, 2, 1)
    fig1 = plt.plot(r, spdm, color='k', linewidth=3,            label='SAR-derived wind speed')    # V_obs
    fig2 = plt.plot(r, V_rankine, color='darkorange',           label='Rankine profile')           # V_rankine
    fig3 = plt.plot(r, V_holland, color='steelblue',            label='Holland profile')           # V_holland
    fig4 = plt.plot(r, V_willoughby_no_smooth, color='orchid',  label='Willoughby - no smoothing') # V_willoughby_no_smooth
    fig5 = plt.plot(r_chavas[:upper_bound], V_chavas[:upper_bound], color='forestgreen',   label='Chavas profile') # V_chavas    
    
    plt.xlabel('Radius (km)')
    plt.ylabel('Wind speed (m/s)')
    plt.legend();plt.grid()

    # Small scale
    plt.subplot(1, 2, 2)
    if Rmax >= 25:
        fig6 = plt.plot(r[Rmax - 25 : Rmax + 25], spdm[Rmax - 25 : Rmax + 25], color='k', linewidth=3, label='SAR-derived wind speed') # V_obs
        fig7 = plt.plot(r[Rmax - 25 : Rmax + 25], V_rankine[Rmax - 25 : Rmax + 25], color='darkorange', label='Rankine profile')       # V_rankine
        fig8 = plt.plot(r[Rmax - 25 : Rmax + 25], V_holland[Rmax - 25 : Rmax + 25], color='steelblue',  label='Holland profile')       # V_holland
        fig9 = plt.plot(r[Rmax - 25 : Rmax + 25], V_willoughby_no_smooth[Rmax - 25 : Rmax + 25], color='orchid',  label='Willoughby - no smoothing') # V_willoughby_no_smooth
        fig10= plt.plot(r_chavas[Rmax_chavas - index25 : Rmax_chavas + index25], V_chavas[Rmax_chavas - index25 : Rmax_chavas + index25], color='forestgreen',   label='Chavas profile') # V_chavas
    else:
        fig6 = plt.plot(r[:50], spdm[:50], color='k', linewidth=3, label='SAR-derived wind speed') # V_obs
        fig7 = plt.plot(r[:50], V_rankine[:50], color='darkorange', label='Rankine profile')       # V_rankine
        fig8 = plt.plot(r[:50], V_holland[:50], color='steelblue',  label='Holland profile')       # V_holland
        fig9 = plt.plot(r[:50], V_willoughby_no_smooth[:50], color='orchid',  label='Willoughby - no smoothing') # V_willoughby_no_smoothing
        fig10= plt.plot(r_chavas[:upper_bound50], V_chavas[:upper_bound50], color='forestgreen',   label='Chavas profile') # V_chavas

    plt.xlabel('Radius (km)')
    plt.ylabel('Wind speed (m/s)')
    plt.legend();plt.grid()
    return True

def save_curves(i, file, ds, r, spdm, INI, FIT):
    # Save path and name
    i        = '{0:03}'.format(i) # convert 1 to '001'
    filename = os.path.basename(os.path.splitext(file)[0])
    filename = 'wProfiles' + i + '_' + filename
    savepath = "/home/arthur/results/windProfiles/" + filename
    
    # Title
    plt.figure(figsize=(25, 19))
    plt.suptitle('N°' + i + ": " + os.path.basename(file), fontsize=14)
    
    # Print TC and spd
    plt.subplot(2, 2, 1)
    plt.pcolormesh(ds['x_coords'], ds['y_coords'], ds['wind_speed']);plt.grid()
    plt.subplot(2, 2, 2)
    radius     = np.arange(501)
    th         = np.arange(361)
    radius, th = np.meshgrid(radius, th)
    ds_r       = np.array(ds['r_polar'])
    ds_th      = np.mod(np.array(ds['theta']) * 180. / np.pi, 360) # convert theta from radians to degrees
    ds_ws      = np.array(ds['wind_speed'])
    # Possible to call griddata() without using meshgrid() before? 
    spd        = griddata((ds_r.flatten(), ds_th.flatten()), ds_ws.flatten(), (radius, th), method='nearest')

    plt.pcolormesh(spd)
    
    # Compute fitted profiles
    V_rankine              = rankine_profile(r, *FIT['Rankine'])
    V_holland              = holland_profile(r, *FIT['Holland'])
    V_willoughby_no_smooth = willoughby_profile_no_smooth(r, *FIT['Willoughby'])
    V_chavas               = FIT['Chavas'][1]         # Different from other profiles: here V_chavas has already been computed before, and is stored in FIT['Chavas'][1]
    r_chavas               = FIT['Chavas'][0] / 1000. # Convert from m to km
    
    Rmax          = np.argmax(spdm[:200]) # center on Rmax
    # Compute indexes to print in the right window (for Chavas profile only)
    upper_bound   = np.int(np.floor(len(r_chavas) * len(spdm) / (FIT['Chavas'][3] / 1000.))) # to not plot the entire Chavas profile (because it goes until r0 which is larger than 500 km in the general case)
    index25       = np.int(np.floor(len(r_chavas) * 25 / (FIT['Chavas'][3] / 1000.)))
    upper_bound50 = np.int(np.floor(len(r_chavas) * 50 / (FIT['Chavas'][3] / 1000.)))
    Rmax_chavas   = np.int(np.floor(len(r_chavas) * Rmax / (FIT['Chavas'][3] / 1000.)))
    
    # Large scale
    plt.subplot(2, 2, 3)
    fig1 = plt.plot(r, spdm, color='k', linewidth=3,            label='SAR-derived wind speed')    # V_obs
    fig2 = plt.plot(r, V_rankine, color='darkorange',           label='Rankine profile')           # V_rankine
    fig3 = plt.plot(r, V_holland, color='steelblue',            label='Holland profile')           # V_holland
    fig4 = plt.plot(r, V_willoughby_no_smooth, color='orchid',  label='Willoughby - no smoothing') # V_willoughby_no_smooth
    fig5 = plt.plot(r_chavas[:upper_bound], V_chavas[:upper_bound], color='forestgreen',   label='Chavas profile') # V_chavas    
    
    plt.xlabel('Radius (km)')
    plt.ylabel('Wind speed (m/s)')
    plt.legend();plt.grid()

    # Small scale
    plt.subplot(2, 2, 4)
    if Rmax >= 25:
        fig6 = plt.plot(r[Rmax - 25 : Rmax + 25], spdm[Rmax - 25 : Rmax + 25], color='k', linewidth=3, label='SAR-derived wind speed') # V_obs
        fig7 = plt.plot(r[Rmax - 25 : Rmax + 25], V_rankine[Rmax - 25 : Rmax + 25], color='darkorange', label='Rankine profile')       # V_rankine
        fig8 = plt.plot(r[Rmax - 25 : Rmax + 25], V_holland[Rmax - 25 : Rmax + 25], color='steelblue',  label='Holland profile')       # V_holland
        fig9 = plt.plot(r[Rmax - 25 : Rmax + 25], V_willoughby_no_smooth[Rmax - 25 : Rmax + 25], color='orchid',  label='Willoughby - no smoothing') # V_willoughby_no_smooth
        fig10= plt.plot(r_chavas[Rmax_chavas - index25 : Rmax_chavas + index25], V_chavas[Rmax_chavas - index25 : Rmax_chavas + index25], color='forestgreen',   label='Chavas profile') # V_chavas
    else:
        fig6 = plt.plot(r[:50], spdm[:50], color='k', linewidth=3, label='SAR-derived wind speed') # V_obs
        fig7 = plt.plot(r[:50], V_rankine[:50], color='darkorange', label='Rankine profile')       # V_rankine
        fig8 = plt.plot(r[:50], V_holland[:50], color='steelblue',  label='Holland profile')       # V_holland
        fig9 = plt.plot(r[:50], V_willoughby_no_smooth[:50], color='orchid',  label='Willoughby - no smoothing') # V_willoughby_no_smoothing
        fig10= plt.plot(r_chavas[:upper_bound50], V_chavas[:upper_bound50], color='forestgreen',   label='Chavas profile') # V_chavas

    plt.xlabel('Radius (km)')
    plt.ylabel('Wind speed (m/s)')
    plt.legend();plt.grid()
    plt.savefig(savepath)
    
    # SAVE INITIAL AND FIT VALUES
    text_file = open(savepath + ".txt", "w")
    text_file.write("========================== FITTING COMMON WIND PROFILES ON SAR DATA ==========================\n" + "N°" + i + ": " + filename + "\n\n\n")
    text_file.write("RANKINE\n") # x, alpha, Vmin, Rmax
    text_file.write("x_ini     = {:.2f}".format(INI['Rankine'][0]) + '   | x_fit     = {:.2f}\n'.format(FIT['Rankine'][0]))
    text_file.write("alpha_ini = {:.2f}".format(INI['Rankine'][1]) + '   | alpha_fit = {:.2f}\n'.format(FIT['Rankine'][1]))
    text_file.write("Vmin_ini  = {:.2f}".format(INI['Rankine'][2]) + '  | Vmin_fit  = {:.2f}\n'.format(FIT['Rankine'][2]))
    text_file.write("Rmax_ini  = {:.2f}".format(INI['Rankine'][3]) + '  | Rmax_fit  = {:.2f}\n\n'.format(FIT['Rankine'][3]))
    text_file.write("HOLLAND\n") # Lat, pn, pc, Vmin, Rmax, Vmax
    text_file.write("Lat      = {:.2f}\n".format(INI['Holland'][0]))
    text_file.write("pn_ini   = {:.2f}".format(INI['Holland'][1] / 100) + ' | pn_fit   = {:.2f}\n'.format(FIT['Holland'][1] / 100))
    text_file.write("pc_ini   = {:.2f}".format(INI['Holland'][2]/ 100) + '  | pc_fit   = {:.2f}\n'.format(FIT['Holland'][2] / 100))
    text_file.write("Vmin_ini = {:.2f}".format(INI['Holland'][3]) + '   | Vmin_fit = {:.2f}\n'.format(FIT['Holland'][3]))
    text_file.write("Rmax_ini = {:.2f}".format(INI['Holland'][4]) + '   | Rmax_fit = {:.2f}\n'.format(FIT['Holland'][4]))
    text_file.write("Vmax_ini = {:.2f}".format(INI['Holland'][5]) + '   | Vmax_fit = {:.2f}\n\n'.format(FIT['Holland'][5]))
    text_file.write("WILLOUGHBY\n") # n, X1, Vmin, Rmax, Vmax
    text_file.write("n_ini    = {:.2f}".format(INI['Willoughby'][0]) + '    | n_fit    = {:.2f}\n'.format(FIT['Willoughby'][0]))
    text_file.write("X1_ini   = {:.2f}".format(INI['Willoughby'][1]) + '   | X1_fit   = {:.2f}\n'.format(FIT['Willoughby'][1]))
    text_file.write("Vmin_ini = {:.2f}".format(INI['Willoughby'][2]) + '   | Vmin_fit = {:.2f}\n'.format(FIT['Willoughby'][2]))
    text_file.write("Rmax_ini = {:.2f}".format(INI['Willoughby'][3]) + '   | Rmax_fit = {:.2f}\n'.format(FIT['Willoughby'][3]))
    text_file.write("Vmax_ini = {:.2f}".format(INI['Willoughby'][4]) + '   | Vmax_fit = {:.2f}\n\n'.format(FIT['Willoughby'][4]))
    text_file.write("CHAVAS\n") # x, alpha, Vmin, Rmax
    text_file.write("Vmax_ini = {:.2f}\n".format(INI['Chavas'][0]))
    text_file.write("Rfit_ini = {:.2f}\n".format(INI['Chavas'][1] / 1000))
    text_file.write("Vfit_ini = {:.2f}\n".format(INI['Chavas'][2]))
    text_file.write('                   | Rmax_fit   = {:.2f}\n'.format(FIT['Chavas'][2] / 1000))
    text_file.write('                   | R0_fit     = {:.2f}\n'.format(FIT['Chavas'][3] / 1000))
    text_file.write('                   | Rmerge_fit = {:.2f}\n'.format(FIT['Chavas'][4] / 1000))
    text_file.write('                   | Vmerge_fit = {:.2f}\n\n'.format(FIT['Chavas'][5]))
    text_file.write('CAVEAT 1: The Willoughby profile is fitted without the polynomial ramp smoothing method.\n')
    text_file.write('CAVEAT 2: The Chavas profile is fitted on the total wind speed and not the tangential wind speed.\n')
    
    return True

#================================= B SENSITIVITY FUNCTIONS =====================================

def initialize_B_sensitivity_experiment(spdm, power_law, rho, Lat, pn, pc, print_values):
    Lat                  = np.float64(Lat)
    R1, R2, V1, V2, Vmin = initialize_radii(spdm, power_law)
    Lat, pn, pc, A, B    = initialize_A_and_B(rho, Lat, pn, pc, Rmax=R1, Vmax=V1)
    if print_values:
        print(
            'B SENSITIVITY - Initial values',
            '\n rho_ini  =', "{:.2f}".format(rho),
            '\n Lat_ini  =', "{:.2f}".format(Lat),
            '\n R1_ini   =', "{:.2f}".format(R1),
            '\n R2_ini   =', "{:.2f}".format(R2),
            '\n V1_ini   =', "{:.2f}".format(V1),
            '\n V2_ini   =', "{:.2f}".format(V2),
            '\n pn_ini   =', "{:.2f}".format(pn),
            '\n pc_ini   =', "{:.2f}".format(pc),
            '\n Vmin_ini =', "{:.2f}".format(Vmin),
            '\n A_ini    =', "{:.2f}".format(A),
            '\n B_ini    =', "{:.2f}".format(B)
        )
    return rho, Lat, R1, R2, V1, V2, pn, pc, Vmin, A, B

def initialize_radii(spdm, power_law):
    '''Given the spdm of a SAR image, returns the radii R1 (= Rmax) and R2, as well as V1 (=Vmax) and V2, and Vmin.
    V2 is computed as a cubic fraction of Vmax. Then R2 is derived accordingly.'''
    p    = power_law # test 2 and 3
    R1   = np.argmax(spdm[:200])
    V1   = spdm[R1]
    Vmin = np.nanmin(spdm)
    V2   = (1 - (V1 ** p / ((V1 + Vmin) ** p))) * V1
    diff = np.abs(spdm - V2)
    R2   = (np.argmin(diff[R1:]) + R1) # We look for R2 to the right of Rmax (R1)
    return R1, R2, V1, V2, Vmin

def initialize_A_and_B(rho, Lat, pn, pc, Rmax, Vmax):
    B     = (Vmax ** 2) * rho * np.e / (pn - pc)
    A     = Rmax ** B
    return Lat, pn, pc, A, B 

def B_sensitivity_holland_profile(r, rho, Lat, pn, pc, Vmin, A, B):
    '''Returns the Holland profile (that will be used separately on the 3 pieces), as a function of A and B (rather than Vmax and Rmax)'''
    fcor  = coriolis(Lat)
    V   = Vmin + np.sqrt(A * B * (pn - pc) * np.exp((-1) * A / (r ** B)) / (rho * r ** B) + (r ** 2 * fcor ** 2) / 4) - (r * fcor / 2)
    return V
    
def B_sensitivity_complete_profile(r, rho, Lat, R1, R2, pn, pc, Vmin, A, B0, B1, B2):
    '''Given some arguments, returns the complete profile, composed of 3 Holland profiles '''
    V              = r * 0. # initialize
    Vi             = B_sensitivity_holland_profile(r, rho, Lat, pn, pc, Vmin, A, B0)
    Vtrans         = B_sensitivity_holland_profile(r, rho, Lat, pn, pc, Vmin, A, B1)
    Vo             = B_sensitivity_holland_profile(r, rho, Lat, pn, pc, Vmin, A, B2)
    r_under_R1     = (r <= R1)
    r_in_between   = (r > R1) & (r < R2)
    r_over_R2      = (r >= R2)
    V[r_under_R1]  = Vi[r_under_R1]
    V[r_in_between]= Vtrans[r_in_between]
    V[r_over_R2]   = Vo[r_over_R2]
    return V

def fit_B_sensitivity_experiment(r, spdm, rho, Lat, R1, R2, V1, V2, pn, pc, Vmin, A, B, print_values):
    '''Note: V1 and V2 are useless but they are contained in INI so they are passed as arguments here.'''
    R1 = min(2 * R1, 200)
    popt, pcov = curve_fit(lambda r, R1, R2, pn, pc, Vmin, A, B0, B1, B2: B_sensitivity_complete_profile(r, rho, Lat, R1, R2, pn, pc, Vmin, A, B0, B1, B2), r, spdm, p0=[R1, R2, pn, pc, Vmin, A, B, B, B], bounds=((0, 0, 1000 * 100, 850 * 100, 0, 0, 0, 0, 0), (200, 500, 1100 * 100, 1000 * 100, 50, 10000, 3, 3, 3))) # Lat, rho, R1, R2 are fixed
    R1, R2, pn, pc, Vmin, A, B0, B1, B2 = popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6], popt[7], popt[8]
    if print_values:
        print(
            'B SENSITIVITY - Fit values',
            '\n R1_fit   =', "{:.2f}".format(R1),
            '\n R2_fit   =', "{:.2f}".format(R2),
            '\n pn_fit   =', "{:.2f}".format(pn),
            '\n pc_fit   =', "{:.2f}".format(pc),
            '\n Vmin_fit =', "{:.2f}".format(Vmin),
            '\n A_fit    =', "{:.2f}".format(A),
            '\n B0_fit   =', "{:.2f}".format(B0),
            '\n B1_fit   =', "{:.2f}".format(B1),
            '\n B2_fit   =', "{:.2f}".format(B2)
        )
    return R1, R2, pn, pc, Vmin, A, B0, B1, B2

def plot_B_sensitivity_experiment(i, file, r, spdm, INI, FIT):
    # Compute fitted profile and fitted parameters
    rho, Lat, _, _, _, _, _, _, _, _, _   = INI['B_sens']
    V_fit                                 = B_sensitivity_complete_profile(r, rho, Lat, *FIT['B_sens'])
    V_holland                             = holland_profile(r, *FIT['Holland'])
    Rmax                                  = np.argmax(spdm[:200]) # to center on Rmax
    _, pn, pc, _, Rmax_fit, Vmax_fit      = FIT['Holland'] # Lat, pn, pc, Vmin, Rmax, Vmax
    _, _, _, A_fit, B_fit                 = initialize_A_and_B(rho, Lat, pn, pc, Rmax_fit, Vmax_fit)
    _, _, _, _, _, A_sen, B0_sen, B1_sen, B2_sen= FIT['B_sens']  # pn, pc, Vmin, A, B0, B1, B2
    
    # Compute RMSEs
    holland_rmse = rmse(V_holland, spdm)
    B_sens_rmse  = rmse(V_fit,     spdm)
    
    # Title
    plt.figure(figsize=(18, 8))
    plt.suptitle('N°' + '%d'%i + ": " + os.path.basename(file), fontsize=14)
    
    # Large scale
    plt.subplot(1, 2, 1)
    fig1 = plt.plot(r, spdm, color='k', linewidth=3, label='SAR-derived wind speed') # V_obs
    fig2 = plt.plot(r, V_holland, color='steelblue', label='Holland profile; A_fit = ' + "{:.2f}".format(A_fit) + '; B_fit = ' + "{:.2f}".format(B_fit) + '; => RMSE = ' + "{:.2f}".format(holland_rmse)) # V_holland
    fig3 = plt.plot(r, V_fit, color='fuchsia',       label='Piecewise Holland; A = ' + "{:.2f}".format(A_sen) + '; B0 = ' + "{:.2f}".format(B0_sen) + ';\n B1 = ' + "{:.2f}".format(B1_sen) + '; B2 = ' + "{:.2f}".format(B2_sen) + ';                                  => RMSE = ' + "{:.2f}".format(B_sens_rmse)) # V_fit  
    
    plt.xlabel('Radius (km)')
    plt.ylabel('Wind speed (m/s)')
    plt.legend();plt.grid()

    # Small scale
    plt.subplot(1, 2, 2)
    if Rmax >= 25:
        fig4 = plt.plot(r[Rmax - 25 : Rmax + 25], spdm[Rmax - 25 : Rmax + 25],  color='k', linewidth=3, label='SAR-derived wind speed') # V_obs
        fig8 = plt.plot(r[Rmax - 25 : Rmax + 25], V_holland[Rmax - 25 : Rmax + 25], color='steelblue',  label='Holland profile')        # V_holland
        fig6 = plt.plot(r[Rmax - 25 : Rmax + 25], V_fit[Rmax - 25 : Rmax + 25], color='fuchsia',        label='V_fit')                  # V_fit
    else:
        fig4 = plt.plot(r[:50], spdm[:50],  color='k', linewidth=3, label='SAR-derived wind speed') # V_obs
        fig8 = plt.plot(r[:50], V_holland[:50], color='steelblue',  label='Holland profile')        # V_holland
        fig6 = plt.plot(r[:50], V_fit[:50], color='fuchsia',        label='V_fit')                  # V_fit

    plt.xlabel('Radius (km)')
    plt.ylabel('Wind speed (m/s)')
    plt.legend();plt.grid()
    return True

def rmse(predictions, targets):
    return np.sqrt(np.mean((predictions - targets) ** 2))


#=====================DEBUG======================+-

def print_spd(ds):
    ### DEFINE (r, theta) GRID
    r     = np.arange(501)
    th    = np.arange(361)
    r, th = np.meshgrid(r, th)
    ds_r  = np.array(ds['r_polar'])

    ds_th = np.mod(np.array(ds['theta']) * 180. / np.pi, 360) # convert theta from radians to degrees
    ds_ws = np.array(ds['wind_speed'])
    # Possible to call griddata() without using meshgrid() before? 
    spd   = griddata((ds_r.flatten(), ds_th.flatten()), ds_ws.flatten(), (r, th), method='nearest')

    plt.figure()
    plt.pcolormesh(spd)
    return True

def print_ds(ds):
    plt.figure(figsize=(18, 8))
    plt.subplot(1, 2, 1)
    plt.pcolormesh(ds['x'], ds['y'], ds['wind_speed'])
    plt.subplot(1, 2, 2)
    plt.pcolormesh(ds['x_coords'], ds['y_coords'], ds['wind_speed']);plt.grid()
    return True


