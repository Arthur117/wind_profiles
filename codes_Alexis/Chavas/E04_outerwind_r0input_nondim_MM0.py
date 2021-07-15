
import numpy as np

import pdb

def E04_outerwind_r0input_nondim_MM0(r0,fcor,Cdvary,C_d,w_cool,Nr):
    
    ## Initialization
    fcor = abs(fcor)
    M0 = .5*fcor*r0**2 #[m2/s] M at outer radius

    drfracr0 = .001
    if ((r0>2500*1000) | (r0<200*1000)):
        drfracr0 = drfracr0/10 
        #extra precision for very large storm to avoid funny bumps near r0 (though rest of solution is stable!)
        #or for tiny storm that requires E04 extend to very small radii to match with ER11

    if Nr > 1/drfracr0:
        Nr = 1/drfracr0    #grid radii must be > 0


    rfracr0_max = 1   #[-] start at r0, move radially inwards
    rfracr0_min = rfracr0_max - (Nr-1)*drfracr0    #[-] inner-most node
    rrfracr0 = np.arange(rfracr0_min,rfracr0_max+drfracr0,drfracr0) #[] r/r0 vector
    MMfracM0 = np.float('NaN')*np.zeros(rrfracr0.size)  #[] M/M0 vector initialized to 1 (M/M0 = 1 at r/r0=1)
    MMfracM0[-1] = 1

    ## First step inwards from r0: d(M/M0)/d(r/r0) = 0 by definition
    rfracr0_temp = rrfracr0[-2] #one step inwards from r0
    #dMfracM0_drfracr0_temp = 0  #[] d(M/M0)/d(r/r0) = 0 at r/r0 = 1
    MfracM0_temp = MMfracM0[-1]
    MMfracM0[-2] = MfracM0_temp
    ##################################################################
    ## Variable C_d: code from Cd_Donelan04.m (function call is slow) ######
    ##Piecewise linear fit parameters estimated from Donelan2004_fit.m
    C_d_lowV = 6.2e-4
    V_thresh1 = 6  #m/s transition from constant to linear increasing
    V_thresh2 = 35.4  #m/s transition from linear increasing to constant
    C_d_highV = 2.35e-3
    linear_slope = (C_d_highV-C_d_lowV)/(V_thresh2-V_thresh1)
    ##################################################################


    ## Integrate inwards from r0 to obtain profile of M/M0 vs. r/r0
    for ii in range(0,np.int(Nr)-2,1):   #first two nodes already done above
        
        ## Calculate C_d varying with V, if desired
        if Cdvary==1:
            
            ##Calculate V at this r/r0 (for variable C_d only)
            V_temp = (M0/r0)*((MfracM0_temp/rfracr0_temp)-rfracr0_temp)
            
            ##Calculate C_d
            if V_temp<=V_thresh1:
                C_d = C_d_lowV
            elif V_temp>V_thresh2:
                C_d = C_d_highV
            else:
                C_d = C_d_lowV + linear_slope*(V_temp-V_thresh1)

        ## Calculate model parameter, gamma
        gam = C_d*fcor*r0/w_cool   #[] non-dimensional model parameter

        ## Update dMfracM0_drfracr0 at next step inwards
        dMfracM0_drfracr0_temp = gam*((MfracM0_temp-rfracr0_temp**2)**2)/(1-rfracr0_temp**2)

        ## Integrate M/M0 radially inwards
        MfracM0_temp = MfracM0_temp - dMfracM0_drfracr0_temp*drfracr0

        ## Update r/r0 to follow M/M0
        rfracr0_temp = rfracr0_temp - drfracr0 #[] move one step inwards

        ## Save updated values
        MMfracM0[MMfracM0.shape[0]-1-ii-2] = MfracM0_temp

    #pdb.set_trace()
    return rrfracr0,MMfracM0


    #pdb.set_trace()
