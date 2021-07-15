
import numpy as np
import sys
import pdb


from ER11E04_nondim_rfitinput import ER11E04_nondim_rfitinput

def gochv(Vmax,Rfit,Lat=15.,rrr=None,vfit=17.):





    
    #=======================================
    # Default Parameters
    #=======================================

    #fcor = 5e-5                    #[s-1] {5e-5} Coriolis parameter at storm center
    omega=2*np.pi/(24*3600)
    fcor=2*omega*np.sin(Lat*np.pi/180)
    # vfm1=5.7
    
    ## Environmental parameters
    ##Outer region
    Cdvary = 1                     #[-] {1} 0 : Outer region Cd = constant (defined on next line) 1 : Outer region Cd = f(V) (empirical 	Donelan et al. 2004)
    Cd = 1.5e-3                #[-] {1.5e-3} ignored if Cdvary = 1 surface momentum exchange (i.e. drag) coefficient
    w_cool = 2/1000                #[ms-1] {2/1000 Chavas et al 2015} radiative-subsidence rate in the rain-free tropics above the boundary layer top

    ##Inner region
    CkCdvary = 1                   #[-] {1} 0 : Inner region Ck/Cd = constant (defined on next line) 1 : Inner region Ck/Cd = f(Vmax) (empirical Chavas et al. 2015)
    CkCd = 1                   #[-] {1} ignored if CkCdvary = 1 ratio of surface exchange coefficients of enthalpy and momentum capped at 1.9 (things get weird >=2)

    ## Eye adjustment
    eye_adj = 0                    #[-] {1} 0 = use ER11 profile in eye 1 = empirical adjustment
    alpha_eye = .15            #[-] {.15 empirical Chavas et al 2015} V/Vm in eye is reduced by factor (r/rm)^alpha_eye ignored if eye_adj=0



    #pdb.set_trace()
    ## Get profile: rfit input

    for cpt, _rfit in enumerate(Rfit):


        rr,VV,rmax,r0,rmerge,Vmerge= ER11E04_nondim_rfitinput(Vmax[cpt],Rfit[cpt],Vfit,fcor,Cdvary,Cd,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye)


        VVV = interp1D(rr,VV,rrr)
        
