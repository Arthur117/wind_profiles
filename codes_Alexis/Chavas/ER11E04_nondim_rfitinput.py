
import numpy as np

import pdb
from ER11 import ER11_radprof_raw, ER11_radprof
from  E04_outerwind_r0input_nondim_MM0 import  E04_outerwind_r0input_nondim_MM0

from shapely.geometry import LineString
from scipy.interpolate import interp1d

import matplotlib.pyplot as plt
 
def ER11E04_nondim_rfitinput(Vmax,rfit,Vfit,fcor,Cdvary,C_d,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye): 

    ## Initialization
    fcor = np.abs(fcor)
    if CkCdvary==1 :
        CkCd_coefquad = 5.5041e-04
        CkCd_coeflin = -0.0259
        CkCd_coefcnst = 0.7627
        CkCd = CkCd_coefquad*Vmax**2 + CkCd_coeflin*Vmax + CkCd_coefcnst
    CkCd = np.min((1.9,CkCd))
    Mfit = rfit*Vfit + .5*fcor*rfit**2
    soln_converged = 0
    
    while soln_converged==0:
        
        rmaxrfit_min = .01
        rmaxrfit_max = 1.
        rmaxrfit_new = (rmaxrfit_max+rmaxrfit_min)/2    #first guess -- in the middle
        rmaxrfit = rmaxrfit_new    #initialize
        drmaxrfit = rmaxrfit_max - rmaxrfit    #initialize
        drmaxrfit_thresh = .0001
        jterN = 0
        while np.abs(drmaxrfit)>=drmaxrfit_thresh: #keep looping til changes in estimate are very small
            jterN = jterN+1
            
            rmax = rmaxrfit_new*rfit
            ## Step 1: Calculate ER11 M/Mm vs. r/rm
            #[~,~,rrfracrm_ER11,MMfracMm_ER11] = ER11_radprof_nondim(Vmax,rmax,fcor,CkCdvary,CkCd)
            drfracrm = .01
            if(rmax>100*1000):
                drfracrm = drfracrm/10. #extra precision for large storm
                
            rfracrm_min = 0.   #[-] start at r=0
            rfracrm_max = 50.    #[-] extend out to many rmaxs
            rrfracrm_ER11 = np.arange(rfracrm_min,rfracrm_max+drfracrm,drfracrm) #[] r/r0 vector
            rr_ER11 = rrfracrm_ER11*rmax
            rmax_or_r0 = 'rmax'
            #pdb.set_trace()
            VV_ER11,dummy = ER11_radprof(Vmax,rmax,rmax_or_r0,fcor,CkCd,rr_ER11)
            #fig=plt.figure()
            #plt.plot(rr_ER11,VV_ER11);
            #plt.grid()
            #plt.savefig('ER11.png')
            #pdb.set_trace()

            if not np.isnan(np.max(VV_ER11)):    #ER11_radprof converged
                Mm = rmax*Vmax + .5*fcor*rmax**2
                MMfracMm_ER11 = (rr_ER11*VV_ER11 + .5*fcor*rr_ER11**2)/Mm
                rmaxr0_min = .01
                rmaxr0_max = .75
                rmaxr0_new = (rmaxr0_max+rmaxr0_min)/2    #first guess -- in the middle
                rmaxr0 = rmaxr0_new    #initialize
                drmaxr0 = rmaxr0_max - rmaxr0    #initialize
                drmaxr0_thresh = .000001
                
                iterN = 0
                while np.abs(drmaxr0)>=drmaxr0_thresh:
                    iterN = iterN + 1
                    ##Calculate E04 M/M0 vs r/r0
                    r0 = rmax/rmaxr0_new   #[m]
                    Nr = 100000

                    
                    rrfracr0_E04,MMfracM0_E04 = E04_outerwind_r0input_nondim_MM0(r0,fcor,Cdvary,C_d,w_cool,Nr)
                                   
                    
                    ##Convert ER11 to M/M0 vs. r/r0 space
                    rrfracr0_ER11 = rrfracrm_ER11*(rmaxr0_new)
                    M0_E04 = .5*fcor*r0**2
                    MMfracM0_ER11 = MMfracMm_ER11*(Mm/M0_E04)
                    
                    #fig=plt.figure()
                    #plt.plot(rrfracr0_E04,MMfracM0_E04)
                    #plt.plot(rrfracr0_ER11,MMfracM0_ER11)
                    #plt.grid();plt.ylim([0,1]);plt.xlim([0,1])
                    #plt.savefig('moment_'+str(iterN)+'.png')
                    
                    l1 = LineString(zip(rrfracr0_E04,MMfracM0_E04))
                    l2 = LineString(zip(rrfracr0_ER11,MMfracM0_ER11))
                    intersection = l1.intersection(l2)
                    #if intersection.wkt == 'GEOMETRYCOLLECTION EMPTY':   ##no intersections r0 too large --> rmaxr0 too small   
                    if intersection.wkt == 'LINESTRING EMPTY':   ##no intersections r0 too large --> rmaxr0 too small
                        drmaxr0 = np.abs(drmaxr0)/2
                    else:			##at least one intersection -- r0 too small --> rmaxr0 too large
                        if intersection.wkt.split(' ')[0]=='POINT':
                            X0,Y0 = intersection.coords[0]
                        elif intersection.wkt.split(' ')[0]=='MULTIPOINT':
                            #X0,Y0 = intersection[0].coords[0]
                            X0 = [];Y0 =[]
                            for ikt in np.arange(len(intersection)):
                                X0.append(intersection[ikt].coords[0][0])
                                Y0.append(intersection[ikt].coords[0][1])
                            X0 = np.array(X0); Y0 = np.array(Y0)

                        #pdb.set_trace()

                        #at least one intersection -- rmaxr0 too large
                        drmaxr0 = -np.abs(drmaxr0)/2
                        rmerger0 = np.mean(X0)
                        MmergeM0 = np.mean(Y0)
                        
	            ###update value of rmaxr0
                    rmaxr0 = rmaxr0_new    #this is the final one
                    rmaxr0_new = rmaxr0_new + drmaxr0

                    #print(iterN,rmaxr0,rmaxr0_new,drmaxr0)
                    #print(X0,Y0)
                    #print(intersection.wkt, rmerger0, MmergeM0)
                    #pdb.set_trace()
                                    ## Calculate some things
                M0 = .5*fcor*r0**2
                Mm = .5*fcor*rmax**2 + rmax*Vmax
                MmM0 = Mm/M0

                ## Define merged solution
                ii_ER11 = np.argwhere((rrfracr0_ER11<rmerger0) & (MMfracM0_ER11<MmergeM0))[:,0]
                ii_E04 = np.argwhere((rrfracr0_E04>=rmerger0) & (MMfracM0_E04>=MmergeM0))[:,0]
                MMfracM0_temp =np.hstack((MMfracM0_ER11[ii_ER11],MMfracM0_E04[ii_E04]))
                rrfracr0_temp =np.hstack((rrfracr0_ER11[ii_ER11],rrfracr0_E04[ii_E04]))
                del ii_ER11
                del ii_E04

                ## Check to see how close solution is to input value of (rfitr0,MfitM0)
                rfitr0 = rfit/r0
                MfitM0 = Mfit/M0

                ## 2020-06-23 fixed bug returning NaN (and wrong solution) if rfit > current r0
                if rfitr0<=1:
                    f = interp1d(rrfracr0_temp,MMfracM0_temp)
                    #print MMfracM0_temp
                    MfitM0_temp = f(rfitr0)
                    MfitM0_err = MfitM0 - MfitM0_temp
                else:    #true rfit exceeds current r0, so doesnt exist in profile!
                    MfitM0_err = 1000000;  #simply need smaller r0 -- set MfitM0_err to any positive number
                #print(rfitr0, MfitM0, MfitM0_err)
                if MfitM0_err>0:    #need smaller rmax (r0)
                    drmaxrfit = np.abs(drmaxrfit)/2.
                else:    #need larger rmax (r0)
                    drmaxrfit = -np.abs(drmaxrfit)/2.
                #pdb.set_trace()

                
            else:    #ER11_radprof did not converge -- convergence fails for low CkCd and high Ro = Vm/(f*rm)
                ##Must reduce rmax (and thus reduce Ro)
                drmaxrfit = -abs(drmaxrfit)/2

            ##update value of rmaxrfit
            rmaxrfit = rmaxrfit_new    #this is the final one
            rmaxrfit_new = rmaxrfit + drmaxrfit

            #print(rmaxrfit, drmaxrfit, rmaxrfit_new)
            #print(jterN)
            #pdb.set_trace()
            
        #pdb.set_trace()   
        ##Check if solution converged
        if not np.isnan(np.max(VV_ER11)):
            soln_converged = 1
        else:
            soln_converged = 0
            CkCd = CkCd + .1
            #print('Adjusting CkCd to find convergence')


            
    #pdb.set_trace()
    ## Finally: Interpolate to a grid
    # drfracr0 = .0001
    # rfracr0_min = 0    #[-] r=0
    # rfracr0_max = 1   #[-] r=r0
    # rrfracr0 = rfracr0_min:drfracr0:rfracr0_max #[] r/r0 vector
    # MMfracM0 = interp1(rrfracr0_temp,MMfracM0_temp,rrfracr0,'pchip',NaN)
    drfracrm = .01 #calculating VV at radii relative to rmax ensures no smoothing near rmax!
    rfracrm_min = 0.    #[-] r=0
    rfracrm_max = r0/rmax   #[-] r=r0
    rrfracrm = np.arange(rfracrm_min,rfracrm_max,drfracrm) #[] r/r0 vector
    f = interp1d(rrfracr0_temp*(r0/rmax),MMfracM0_temp*(M0/Mm))
    MMfracMm = f(rrfracrm)
    rrfracr0 = rrfracrm*rmax/r0    #save this as output
    MMfracM0 = MMfracMm*Mm/M0;

    ## Calculate dimensional wind speed and radii
    # VV = (M0/r0)*((MMfracM0./rrfracr0)-rrfracr0)  #[ms-1]
    # rr = rrfracr0*r0   #[m]
    # rmerge = rmerger0*r0
    # Vmerge = (M0/r0)*((MmergeM0./rmerger0)-rmerger0)  #[ms-1]
    VV = (Mm/rmax)*(MMfracMm/rrfracrm)-.5*fcor*rmax*rrfracrm  #[ms-1]
    rr = rrfracrm*rmax   #[m]

    ## Make sure V=0 at r=0
    VV[rr==0] = 0

    rmerge = rmerger0*r0
    Vmerge = (M0/r0)*((MmergeM0/rmerger0)-rmerger0)  #[ms-1]
    ## Adjust profile in eye, if desired
    #if(eye_adj==1)
    #    r_eye_outer = rmax
    #    V_eye_outer = Vmax
    #    [VV] = radprof_eyeadj(rr,VV,alpha_eye,r_eye_outer,V_eye_outer)
    #
    #    sprintf('EYE ADJUSTMENT: eye alpha = #3.2f',alpha_eye)
    return rr,VV,rmax,r0,rmerge,Vmerge

    pdb.set_trace()
