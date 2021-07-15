#/usr/bin/env python
import numpy as np
import sys
import copy
from shapely.geometry import LineString
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math
import pdb

#######################################################################
# NOTES FOR USER:
# Parameter units listed in []
# Characteristic values listed in {}
#######################################################################


def glochv(Rfit,Vfit,Az,Vmax,Vfm,Dfm,Lat,Rlim):
	# I) Vmax distribution.

	# Vector decomposition: TC rotational wind component (without translation)
	Vmax=Vmax-0.5*Vfm # withdraw the transport component (Chavas2015)
	ux=Vmax*np.sin(Az*np.pi/180) # zonal wind component
	uy=Vmax*np.cos(Az*np.pi/180) # meridional wind component

	# Translation wind influence: Geostrophic component
	ux_g=0.5*Vfm*np.sin(Dfm*np.pi/180)
	uy_g=0.5*Vfm*np.cos(Dfm*np.pi/180)

	# Total wind component
	u=ux+ux_g
	v=uy+uy_g 
	vmax_azimuth=np.sqrt(u**2+v**2) # azimuthal vmax 
	#--------------------------
	# II) Azimuthal Individual Chavas computation

	Z=Az.size # Number of points.
	# Wind field Parameters
	WW_cell=[] # Store the different wind values of individual transects. 
	RR_cell=[] # Store the radius vector of the different individual transects.

	# Parameters for scaling issues.
	D=np.zeros((Z)) # Store the max extent of every transect.
	deltax=np.zeros((Z)) # Store the radial interval of every transect.
	N=np.zeros((Z)) # length of every transect.

	#pdb.set_trace()
	# Computation of Chavas transect, over all the azimuthal sampling.
	for k in range(0,Z):
		print(k)
		radi=Rfit[k] # select azimuthal value of rfit.
		Vm=vmax_azimuth[k]# select the following vmax.
		pdb.set_trace()
		rr,VV=chv(radi*1000,Vfit,Vm,Lat) # Invoke Chavas power
		pdb.set_trace()
		rr=rr/1000 # set in km (more convenient units)
		VV=VV[rr<Rlim] 
		rr=rr[rr<Rlim]
		WW_cell.append(VV)
		RR_cell.append(rr)
		N[k]=rr.size # number of points for each transects
		D[k]=rr.max() # take the maximum distance for each transects
		deltax[k]=D[k]/N[k] # interval dx for each transects

	pdb.set_trace()
	###### SCALING PARAMETER
	############## wanted parameters
	Df=Rlim  #max(D) the max distance among all transects
	Delta=1  #min(deltax) the best resolution
	#############

	xq=np.arange(Df) # the new spatial scaling to apply uniformaly for every transects
	xq=xq[0:Df:Delta]

	# III) interpolation of each transect, to obtain homogeneous transect
	# each cells must be same interval
	WWF=np.zeros((xq.size,Z)) # Final Wind matrix
	RRF=np.zeros((xq.size,Z)) # Final radius matrix
	WW_cell=np.array(WW_cell)
	RR_cell=np.array(RR_cell)

	for k in range(0,Z):
		WWF[:,k]=np.interp(xq,np.array(RR_cell[k]),np.array(WW_cell[k]))

	# IV) Cartesian Grid

	dim=xq.size # radial extent
	dir1=np.ones((dim,Z),float) # azimuth matrix

	for j in range(0,dim):
		dir1[j,:]=Az
	
	for j in range(0,Z):
		RRF[:,j]=xq
	

	# New Cartesian grid 
	xf=RRF*np.cos(dir1*np.pi/180)    
	yf=RRF*np.sin(dir1*np.pi/180)

	return RRF, WWF, xf,yf

def chv(Rfit,Vfit,Vmax,Lat):

	## Storm parameters
	Vmax = 39                      #[ms-1] {50} maximum azimuthal-mean wind speed
	rfit = 93*1000                #[m] {300*1000} a wind radius
	Vfit = 17                      #[m] {12} wind speed at rfit
	lat=12.3
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
	    
	## Get profile: rfit input

	rr,VV,rmax,r0,rmerge,Vmerge= ER11E04_nondim_rfitinput(Vmax,Rfit,Vfit,fcor,Cdvary,Cd,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye)

	return rr, VV


def E04_outerwind_r0input_nondim_MM0(r0,fcor,Cdvary,C_d,w_cool,Nr):

    ## Initialization
	fcor = abs(fcor)
	M0 = .5*fcor*r0**2 #[m2/s] M at outer radius

	drfracr0 = .001
	if ((r0>2500*1000)|(r0<200*1000)):
		drfracr0 = drfracr0/10 #extra precision for very large storm to avoid funny bumps near r0 (though rest of solution is stable!) #or for tiny storm that requires E04 extend to very small radii to match with ER11
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

	return rrfracr0,MMfracM0 

def ER11_radprof_raw(Vmax,r_in,rmax_or_r0,fcor,CkCd,rr_ER11):
    fcor = np.abs(fcor)
    if rmax_or_r0 == 'rmax':
       rmax = r_in
    else:
      print ('rmax_or_r0 must be set to"rmax"')
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
      print ('rmax_or_r0 must be set to"rmax"')

    #print('figure de test')
    #from matplotlib import pyplot as plt
    #fig=plt.figure()
    #plt.plot(rr_ER11/1000.,V_ER11)
    #plt.grid();plt.xlim([0,500])
    #plt.savefig('test.png')
    #pdb.set_trace()
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
			drin_temp = r_in_save-f(0)

	return V_ER11,r_out

def ER11E04_nondim_r0input(Vmax,r0,fcor,Cdvary,C_d,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye):


    ## Initialization
	fcor = np.abs(fcor)
	if CkCdvary==1 :
		CkCd_coefquad = 5.5041e-04
		CkCd_coeflin = -0.0259
		CkCd_coefcnst = 0.7627
		CkCd = CkCd_coefquad*Vmax**2 + CkCd_coeflin*Vmax + CkCd_coefcnst
    ## 'Ck/Cd is capped at 1.9 and has been set to this value. If CkCdvary=1, 
    ## then Vmax is much greater than the range of data used to estimate 
    ## CkCd as a function of Vmax -- here be dragons!')
	CkCd = np.min((1.9,CkCd))

	## Step 1: Calculate E04 M/M0 vs. r/r0
	Nr = 100000
	rrfracr0_E04,MMfracM0_E04 = E04_outerwind_r0input_nondim_MM0(r0,fcor,Cdvary,C_d,w_cool,Nr)

	M0_E04 = .5*fcor*r0**2

	##Step 2: Converge rmaxr0 geometrically until ER11 M/M0 has tangent point with E04 M/M0

	count = 0
	soln_converged = 0
	while soln_converged==0:
		count += 1
		print (count)
		##Break up interval into 3 points, take 2 between which intersection vanishes, repeat til converges
		rmaxr0_min = .001
		rmaxr0_max = .75
		rmaxr0_new = (rmaxr0_max+rmaxr0_min)/2.    #first guess -- in the middle
		rmaxr0 = rmaxr0_new    #initialize
		drmaxr0 = rmaxr0_max - rmaxr0    #initialize
		drmaxr0_thresh = .000001
		iterN = 0
		rfracrm_min = 0.   #[-] start at r=0
		rfracrm_max = 50.    #[-] extend out to many rmaxs
		while np.abs(drmaxr0)>=drmaxr0_thresh:
			iterN = iterN + 1
			##Calculate ER11 M/Mm vs r/rm
			rmax = rmaxr0_new*r0   #[m]
			print (rmax)
			#     [~,~,rrfracrm_ER11,MMfracMm_ER11] =
			#     ER11_radprof_nondim(Vmax,rmax,fcor,CkCd) #FAILS FOR LOW CK/CD NOT
			#     SURE WHY
			drfracrm = .01
			if rmax>100.*1000 :
				drfracrm = drfracrm/10. #extra precision for large storm

				rrfracrm_ER11 = np.arange(rfracrm_min,rfracrm_max+drfracrm,drfracrm) #[] r/r0 vector
				rr_ER11 = rrfracrm_ER11*rmax
				rmax_or_r0 = 'rmax'
				VV_ER11,dummy = ER11_radprof(Vmax,rmax,rmax_or_r0,fcor,CkCd,rr_ER11)

			if not np.isnan(np.max(VV_ER11)):    #ER11_radprof converged
				rrfracr0_ER11 = rr_ER11/r0
				MMfracM0_ER11 = (rr_ER11*VV_ER11 + .5*fcor*rr_ER11**2)/M0_E04
				l1 = LineString(zip(rrfracr0_E04,MMfracM0_E04))
				l2 = LineString(zip(rrfracr0_ER11,MMfracM0_ER11))
				intersection = l1.intersection(l2)
				if intersection.wkt == 'GEOMETRYCOLLECTION EMPTY':    #no intersections -- rmaxr0 too small
					drmaxr0 = np.abs(drmaxr0)/2
				else:
					if intersection.wkt.split(' ')[0]=='POINT':
						X0,Y0 = intersection.coords[0]
					elif intersection.wkt.split(' ')[0]=='MULTIPOINT':
						X0,Y0 = intersection[0].coords[0]
					#at least one intersection -- rmaxr0 too large
						drmaxr0 = -np.abs(drmaxr0)/2
						rmerger0 = np.mean(X0)
						MmergeM0 = np.mean(Y0)
			else:   #ER11_radprof did not converge -- convergence fails for low CkCd and high Ro = Vm/(f*rm)
				    ##Must reduce rmax (and thus reduce Ro)
					drmaxr0 = -abs(drmaxr0)/2
			
			## update value of rmaxr0
			rmaxr0 = rmaxr0_new    #this is the final one
			rmaxr0_new = rmaxr0_new + drmaxr0


		## Check if solution converged
		if ((not np.isnan(np.max(VV_ER11))) and ('rmerger0' in locals())):
			soln_converged = 1
		else:
			soln_converged = 0
			CkCd = CkCd + .1
			print('Adjusting CkCd to find convergence')

	## Calculate some things
	M0 = .5*fcor*r0**2
	Mm = .5*fcor*rmax**2 + rmax*Vmax
	MmM0 = Mm/M0

	## Finally: Interpolate to a grid
	ii_ER11 = np.argwhere((rrfracr0_ER11<rmerger0) & (MMfracM0_ER11<MmergeM0))[:,0]
	ii_E04 = np.argwhere((rrfracr0_E04>=rmerger0) & (MMfracM0_E04>=MmergeM0))[:,0]
	MMfracM0_temp =np.hstack((MMfracM0_ER11[ii_ER11],MMfracM0_E04[ii_E04]))
	rrfracr0_temp =np.hstack((rrfracr0_ER11[ii_ER11],rrfracr0_E04[ii_E04]))
	del ii_ER11
	del ii_E04

	# drfracr0 = .0001
	# rfracr0_min = 0    #[-] r=0
	# rfracr0_max = 1   #[-] r=r0
	# rrfracr0 = rfracr0_min:drfracr0:rfracr0_max #[] r/r0 vector
	# MMfracM0 = interp1(rrfracr0_temp,MMfracM0_temp,rrfracr0,'pchip',NaN)
	drfracrm = .01 #calculating VV at radii relative to rmax ensures no smoothing near rmax!
	rfracrm_min = 0    #[-] r=0
	rfracrm_max = r0/rmax   #[-] r=r0
	rrfracrm = np.arange(rfracrm_min,rfracrm_max,drfracrm) #[] r/r0 vector
	f = interp1d(rrfracr0_temp*(r0/rmax),MMfracM0_temp*(M0/Mm))
	MMfracMm = f(rrfracrm)

	rrfracr0 = rrfracrm*rmax/r0    #save this as output
	MMfracM0 = MMfracMm*Mm/M0

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
	#	r_eye_outer = rmax
	#	V_eye_outer = Vmax
	#	[VV] = radprof_eyeadj(rr,VV,alpha_eye,r_eye_outer,V_eye_outer)

	#    sprintf('EYE ADJUSTMENT: eye alpha = #3.2f',alpha_eye)
	return rr,VV,rmerge,Vmerge,rmax

def ER11E04_nondim_rfitinput(Vmax,Rfit,Vfit,fcor,Cdvary,C_d,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye): 
	## Initialization
	fcor = np.abs(fcor)
	if CkCdvary==1 :
		CkCd_coefquad = 5.5041e-04
		CkCd_coeflin = -0.0259
		CkCd_coefcnst = 0.7627
		CkCd = CkCd_coefquad*Vmax**2 + CkCd_coeflin*Vmax + CkCd_coefcnst

	CkCd = np.min((1.9,CkCd))
	Mfit = Rfit*Vfit + .5*fcor*Rfit**2

	soln_converged = 0
	while soln_converged==0:
                
		rmaxrfit_min = .01
		rmaxrfit_max = 1.
		rmaxrfit_new = (rmaxrfit_max+rmaxrfit_min)/2    #first guess -- in the middle
		rmaxrfit = rmaxrfit_new    #initialize
		drmaxrfit = rmaxrfit_max - rmaxrfit    #initialize
		drmaxrfit_thresh = .0001
                #print(rmaxrfit_new*Rfit)
		while np.abs(drmaxrfit)>=drmaxrfit_thresh: #keep looping til changes in estimate are very small
			rmax = rmaxrfit_new*Rfit
			## Step 1: Calculate ER11 M/Mm vs. r/rm
			#[~,~,rrfracrm_ER11,MMfracMm_ER11] = ER11_radprof_nondim(Vmax,rmax,fcor,CkCdvary,CkCd)
			drfracrm = .01; #print(Rfit)
			if(rmax>100*1000):
				drfracrm = drfracrm/10. #extra precision for large storm

			rfracrm_min = 0.   #[-] start at r=0
			rfracrm_max = 50.    #[-] extend out to many rmaxs
			rrfracrm_ER11 = np.arange(rfracrm_min,rfracrm_max+drfracrm,drfracrm) #[] r/r0 vector
			rr_ER11 = rrfracrm_ER11*rmax
			rmax_or_r0 = 'rmax'
			VV_ER11,dummy = ER11_radprof(Vmax,rmax,rmax_or_r0,fcor,CkCd,rr_ER11); print(VV_ER11)
                        
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
					rrfracr0_E04,MMfracM0_E04 = E04_outerwind_r0input_nondim_MM0(r0,fcor,Cdvary,C_d,w_cool,Nr); print('toto'); print(np.abs(drmaxr0),drmaxr0_thresh)
					##Convert ER11 to M/M0 vs. r/r0 space
					rrfracr0_ER11 = rrfracrm_ER11*(rmaxr0_new)
					M0_E04 = .5*fcor*r0**2
					MMfracM0_ER11 = MMfracMm_ER11*(Mm/M0_E04)
					l1 = LineString(zip(rrfracr0_E04,MMfracM0_E04))
					l2 = LineString(zip(rrfracr0_ER11,MMfracM0_ER11))
					intersection = l1.intersection(l2)

					if intersection.wkt == 'GEOMETRYCOLLECTION EMPTY':   ##no intersections r0 too large --> rmaxr0 too small 
						drmaxr0 = np.abs(drmaxr0)/2
					else:                   ##at least one intersection -- r0 too small --> rmaxr0 too large
						if intersection.wkt.split(' ')[0]=='POINT':
							X0,Y0 = intersection.coords[0]
							drmaxr0 = -np.abs(drmaxr0)/2
							rmerger0 = np.mean(X0)
							MmergeM0 = np.mean(Y0)
						elif intersection.wkt.split(' ')[0]=='MULTIPOINT':
							X0,Y0 = intersection[0].coords[0]  #at least one intersection -- rmaxr0 too large
							drmaxr0 = -np.abs(drmaxr0)/2
							rmerger0 = np.mean(X0)
							MmergeM0 = np.mean(Y0)
			        ###update value of rmaxr0
					rmaxr0 = rmaxr0_new    #this is the final one
					rmaxr0_new = rmaxr0_new + drmaxr0
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
				rfitr0 = Rfit/r0
				MfitM0 = Mfit/M0
				f = interp1d(rrfracr0_temp,MMfracM0_temp)
					#print MMfracM0_temp
				MfitM0_temp = f(rfitr0)
				MfitM0_err = MfitM0 - MfitM0_temp

				if MfitM0_err>0:    #need smaller rmax (r0)
					drmaxrfit = np.abs(drmaxrfit)/2.
				else:    #need larger rmax (r0)
					drmaxrfit = -np.abs(drmaxrfit)/2.

			else:    #ER11_radprof did not converge -- convergence fails for low CkCd and high Ro = Vm/(f*rm)
		##Must reduce rmax (and thus reduce Ro)
				drmaxrfit = -abs(drmaxrfit)/2

		##update value of rmaxrfit
			rmaxrfit = rmaxrfit_new    #this is the final one
			rmaxrfit_new = rmaxrfit + drmaxrfit
		##Check if solution converged
		if not np.isnan(np.max(VV_ER11)):
		   soln_converged = 1
		else:
		   soln_converged = 0
		   CkCd = CkCd + .1
		   print('Adjusting CkCd to find convergence')

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
	MMfracM0 = MMfracMm*Mm/M0

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




