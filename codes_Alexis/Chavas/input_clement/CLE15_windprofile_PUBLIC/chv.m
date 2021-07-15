%CLE15_plot_rfitinput.m -- Plot profile from Chavas et al. (2015), (rfit,Vfit) in
% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 2015-05-12; Last revision:

%------------- BEGIN CODE --------------
function [rr,VV]=chv(rfit,Vfit,Vmax,lat)
addpath(genpath('mfiles/'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Storm parameters
% Vmax = 39;                      %[ms-1] {50}; maximum azimuthal-mean wind speed
% rfit = 93*1000;                %[m] {300*1000}; a wind radius
% Vfit = 17;                      %[m] {12}; wind speed at rfit
% lat=12.3;
%fcor = 5e-5;                    %[s-1] {5e-5}; Coriolis parameter at storm center
omega=2*pi/(24*3600);
fcor=2*omega.*sin(lat*pi/180);
% vfm1=5.7;
%% Environmental parameters
%%Outer region
Cdvary = 1;                     %[-] {1}; 0 : Outer region Cd = constant (defined on next line); 1 : Outer region Cd = f(V) (empirical Donelan et al. 2004)
    Cd = 1.5e-3;                %[-] {1.5e-3}; ignored if Cdvary = 1; surface momentum exchange (i.e. drag) coefficient
w_cool = 2/1000;                %[ms-1] {2/1000; Chavas et al 2015}; radiative-subsidence rate in the rain-free tropics above the boundary layer top

%%Inner region
CkCdvary = 1;                   %[-] {1}; 0 : Inner region Ck/Cd = constant (defined on next line); 1 : Inner region Ck/Cd = f(Vmax) (empirical Chavas et al. 2015)
    CkCd = 1;                   %[-] {1}; ignored if CkCdvary = 1; ratio of surface exchange coefficients of enthalpy and momentum; capped at 1.9 (things get weird >=2)

%% Eye adjustment
eye_adj = 0;                    %[-] {1}; 0 = use ER11 profile in eye; 1 = empirical adjustment
    alpha_eye = .15;            %[-] {.15; empirical Chavas et al 2015}; V/Vm in eye is reduced by factor (r/rm)^alpha_eye; ignored if eye_adj=0
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Get profile: rfit input
% [rr,VV,rmax,r0,rmerge,Vmerge,rrfracr0,MMfracM0,rmaxr0,MmM0,rmerger0,MmergeM0] = ...
%    ER11E04_nondim_rfitinput(Vmax,rfit,Vfit,fcor,Cdvary,Cd,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye);
[rr,VV,rmax,r0,rmerge,Vmerge] = ...
    ER11E04_nondim_rfitinput(Vmax,rfit,Vfit,fcor,Cdvary,Cd,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING
%% Default options -- as desired
% set(0,'defaultaxesfontsize',18,'defaultaxesfontweight','normal',...
%                 'defaultlinelinewidth',4,'DefaultAxesFontName','Helvetica')

%% Initializaiton
clear dir1  rr1 dd1 U1 U2 VVV VV1 VVV10 RR x1 y2 beta1 WW 



end
