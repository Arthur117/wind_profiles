%CLE15_plot_rmaxinput.m -- Plot profile from Chavas et al. (2015), rmax in
%
% Outputs:
%   rr [m] - vector of radii
%   VV [ms-1] - vector of wind speeds at rr
%   r0 [m] - outer radius
%   rmerge [m] - radius of merge point
%   Vmerge [ms-1] - wind speed of merge point

% Author: Dan Chavas
% CEE Dept, Princeton University
% email: drchavas@gmail.com
% Website: -
% 2015-05-12; Last revision:

%------------- BEGIN CODE --------------


addpath(genpath('mfiles/'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTES FOR USER:
% Parameter units listed in []
% Characteristic values listed in {}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Storm parameters
% Vmax = 35;                      %[ms-1] {50}; maximum azimuthal-mean wind speed
% rmax = 20*1000;                 %[m] {40*1000}; radius of Vmax
% fcor = 5e-5;                    %[s-1] {5e-5}; Coriolis parameter at storm center
omega=4*pi/(24*3600);
fcor=2*omega.*sin(lat*pi/180);
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


%% Get profile: rmax input
% [rr,VV,r0,rmerge,Vmerge,rrfracr0,MMfracM0,rmaxr0,MmM0,rmerger0,MmergeM0] = ...
%     ER11E04_nondim_rmaxinput(Vmax,rmax,fcor,Cdvary,Cd,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye);
[rr,VV,r0,rmerge,Vmerge] = ...
    ER11E04_nondim_rmaxinput(Vmax,rmax,fcor,Cdvary,Cd,w_cool,CkCdvary,CkCd,eye_adj,alpha_eye);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOTTING
%% Default options -- as desired
set(0,'defaultaxesfontsize',18,'defaultaxesfontweight','normal',...
                'defaultlinelinewidth',4,'DefaultAxesFontName','Helvetica')

%% Initializaiton
hh = figure(1);
clf(hh)

%%Position/size
set(hh,'Units','centimeters');
hpos = [0 0 30 30];
set(hh,'Position',hpos);
set(hh,'PaperUnits','centimeters');
set(hh,'PaperPosition',hpos);
set(hh,'PaperSize',hpos(3:4));
set(gca,'position',[0.12    0.1    0.84    0.81]);

figure(1)
clear hpl input_legend
hold off
hpl(1) = plot(rr/1000,VV,'b','LineWidth',3);
input_legend{1} = 'Model';
hold on
hpl(2) = plot(rmax/1000,Vmax,'bx','MarkerSize',20);
input_legend{2} = 'r_{max}';
hpl(3) = plot(r0/1000,0,'r.','MarkerSize',40);
input_legend{3} = 'r_0';
%plot(rmerge/1000,Vmerge,'bx','MarkerSize',15);
%input_legend{2} = 'r0';

axis([0 1.05*r0/1000 0 1.1*Vmax])
legend(hpl,input_legend,'Location','NorthEast'); legend boxoff

%%Save plot
plot_filename = sprintf('CLE15_plot_rmaxinput.pdf');
saveas(gcf,plot_filename,'pdf')

