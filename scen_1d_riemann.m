% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem
% set(0,'DefaultFigureWindowStyle','docked') %'normal'
% close all; clear;
ps = sph_scenario();

%% general parameter
ps.Ntot      = 500; %500
ps.equalmass = true;
% ps.dt        = 1e-5;
% ps.dtfactor  = 0.1;
ps.tend      = 2.5; 
% ps.tpause   =[1.3,2];%0.9:0.02:2;%1.0:0.2:3; $and save a figure
ps.eta     = 1.2;     
ps.set_kernel('Gauss');

ps.scheme  = 'm';

ps.h_const   = false;
%experimental settings:
ps.compOmegaj = true;
ps.normalizeOmega = false;
ps.exp_settings.tweakmass = false;


%IO
ps.plot_dt = 5e-2;  
ps.plot_quantity = 'vpd';  %p
%  ps.fixaxes.v = [-0.004, 0.004];
%  ps.fixaxes.p =  [-0.1e-2, 0.6e-2];
%  ps.fixaxes.d = [1-1e-4 ,1+6e-4 ];
% ps.fixaxes.e = [1.49 ,1.51];%1e-5 ];
ps.plotconfig.fixaxes.f = [-1e-3 ,1e-3];%1e-5 ];
ps.plotconfig.latexplot = true;
% ps.plotconfig.drawfluxes = false;
ps.plotconfig.transpose  = false;
% ps.plotconfig.figuresize = [3,3,8,10];
ps.plotconfig.figurename = '';%'cancelation_v';
% ps.Neval = 300;
%output
ps.save_dt = 5e-2;
ps.write_data = false;
ps.output_name ='data/riemann_boun_mass';
ps.save_as_movie = false;
ps.movie_name = 'reflection';
 %% material parameter
 ps.EOS     = 'ISO';

rho0 = 1;     % density
c0   = 1.0; 
e0=0;
MG_Gamma = 2; %(alu)
MG_S     = 1.338;
% ps.art_diss_para.alpha_mass = 0.3;
% ps.art_diss_para.beta_mass = 0.4;
% ps.art_diss_para.alpha_viscosity = 0.2;
% ps.art_diss_para.beta_viscosity = 0.5;
%  ps.art_diss_para.alpha_energy = 0;
%  ps.art_diss_para.beta_energy = 0;
%% domain         
ps.Omega = [-1.5, 1.5]; 

%% active particles
xl = -1.3;
xm = 1.0;
xr = 1.3;
%left
omega_geo = [xl,-xm]; %x
v0   = 0.2;
ps.add_geometry(omega_geo, rho0, v0, c0,e0, MG_Gamma, MG_S)
%% set BC
% ps.add_bc('nrc',-1.2,-1);

%middle
omega_geo = [-xm,xm]; %x
v0   = 0.;
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)
% damping_area = [0.7,1.05];
 ps.add_bc('nrc',0,1);



% right
% omega_geo = [xm,xr]; %x
% v0   = 0.005;
% ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)
% set BC
% ps.add_bc('nrc',xr,1.5);

% right %bigger particles:
% omega_geo = [0.6,0.7];
% v0   = 0.00;
% Nfactor = 1;
% ps.add_geometry(omega_geo, rho0, v0, c0,Nfactor)

%experinamtal: damping area:
%ps.bc(1).damping_area=[0.4;0.5];

%% -----------------------------------------------------
% generate particles
ps.create_geometry;
%% disp data:
%dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);

%obj_particles.IO.plot_data(obj_particles,'dpv');

%% start simulation
start_simulation(obj_particles)