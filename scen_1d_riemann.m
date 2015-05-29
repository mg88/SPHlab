% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem
set(0,'DefaultFigureWindowStyle','docked') %'normal'
close all; clear;
ps = sph_scenario();

%% general parameter
ps.Ntot      = 200;
ps.equalmass = false;
% ps.dt        = 1e-5;
%  ps.dtfactor  = 0.02;
ps.tend      = 2.8; 
% ps.tpause   = 0.9;%0.9:0.02:2;%1.0:0.2:3;
ps.eta     = 1.2;     
ps.set_kernel('Wendland');

ps.scheme  = 'm';
ps.EOS     = 'ISO';
e0 = 0;
% ps.EOS     = 'IdealGas53'; 
% e0 = 1.5;

ps.h_const = false;

%experimental settings:
% ps.exp_settings.tweakmass = true;


%IO
ps.plot_dt = 5e-2;  
ps.save_as_movie = false;
ps.movie_name = 'out2';
ps.plot_quantity = 'pv';  %p
ps.fixaxes.v = [-0.005, 0.005];

% ps.fixaxes.p = [0.995, 1.005];
ps.fixaxes.p = [-3e-3, 3e-3];
ps.fixaxes.d = [1-2e-3 ,1+2e-3 ];
% ps.fixaxes.e = [1.49 ,1.51];%1e-5 ];
ps.fixaxes.f = [-1e-3 ,1e-3];%1e-5 ];

% ps.Neval = 300;
%output
ps.save_dt = 5e-2;
ps.write_data = false;
ps.output_name ='data/riemann_boun_mass';
 %% material parameter
rho0 = 1;     % density
c0   = 1.0; 

% ps.art_diss_para.alpha_mass = 0.5;
% ps.art_diss_para.alpha_energy = 0.5;
% ps.art_diss_para.beta_mass = 0;
% ps.art_diss_para.beta_viscosity = 0;
%% domain         
ps.Omega = [-2.7, 2.7]; 

%% active particles
xl=-1.5;
xm=0.3;
xr=1.5;
%left
omega_geo = [xl,-xm]; %x
v0   = -0.00;
ps.add_geometry(omega_geo, rho0, v0, c0,e0)
% set BC
ps.add_bc('nrp',xl,0,-1);

%middle
omega_geo = [-xm,xm]; %x
v0   = 0.004;
ps.add_geometry(omega_geo, rho0, v0, c0,e0)

% right
omega_geo = [xm,xr]; %x
v0   = 0.00;
ps.add_geometry(omega_geo, rho0, v0, c0,e0)
%% set BC
ps.add_bc('nrc',xr,0,1);
% ps.add_bc_noflow(xr,0,1);

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