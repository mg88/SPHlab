% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d shock tube problem - paper of zisis2014
% close all; clear;
ps = SPHlab_scenario();


%% general parameter
ps.Ntot    = 200;
ps.equalmass = true;
ps.tend    = 0.4;   
% ps.tpause  = 0; 
% ps.dtfactor = 0.01;
ps.eta     = 1.2;     
ps.set_kernel('Gauss');

ps.scheme     = 'm';
ps.h_const    = false;
ps.compOmegaj = true;
ps.normalizeOmega = false;

%IO
ps.plot_dt = 1e-2;  
ps.save_as_movie = false;
ps.plot_quantity = 'dpv';
ps.fixaxes.d = [0.2,1.8];
ps.fixaxes.v = [-0.1,1.2];
ps.fixaxes.p = [-0.1,0.8];
ps.fixaxes.e = [1.5,3];
ps.plotconfig.latexplot = true;
%ps.Neval = ps.Ntot ;
% ps.Nss = 4;

 %% material parameter
ps.EOS = 'ISO';

% ps.art_diss_para.alpha_mass      = 0.3;%0.3;
% ps.art_diss_para.beta_mass       = 0;
% ps.art_diss_para.alpha_viscosity = 1;
% ps.art_diss_para.beta_viscosity  = 2;
% ps.art_diss_para.alpha_energy    = 1;%0.5;
% ps.art_diss_para.beta_energy     = 0;


%% domain         
ps.Omega = [-0.2, 1.6]; 

%% active particles

%% test 1
%1
omega_geo = [0.0,0.6];
v0   = 1;
rho0 = 1;     % density
c0   = 1; %ma=1 (Ma=vimp/c0)
ps.add_geometry(omega_geo, rho0, v0, c0)

%2
omega_geo = [0.6,0.8];
v0   = 0;
rho0 = 1;     % density
c0   = 1;  %ma=1
ps.add_geometry(omega_geo, rho0, v0, c0)

%3
omega_geo = [0.8,1.2];
v0   = 0;
rho0 = 0.25;     % density
c0   = 0.5;  %ma=2
ps.add_geometry(omega_geo, rho0, v0, c0)

%4
omega_geo = [1.2,1.4];
v0   = 0;
rho0 = 1;     % density
c0   = 1;  %ma=1
ps.add_geometry(omega_geo, rho0, v0, c0)


%% BC:
% ps.add_bc('noflow',0,-1); %left
% ps.add_bc('noflow',1.4,1); %right

%% -----------------------------------------------------
% generate particles
ps.create_geometry;
% dispdata(ps);

%% create particle class
obj_particles = SPHlab_particles(ps);
%% start simulation
start_simulation(obj_particles)

% exact_shocktube