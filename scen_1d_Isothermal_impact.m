% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d isothermal impact - paper of zisis2014

close all; clear; clc;
ps = sph_scenario();


%% general parameter
ps.Ntot    = 500;
ps.equalmass = true;
ps.tend    = 0.4;    
 
ps.eta     = 1.2;     
ps.eta2    = 2;
ps.kernel  = 'Wendland';
ps.scheme   = 'v';
ps.h_const  = false;

%IO
ps.plot_dt = 2e-2;  
ps.save_as_movie = false;
ps.plot_quantity = 'dpv';

 %% material parameter


%% domain         
ps.Omega = [-0.05, 1.8]; 

%% active particles
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


%% -----------------------------------------------------
% generate particles
ps.create_geometry;
%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%obj_particles.showInitial()
%% start simulation
start_simulation(obj_particles)
