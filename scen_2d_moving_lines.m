% SPH for HVI  - Markus Ganser - TU/e - 2015
%two moving lines colliding 

close all; clear; clc;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 1000;
ps.tend    = 0.1;

ps.eta     = 1.2;     
ps.scheme  = 'v';

%IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plot_quantity = '';
ps.fixaxes.p = [-1,1];


%% material parameter
rho0 = 1;     %relativ density
c0   = 20;

%% domain
ps.Omega = [0, 1;  %x
            -0.3, 0.3]; %y

%line 1
omega_geo = [0.2,0.5;
             0.02,0.07];
v0    = [0,-1];
ps.add_geometry(omega_geo, rho0, v0, c0)

%line 2
omega_geo = [0.4,0.8;
             -0.05,0];
v0   = [0,0];
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


