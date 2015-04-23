% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

%read input?
ps.read_data  = false;
ps.input_name = 'test';
%% general parameter
ps.Ntot      = 100;
ps.equalmass = false;
%ps.dt        = 1e-3;
ps.tend      = 1.75; 

ps.eta     = 1.2;     
ps.eta2    = 2;
ps.kernel  = 'Wendland';
ps.scheme  = 'm';
ps.h_const = false;

%IO
ps.plot_dt = 5e-2;  
ps.save_as_movie = false;
ps.movie_name = 'out2';
ps.plot_quantity = 'pv';
ps.fixaxes.v = [-0.003, 0.005];
ps.fixaxes.p = [-0.005 , 0.005];
%ps.fixaxes.d = [1-1e-2 ,1+1e-2 ];


 %% material parameter
rho0 = 1;     % density
c0   = 1.0; 

%% domain         
ps.Omega = [-0.2, 1.2]; 

%% active particles

%left
omega_geo = [0,0.1]; %x
v0   = 0.004;
ps.add_geometry(omega_geo, rho0, v0, c0)

% middle
omega_geo = [omega_geo(2),0.7]; %x
v0   = 0;
ps.add_geometry(omega_geo, rho0, v0, c0)

% right %bigger particles:
% omega_geo = [0.6,0.7];
% v0   = 0.00;
% Nfactor = 1;
% ps.add_geometry(omega_geo, rho0, v0, c0,Nfactor)

%% set BC

p1 = 0.7;
outer_normal = 1;
ps.add_bc_nr(p1,0,outer_normal);
%ps.bc(1).damping_area=[0.4;0.5];

%% -----------------------------------------------------
% generate particles
ps.create_geometry;
%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

