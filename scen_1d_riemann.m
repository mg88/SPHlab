% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.read_data = false;
ps.input_name = 'test';

ps.dx      = 2e-3;
ps.dtfactor   = 0.4;  
ps.tend    = 0.11;    
 
ps.eta     = 1.2;     
ps.eta2    = 2;
ps.kernel  = 'Wendland';
ps.scheme   = 'm';
ps.h_const  = true;

%IO
ps.plot_dt = 1e-3;  
ps.save_as_movie = false;
ps.plotstyle = 'dpv';

 %% material parameter
rho0 = 1;     % density
c0   = 10.0; 

%% domain         
ps.Omega = [-0.2, 1]; 
Vp = ps.dx; %volume per particle

% non-reflecting bc
%ps.omega_nr_bc=[0.25;0.6];

%% active particles

%left
leftpoint = 0.0;
rightpoint= -0.02;
v0   = 1;
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0, false)

% right
leftpoint = leftpoint+ps.dx;
rightpoint= 0.6;
v0   = 0;
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0, false)

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

