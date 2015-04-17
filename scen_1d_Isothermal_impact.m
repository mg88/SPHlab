% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.read_data = false;
ps.input_name = 'test';

ps.dx      = 2e-3;
ps.tend    = 0.4;    
 
ps.eta     = 1.2;     
ps.eta2    = 2;
ps.kernel  = 'Wendland';
ps.scheme   = 'v';
ps.h_const  = false;

%IO
ps.plot_dt = 2e-2;  
ps.save_as_movie = false;
ps.plotstyle = 'dpv';

 %% material parameter


%% domain         
ps.Omega = [-0.05, 1.8]; 
Vp = ps.dx; %volume per particle

%% active particles

%1
leftpoint = 0.0;
rightpoint= 0.6-ps.dx;
v0   = 1;
rho0 = 1;     % density
c0   = 1; %ma=1 (Ma=vimp/c0)
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0)

%2
leftpoint = 0.6;
rightpoint= 0.8-ps.dx;
v0   = 0;
rho0 = 1;     % density
c0   = 1;  %ma=1
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0)

%3
leftpoint = 0.8;
rightpoint= 1.2-ps.dx;
v0   = 0;
rho0 = 0.25;     % density
c0   = 0.5;  %ma=2
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0)


%4
leftpoint = 1.2;
rightpoint= 1.4-ps.dx;
v0   = 0;
rho0 = 1;     % density
c0   = 1;  %ma=1
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0)

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

%Elapsed time is 242.124032 seconds.
