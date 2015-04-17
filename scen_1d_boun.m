% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d boundaries on both sides

close all; clear; clc;
ps = sph_scenario();

ps.dx      = 1e-3;
ps.tend    = 1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 6e-2;   
ps.save_as_movie = false;
ps.plotstyle = 'vd';

 %% material parameter
rho0 = 1;     %relativ density
c0   = 30; 


%% domain         
ps.Omega = [0,1.1]; 
Vp = ps.dx; %volume per particle

%% active particles
leftpoint = 0.18;
rightpoint= 0.88;
v0   = 1;
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0)

%% BC:
p1 = 0.05;
n = -1;
ps.add_bc_noflow(p1,0,n);

p1  = 0.95;
n=1;
ps.add_bc_noflow(p1,0,n);

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

