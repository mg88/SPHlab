% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d boundaries on both sides

close all; clear; clc;
ps = sph_scenario();

ps.dx      = 1e-3;
ps.tend    = 1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-3;   
ps.save_as_movie = false;
ps.plotstyle = 'vd';

 %% material parameter
rho0 = 1;     %relativ density
c0   = 10; 


%% domain         
ps.Omega = [0,1.1]; 
Vp = ps.dx; %volume per particle

%% active particles
leftpoint = 0.08;
rightpoint= 0.88;
v0   = -1;
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0, false)

%% boundary:
bounleft1  = 0;
bounright1 = 0.05;
bounleft2  = 0.95;
bounright2 = 1;
v0=0;
I = add_line1d(ps,[bounleft1;bounleft2]...
                 ,[bounright1;bounright2]);
addproperties(ps, I, Vp, rho0, v0,c0, true)

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

