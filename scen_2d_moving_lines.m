% SPH for HVI  - Markus Ganser - TU/e - 2015
%two moving lines colliding 

close all; clear; clc;
ps = sph_scenario();

ps.dx      = 1e-2;
ps.tend    = 0.5;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';


%% material parameter
rho0 = 1;     %relativ density
c0   = 50;

%% domain
ps.Omega = [0, 1;  %x
            0, 1]; %y
Vp = ps.dx*ps.dx; %volume per particle

%line 1
startpoint = [0.22,0.49];
endpoint   = [0.55,0.49];
v0   = [0,-1];
layer = 4;
I = add_line2d(ps,startpoint,endpoint,layer);   
addproperties(ps, I, Vp, rho0, v0,c0, false)

%line 2
startpoint = [0.42,0.41];
endpoint   = [0.8,0.41];
v0   = [0,0];
I = add_line2d(ps,startpoint,endpoint,layer);
addproperties(ps, I, Vp, rho0, v0,c0, false)

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

