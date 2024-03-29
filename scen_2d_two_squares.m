% SPH for HVI  - Markus Ganser - TU/e - 2015
% two squares with surface tension

error('outdated!')

close all; clear; clc;
ps = sph_scenario();

ps.dx      = 1e-2;
ps.tend    = 5;    
ps.eta     = 1.2;     

%IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plot_quantity = 'p';

%% material parameter
rho0 = 1;     %relativ density
c0   = 50;
ps.beta = 0.1;    %for surface tension

%general
ps.Omega = [0, 1;
            0, 1];  
Vp = ps.dx*ps.dx; %volume per particle

lowerleftcorner1 = [0.45,0.35];
upperrightcorner1= [ 0.55,0.6];
lowerleftcorner2 = [0.35,0.35];
upperrightcorner2= [ 0.45,0.45];
v0 = [0,0];
I = add_rectangle2d(ps,[lowerleftcorner1;lowerleftcorner2],...
                 [upperrightcorner1;upperrightcorner2]);
addproperties(ps, I, Vp, rho0, v0,c0)

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

