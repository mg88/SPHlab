% SPH for HVI  - Markus Ganser - TU/e - 2015
% two squares with surface tension

close all; clear; clc;
ps = sph_scenario();

ps.dim     = 2;
ps.dx      = 1e-2;
ps.tend    = 5;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;

%% material parameter
ps.rho0 = 1;     %relativ density
ps.c0   = 50;
ps.beta = 0.1;    %for surface tension

%% domain
ps.Omega = [0, 1;
            0, 1];  

lowerleftcorner1 = [0.45,0.35];
upperrightcorner1= [ 0.55,0.6];
lowerleftcorner2 = [0.35,0.35];
upperrightcorner2= [ 0.45,0.45];
v0=[0,0];
noise = 0;
I_new = add_rectangle2d(ps.obj_geo,[lowerleftcorner1;lowerleftcorner2],...
                 [upperrightcorner1;upperrightcorner2],v0,ps.dx,ps.dx,noise);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

