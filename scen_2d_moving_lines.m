% SPH for HVI  - Markus Ganser - TU/e - 2015
%two moving lines colliding 

close all; clear; clc;
ps = sph_scenario();

ps.dim     = 2;
ps.dx      = 1e-2;
ps.tend    = 0.5;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';


%% material parameter
ps.rho0 = 1;     %relativ density
ps.c0   = 100;

%% domain
ps.Omega = [1,1];  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]

%line 1
startpoint = [0.22,0.49];
endpoint   = [0.55,0.49];
v0   = [0,-0.5];
noise = 0;
layer = 2;
I_new = add_line2d(ps.obj_geo,startpoint,endpoint,...
    v0,ps.dx,layer,noise);   
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];

%line 2
startpoint = [0.42,0.41];
endpoint   = [0.8,0.41];
v0   = [0,0];
noise = 0;
layer = 2;
I_new = add_line2d(ps.obj_geo,startpoint,endpoint,...
    v0 ,ps.dx,layer,noise);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

