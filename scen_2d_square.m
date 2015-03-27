% SPH for HVI  - Markus Ganser - TU/e - 2015
% one square with an a singularity

close all; clear; clc;
ps = sph_scenario();

ps.dim     = 2;
ps.dx      = 1e-2;
ps.tend    = 2e-1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';

%% material parameter
ps.rho0 = 1.0;     %relativ density
ps.c0   = 10;

%% domain
ps.Omega = [1,1]; 

lowerleftcorner1 = [0.3,0.3];
upperrightcorner1= [ 0.7,0.7];
v0    =[0,0];
noise = 0;
I_new = add_rectangle2d(ps.obj_geo,lowerleftcorner1,...
                 upperrightcorner1,v0,ps.dx,ps.dx,noise);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];                            

%set rho0
N        = size(ps.Iin,1);
ps.rho0  = ones(N,1);
ps.rho0(floor(N/2)-100) = 2;  % peak in the center

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

