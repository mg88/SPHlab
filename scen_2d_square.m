% SPH for HVI  - Markus Ganser - TU/e - 2015
% one square with an a singularity

close all; clear; clc;
ps = sph_scenario();

ps.dx      = 1e-2;
ps.tend    = 2e-1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';

%% material parameter
rho0 = 1.0;  
c0   = 10;

%% domain
ps.Omega = [0,1;
            0,1]; 
Vp = ps.dx*ps.dx; %volume per particle
 
        
lowerleftcorner1 = [0.3,0.3];
upperrightcorner1= [ 0.7,0.7];
v0    =[0,0];
I = add_rectangle2d(ps,lowerleftcorner1,...
                 upperrightcorner1);
addproperties(ps, I, Vp, rho0, v0,c0, false)
             
%set rho0
N        = size(ps.Iin,1);
ps.rho0j  = ones(N,1);
ps.rho0j(floor(N/2)-100) = 2;  % peak in the center

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

