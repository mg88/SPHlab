% SPH for HVI  - Markus Ganser - TU/e - 2015
% one square with an a singularity

close all; clear; clc;
ps = sph_scenario();

ps.dx      = 2e-2;
ps.tend    = 0.1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-3;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';
ps.fixaxes.p = [-0.5,0.5];

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
addproperties(ps, I, Vp, rho0, v0,c0)
             
%set rho0
N        = size(ps.Iin,1);
ps.rhoj(floor(N/2-10)) = 1.5;  % peak in the center

%% set BC
y  = max(ps.Xj(ps.Iin,1));
kb = abs(ps.Xj(ps.Iin,1) - y) < ps.dx/2;
             
ps.mirrorParticlesj = kb;

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

