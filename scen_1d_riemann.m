% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.read_data = false;
ps.input_name = 'test';

ps.dx      = 2e-2;
ps.dtfactor   = 0.4;  
ps.tend    = 0.1;    
 
ps.eta     = 1.2;     
ps.eta2    = 2;
ps.kernel  = 'M4';
ps.scheme   = 'v';

%IO
ps.plot_dt = 2e-3;  
ps.save_as_movie = false;
ps.plotstyle = 'vp';

 %% material parameter
rho0 = 10;     % density
c0   = 10; 

%% domain         
ps.Omega = [-0.2, 0.6]; 
Vp = ps.dx; %volume per particle

%% active particles

%left
leftpoint = 0.22;
rightpoint= -0.2;
v0   = 2;
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0, false)

% right
leftpoint = leftpoint+ps.dx;
rightpoint= 0.4;
v0   = 0;
I = add_line1d(ps,leftpoint,rightpoint);
addproperties(ps, I, Vp, rho0, v0,c0, false)

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

