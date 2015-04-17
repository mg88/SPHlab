% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.read_data  = false;
ps.input_name = 'test';

ps.dx       = 4e-3;
ps.tend     = 1.4;    
 
ps.eta     = 1.2;     
ps.eta2    = 2;
ps.kernel  = 'Wendland';
ps.scheme  = 'm';
ps.h_const = false;

%IO
ps.plot_dt = 1e-2;  
ps.save_as_movie = false;
ps.movie_name = 'out2';
ps.plotstyle = 'pv';
ps.fixaxes.v = [-0.003, 0.003];
ps.fixaxes.p = [-0.005 , 0.005];

 %% material parameter
rho0 = 1;     % density
c0   = 1.0; 

%% domain         
ps.Omega = [-0.2, 1.2]; 
Vp = ps.dx; %volume per particle

%% active particles

%left
leftpoint = 0.0;
rightpoint= 0.3;
v0   = 0.004;
I = ps.add_line1d(leftpoint,rightpoint);
ps.addproperties(I, Vp, rho0, v0,c0)

% middle
leftpoint = max(ps.Xj)+ps.dx;
rightpoint= 0.9;
v0   = 0;
I = ps.add_line1d(leftpoint,rightpoint);
ps.addproperties(I, Vp, rho0, v0,c0)

% % right
% leftpoint = max(ps.Xj)+ps.dx;
% rightpoint= leftpoint+0.1;
% v0   = -0.01;
% I = ps.add_line1d(leftpoint,rightpoint);
% ps.addproperties(I, Vp, rho0, v0,c0)


%% set BC

mirrorParticlesj = false(size(ps.Xj));
mirrorParticlesj(end) = true;    %last particle is mirror particle
outer_normal = 1;
ps.add_bc_nr(mirrorParticlesj,outer_normal);
%ps.bc(1).damping_area=[0.4;0.5];

% %% -----------------------------------------------------
%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

