% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.read_data  = false;
ps.input_name = 'test';

ps.dx       = 2e-3;
ps.dtfactor = 0.4;  
ps.tend     = 0.3;    
 
ps.eta     = 1.2;     
ps.eta2    = 2;
ps.kernel  = 'Wendland';
ps.scheme  = 'm';
ps.h_const = false;

%IO
ps.plot_dt = 1e-2;  
ps.save_as_movie = false;
ps.plotstyle = 'pdvf';
ps.fixaxes.v = [-0.001, 0.005];
ps.fixaxes.p = [-2e-2 , 2e-2];

 %% material parameter
rho0 = 1;     % density
c0   = 5.0; 

%% domain         
ps.Omega = [-0.2, 1.2]; 
Vp = ps.dx; %volume per particle

%% active particles

%left
leftpoint = 0.0;
rightpoint= 0.2;
v0   = 0.005;
I = ps.add_line1d(leftpoint,rightpoint);
ps.addproperties(I, Vp, rho0, v0,c0)

% middle
leftpoint = max(ps.Xj)+ps.dx;
rightpoint= 1;
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
%ps.damping_area=[0.4;0.5];
kb = zeros(size(ps.Xj));
kb(end)=1; %last particle is mirror particle
ps.mirrorParticlesj = logical(kb); 

% %% -----------------------------------------------------
%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

