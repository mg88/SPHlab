% SPH for HVI  - Markus Ganser - TU/e - 2015
% 2d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.dx      = 4e-3;
ps.tend    = 0.1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-3;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';
ps.fixaxes.p = [-0.5,0.5];

 %% material parameter
rho0 = 1;   
c0   = 10;

%% domain         
ps.Omega = [0.1,0.6  ;  %x
            0,0.5]; %y
Vp = ps.dx*ps.dx; %volume per particle
  
%% active particles
l  = 0.2;
r  = 0.5;
h1 = 0.1;
h2 = 0.4;
interface=0.35;
%left
lowerleftcorner1 = [l,h1];
upperrightcorner1= [ interface,h2];
v0    = [0.1,0];
I = add_rectangle2d(ps,lowerleftcorner1, upperrightcorner1);
addproperties(ps, I, Vp, rho0, v0, c0)
             
% right
lowerleftcorner1 = [max(ps.Xj(:,1))+ps.dx,h1];
upperrightcorner1= [ r,h2];
v0    = [0,0];
I = add_rectangle2d(ps,lowerleftcorner1, upperrightcorner1);
addproperties(ps, I, Vp, rho0, v0, c0)
%% boundary
%top
% startpoint = [0,h1-ps.dx];
% endpoint   = [1,h1-ps.dx];
% v0   = [0,0];
% layer = 1;
% I = add_line2d(ps,startpoint,endpoint,layer);   
%addproperties(ps, I, Vp, rho0, v0,c0, false)

%  %bottom
% startpoint = [0,h2+ps.dx];
% endpoint   = [1,h2+ps.dx];
% v0   = [0,0];
% layer = 1;
% I = add_line2d(ps,startpoint,endpoint,layer);   
%addproperties(ps, I, Vp, rho0, v0,c0, false)

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

