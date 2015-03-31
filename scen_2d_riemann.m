% SPH for HVI  - Markus Ganser - TU/e - 2015
% 2d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.dim     = 1;
ps.dx      = 1e-2;
ps.tend    = 0.1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-3;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';

 %% material parameter
ps.rho0 = 1;   
ps.c0   = 10;

%% domain         
ps.Omega = [0,1  ;  %x
            0,0.5]; %y

%% active particles
l  = 0.2;
r  = 0.7;
h1 = 0.1;
h2 = 0.4;
interface=0.35;
%left
lowerleftcorner1 = [l,h1];
upperrightcorner1= [ interface,h2];
v0    =[3,0];
noise = 0;
I_new = add_rectangle2d(ps,lowerleftcorner1,...
                 upperrightcorner1,v0,ps.dx,ps.dx,noise);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];                            
% right
lowerleftcorner1 = [interface+ps.dx,h1];
upperrightcorner1= [ r,h2];
v0    = [0,0];
I_new = add_rectangle2d(ps,lowerleftcorner1,...
                 upperrightcorner1,v0,ps.dx,ps.dx,noise);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];                            
%% boundary
%top
% startpoint = [0,h1-ps.dx];
% endpoint   = [1,h1-ps.dx];
% v0   = [0,0];
% noise = 0;
% layer = 1;
% I_new = add_line2d(ps.obj_geo,startpoint,endpoint,...
%     v0,ps.dx,layer,noise);   
% ps.Iboun=[ps.Iboun;I_new];
% ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];
% 
%  %bottom
% startpoint = [0,h2+ps.dx];
% endpoint   = [1,h2+ps.dx];
% v0   = [0,0];
% noise = 0;
% layer = 1;
% I_new = add_line2d(ps.obj_geo,startpoint,endpoint,...
%     v0,ps.dx,layer,noise);   
% ps.Iboun=[ps.Iboun;I_new];
% ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];


%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

