% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear; clc;
ps = sph_scenario();

ps.read_file = false;
ps.input_name = 'test';

ps.dim     = 2;
ps.dx      = 4e-2;
ps.omega   = 0.5;
ps.tend    = 0.5;    
ps.eta     = 1.3;     
ps.eta2    = 2;  
ps.kernel  = 'M4';

% IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plotstyle = 'scatter';

%% material parameter
ps.rho0 = 1;   
ps.c0   = 10; 

%% domain
ps.Omega = [-0.3,0.7;  %x
            -0.5,0.9]; %y  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]

%layer1

lowerleftcorner1 = [0.0,  -0.25];
upperrightcorner1= [0.05, 0.25];
v0 = [0,0];
noise = 0;
I_new = add_rectangle2d(ps.obj_geo,lowerleftcorner1,...
                 upperrightcorner1,v0,ps.dx,ps.dx,noise);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];

%             %layer 2
%             lowerleftcorner1 = [0.85,0.1];
%             upperrightcorner1= [0.9,1.9];
%             v0=[0,0];
%             noise = 0;
%             I_new = add_rectangle2d(ps.obj_geo,lowerleftcorner1,...
%                              upperrightcorner1,v0,ps.dx,ps.dx,noise);
%             ps.Iin=[ps.Iin;I_new];

%projectile
lowerleftcorner1 = [-0.2 , -0.02] ;
upperrightcorner1= [-0.02, 0.02] ;
v0=[1,0];
noise = 0;
I_new = add_rectangle2d(ps.obj_geo,lowerleftcorner1,...
                 upperrightcorner1,v0,ps.dx,ps.dx,noise);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

