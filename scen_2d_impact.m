% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear; clc;
ps = sph_scenario();

ps.read_file = false;
ps.input_name = 'test';

ps.dim     = 2;
ps.dx      = 1e-2;
ps.dt      = 1e-5;   
ps.tend    = 0.3;    
ps.eta     = 1.3;     
ps.eta2    = 2;  

% IO
ps.plot_dt = 500*ps.dt;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';

%% material parameter
ps.rho0 = 1;     %relativ density
ps.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )     K    = 73e9; bulk modulus 
ps.beta = 0;    %for surface tension
ps.mu   = 1;    %for dissipation

%% domain
ps.Omega = [1,2];  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]

%layer1

%% todo: geht nicht wenn y={0.5,1,5) - wieso?
lowerleftcorner1 = [0.3,  0.75];
upperrightcorner1= [0.35, 1.25];
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
lowerleftcorner1 = [0.15 , 0.98];
upperrightcorner1= [0.28, 1.02];
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

