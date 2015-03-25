% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d boundaries on both sides

close all; clear; clc;
ps = sph_scenario();

ps.dim     = 1;
ps.dx      = 1e-3;
ps.dt      = 0.0001;   
ps.tend    = 1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 100*ps.dt;   
ps.save_as_movie = false;
ps.plotstyle = 'vd';

 %% material parameter
ps.rho0 = 1;     %relativ density
ps.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )
ps.beta = 0;    %for surface tension
ps.mu   = 0.5;%3.5;    %for dissipation

%% domain         
ps.Omega = 1; 
ps.obj_geo = sph_geometry();

%% active particles
leftpoint = 0.08;
rightpoint= 0.88;
v0   = -1;
I_new = add_line1d(ps.obj_geo,leftpoint,rightpoint,v0,ps.dx);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];            

%% boundary:
bounleft1  = 0;
bounright1 = 0.05;
bounleft2  = 0.95;
bounright2 = 1;
v0=0;
I_new=add_line1d(ps.obj_geo,[bounleft1;bounleft2]...
                 ,[bounright1;bounright2],v0,ps.dx);
ps.Iboun=[ps.Iboun;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];                         

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

