% SPH for HVI  - Markus Ganser - TU/e - 2015
% one square with an a singularity

close all; clear; clc;
ps = sph_scenario();

ps.dim     = 2;
ps.dx      = 5e-3;
ps.dt      = 5e-5;   
ps.tend    = 2e-1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 100*ps.dt;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';

%% material parameter
ps.rho0 = 1.0;     %relativ density
ps.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )     K    = 73e9; bulk modulus 
ps.beta = 0;    %for surface tension
ps.mu   = 1;    %for dissipation

%% domain
ps.Omega = [1,1]; 

lowerleftcorner1 = [0.3,0.3];
upperrightcorner1= [ 0.7,0.7];
v0    =[0,0];
noise = 0;
I_new = add_rectangle2d(ps.obj_geo,lowerleftcorner1,...
                 upperrightcorner1,v0,ps.dx,ps.dx,noise);
ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];                            

%set rho0
N        =size(ps.Iin,1);
ps.rho0 =ones(N,1);
ps.rho0(floor(N/2)-100)=1.0001;  % peak in the center

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

