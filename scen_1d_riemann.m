% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.read_file = false;
ps.input_name = 'test';

ps.dim     = 1;
ps.dx      = 2e-4;
ps.omega   = 0.4;  
ps.tend    = 0.3;    
 
ps.eta     = 1.2;     
ps.eta2    = 2;
ps.kernel  = 'M4';

%IO
ps.plot_dt = 1e-4;  
ps.save_as_movie = false;
ps.plotstyle = 'vp';

 %% material parameter
ps.rho0 = 10;     % density
Ca      = 3.7;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )
ps.c0   = Ca^(-0.5);    % speed of sound  
ps.c0   = 10; 
ps.beta = 0;    %for surface tension
ps.mu   = 0;    %for dissipation

%% domain         
ps.Omega = 0.5;  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]
ps.obj_geo = sph_geometry();


%% active particles

%left
leftpoint = 0.22;
rightpoint= 0.1;
v0   = 2;
I_new = add_line1d(ps.obj_geo,leftpoint,rightpoint,v0,ps.dx);

ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ]; 

% right
leftpoint = leftpoint+ps.dx;
rightpoint= 0.4;
v0   = 0;
I_new = add_line1d(ps.obj_geo,leftpoint,rightpoint,v0,ps.dx);

ps.Iin=[ps.Iin;I_new];
ps.Imaterial = [ps.Imaterial; [I_new(1) I_new(end)] ];     

%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

