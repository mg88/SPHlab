% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear;
ps = sph_scenario();

ps.write_data = false;
ps.read_data  = false;
ps.input_name = 'test';

ps.dx       = 8e-3;
ps.dtfactor = 0.5;
ps.tend     = 0.3;    
ps.eta      = 1.3;     
ps.eta2     = 2;  
ps.kernel   = 'M4';
ps.scheme   = 'c';
% IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plotstyle = 'patches';

%% material parameter
rho0 = 1;   
c0   = 2; 

%general
ps.Omega = [-0.5,0.7;  %x
            -0.5,0.5]; %y  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]
Vp = ps.dx*ps.dx; %volume per particle
  
%layer1
lowerleftcorner1 = [0.0,  -0.25];
upperrightcorner1= [0.05, 0.25];
v0 = [0,0];
I = add_rectangle2d(ps,lowerleftcorner1, upperrightcorner1);
addproperties(ps, I, Vp, rho0, v0,c0, false)


  % layer 2
lowerleftcorner1 = [0.2,  -0.25];
upperrightcorner1= [0.25, 0.25];
I = add_rectangle2d(ps,lowerleftcorner1, upperrightcorner1);
addproperties(ps, I, Vp, rho0, v0,c0, false)

%projectile
lowerleftcorner1 = [-0.1 , -0.01] ;
upperrightcorner1= [-0.02, 0.01] ;
c0 = 5;
v0 = [3,0];
I = add_rectangle2d(ps,lowerleftcorner1, upperrightcorner1);
addproperties(ps, I, Vp, rho0, v0,c0, false)


%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

