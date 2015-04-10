% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear;
ps = sph_scenario();

ps.write_data = false;
ps.read_data  = false;
ps.input_name = 'test';

ps.dx       = 4e-3;
ps.dtfactor = 0.5;
ps.tend     = 0.05;    
ps.eta      = 1.2;     
ps.eta2     = 2;  
ps.kernel   = 'M4';
ps.scheme   = 'v';
ps.h_const  = false;
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


%   % layer 2
% lowerleftcorner1 = [0.2,  -0.25];
% upperrightcorner1= [0.25, 0.25];
% I = add_rectangle2d(ps,lowerleftcorner1, upperrightcorner1);
% addproperties(ps, I, Vp, rho0, v0,c0, false)

%projectile
lowerleftcorner1 = [-0.1 , -0.01] ;
upperrightcorner1= [-0.02, 0.01] ;
c0 = 20;
v0 = [5,0];

I = add_rectangle2d(ps,lowerleftcorner1, upperrightcorner1);
addproperties(ps, I, Vp, rho0, v0,c0, false)


%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

% t=0.000104s (1 iter/sec ,Ncon = 42526)
% t=0.003952s (38 iter/sec ,Ncon = 42743)
% t=0.0076962s (37 iter/sec ,Ncon = 43197)
% t=0.010194s (25 iter/sec ,Ncon = 43691)
% t=0.013536s (33 iter/sec ,Ncon = 44337)
% t=0.016803s (32 iter/sec ,Ncon = 45206)
% t=0.020016s (31 iter/sec ,Ncon = 45814)
% t=0.021754s (17 iter/sec ,Ncon = 45909)
% t=0.024511s (26 iter/sec ,Ncon = 45945)
% t=0.027194s (25 iter/sec ,Ncon = 45945)
% t=0.030013s (26 iter/sec ,Ncon = 45945)
% t=0.031594s (15 iter/sec ,Ncon = 45945)
% t=0.034295s (25 iter/sec ,Ncon = 45945)
% t=0.036874s (24 iter/sec ,Ncon = 45945)
% t=0.039565s (25 iter/sec ,Ncon = 45945)

%Elapsed time is 20.941137 seconds.

