% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear;
ps = sph_scenario();

ps.write_data = false;
ps.read_data  = false;
ps.input_name = 'test';

ps.dx       = 4e-3;
ps.dtfactor = 0.5;
ps.tend     = 0.02;    
ps.eta      = 1.2;     
ps.eta2     = 2;  
ps.kernel   = 'Wendland';
ps.scheme   = 'm';
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


% t=9.6e-05s (1 iter/sec ,Ncon = 45374)
% t=0.00288s (30 iter/sec ,Ncon = 45452)
% t=0.0056077s (30 iter/sec ,Ncon = 45650)
% create new cell-structure
% t=0.0075279s (21 iter/sec ,Ncon = 57161)
% t=0.0098751s (28 iter/sec ,Ncon = 57645)
% create new cell-structure
% t=0.010409s (8 iter/sec ,Ncon = 64953)
% create new cell-structure
% t=0.011153s (11 iter/sec ,Ncon = 80482)
% create new cell-structure
% t=0.012303s (15 iter/sec ,Ncon = 95553)
% t=0.013888s (19 iter/sec ,Ncon = 95783)
% create new cell-structure
% t=0.014075s (3 iter/sec ,Ncon = 116175)
% t=0.015324s (19 iter/sec ,Ncon = 116428)
% t=0.016557s (19 iter/sec ,Ncon = 116789)
% t=0.017778s (19 iter/sec ,Ncon = 116827)
% t=0.019037s (19 iter/sec ,Ncon = 117102)
% Elapsed time is 24.061831 seconds.