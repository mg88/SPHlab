% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 1000;
ps.equalmass = false;

ps.dtfactor = 0.5;
ps.tend     = 0.5;    
ps.eta      = 1.2;     
ps.set_kernel('Wendland');
ps.scheme   = 'm';
ps.h_const  = false;
% IO
ps.plot_dt = 1e-2;   
ps.save_as_movie = false;
ps.plot_quantity = '';%vpd';
ps.plot_style.p = 'patches';
ps.plot_style.d = 'trisurf';
ps.fixaxes.p = [-0.5,0.5];
%output
ps.save_dt = 5e-2;
ps.write_data = true;
ps.output_name ='data/impact';

%% material parameter
rho0 = 1;   
c0   = 2; 

%general
ps.Omega = [-0.5,0.7;  %x
            -0.5,0.5]; %y  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]

%layer1
omega_geo = [0,0.05;     %x
             -0.25,0.25];%y 
v0 = [0,0];
ps.add_geometry(omega_geo, rho0, v0, c0)


%   % layer 2
% lowerleftcorner1 = [0.2,  -0.25];
% upperrightcorner1= [0.25, 0.25];
% I = add_rectangle2d(ps,lowerleftcorner1, upperrightcorner1);
% addproperties(ps, I, Vp, rho0, v0,c0)

%projectile
omega_geo = [-0.05, -0.02;
             -0.005,0.005];
c0 = 2;
v0 = [2,0];
ps.add_geometry(omega_geo, rho0, v0, c0)


%% -----------------------------------------------------
% generate particles
ps.create_geometry;
%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%obj_particles.showInitial()
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

%rel dissipation:
%hconst,v: 0.014
%hconst,m: 0.02
%hvariable,v, F(dWij): 0.0081
%hvariable,m, F(dWij): 0.00088
%hvariable,v, F(dWij_hi,dWij_hj): 0.001
%hvariable,m, F(dWij_hi,dW_hj): 0.01