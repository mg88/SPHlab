% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 5000;
ps.equalmass = false;

ps.dtfactor = 0.5;
ps.tend     = 1e-5;    
ps.eta      = 1.2;     
ps.set_kernel('Wendland');
ps.EOS     = 'Water';
ps.scheme   = 'm';
ps.h_const  = false;
ps.exp_settings.tweakmass = true;

% IO
ps.plot_dt = 1e-6;   
ps.save_as_movie = false;
ps.plot_quantity = 'e';%vpd';
ps.plot_style.p = 'patches';
ps.plot_style.d = 'trisurf';
%ps.fixaxes.p = [-0.5,0.5];
%output
ps.save_dt = 5e-4;
ps.write_data = false;
ps.output_name ='data/impact';

%% material parameter

%general
ps.Omega = [-0.1,0.2;  %x
            -0.1,0.1]; %y  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]

%data from LimeSPH-file        

omega_geo = [0, 20e-3;     %x
             -75e-3,75e-3];%y 
           %  -225e-3,225e-3];%y 

%% set BC
p1 = [0,omega_geo(2,1)];
p2 = [1,omega_geo(2,1)];
outer_normal = [0,-1];
%ps.add_bc_nr(p1,p2,outer_normal);
p1 = [0,omega_geo(2,2)];
p2 = [1,omega_geo(2,2)];
outer_normal = [0,1];
%ps.add_bc_nr(p1,p1,outer_normal);

rho0 = 2785;   
c0   = 5330; 
v0 = [0,0];
ps.add_geometry(omega_geo, rho0, v0, c0)

   

%protectile
omega_geo = [-5e-3-2e-3, -2e-3;
             -2.5e-3, 2.5e-3];
rho0 = 2785;
c0 = 5330;
v0 = [4500,0];
ps.add_geometry(omega_geo, rho0, v0, c0)


%% -----------------------------------------------------
% generate particles
ps.create_geometry;
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

%rel dissipation:
%hconst,v: 0.014
%hconst,m: 0.02
%hvariable,v, F(dWij): 0.0081
%hvariable,m, F(dWij): 0.00088
%hvariable,v, F(dWij_hi,dWij_hj): 0.001
%hvariable,m, F(dWij_hi,dW_hj): 0.01