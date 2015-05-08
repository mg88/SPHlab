% SPH for HVI  - Markus Ganser - TU/e - 2015
% 2d riemann problem

close all; clear; clc;
ps = sph_scenario();
%% general parameter
%ps.Ntot    = 2000;
ps.Ntot    = [319,1680];
ps.tend    = 1.6;   
%ps.dtfactor  =0.1;

ps.eta     = 1.2;     
ps.set_kernel('Wendland');

%IO
ps.plot_dt = 5e-2;  
ps.save_as_movie = false;
ps.plot_quantity = 'xvpe';
ps.plot_style.p = 'patches';% 'plot3';
ps.plot_style.d = 'trisurf';
ps.fixaxes.p = [-0.005,0.005];
%ps.fixaxes.d = [1-1e-2,1+1e-2];
%output
ps.save_dt = 5e-2;
ps.write_data = false;
ps.output_name ='data/riemann2d_boun';

 %% material parameter
rho0 = 1;   
c0   = 1;

%% domain         
ps.Omega = [-0.2,2.2  ;  %x
            -0.3,0.3]; %y
  
%% active particles
l  = 0;
interface=0.15;
r  = 1;
%r  = 1.85;

h1 = -0.2;
h2 = 0.2;
%left
omega_geo = [l,interface; %x
             h1,h2]; %y
v0    = [0.004,0];
ps.add_geometry(omega_geo, rho0, v0, c0)
             
% right
omega_geo = [interface,r; %x
             h1,h2]; %y
v0    = [0,0];
e0  = 0;
ps.add_geometry(omega_geo, rho0, v0, c0, e0)
   
%% set BC
%no-reflecting bc right
p1 = [r,0];
p2 = [r,1];
outer_normal = [1,0];
ps.add_bc_nr(p1,p2,outer_normal);

%no-flow on bottom
p1 =  [0,h1];
p2 =  [1,h1];
outer_normal = [0,-1];
ps.add_bc_noflow(p1,p2,outer_normal);
%no-flow on top
p1 =  [0,h2];
p2 =  [1,h2];
outer_normal = [0,1];
ps.add_bc_noflow(p1,p2,outer_normal);

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

