% SPH for HVI  - Markus Ganser - TU/e - 2015
% 2d riemann problem

close all; clear; clc;
ps = sph_scenario();
%% general parameter
ps.Ntot    = 500;
% ps.Ntot    = [319,1680];
ps.tend    = 1.8;   
%ps.dtfactor  =0.1;

ps.eta     = 1.7;     
ps.set_kernel('M4');

%IO
ps.plot_dt = 5e-2;  
ps.save_as_movie = false;
ps.plot_quantity = 'vd';%'pv';
ps.plot_style.p = 'trisurf';% 'plot3';patches
% ps.plot_style.d = 'trisurf';
ps.fixaxes.p = [-0.5,0.5];
%ps.fixaxes.d = [1-1e-2,1+1e-2];
%output
ps.save_dt = 5e-2;
ps.write_data = false;
ps.output_name ='data/riemann2d_boun';

 %% material parameter
rho0 = 1;   
c0   = 1;
ps.EOS     = 'ISO';

%ps.art_diss_para.alpha_mass = 0;
%ps.art_diss_para.alpha_energy = 0;
%ps.art_diss_para.alpha_viscosity = 0;
%ps.art_diss_para.beta_viscosity = 0;

%% domain         
ps.Omega = [-0.5,2.5  ;  %x
            -0.5,0.5]; %y
  
%% active particles
l  = 0;
interface=0.15;
r  = 1;
%r  = 1.85;

h1 = -0.1;
h2 = 0.1;
%left
omega_geo = [l,interface; %x
             h1,h2]; %y
v0    = [0.2,0];
ps.add_geometry(omega_geo, rho0, v0, c0)
             
% right
omega_geo = [interface,r; %x
             h1,h2]; %y
v0    = [0,0];
e0  = 0;
ps.add_geometry(omega_geo, rho0, v0, c0, e0)
   
%% set BC
%no-reflecting bc right
bp = [r,0];
outer_normal = [1,0];
ps.add_bc('nrm',bp,outer_normal);

%no-flow on bottom
bp =  [0,h1];
outer_normal = [0,-1];
ps.add_bc('noflow',bp,outer_normal);
%no-flow on top
bp =  [0,h2];
outer_normal = [0,1];
ps.add_bc('noflow',bp,outer_normal);

%% -----------------------------------------------------
% generate particles
ps.create_geometry;

% dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%obj_particles.showInitial()
%% start simulation
start_simulation(obj_particles)

