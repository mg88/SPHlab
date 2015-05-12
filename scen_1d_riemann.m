% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

%% general parameter
ps.Ntot      = 80;
ps.equalmass = false;
%ps.dt        = 1e-3;
%ps.dtfactor  = 0.01;
ps.tend      = 1.9; 

ps.eta     = 1.2;     
ps.set_kernel('Wendland');

ps.scheme  = 'm';
ps.EOS     = 'MG';
ps.h_const = false;

%experimental settings:
ps.exp_settings.tweakmass = false;


%IO
ps.plot_dt = 5e-2;  
ps.save_as_movie = false;
ps.movie_name = 'out2';
ps.plot_quantity = 'v';  %p
ps.fixaxes.v = [-0.005, 0.005];
ps.fixaxes.p = [-0.005, 0.005];
ps.fixaxes.d = [1-1e-1 ,1+1e-1 ];
ps.Neval = 100;
%output
ps.save_dt = 5e-2;
ps.write_data = false;
ps.output_name ='data/riemann_boun_mass';
 %% material parameter
rho0 = 1;     % density
c0   = 1.0; 

ps.art_diss_para.alpha_energy = 0;
%% domain         
ps.Omega = [-1.5, 1.5]; 

%% active particles
l=-1;
m=0.2;
r=1;
%left
omega_geo = [l,-m]; %x
v0   = 0;
ps.add_geometry(omega_geo, rho0, v0, c0)


%middle
omega_geo = [-m,m]; %x
v0   = 0.004;
ps.add_geometry(omega_geo, rho0, v0, c0)

% right
omega_geo = [m,r]; %x
v0   = 0;
ps.add_geometry(omega_geo, rho0, v0, c0)

% right %bigger particles:
% omega_geo = [0.6,0.7];
% v0   = 0.00;
% Nfactor = 1;
% ps.add_geometry(omega_geo, rho0, v0, c0,Nfactor)

%% set BC
%left
p1 = l;
outer_normal = -1;
ps.add_bc_nr(p1,0,outer_normal);

%right
p1 = r;
outer_normal = 1;
ps.add_bc_nr(p1,0,outer_normal);

%experinamtal: damping area:
%ps.bc(1).damping_area=[0.4;0.5];

%% -----------------------------------------------------
% generate particles
ps.create_geometry;
%% disp data:
%dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);

%obj_particles.IO.plot_data(obj_particles,'dpv');

%% start simulation
start_simulation(obj_particles)

