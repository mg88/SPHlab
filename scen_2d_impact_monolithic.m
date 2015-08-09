% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario like in Zisis2014 paper 

close all; clear;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 800; %40000 %1600
ps.equalmass = false;

ps.dtfactor = 0.5;
% ps.dt       = 1e-9;
ps.tend     = 4e-7;%0.45e-6;   
%  ps.tpause   = 0;
ps.eta      = 1.2;     
ps.set_kernel('M4');
ps.scheme   = 'm';
ps.compOmegaj  = true;
ps.h_const  = false;

% IO
ps.plot_dt = 1e-8;   
ps.save_as_movie = false;
ps.plot_quantity = 'p';%vpd';
ps.plot_style.p = 'patches';
% ps.plot_style.d = 'trisurf';
ps.plot_style.e = 'patches';

ps.fixaxes.p = [0,1e11];
%output
ps.save_dt = 1e-8;
ps.write_data = false;
ps.output_name ='data/impactMG_monolithic';
ps.plotconfig.latexplot = true;
%% material parameter
ps.EOS     = 'MG';
e0 = 0;


% ps.art_diss_para.alpha_mass = 0.5;
% ps.art_diss_para.beta_mass = 0;
% ps.art_diss_para.alpha_viscosity = 1;
% ps.art_diss_para.beta_viscosity = 2;
% ps.art_diss_para.alpha_energy = 1;
% ps.art_diss_para.beta_energy = 0;

%general
ps.Omega = [-0.005,0.005;  %x
            -0.007,0.007]; %y  

%% plate
omega_geo = [0, 3e-3;     %x    
              -5e-3,5e-3];%y 

omega_geo_cut = [0, 3e-3;     %x    
             -5e-3,5e-3];%y 


% set BC
bp = [0,omega_geo_cut(2,1)]; %bottom
outer_normal = [0,-1];
ps.add_bc('nrc',bp,outer_normal);
bp = [0,omega_geo_cut(2,2)]; % top
outer_normal = [0,1];
ps.add_bc('nrc',bp,outer_normal);

% rho0 = 2700;   
% c0   = 6320; %m/s
% MG_Gamma = 2;
% MG_S     = 1.338;
MG_Gamma = 1.7;
MG_S     = 1.5;
rho0 = 2710;   
c0   = 5300; 
v0 = [0,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)

   

%% protectile
shift = 1e-3;
omega_geo = [-1e-3-shift,0-shift;
             -2.5e-3, 2.5e-3];
% rho0 = 2700;
% c0 = 6320;
v0 = [1.5*c0,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)


%% use symmetry 
% ps.Ntot=ps.Ntot/2;
% for i = 1:(ps.iter_geo-1)
%     ps.geo(i).omega_geo(2,1) = 0;
% end
% p1 = [0,0];
% outer_normal = [0,-1];
% ps.add_bc('noflow',bp,outer_normal);
% ps.Omega(2,1)= -0.02;

%% -----------------------------------------------------
% generate particles
ps.create_geometry;
%dispdata(ps);
size(ps.Xj)
%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
tttic = tic;
start_simulation(obj_particles)
toc(tttic)
% if isempty(ps.plot_quantity)
%     obj_particles.IO.plot_data(obj_particles,'pve')
% end