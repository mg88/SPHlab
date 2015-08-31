% SPH for HVI  - Markus Ganser - TU/e - 2015
% Impact scenario 2d - Al-Al

% close all; clear;
ps = SPHlab_scenario();

%% general parameter
ps.Ntot    = 4000; 
ps.equalmass = false;
               
% ps.dtfactor = 0.5;
ps.dt       = 2.5e-9;
ps.tend     = 5e-5;   
%  ps.tpause   = 0;
ps.eta      = 1.2;     
ps.set_kernel('M4');
ps.scheme   = 'm';
ps.compOmegaj  = true;
ps.h_const  = false;

% IO
ps.plot_dt = 1e-7;   
ps.save_as_movie = false;
ps.plot_quantity = 'x';%vpd';
ps.plot_style.p = 'patches';
% ps.plot_style.d = 'trisurf';
%ps.fixaxes.p = [-0.5,0.5];
%output
ps.save_dt = 2e-9;
ps.write_data = true;
ps.output_name ='data/impactMG_40000_bc_hexa_b';
ps.plotconfig.latexplot = false;
%% material parameter
ps.EOS     = 'MG';
e0 = 0;
MG_Gamma = 2;
MG_S     = 1.338;

%general
ps.Omega = [-0.02,0.05;  %x 0.06
            -0.11,0.11]; %y  

%% plate
omega_geo = [0, 5e-3;     %x    
              -100e-3,100e-3];%y 

omega_geo_cut = [0, 5e-3;     %x    
             -20e-3,20e-3];%y 


% set BC
bp = [0,omega_geo_cut(2,1)]; %bottom
outer_normal = [0,-1];
ps.add_bc('cut',bp,outer_normal);
bp = [0,omega_geo_cut(2,2)]; % top
outer_normal = [0,1];
ps.add_bc('cut',bp,outer_normal);

rho0 = 2790;   
c0   = 5330; 
v0 = [0,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)

   

%% protectile
r = 1e-3;
center = [-2e-3 - r, 0];

omega_geo = [r, center];
rho0 = 2790;
c0 = 5330;
v0 = [4360,0];
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

%% create particle class
obj_particles = SPHlab_particles(ps);
%% start simulation
start_simulation(obj_particles)

if isempty(ps.plot_quantity)
    obj_particles.IO.plot_data(obj_particles,'pve')
end