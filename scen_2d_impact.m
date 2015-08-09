% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

% close all; clear;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 50000; %40000
ps.equalmass = false;

% ps.dtfactor = 0.5;
ps.dt       = 2e-9;
ps.tend     = 5e-5;   
%  ps.tpause   = 0;
ps.eta      = 1.2;     
ps.set_kernel('M4');
ps.scheme   = 'm';
ps.compOmegaj  = true;
ps.h_const  = false;

% IO
ps.plot_dt = 1e-6;   
ps.save_as_movie = false;
ps.plot_quantity = 'x';%vpd';
ps.plot_style.p = 'patches';
% ps.plot_style.d = 'trisurf';
%ps.fixaxes.p = [-0.5,0.5];
%output
ps.save_dt = 1e-8;
ps.write_data = true;
ps.output_name ='data/impactMG_50001';
ps.plotconfig.latexplot = false;
%% material parameter
ps.EOS     = 'MG';
e0 = 0;
MG_Gamma = 2;
MG_S     = 1.338;

%general
ps.Omega = [-0.10,0.20;  %x
            -0.225,0.225]; %y  

%% plate
omega_geo = [0, 20e-3;     %x    
              -200e-3,200e-3];%y 

omega_geo_cut = [0, 20e-3;     %x    
             -70e-3,70e-3];%y 


% set BC
bp = [0,omega_geo_cut(2,1)]; %bottom
outer_normal = [0,-1];
% ps.add_bc('cut',bp,outer_normal);
bp = [0,omega_geo_cut(2,2)]; % top
outer_normal = [0,1];
% ps.add_bc('cut',bp,outer_normal);

rho0 = 2785;   
c0   = 5330; 
v0 = [0,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)

   

%% protectile
shift = -2e-3;
omega_geo = [-10e-3+shift, shift;
             -5e-3, 5e-3];
rho0 = 2785;
c0 = 5330;
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

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

if isempty(ps.plot_quantity)
    obj_particles.IO.plot_data(obj_particles,'pve')
end