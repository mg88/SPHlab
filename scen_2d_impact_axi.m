% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

warning('this version is not working right now!')
keyboard

close all; clear;
ps = SPHlab_scenario();

%% general parameter
ps.Ntot    = 200; %40000
ps.equalmass = false;

ps.dtfactor = 0.2;
% ps.dt       = 2e-9;
ps.tend     = 7e-5;   
% ps.tpause   = (1:1:10)*1e-7;
ps.eta      = 1.2;     
ps.set_kernel('M3');
ps.scheme   = 'a';
ps.compOmegaj  = false;
ps.h_const  = true;

% IO
ps.plot_dt = 1e-8;   
ps.save_as_movie = false;
ps.plot_quantity = 'xpfd';%vpd';
% ps.plot_style.x = 'ringcloud';
ps.plot_style.p = 'patches';
ps.plot_style.d = 'patches';
% ps.plot_style.d = 'trisurf';
%ps.fixaxes.p = [-0.5,0.5];
%output
ps.save_dt = 1e-7;
ps.write_data = false;
ps.output_name ='data/impactISO_org3';
ps.plotconfig.latexplot=false;


%% material parameter
ps.EOS     = 'ISO';
e0 = 0;
MG_Gamma = 2;
MG_S     = 1.338;

ps.art_diss_para.alpha_mass      = 0;
ps.art_diss_para.beta_mass       = 0;
ps.art_diss_para.alpha_viscosity = 0;
ps.art_diss_para.beta_viscosity  = 0;
ps.art_diss_para.alpha_energy    = 0;
ps.art_diss_para.beta_energy     = 0;

%general
ps.Omega = [-0.1,0.15;  %x
            -0.25,0.25]; %y  

%% plate
omega_geo = [0, 20e-3;     %x    
              -225e-3,225e-3];%y 

omega_geo_cut = [0, 20e-3;     %x    
             -60e-3,60e-3];%y 


% set BC
bp = [0,omega_geo_cut(2,2)]; % top
outer_normal = [0,1];
% ps.add_bc('nrc',bp,outer_normal);

rho0 = 2785;   
c0   = 5330; 
v0 = [0,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)
    
   

%% protectile
omega_geo = [-5e-3-2e-3, -2e-3;
             -2.5e-3, 2.5e-3];
rho0 = 2785;
c0 = 5330;
v0 = [0.*c0,0];
% ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)

%% use symmetry 
ps.Ntot = ps.Ntot/2;
for i = 1:(ps.iter_geo-1)
    ps.geo(i).omega_geo(2,1) = 0;
end
bp=[0,0];
outer_normal = [0,-1];
ps.add_bc('noflow',bp,outer_normal);
ps.Omega(2,1)= -0.2;

%% -----------------------------------------------------
% generate particles
ps.create_geometry;


%dispdata(ps);

%% create particle class
obj_particles = SPHlab_particles(ps);
%% start simulation
start_simulation(obj_particles)

if isempty(ps.plot_quantity)
    obj_particles.IO.plot_data(obj_particles,'pve')
end