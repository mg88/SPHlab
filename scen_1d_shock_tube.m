% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d shock tube problem - paper of zisis2014
%close all; clear;
ps = sph_scenario();


%% general parameter
ps.Ntot    = 400;
ps.equalmass = true;
ps.tend    = 0.25;   
% ps.tpause  = 0; 
ps.eta     = 1.2;     
ps.set_kernel('Gauss');

ps.scheme     = 'm';
ps.h_const    = false;
ps.compOmegaj = true;
ps.normalizeOmega = false;

%IO
ps.plot_dt = 1e-2;  
ps.save_as_movie = false;
ps.plot_quantity = 'dvpe';
ps.fixaxes.d = [0,1.2];
ps.fixaxes.v = [-0.1,1.1];
ps.fixaxes.p = [-0.1,1.1];
ps.fixaxes.e = [1.5,3];
ps.plotconfig.latexplot = true;
%ps.Neval = ps.Ntot ;
% ps.Nss = 4;

 %% material parameter

%% domain         
ps.Omega = [-0.2, 1.4]; 

%% active particles

%% test 1
% %1
% omega_geo = [0.0,0.6];
% v0   = 1;
% rho0 = 1;     % density
% c0   = 1; %ma=1 (Ma=vimp/c0)
% ps.add_geometry(omega_geo, rho0, v0, c0)
% 
% %2
% omega_geo = [0.6,0.8];
% v0   = 0;
% rho0 = 1;     % density
% c0   = 1;  %ma=1
% ps.add_geometry(omega_geo, rho0, v0, c0)
% 
% %3
% omega_geo = [0.8,1.2];
% v0   = 0;
% rho0 = 0.25;     % density
% c0   = 0.5;  %ma=2
% ps.add_geometry(omega_geo, rho0, v0, c0)
% 
% %4
% omega_geo = [1.2,1.4];
% v0   = 0;
% rho0 = 1;     % density
% c0   = 1;  %ma=1
% ps.add_geometry(omega_geo, rho0, v0, c0)


%% test2
% ps.dt  = 1e-4;
ps.EOS = 'IdealGas14';
Gamma  = 1.4;
% 
% ps.art_diss_para.alpha_mass      = 0.3;%0.3;
% ps.art_diss_para.beta_mass       = 0;
% ps.art_diss_para.alpha_viscosity = 1;
% ps.art_diss_para.beta_viscosity  = 2;
% ps.art_diss_para.alpha_energy    = 0.3;%0.5;
% ps.art_diss_para.beta_energy     = 0;
% 

xl = 0;
x0 = 0.5;
xr = 1.2;
%1 (p=1)
omega_geo = [xl,x0];
vL   = 0;
rhoL = 1;     % density
c0   = 1.1832;  %c0=sqrt(Gamma*p/q)
eL = 2.5;
ps.add_geometry(omega_geo, rhoL, vL, c0,eL)

%2 (p=0.1)
omega_geo = [x0,xr];
vR   = 0;
rhoR = 0.125;     % density
c0   = 1.0583; 
eR = 2;

ps.add_geometry(omega_geo, rhoR, vR, c0,eR)

% %2
% omega_geo = [0.7,1];
% v0   = 0.1;
% rho0 = 1;     % density
% c0   = sqrt(5/4); 
% e0 = 1.5;
% ps.add_geometry(omega_geo, rho0, v0, c0,e0)


%exact solution:
pL = (Gamma-1)*eL*rhoL;
pR = (Gamma-1)*eR*rhoR;

ps.exact_sol = exact_rarefactionShock(rhoL,pL,vL,rhoR,pR,vR,Gamma,x0);

%% BC:
ps.add_bc('noflow',xl,-1); %left
ps.add_bc('nrc',0.8,1); %right

%% -----------------------------------------------------
% generate particles
ps.create_geometry;
% dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)
