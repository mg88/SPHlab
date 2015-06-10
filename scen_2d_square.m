% SPH for HVI  - Markus Ganser - TU/e - 2015
% one square with an a singularity

close all; clear; clc;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 1000;
ps.dtfactor= 0.4;

ps.tend    = 0.08;    
ps.eta     = 1.2;     

%IO
ps.plot_dt = 1e-3;   
ps.save_as_movie = false;
ps.movie_name = 'test_nrbc1';

ps.plot_quantity = 'vpde';
ps.plot_style.p = 'patches';
ps.plot_style.d = 'trisurf';
ps.fixaxes.p = [-0.1,0.1];
ps.fixaxes.d = [1-1e-3,1+1e-3];


%% material parameter
rho0 = 1.0;  
c0   = 10;
MG_Gamma = 2; %(alu)
MG_S     = 1.338;
ps.EOS     = 'MG';

%% domain
ps.Omega = [-0.1,0.6;
            -0.1,1.1]; 
% 

l= 0.1;
r= 0.4;
h1= 0.4;
h2= 0.9;
omega_geo = [l,r;
             h1,h2];        
v0 = [0,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)
             
%% set BC
%no-reflecting bc right
p1 = [r,0];
p2 = [r,1];
outer_normal =[1,0];
ps.add_bc('nrc',p1,p2,outer_normal);
%no-reflecting bc left
p1 = [l,0];
p2 = [l,1]; 
outer_normal =[-1,0];
% ps.add_bc('nrc',p1,p2,outer_normal);

%no-reflecting bc bottom
p1 =  [0,h1];
p2 =  [1,h1];
outer_normal =[0,-1];
ps.add_bc('nrc',p1,p2,outer_normal);

%no-reflecting bc top
p1 =  [0,h2];
p2 =  [1,h2];
outer_normal =[0,1];
ps.add_bc('nrc',p1,p2,outer_normal);


%% -----------------------------------------------------
% generate particles
ps.create_geometry;

%% 
%set rho0
N        = size(ps.Iin,1);
ps.rhoj(floor(N/2)) = 1.6;  % peak in the center


% dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%obj_particles.showInitial()
%% start simulation
start_simulation(obj_particles)
