% SPH for HVI  - Markus Ganser - TU/e - 2015
% one square with an a singularity

close all; clear; clc;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 500;
ps.dtfactor= 0.4;
% ps.tpause  = 0;%(1:10)*2e-5;
ps.tend    = 8e-4;    
ps.eta     = 1.2;     

%IO
ps.plot_dt = 1e-6;   
ps.save_as_movie = true;
ps.movie_name = 'test_nrbc_square';

ps.plot_quantity = 'vdp';
ps.plot_style.p = 'trisurf';
ps.plot_style.d = 'trisurf';
ps.plot_style.e = 'trisurf';

% ps.fixaxes.p = 1e9*[-0.5,0.5];
%  ps.fixaxes.d = [2980,3050];

ps.plotconfig.latexplot = true;

%% material parameter
rho0 = 3000.0;  
c0   = 5000;
e0=0;
MG_Gamma = 2; %(alu)
MG_S     = 1.338;
ps.EOS     = 'MG';

%% domain
ps.Omega = [-0.3,0.3;
            -0.6,0.6]; 
% 

l= -0.15;
r= 0.15;
h1= -0.4;
h2= 0.4;
omega_geo = [l,r;
             h1,h2];        
v0 = [0,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)
             
%% set BC
%no-reflecting bc right
bp = [r,0];
outer_normal =[1,0];
% ps.add_bc('nrc',bp,outer_normal);
%no-reflecting bc left
bp = [l,0];
outer_normal =[-1,0];
ps.add_bc('nrc',bp,outer_normal);

%no-reflecting bc bottom
bp =  [0,h1];
outer_normal =[0,-1];
ps.add_bc('nrc',bp,outer_normal);

%no-reflecting bc top
bp =  [0,h2];
outer_normal =[0,1];
% ps.add_bc('nrc',bp,outer_normal);


%% -----------------------------------------------------
% generate particles
ps.create_geometry;

%% 
%set rho0
% N        = size(ps.Iin,1);
% ps.rhoj(floor(N/2)) = ps.rhoj(floor(N/2))*3;  % peak in the center

%set rho0
x1 = -0.02;
x2 = 0.02;
y1 = -0.02;
y2 = 0.02;
k=(ps.Xj(:,1) > x1).*(ps.Xj(:,1) < x2) .* (ps.Xj(:,2)>y1) .*(ps.Xj(:,2)<y2);
ps.rhoj(logical(k)) = ps.rhoj(logical(k))*1.2;  % peak in the center



% dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%obj_particles.showInitial()
%% start simulation
start_simulation(obj_particles)
