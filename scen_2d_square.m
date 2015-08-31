% SPH for HVI  - Markus Ganser - TU/e - 2015

% A alumnium square with a density singularity in the center. Test scenario
% for non-reflecting BC.

% close all; clear; clc;
ps = SPHlab_scenario();

%% general parameter
ps.Ntot    = 1500;
ps.tpause  = 10e-4;%[1e-5,6e-5,12e-5];
ps.tend    = 4e-4;    
ps.eta     = 1.2;     

%IO
ps.plot_dt = 1e-6;   
ps.save_as_movie = true;
ps.movie_name = 'test_nrbc_square2';

ps.plot_quantity = 'dvp';
ps.plot_style.p = 'patches';
ps.plot_style.d = 'patches';
ps.plot_style.e = 'trisurf';

ps.fixaxes.p = 1e9*[-0.5,0.5];
ps.fixaxes.d = [2784,2810];

ps.plotconfig.latexplot = true;

%% material parameter
rho0 = 2790.0;  
c0   = 5330;
e0=0;
MG_Gamma = 2; %(alu)
MG_S     = 1.338;
ps.EOS     = 'MG';

%% domain
ps.Omega = [-0.35,0.25;
            -0.65,0.55]; 
% 

l= -0.15;
r= 0.15;
h1= -0.405;
h2= 0.395;
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
% ps.rhoj(logical(k)) = ps.rhoj(logical(k))*1.2;  % peak in the center
% 

%set rho0
r=0.02;
% k=(sum((ps.Xj.^2),2).^0.5)<r;
ps.rhoj(logical(k)) = ps.rhoj(logical(k))*1.2;  % peak in the center



% dispdata(ps);

%% create particle class
obj_particles = SPHlab_particles(ps);
%obj_particles.showInitial()
%% start simulation
start_simulation(obj_particles)
