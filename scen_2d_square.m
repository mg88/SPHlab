% SPH for HVI  - Markus Ganser - TU/e - 2015
% one square with an a singularity

close all; clear; clc;
ps = sph_scenario();

ps.dx      = 1e-2;
ps.dtfactor= 0.4;
ps.tend    = 0.1;    
ps.eta     = 1.2;     
ps.eta2    = 2;  

%IO
ps.plot_dt = 1e-3;   
ps.save_as_movie = false;
ps.plot_param = 'vp';
ps.plot_style.p = 'patches';
ps.fixaxes.p = [-0.1,0.1];

%% material parameter
rho0 = 1.0;  
c0   = 10;

%% domain
ps.Omega = [0,0.5;
            0,1]; 
Vp = ps.dx*ps.dx; %volume per particle
 
        
lowerleftcorner1 = [0.1,0.1];
upperrightcorner1= [ 0.4,0.9];
v0    = [0,0];
I = add_rectangle2d(ps,lowerleftcorner1,...
                 upperrightcorner1);
addproperties(ps, I, Vp, rho0, v0,c0)
             
%set rho0
N        = size(ps.Iin,1);
ps.rhoj(floor(N/2-20)) = 1.1;  % peak in the center

%% set BC
%no-reflecting bc right
y  = max(ps.Xj(ps.Iin,1));
mirrorParticlesj = abs(ps.Xj(ps.Iin,1) - y) < ps.dx/2;    
outer_normal =[1,0];
ps.add_bc_nr(mirrorParticlesj,outer_normal);
%no-reflecting bc left
y  = min(ps.Xj(ps.Iin,1));
mirrorParticlesj = abs(ps.Xj(ps.Iin,1) - y) < ps.dx/2;    
outer_normal =[-1,0];
ps.add_bc_nr(mirrorParticlesj,outer_normal);

%no-reflecting bc bottom
y  = min(ps.Xj(ps.Iin,2));
mirrorParticlesj = abs(ps.Xj(ps.Iin,2) - y) < ps.dx/2;    
outer_normal =[0,-1];
ps.add_bc_nr(mirrorParticlesj,outer_normal);

%no-reflecting bc top
y  = max(ps.Xj(ps.Iin,2));
mirrorParticlesj = abs(ps.Xj(ps.Iin,2) - y) < ps.dx/2;    
outer_normal =[0,1];
ps.add_bc_nr(mirrorParticlesj,outer_normal);


%% -----------------------------------------------------

%% disp data:
dispdata(ps);

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

