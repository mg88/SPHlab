% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

% close all; clear;
ps = SPHlab_scenario();

%% general parameter
ps.Ntot    = 2000;
ps.equalmass = false;

% ps.dt       = 2e-10;
ps.tend     = 1e-6;  
% ps.tpause   = [4e-6];
ps.eta      = 1.2;     
ps.set_kernel('M4');
ps.EOS     = 'MG';
ps.scheme   = 'n';
ps.compOmegaj  = true;
ps.h_const  = false;

% IO
ps.plot_dt = 2e-7;   
ps.save_as_movie = false;
ps.movie_name = 'impact_laminateV2';
ps.plot_quantity = 'x';%vpd';
ps.plot_style.p = 'patches';
ps.plot_style.d = 'trisurf';
ps.fixaxes.p = [-5e8,5e8];
%output
ps.save_dt = 1e-8;
ps.write_data = true;
ps.output_name ='data/impact_laminate_test';

%% material parameter

%general
% ps.Omega = [-0.02,0.05;  %x 0.06
%             -0.025,0.025]; %y  
ps.Omega = [-0.01,0.015;  %x 0.06
    -0.02,0.02]; %y          

%% plate
%Py= [-20e-3,20e-3];
Py= [-15e-3,15e-3];

n_layer =5;
Px= [0, 5e-3];
x=Px(1);
dxAlu=1e-3;
dxEpoxy=1e-3;
v0 = [0,0];
for i=2:n_layer+1
    if mod(i,2)==0 %alu2004
        omega_geo=[x,x+dxAlu ;
           Py];
        x = x+dxAlu;
        rho0 = 2790;   
        c0   = 5330; 
        e0   = 0;
        MG_Gamma = 2;
        MG_S     = 1.338;
        ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)
    else  %mat2 (Stycast Epoxy)
        omega_geo=[x,x+dxEpoxy ;
           Py];
       x=x+dxEpoxy;
        rho0 = 1140;   
        c0   = 2580; 
        e0   = 0;
        MG_Gamma = 0.67; %gruneisen parameter
        MG_S     = 1.47;
        ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)        
    end
end


 % set BC
bp = [0,omega_geo(2,1)]; %bottom
outer_normal = [0,-1];
ps.add_bc('nrc',bp,outer_normal);
bp = [0,omega_geo(2,2)]; % top
outer_normal = [0,1];
ps.add_bc('nrc',bp,outer_normal); 

%% protectile
r = 1e-3;
center = [-2e-3 - r, 0];

omega_geo = [r, center];
rho0 = 2790;
c0 = 5330;
v0 = [4360,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)



% %% use symmetry 
% ps.Ntot=ps.Ntot/2;
% for i = 1:(ps.iter_geo-1)
%     ps.geo(i).omega_geo(2,1) = 0;
% end
% bp = [0,0];
% outer_normal = [0,-1];
% ps.add_bc('noflow',bp,outer_normal);
% ps.Omega(2,1)= -0.02;

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