% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 20000;
ps.equalmass = false;

ps.dtfactor = 0.5;
ps.dt       = 1e-8;
ps.tend     = 2e-4;   
% ps.tpause   = 0;
ps.eta      = 1.2;     
ps.set_kernel('M4');
ps.EOS     = 'MG';
ps.scheme   = 'm';
ps.compOmegaj  = false;
ps.h_const  = false;

% IO
ps.plot_dt = 1e-6;   
ps.save_as_movie = true;
ps.movie_name = 'impact_laminate';
ps.plot_quantity = 'x';%vpd';
ps.plot_style.p = 'patches';
ps.plot_style.d = 'trisurf';
%ps.fixaxes.p = [-0.5,0.5];
%output
ps.save_dt = 1e-6;
ps.write_data = false;
ps.output_name ='data/impactLayer_bc';

%% material parameter

%general
ps.Omega = [-0.1,0.2;  %x
            -0.1,0.1]; %y  

%data from LimeSPH-file        

%% plate
Py= [-70e-3,70e-3];
n_layer =3;
Px= [0, 20e-3];
x=Px(1);
dx=(Px(2)-Px(1))/n_layer;
v0 = [0,0];
for i=2:n_layer+1
    omega_geo=[x,x+dx ;
           Py];
    if mod(i,2)==0 %mat1 (Alu series 1000 Libersky 1993)
        rho0 = 2710;   
        c0   = 5300; 
        e0   = 0;
        MG_Gamma = 1.7; %gruneisen parameter
        MG_S     = 1.5;
        ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)
    else  %mat2 (Stycast Epoxy)
        rho0 = 2300;   
        c0   = 1800; 
        e0   = 0;
        MG_Gamma = 0.67; %gruneisen parameter
        MG_S     = 1.63;
        ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)        
    end
    x=x+dx;
end


 % set BC
p1 = [0,omega_geo(2,1)]; %bottom
p2 = [1,omega_geo(2,1)];
outer_normal = [0,-1];
ps.add_bc('nrc',p1,p2,outer_normal);
p1 = [0,omega_geo(2,2)]; % top
p2 = [1,omega_geo(2,2)];
outer_normal = [0,1];
ps.add_bc('nrc',p1,p2,outer_normal); 

%% protectile %alu series 2024
omega_geo = [-5e-3-2e-3, -2e-3;
             -2.5e-3, 2.5e-3];
rho0 = 2790;
c0 = 5330;
e0 = 0;
MG_Gamma = 2;
MG_S     = 1.338;
v0 = [4500,0];
ps.add_geometry(omega_geo, rho0, v0, c0, e0, MG_Gamma, MG_S)


%% use symmetry 
ps.Ntot=ps.Ntot/2;
for i = 1:(ps.iter_geo-1)
    ps.geo(i).omega_geo(2,1) = 0;
end
p1 = [0,0];
p2 = [1,0];
outer_normal = [0,-1];
ps.add_bc('noflow',p1,p2,outer_normal);
ps.Omega(2,1)= -0.02;

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