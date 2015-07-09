% SPH for HVI  - Markus Ganser - TU/e - 2015
% impact scenario

close all; clear;
ps = sph_scenario();

%% general parameter
ps.Ntot    = 20000; %40000
ps.equalmass = false;

ps.dtfactor = 0.5;
ps.dt       = 2e-8;
ps.tend     = 5.5e-5;   
%  ps.tpause   = 0;
ps.eta      = 1.2;     
ps.set_kernel('M4');
ps.scheme   = 'n';
ps.compOmegaj  = true;
ps.h_const  = false;

% IO
ps.plot_dt = 1e-6;   
ps.save_as_movie = false;
ps.plot_quantity = 'xpce';%vpd';
ps.plot_style.p = 'patches';
% ps.plot_style.d = 'trisurf';
%ps.fixaxes.p = [-0.5,0.5];
%output
ps.save_dt = 1e-7;
ps.write_data = true;
ps.output_name ='data/impactMG_laminate2_nobc';
ps.plotconfig.latexplot = false;
%% material parameter
ps.EOS     = 'MG';

%general
ps.Omega = [-0.10,0.20;  %x
            -0.15,0.15]; %y  

%% plate
Py= [-225e-3,225e-3];
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


omega_geo_cut = [0, 20e-3;     %x    
             -100e-3,100e-3];%y 

% set BC
bp = [0,omega_geo_cut(2,1)]; %bottom
outer_normal = [0,-1];
ps.add_bc('cut',bp,outer_normal);
bp = [0,omega_geo_cut(2,2)]; % top
outer_normal = [0,1];
ps.add_bc('cut',bp,outer_normal);


%% protectile
shift = -2e-3;
omega_geo = [-10e-3+shift, shift;
             -5e-3, 5e-3];
rho0 = 2710;   
c0   = 5300; 
e0   = 0;
MG_Gamma = 1.7; %gruneisen parameter
MG_S     = 1.5;
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