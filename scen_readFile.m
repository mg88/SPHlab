% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.read_file = true;
ps.input_name = 'test';


%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

            %% in
            %Gamma,
            %Gmod,
            %S, Y0, 
            %c, c0,
            %e,
            %m,
            %p,
            %phi,
            %rho, rho0, 
            % tauXX, tauXY, tauYY,
            % u, v, x,y