% SPH for HVI  - Markus Ganser - TU/e - 2015
% 1d riemann problem

close all; clear; clc;
ps = sph_scenario();

ps.read_data = true;
ps.input_name = 'in-hvi.h5';

ps.write_data = false;

ps.kernel = 'M4';
ps.Omega = [ -0.1, 0.4;
            -0.3, 0.3];
%    --scheme n \
%    --weakly-compressible 0 \
    ps.eta  = 1.2;
    ps.eta2  = 2;

    %--alphamom 10 \
    %--alphamass 0.5 \
    %--alphaen 0.0 \
    %--load $h5name \
    ps.plot_dt = 1e-4;
%     --plot-every 100 \
%    --out $outname \
%     --ids @@ \
%     --dt 1e-6 \
    ps.tend = 0.001;
    ps.dtfactor = 0.5;
%     --help \
%     --verbose

%IO
ps.plot_dt = 1e-6;  
ps.save_as_movie = false;
ps.plotstyle = 'patches';

%% create particle class
obj_particles = sph_particles(ps);
%% start simulation
start_simulation(obj_particles)

            %% in
            %Gamma,
            %Gmod, Shear modulus
            %S,  S-EOS
            %Y0,  Yield strngth
            
            %c, c0,
            %e,
            %m,
            %p,
            %phi,  %surface tension?
            %rho, rho0, 
            % tauXX, tauXY, tauYY,
            % u, v, x,y