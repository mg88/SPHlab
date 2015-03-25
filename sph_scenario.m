classdef sph_scenario < handle
    % SPH for HVI  - Markus Ganser - TU/e - 2015
    % scenario class - defines several testscenario
    
    properties

        %% simulation parameter
        dim         % 1 or 2
        dx          % initial particle distance
        dt          % timestep (fixed)
        tend        % Simulation time
        eta         % h=eta*dx
        eta2        % eta2*h is the cutoff radius 

        %% geometry
        Omega       % domain
        obj_geo     % geometry class
        
        %% indices domain
        Iboun
        Iin  
        Imaterial  % [1st indice of n-st material, last indice of n-st material]^n
       
        %% material parameter
        rho0     % relative density
        Ca       % Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )
        beta     % for surface tension
        mu       % for dissipation
        
        % 
        g_ext    %gravity

        %% IO
        % input
        read_file
        input_name
        % output
        
        
        % some plotting properties
        plot_dt         % plotting timestep
        plotstyle       % 1D: p-position, v-velocity, p-pressure, d-density, f-forces
                        % 2D: scatter | trisurf | patches  
        %movie settings                        
        save_as_movie
        movie_name
    end
    
    methods
        % Constructor
        function obj = sph_scenario()
           obj.Iboun = [];
           obj.Iin   = [];
           obj.Imaterial = [];
           obj.g_ext = [0,0];
           obj.save_as_movie = false;
           obj.movie_name = 'out';
           obj.obj_geo = sph_geometry();
           obj.plotstyle = 'scatter';
           obj.read_file = false;
        end       
        
        function checkIfAlreadySet(obj)
            if ~isempty(obj.Iin)
                error('A scenario is already set!');
            end
        end
        
        function name = get_movie_name(obj)
            movie_dir='movies/'; %ToDo: anders machen
            movie_format='.avi';
            name=[movie_dir,obj.movie_name,movie_format];
        end
          
        function dispdata(obj)
           disp(obj)             
        end
        
    end    
end

