classdef sph_scenarios < handle
    % SPH for HVI  - Markus Ganser - TU/e - 2015
    % scenario class - defines several testscenario
    
    properties
        %% simulation parameter
        dim         % 1 or 2
        dx          % initial particle distance
        dt          % timestep (fixed)
        tend        % Simulation time
        plot_dt     % plotting timestep
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
        rho0     % relativ density
        Ca       % Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )
        beta     % for surface tension
        mu       % for dissipation
        
        % 
        g_ext    %gravity
        
        % some plotting properties
        plotstyle       % 1D: p-position, v-velocity, p-pressure, d-density, f-forces
                        % 2D: scatter | trisurf | patches  
        save_as_movie
        movie_name
    end
    
    methods
        % Constructor
        function obj = sph_scenarios()
           obj.Iboun = [];
           obj.Iin   = [];
           obj.Imaterial = [];
           obj.g_ext = [0,0];
           obj.save_as_movie = false;
           obj.movie_name = 'out';
           obj.obj_geo = sph_geometry();
           obj.plotstyle = 'scatter';
        end       
        
        function checkIfAlreadySet(obj)
            if ~isempty(obj.Iin)
                error('A scenario is already set!');
            end
        end
        
        function name = get_movie_name(obj)
            movie_dir='movies/';
            movie_format='.avi';
            name=[movie_dir,obj.movie_name,movie_format];
        end
          
        function dispdata(obj)
           disp(obj)             
        end
        
        %% test scenarios
        function set_1d_boun(obj)
            checkIfAlreadySet(obj);
            obj.dim     = 1;
            obj.dx      = 1e-3;
            obj.dt      = 0.0001;   
            obj.tend    = 1;    
            obj.plot_dt = 100*obj.dt;   
            obj.eta     = 1.2;     
            obj.eta2    = 2;  
            obj.save_as_movie = false;

             %% material parameter
            obj.rho0 = 1;     %relativ density
            obj.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )
            obj.beta = 0;    %for surface tension
            obj.mu   = 0.5;%3.5;    %for dissipation
            
            %% domain         
            obj.Omega = 1;  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]
            obj.obj_geo = sph_geometry();
           
            
            %% active particles
            leftpoint = 0.08;
            rightpoint= 0.88;
            v0   = -1;
            I_new = add_line1d(obj.obj_geo,leftpoint,rightpoint,v0,obj.dx);
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];            

            %% boundary:
            bounleft1  = 0;
            bounright1 = 0.05;
            bounleft2  = 0.95;
            bounright2 = 1;
            v0=0;
            I_new=add_line1d(obj.obj_geo,[bounleft1;bounleft2]...
                             ,[bounright1;bounright2],v0,obj.dx);
            obj.Iboun=[obj.Iboun;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];                         
          %%
            dispdata(obj);                         
        end

        function set_1d_riemann(obj)
            checkIfAlreadySet(obj);
            obj.dim     = 1;
            obj.dx      = 5e-4;
            obj.dt      = 1e-5;   
            obj.tend    = 0.3;    
            obj.plot_dt = 10*obj.dt;   
            obj.eta     = 1.2;     
            obj.eta2    = 2;  
            obj.save_as_movie = false;
            obj.plotstyle = 'vd';

             %% material parameter
            obj.rho0 = 1;     %relativ density
            obj.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )
            obj.beta = 0;    %for surface tension
            obj.mu   = 0.5;%3.5;    %for dissipation
            
            %% domain         
            obj.Omega = 1;  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]
            obj.obj_geo = sph_geometry();
           
            
            %% active particles
            
            %left
            leftpoint = 0.22;
            rightpoint= 0.1;
            v0   = 1;
            I_new = add_line1d(obj.obj_geo,leftpoint,rightpoint,v0,obj.dx);
            
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ]; 
            
            % right
            leftpoint = leftpoint+obj.dx;
            rightpoint= 0.7;
            v0   = 0;
            I_new = add_line1d(obj.obj_geo,leftpoint,rightpoint,v0,obj.dx);
            
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];     
            
                       
          %%
            dispdata(obj);                         
        end
  
        function set_2d_riemann(obj)
            checkIfAlreadySet(obj);
            obj.dim     = 1;
            obj.dx      = 1e-2;
            obj.dt      = 1e-5;   
            obj.tend    = 0.1;    
            obj.plot_dt = 10*obj.dt;   
            obj.eta     = 1.2;     
            obj.eta2    = 2;  
            obj.save_as_movie = false;
            obj.plotstyle = 'patches';

             %% material parameter
            obj.rho0 = 1;     %relativ density
            obj.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )
            obj.beta = 0;    %for surface tension
            obj.mu   = 3;%3.5;    %for dissipation
            
            %% domain         
            obj.Omega = [1,0.5]; 
            obj.obj_geo = sph_geometry();
           
            
            %% active particles
            l  = 0.2;
            r  = 0.7;
            h1 = 0.1;
            h2 = 0.4;
            interface=0.35;
            %left
            lowerleftcorner1 = [l,h1];
            upperrightcorner1= [ interface,h2];
            v0    =[3,0];
            noise = 0;
            I_new = add_rectangle2d(obj.obj_geo,lowerleftcorner1,...
                             upperrightcorner1,v0,obj.dx,obj.dx,noise);
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];                            
            % right
            lowerleftcorner1 = [interface+obj.dx,h1];
            upperrightcorner1= [ r,h2];
            v0    =[0,0];
            noise = 0;
            I_new = add_rectangle2d(obj.obj_geo,lowerleftcorner1,...
                             upperrightcorner1,v0,obj.dx,obj.dx,noise);
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];                            
            %% boundary
            %top
            startpoint = [0,h1-obj.dx];
            endpoint   = [1,h1-obj.dx];
            v0   = [0,0];
            noise = 0;
            layer = 1;
            I_new = add_line2d(obj.obj_geo,startpoint,endpoint,...
                v0,obj.dx,layer,noise);   
            obj.Iboun=[obj.Iboun;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];
            
             %bottom
            startpoint = [0,h2+obj.dx];
            endpoint   = [1,h2+obj.dx];
            v0   = [0,0];
            noise = 0;
            layer = 1;
            I_new = add_line2d(obj.obj_geo,startpoint,endpoint,...
                v0,obj.dx,layer,noise);   
            obj.Iboun=[obj.Iboun;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];
                      
            
          %%
            dispdata(obj);                         
        end        
        
        function set_2d_moving_lines(obj)
            checkIfAlreadySet(obj);           
            obj.dim     = 2;
            obj.dx      = 1e-2;
            obj.dt      = 2e-5;   
            obj.tend    = 0.5;    
            obj.plot_dt = 100*obj.dt;   
            obj.eta     = 1.2;     
            obj.eta2    = 2;  
            obj.save_as_movie = false;
            obj.plotstyle = 'patches';

            %% material parameter
            obj.rho0 = 1;     %relativ density
            obj.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )     K    = 73e9; bulk modulus 
            obj.beta = 0;    %for surface tension
            obj.mu   = 0;    %for dissipation
             
            %% domain
            obj.Omega = [1,1];  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]
 
            %line 1
            startpoint = [0.22,0.49];
            endpoint   = [0.55,0.49];
            v0   = [0,-0.5];
            noise = 0;
            layer = 2;
            I_new = add_line2d(obj.obj_geo,startpoint,endpoint,...
                v0,obj.dx,layer,noise);   
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];
            
            %line 2
            startpoint = [0.42,0.41];
            endpoint   = [0.8,0.41];
            v0   = [0,0];
            noise = 0;
            layer = 2;
            I_new = add_line2d(obj.obj_geo,startpoint,endpoint,...
                v0 ,obj.dx,layer,noise);
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];
            
            %%
            dispdata(obj);
        end
    
        function set_2d_impact(obj)
            checkIfAlreadySet(obj);
            obj.dim     = 2;
            obj.dx      = 5e-3;
            obj.dt      = 1e-5;   
            obj.tend    = 0.3;    
            obj.plot_dt = 500*obj.dt;   
            obj.eta     = 1.3;     
            obj.eta2    = 2;  
            obj.save_as_movie = false;
            obj.plotstyle = 'patches';

            %% material parameter
            obj.rho0 = 1;     %relativ density
            obj.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )     K    = 73e9; bulk modulus 
            obj.beta = 0;    %for surface tension
            obj.mu   = 1;    %for dissipation
             
            %% domain
            obj.Omega = [1,2];  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]
 
            %layer1
            lowerleftcorner1 = [0.3,0.1];
            upperrightcorner1= [0.35,1.9];
            v0=[0,0];
            noise = 0;
            I_new = add_rectangle2d(obj.obj_geo,lowerleftcorner1,...
                             upperrightcorner1,v0,obj.dx,obj.dx,noise);
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];

%             %layer 2
%             lowerleftcorner1 = [0.85,0.1];
%             upperrightcorner1= [0.9,1.9];
%             v0=[0,0];
%             noise = 0;
%             I_new = add_rectangle2d(obj.obj_geo,lowerleftcorner1,...
%                              upperrightcorner1,v0,obj.dx,obj.dx,noise);
%             obj.Iin=[obj.Iin;I_new];
            
            %projectile
            lowerleftcorner1 = [0.15 , 0.98];
            upperrightcorner1= [0.28, 1.02];
            v0=[1,0];
            noise = 0;
            I_new = add_rectangle2d(obj.obj_geo,lowerleftcorner1,...
                             upperrightcorner1,v0,obj.dx,obj.dx,noise);
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];

            %%
            %%
            dispdata(obj);
        end
        
        function set_2d_two_squares(obj)
            checkIfAlreadySet(obj);
            obj.dim     = 2;
            obj.dx      = 1e-2;
            obj.dt      = 2e-4;   
            obj.tend    = 5;    
            obj.plot_dt = 100*obj.dt;   
            obj.eta     = 1.2;     
            obj.eta2    = 2;  
            obj.save_as_movie = false;
            
            %% material parameter
            obj.rho0 = 1;     %relativ density
            obj.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )     K    = 73e9; bulk modulus 
            obj.beta = 0.1;    %for surface tension
            obj.mu   = 10;    %for dissipation
             
            %% domain
            obj.Omega = [1,1];  %upper right corner  lower-left ist zero: [0,Omega(1)]x[0,Omega(2)]

            lowerleftcorner1 = [0.45,0.35];
            upperrightcorner1= [ 0.55,0.6];
            lowerleftcorner2 = [0.35,0.35];
            upperrightcorner2= [ 0.45,0.45];
            v0=[0,0];
            noise = 0;
            I_new = add_rectangle2d(obj.obj_geo,[lowerleftcorner1;lowerleftcorner2],...
                             [upperrightcorner1;upperrightcorner2],v0,obj.dx,obj.dx,noise);
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];
            dispdata(obj);
        end
        
        function set_2d_square(obj)
            checkIfAlreadySet(obj);
            obj.dim     = 2;
            obj.dx      = 5e-3;
            obj.dt      = 5e-5;   
            obj.tend    = 2e-1;    
            obj.plot_dt = 100*obj.dt;   
            obj.eta     = 1.2;     
            obj.eta2    = 2;  
            obj.save_as_movie = false;
            obj.plotstyle = 'patches';
            %% material parameter
            obj.rho0 = 1.0;     %relativ density
            obj.Ca   = 1.7e-2;   %Cauchy number (dependend on the buld modulus Ca=rho(x,0)*u0^2/K(x) )     K    = 73e9; bulk modulus 
            obj.beta = 0;    %for surface tension
            obj.mu   = 1;    %for dissipation
             
            %% domain
            obj.Omega = [1,1]; 

            lowerleftcorner1 = [0.3,0.3];
            upperrightcorner1= [ 0.7,0.7];
            v0    =[0,0];
            noise = 0;
            I_new = add_rectangle2d(obj.obj_geo,lowerleftcorner1,...
                             upperrightcorner1,v0,obj.dx,obj.dx,noise);
            obj.Iin=[obj.Iin;I_new];
            obj.Imaterial = [obj.Imaterial; [I_new(1) I_new(end)] ];                            
            
            %set rho0
            N        =size(obj.Iin,1);
            obj.rho0 =ones(N,1);
            obj.rho0(floor(N/2)-100)=1.0001;  % peak in the center
                      
            dispdata(obj);
        end
    end
    
end

