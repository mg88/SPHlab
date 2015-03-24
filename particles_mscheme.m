classdef particles_mscheme < particles
    %PARTICLES_MSCHEME Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        % Constructor
        function obj = particles_mscheme(obj_geo,rho0,h,Iin,Iboun)
             obj = obj@particles(obj_geo,rho0,h,Iin,Iboun);
        end
        
        function test(obj)
           disp (obj.N)
        end
        
        function comp_Fdiss(obj,mu)      %page415-violeau | particle-friction   
            obj.F_diss=zeros(obj.N,obj.dim);
            disp('nichts');
        end
        
        
    end
    
end

