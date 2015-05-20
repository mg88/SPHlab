%exact solution of a Riemannshock after Clavijo2013
classdef exact_rarefactionShock < handle
   
    properties
       xj
       x0
       vj
       rhoj
       ej
       pj
        
       Gamma
       Vhead
       Vtail
       Vcontact
       Vshock
       p1
       p3
       p4
       p6
       v1
       v3
       v4
       v6
       rho1
       rho3
       rho4
       rho6
       a1
    end
        
         %%   
    methods
       
        %% constructor
        function obj = exact_rarefactionShock(rhoL,pL,vL,rhoR,pR,vR,Gamma,x0)
            % pL=1;
% pR=0.1;
% vL=0;
% vR=0;
% rhoL=1;
% rhoR=0.15;
% Gamma = 1.4;
% 
% t=0.25;
% x=linspace(0,1,100);
% x0=0.5;


            obj.Gamma=Gamma;
            obj.x0=x0;
            %region 1 and 6
            obj.v1=vL;
            obj.v6=vR;
            obj.p1=pL;
            obj.p6=pR;
            obj.rho1=rhoL;
            obj.rho6=rhoR;

            aL = sqrt(Gamma*pL/rhoL);
            aR = sqrt(Gamma*pR/rhoR);

            obj.a1=aL;
            a6=aR;
            A6=2/(obj.rho6*(Gamma+1));
            B6=obj.p6*(Gamma-1)/(Gamma+1);

            error =@(ps_guess) (ps_guess-obj.p6)*sqrt(A6/(ps_guess+B6))...
                     + 2*obj.a1/(Gamma-1)*((ps_guess/obj.p1)^((Gamma-1)/(2*Gamma)) - 1) + obj.v6 -obj.v1;


            options = optimoptions('fsolve','Display','off');  % Turn off display
            pstar=fsolve(error,1.0662,options);
            
            %region 3 and 4
            obj.p3=pstar;
            obj.p4=pstar;

            obj.v3= obj.v1 - (2*obj.a1)/(Gamma-1) * ( ( obj.p3/obj.p1)^((Gamma-1)/(2*Gamma)) -1 );
            obj.v4= obj.v6 + (obj.p4-obj.p6)*sqrt(A6/(obj.p4+B6));

            obj.rho3 = obj.rho1 * (obj.p3/obj.p1)^(1/Gamma);
            obj.rho4 = obj.rho6*((obj.p6*(Gamma-1)+obj.p4*(Gamma+1))/(obj.p4*(Gamma-1)+obj.p6*(Gamma+1)));



            a3 = sqrt(Gamma*obj.p3/obj.rho3);

            obj.Vhead = obj.v1-obj.a1;
            obj.Vtail = obj.v3-a3;
            obj.Vcontact = obj.v3;
            obj.Vshock  = obj.v6+a6*sqrt((Gamma+1)*obj.p4/(2*Gamma*obj.p6) + (Gamma-1)/2*Gamma);
%             if (obj.v3-obj.v4)> 1e-10
%                 keyboard
%             end            
        end
        
        function compute (obj,t,x)
            if nargin > 2
                obj.xj = x;
            else 
                x=obj.xj;
            end
             obj.vj=x*0;
             obj.rhoj=x*0;
             obj.pj=x*0;   

            %region 1:
            I1 = x-obj.x0<t*obj.Vhead;
            obj.vj(I1) = obj.v1;
            obj.pj(I1) = obj.p1;
            obj.rhoj(I1) = obj.rho1;



            %region 2:
            I2 = ((t*obj.Vhead < (x-obj.x0))+((x-obj.x0)<t*obj.Vtail))==2;
            brackets = ( 2/(obj.Gamma+1) + (obj.Gamma-1)/(obj.a1*(obj.Gamma+1)) * (obj.v1-(x(I2)-obj.x0)/t));
            obj.rhoj(I2) = obj.rho1 * brackets.^(2/(obj.Gamma-1));
            obj.pj(I2)   = obj.p1*brackets.^(2*obj.Gamma/(obj.Gamma-1));
            obj.vj(I2)   = 2/(obj.Gamma+1) * (obj.a1+0.5*(obj.Gamma-1)*obj.v1+(x(I2)-obj.x0)/t);

            %region 3:
            I3 = ((t*obj.Vtail < (x-obj.x0)) + ((x-obj.x0)<t*obj.Vcontact)) == 2; 
            obj.pj(I3)   = obj.p3;
            obj.rhoj(I3) = obj.rho3;
            obj.vj(I3)   = obj.v3;

            %region 4;
            I4 = ((t*obj.Vcontact < (x-obj.x0)) + ((x-obj.x0)<t*obj.Vshock)) == 2;
            obj.pj(I4)   = obj.p4;
            obj.rhoj(I4) = obj.rho4;
            obj.vj(I4)   = obj.v4;

            %region 6;
            I6 = t*obj.Vshock < x-obj.x0;
            obj.pj(I6)  = obj.p6;
            obj.rhoj(I6)= obj.rho6;
            obj.vj(I6)  = obj.v6;

            %ernergy
            obj.ej =obj.pj./obj.rhoj/(obj.Gamma-1);            
        end
        
        function plot(obj)
            subplot(2,2,1)
            plot(obj.xj,obj.rhoj);
            subplot(2,2,2)
            plot(obj.xj,obj.pj)
            subplot(2,2,3)
            plot(obj.xj,obj.vj)
            subplot(2,2,4)
            plot(obj.xj,obj.ej);
        end
    end

end