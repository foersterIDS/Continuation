classdef NL_Element_DUF < NL_Element
    %% Duffing element
    
    properties
        lambda
    end
    
    methods
        function obj = NL_Element_DUF(lambda)
            % fnl = lambda*x^3
            obj = obj@NL_Element(1,1); % Super class constructor
            obj.lambda = lambda;
        end
        
        function fnl = getFnl(obj,eqn)
            %% gets function of the nonlinear force-element fnl(x)
            fnl = @(x) obj.lambda*x.^3;
        end
        
        function Jnl = getJnl(obj,eqn,var)
            %% gets function of the Jacobi-matrix-element of the nonlinear force
            Jnl = @(x) 3*obj.lambda*x.^2;
        end
        
        function mufnl = getMuFnl(obj,eqn)
            %% gets function of the mean of the nonlinear force E{fnl(x)}=fun(mux,Kxx)
            mufnl = @(mux,Kxx) obj.lambda*mux*(3*Kxx+mux^2);
        end
        
        function Bnl = getBnl(obj,eqn,var)
            %% gets function of the equivalent linear matrix E{Jnl(x)}=fun(mux,Kxx)
            Bnl = @(mux,Kxx) 3*obj.lambda*(Kxx+mux^2);
        end
            
    end
end