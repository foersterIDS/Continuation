classdef NL_Element_ECS < NL_Element
    %% Duffing element
    
    properties
        mu
        fn
        kmu
        chi
    end
    
    methods
        function obj = NL_Element_ECS(mu,fn,kmu)
            % fnl = nu*x^2*dxdt
            obj = obj@NL_Element(2,2); % Super class constructor
            obj.mu = mu;
            obj.fn = fn;
            obj.kmu = kmu;
            obj.chi = mu*fn/kmu;
        end
        
        function fnl = getFnl(obj,equ)
            %% gets function of the nonlinear force-element fnl(x)
            if equ==1
                fnl = @(x) obj.kmu*x(2)*(heaviside(x(2)+obj.chi)-heaviside(x(2)-obj.chi))+obj.kmu*obj.chi*(heaviside(x(2)-obj.chi)*heaviside(x(1))-heaviside(-x(2)-obj.chi)*heaviside(-x(1)));
            elseif equ==2
                fnl = @(x) -x(1)*(heaviside(x(2)+obj.chi)-heaviside(x(2)-obj.chi)+heaviside(x(2)-obj.chi)*heaviside(-x(1))+heaviside(-x(2)-obj.chi)*heaviside(x(1)));
            else
                error('max. number of equations is %d',obj.ne);
            end
        end
        
        function Jnl = getJnl(obj,equ,var)
            %% gets function of the Jacobi-matrix-element of the nonlinear force
            if equ==1
                if var==1
                    Jnl = @(x) obj.chi*obj.kmu*dirac(x(1)).*(heaviside(-x(2)-obj.chi)+heaviside(x(2)-obj.chi));
                elseif var==2
                    Jnl = @(x) obj.kmu*x(2)*(-dirac(-obj.chi+x(2))+dirac(obj.chi+x(2)))+obj.chi*obj.kmu*(dirac(obj.chi+x(2))*heaviside(-x(1))+dirac(-obj.chi+x(2))*heaviside(x(1)))+obj.kmu*(-heaviside(-obj.chi+x(2))+heaviside(obj.chi+x(2)));
                else
                    error('not a non-zeros derivative');
                end
            elseif equ==2
                if var==1
                    Jnl = @(x) -(heaviside(x(1))*heaviside(-obj.chi-x(2))-heaviside(-obj.chi+x(2))+heaviside(-x(1))*heaviside(-obj.chi+x(2))+x(1)*(dirac(x(1))*heaviside(-obj.chi-x(2))-dirac(x(1))*heaviside(-obj.chi+x(2)))+heaviside(obj.chi+x(2)));
                elseif var==2
                    Jnl = @(x) -x(1)*(-dirac(-obj.chi+x(2))+dirac(obj.chi+x(2))+dirac(-obj.chi+x(2))*heaviside(-x(1))-dirac(obj.chi+x(2))*heaviside(x(1)));
                else
                    error('not a non-zeros derivative');
                end
            else
                error('max. number of equations is %d',obj.ne);
            end
        end
        
        function mufnl = getMuFnl(obj,equ)
            %% gets function of the mean of the nonlinear force E{fnl(x)}=fun(mux,Kxx)
            %   Hier wird die Annahme verwendet, dass dqdt und xi
            %   mittelwertfrei sind!
            if equ==1
                mufnl = @(mux,Kxx) 0;
            elseif equ==2
                mufnl = @(mux,Kxx) 0;
            else
                error('max. number of equations is %d',obj.ne);
            end
        end
        
        function Bnl = getBnl(obj,equ,var)
            %% gets function of the equivalent linear matrix E{Jnl(x)}=fun(mux,Kxx)
            %   Hier wird die Annahme verwendet, dass dqdt und xi
            %   mittelwertfrei sind!
            if equ==1
                if var==1
                    Bnl = @(mux,Kxx) obj.kmu*obj.chi/(sqrt(real(Kxx(1,1)))*sqrt(2*pi))*erfc(full(obj.chi/(sqrt(real(Kxx(2,2)))*sqrt(2*(1-(real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))^2)))));
                elseif var==2
                    Bnl = @(mux,Kxx) obj.kmu*erf(full(obj.chi/(sqrt(real(Kxx(2,2)))*sqrt(2))))-obj.kmu*obj.chi/(sqrt(real(Kxx(2,2)))*sqrt(2*pi))*exp(-obj.chi^2/(2*sqrt(real(Kxx(2,2)))^2))*erfc(full((real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))*obj.chi/(sqrt(real(Kxx(2,2)))*sqrt(2*(1-(real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))^2)))));
                else
                    error('not a non-zeros derivative');
                end
            elseif equ==2
                if var==1
                    Bnl = @(mux,Kxx) -(1/2*(1+erf(full(obj.chi/(sqrt(real(Kxx(2,2)))*sqrt(2)))))-1/sqrt(pi)*integral(@(v) exp(-v.^2).*erf(full((real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))*v./sqrt(1-(real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))^2))),obj.chi/(sqrt(real(Kxx(2,2)))*sqrt(2)),inf));
                elseif var==2
                    Bnl = @(mux,Kxx) -(-(real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))*obj.chi*sqrt(real(Kxx(1,1)))/(sqrt(real(Kxx(2,2)))^2*sqrt(2*pi))*exp(-obj.chi^2/(2*sqrt(real(Kxx(2,2)))^2))*(1+erf(full((real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))*obj.chi/(sqrt(real(Kxx(2,2)))*sqrt(2*(1-(real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))^2))))))-sqrt(real(Kxx(1,1)))*sqrt(1-(real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))^2)/(sqrt(real(Kxx(2,2)))*pi)*exp(-obj.chi^2/(2*(1-(real(Kxx(1,2))/sqrt(real(Kxx(1,1))*real(Kxx(2,2))))^2)*sqrt(real(Kxx(2,2)))^2)));
                else
                    error('not a non-zeros derivative');
                end
            else
                error('max. number of equations is %d',obj.ne);
            end
        end
            
    end
end