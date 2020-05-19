classdef NonLinearSystem_1O < handle
    
    properties
        %
        % A*dxdt+B*x+fnl(x)=SM*fex
        %
        name = 'NLS1O';
        nx
        nex
        A
        invA
        B
        fnl
        SM
        fex
        opt
        LSp % Linear part
        LS % Linearized System
        nlsh % homogene differential equation
        isMeanFree
        dim_info % 0: mech. dof; 1: velo. mech. dof; 2: add. state space dim.; 3: filter dim.
    end
    
    methods
        function obj = NonLinearSystem_1O(A,B,fnl,SM,fex,dim_info)
            %% Konstruct Nonlinear System: A*x' + B*x + fnl(x) = SMx*fex
            % A:        Matrix (nx x nx)
            % B:        Matrix (nx x nx)
            % fnl:      NonLinearForce (nx x 1)
            % SM:       Matrix (nx x nex)
            % fex:      RandomProcess (nex x 1)
            % dim_info: Information about dimensions (NaN: unknown; 0: mech. dof; 1: velo. mech. dof; 2: add. state space dim.; 3: filter dim.)
            %% error management:
            n0 = length(A(:,1));
            dim_info_temp = NaN(fnl.n,1);
            if nargin<=5
                dim_info = NaN(fnl.n,1);
            end
            dim_info_temp(1:length(dim_info(:))) = dim_info(:);
            if fnl.n<n0
                error('fnl must be at least (n x 1) is A is (n x n).');
            elseif fnl.n>n0
                temp = eye(fnl.n);
                temp(1:length(A(:,1)),1:length(A(1,:))) = A;
                A = temp;
                temp = zeros(fnl.n);
                temp(1:length(B(:,1)),1:length(B(1,:))) = B;
                B = temp;
                temp = zeros(fnl.n,length(SM(1,:)));
                temp(1:length(SM(:,1)),1:length(SM(1,:))) = SM;
                SM = temp;
                dim_info_temp((n0+1):fnl.n) = 2;
            end
            dim_info = dim_info_temp;
            sizeCheck = [size(A),size(B),length(SM(:,1))];
            if ~prod(sizeCheck==fnl.n) || ~(fex.n==length(SM(1,:)))
                error('Input arguments must have consistent sizes.');
            end
            %% initialize object:
            obj.nx = fnl.n;
            obj.nex = length(SM(1,:));
            obj.A = A;
            obj.invA = inv(A);
            obj.B = B;
            obj.fnl = fnl;
            obj.SM = SM;
            obj.fex = fex;
            obj.opt = optimoptions('fsolve','display','off');
            obj.isMeanFree = false; % can be set true by obj.isMeanFree = true
            obj.dim_info = dim_info;
            %% Extend system in case of non-GWN-excitation:
            if isa(fex,'FilteredProcess')
                obj = fex.applyExcitation(obj);
            end
        end
        
        function obj = setMeanFree(obj,state)
            % set isMeanFree true or false
            obj.isMeanFree = (state==true);
            if ~isempty(obj.nlsh)
                obj.nlsh.isMeanFree = (state==true);
            end
        end
        
        function mux = getMeanX(obj,Kxx)
            if obj.isMeanFree
                mux = zeros(obj.nx,1);
            else
                mfex = obj.SM*obj.fex.getMean();
                R = @(mux) obj.B*mux+obj.fnl.getMuFnl(mux,Kxx)-mfex;
                warning off;
                mux0 = obj.B\mfex;
                warning on;
                mux0(isnan(mux0)) = 0;
                if sum(R(mux0).^2)==0
                    mux = zeros(size(mux0));
                else
                    mux = fsolve(R,mux0,obj.opt);
                end
            end
        end
        
        function Bnl = getBnl(obj,mux,Kxx)
            Bnl = obj.fnl.getBnl(mux,Kxx);
        end
        
        function [fmcs,Gmcs] = getSysMCS(obj,filter)
            if nargin==1
                fmcs = @(t,x) -obj.invA*(obj.B*x+obj.fnl.getFnl(x));
                Gmcs = @(t,x) obj.invA*obj.SM*sqrt(obj.fex.getPSD());
            else
                NLSz = filter.merge(obj);
                [fmcs,Gmcs] = NLSz.getSysMCS();
            end
        end
        
        function [LSp] = getLinearPart(obj)
            %% Linear part, NOT linearized System!
            if isempty(obj.LSp)
                LSp = LinearSystem_1O(obj.A,obj.B,obj.SM,obj.fex.Sxx,obj.invA);
                obj.LSp = LSp;
            else
                LSp = obj.LSp;
            end
        end
        
        function [LS] = getLinearizedSystem(obj,mux,Kxx)
            %% Linearized System
            if isempty(obj.LS)
                LS = LinearSystem_1O(obj.A,obj.B+obj.getBnl(mux,Kxx),obj.SM,obj.fex.Sxx,obj.invA);
                obj.LS = LS;
            else
                obj.LS.updateB(obj.B+obj.fnl.getBnl(mux,Kxx));
                LS = obj.LS;
            end
        end
        
        function [LS] = getLinearizedSystemFromVars(obj,vars)
            %% Linearized System
            [Kxx_r,mux_r] = varsToCov(obj,vars);
            LS = obj.getLinearizedSystem(mux_r,Kxx_r);
        end
    end
    
    methods(Static)
        function [A,B,SMx,dim_info] = getAB(M,C,K,SMq)
            %% Transform Linear System
            %   M*q'' + C*q' + K*q = SMq*fex
            %  To
            %   A*x' + B*x = SMx*fex
            Ndof = length(M(:,1));
            if nargin<4
                SMq = eye(Ndof);
            end
            Nex = length(SMq(1,:));
            test = [size(M),size(C),size(K),length(SMq(:,1))];
            if ~prod(test==Ndof)
                error('Matrices must have the same size');
            end
            A = [eye(Ndof),zeros(Ndof);zeros(Ndof),M];
            B = [zeros(Ndof),-eye(Ndof);K,C];
            SMx = [zeros(Ndof,Nex);SMq];
            dim_info = [zeros(Ndof,1);ones(Ndof,1)];
        end
    end
end