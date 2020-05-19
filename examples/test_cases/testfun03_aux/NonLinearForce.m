classdef NonLinearForce < handle
    
    properties
        %% private properties
        n % number of elements
        fnl
        Jnl
        fnl_info
        Jnl_info
        Cov_info
        Mux_info
    end
    
    methods
        %% public methods
        
        function obj = NonLinearForce(n)
            %% constructor
            obj.n = n;
            obj.fnl = cell(n,1);
            obj.Jnl = cell(n,n);
            obj.Cov_info = [];%(((1:n)-1)*n+(1:n))';
            obj.Mux_info = [];
        end
        
        function obj = addElement(obj,NLE,eqns,vars)
            %% add NL_Element
            neq = length(eqns);
            nva = length(vars);
            for j=1:nva
                for k=j:nva
                    obj.Cov_info = [obj.Cov_info;(vars(k)-1)*obj.n+vars(j)];
                end
                obj.Mux_info = [obj.Mux_info;vars(:)];
            end
            obj.Cov_info = unique(obj.Cov_info);
            obj.Mux_info = unique(obj.Mux_info);
            for i=1:neq
                obj.addFnli(eqns(i),vars,NLE.getFnl(i),NLE.getMuFnl(i));
                for j=1:nva
                    obj.addJnlij([eqns(i),vars(j)],vars,NLE.getJnl(i,j),NLE.getBnl(i,j));
                end
            end
        end
        
        function fnl = getFnl(obj,x)
            %% gets nonlinear force vector fnl(x)
            nfnl = length(obj.fnl_info);
            vs = zeros(nfnl,1);
            for k=1:nfnl
                i = obj.fnl_info(k);
                for j=1:length(obj.fnl(i,1,:))
                    if ~isempty(obj.fnl{i,1,j})
                        dim = obj.fnl{i,1,j}{1};
                        vs(k) = vs(k)+obj.fnl{i,1,j}{2}(x(dim));
                    else
                        break;
                    end
                end
            end
            fnl = sparse(obj.fnl_info,ones(nfnl,1),vs,obj.n,1);
        end
        
        function Jnl = getJnl(obj,x)
            %% gets Jacobi-matrix of the nonlinear force vector
            % Jnl(x)=dfnl/dx
            nJnl = length(obj.Jnl_info(:,1));
            is = zeros(nJnl,1);
            js = zeros(nJnl,1);
            vs = zeros(nJnl,1);
            for k=1:nJnl
                i = obj.Jnl_info(k,1);
                j = obj.Jnl_info(k,2);
                is(k) = i;
                js(k) = j;
                for l=1:length(obj.Jnl(i,j,:))
                    if ~isempty(obj.Jnl{i,j,l})
                        dim = obj.Jnl{i,j,l}{1};
                        vs(k) = vs(k)+obj.Jnl{i,j,l}{2}(x(dim));
                    else
                        break;
                    end
                end
            end
            Jnl = sparse(is,js,vs,obj.n,obj.n);
        end
        
        function mufnl = getMuFnl(obj,mux,Kxx)
            %% gets mean of the nonlinear force E{fnl(x)}=fun(mux,Kxx)
            nfnl = length(obj.fnl_info);
            vs = zeros(nfnl,1);
            for k=1:nfnl
                i = obj.fnl_info(k);
                for j=1:length(obj.fnl(i,1,:))
                    if ~isempty(obj.fnl{i,1,j})
                        dim = obj.fnl{i,1,j}{1};
                        vs(k) = vs(k)+obj.fnl{i,1,j}{3}(mux(dim),Kxx(dim,dim));
                        % TODO: Pruefen, ob 'obj.fnl{i,1,j}{3}' belegt ist.
                        % Sonst 'obj.fnl{i,1,j}{2}' numerisch integrieren?
                    else
                        break;
                    end
                end
            end
            mufnl = sparse(obj.fnl_info,ones(nfnl,1),vs,obj.n,1);
        end
        
        function Bnl = getBnl(obj,mux,Kxx)
            %% gets equivalent linear matrix E{Jnl(x)}=fun(mux,Kxx)
            nJnl = length(obj.Jnl_info(:,1));
            is = zeros(nJnl,1);
            js = zeros(nJnl,1);
            vs = zeros(nJnl,1);
            for k=1:nJnl
                i = obj.Jnl_info(k,1);
                j = obj.Jnl_info(k,2);
                is(k) = i;
                js(k) = j;
                for l=1:length(obj.Jnl(i,j,:))
                    if ~isempty(obj.Jnl{i,j,l})
                        dim = obj.Jnl{i,j,l}{1};
                        vs(k) = vs(k)+obj.Jnl{i,j,l}{3}(full(real(mux(dim))),full(real(Kxx(dim,dim))));
                        % TODO: Pruefen, ob 'obj.Jnl{i,j,l}{3}' belegt ist.
                        % Sonst 'obj.Jnl{i,j,l}{2}' numerisch integrieren?
                    else
                        break;
                    end
                end
            end
            Bnl = sparse(is,js,vs,obj.n,obj.n);
        end
        
        function nlf = expand(obj,n)
            %% expands force to n dimensions and returns new object
            nlf = NonLinearForce(n);
            f = cell(n,1,length(obj.fnl(1,1,:)));
            f(1:obj.n,1,:) = obj.fnl;
            J = cell(n,n,length(obj.Jnl(1,1,:)));
            J(1:obj.n,1:obj.n,:) = obj.Jnl;
            nlf.fnl = f;
            nlf.Jnl = J;
            nlf.fnl_info = obj.fnl_info;
            nlf.Jnl_info = obj.Jnl_info;
            nlf.Cov_info = (ceil(obj.Cov_info/obj.n)-1)*n+obj.Cov_info-floor(obj.Cov_info/obj.n-10^-8)*obj.n;
            nlf.Mux_info = obj.Mux_info;
        end
        
        function nlf = convertTo1O(obj)
            %% expands force to n dimensions and returns new object
            nlf = NonLinearForce(2*obj.n);
            f = cell(2*obj.n,1,length(obj.fnl(1,1,:)));
            f(obj.n+(1:obj.n),1,:) = obj.fnl;
            J = cell(2*obj.n,2*obj.n,length(obj.Jnl(1,1,:)));
            J(obj.n+(1:obj.n),1:obj.n,:) = obj.Jnl;
            nlf.fnl = f;
            nlf.Jnl = J;
            nlf.fnl_info = obj.fnl_info;
            nlf.Jnl_info = obj.Jnl_info;
%             nlf.Cov_info = (ceil(obj.Cov_info/obj.n)-1)*n+obj.Cov_info-floor(obj.Cov_info/obj.n-10^-8)*obj.n;
error('Das funktioniert noch nicht');
            nlf.Mux_info = obj.Mux_info;
        end
    end
    
    methods (Access = private)
        %% private methods
        
        function obj = addFnli(obj,i,vars,fnli,mufnli)
            %% adds element to the nonlinear force
            if ~((i>=1) && (i<=obj.n))
                error('Index out of bounds');
            end
            if nargin>4
                fnli = {vars,fnli,mufnli};
            else
                fnli = {vars,fnli};
            end
            j = 1;
            try
                while ~isempty(obj.fnl{i,1,j})
                    j=j+1;
                end
            catch
                % k neu belegen
            end
            obj.fnl{i,1,j} = fnli;
            obj.fnl_info = unique([obj.fnl_info;i]);
        end
        
        function obj = addJnlij(obj,ij,vars,Jnlij,Blinij)
            %% adds element to the Jacobi-matrix
            i = ij(1);
            j = ij(2);
            if ~((i>=1) && (i<=obj.n) && (j>=1) && (j<=obj.n))
                error('Index out of bounds');
            end
            if nargin>4
                Jnlij = {vars,Jnlij,Blinij};
            else
                Jnlij = {vars,Jnlij};
            end
            k = 1;
            try
                while ~isempty(obj.Jnl{i,j,k})
                    k=k+1;
                end
            catch
                % k neu belegen
            end
            obj.Jnl{i,j,k} = Jnlij;
            if isempty(obj.Jnl_info) || (sum(obj.Jnl_info(obj.Jnl_info(:,1)==i,2)==j)==0)
                obj.Jnl_info = [obj.Jnl_info;i,j];
            end
        end
        
    end
end