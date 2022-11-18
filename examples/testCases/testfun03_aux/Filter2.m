classdef Filter2 < handle
    % Filter mit DGL 2. Ordnung
    
    properties
        Nex;
        ny;
        D;
        Om;
        sigw;
        a;
        b;
        lsy;
        name = 'filter2';
    end
    
    methods
        function obj = Filter2(Nex,D,Om,sigw,ab)
            % Filter2: Filter mit DGL 2. Ordnung
            %   Nex: number of filter processes
            %   D: modal damping vector of the filter
            %   Om: dominant frequencies
            %   sigw: STD of output processes
            %   ab (optional): H(s)=(ab(1)+ab(2)*s)/(s^2+2*D*Om*s+Om^2)
            obj.Nex = Nex;
            obj.ny = 2*Nex;
            obj.D = D(:);
            obj.Om = Om(:);
            if nargin>=5
                sz = size(ab);
                if sz(2)~=2 || (sz(1)~=1 && sz(1)~=Nex)
                    error('ab must be (1 x 2) or (Nex x 2)');
                else
                    obj.a = ab(:,1);
                    obj.b = ab(:,2);
                end
            else
                obj.a = 1;
                obj.b = 0;
            end
            if length(sigw)==1 || length(sigw)==Nex
                obj.sigw = sigw(:);
            elseif mod(Nex,length(sigw))==0
                ns = Nex/length(sigw);
                obj.sigw = kron(ones(ns,1),sigw(:));
            else
                error('Cannot handle STD input.');
            end
            Ay = eye(2*Nex);
            By = [zeros(Nex),-eye(Nex);diag(Om.^2).*eye(Nex),diag(2*D.*Om).*eye(Nex)];
            SMy = [zeros(Nex);eye(Nex)];
            Sy = diag((4*obj.D.*obj.Om.*obj.sigw.^2)./(obj.a.^2./obj.Om.^2+obj.b.^2)).*eye(obj.Nex);
            obj.lsy = LinearSystem_1O(Ay,By,SMy,Sy);
        end
        
        function obj = updateOm(obj,om)
            obj.Om = om;
            By = [zeros(obj.Nex),-eye(obj.Nex);diag(obj.Om.^2).*eye(obj.Nex),diag(2*obj.D.*obj.Om).*eye(obj.Nex)];
            Sy = diag((4*obj.D.*obj.Om.*obj.sigw.^2)./(obj.a.^2./obj.Om.^2+obj.b.^2)).*eye(obj.Nex);
            obj.lsy.updateB(By);
            obj.lsy.SFF = Sy;
        end
        
        function syz = merge(obj,syx)
            % merge: merges filter with a LinearSystem_1O/NonLinearSystem_1O
            %   syx: LinearSystem_1O/NonLinearSystem_1O
            if strcmp(syx.name,'LS1O')
                nx = syx.Ndof;
                Nexx = syx.Nex;
                if Nexx~=obj.Nex
                    error('System input does not match filter output!');
                end
                Az = [syx.A,zeros(nx,2*obj.Nex);
                      zeros(2*obj.Nex,nx),eye(2*obj.Nex)];
                Bz = [syx.B,-syx.SM.*(ones(nx,1)*(obj.a'.*ones(1,Nexx))),-syx.SM.*(ones(nx,1)*(obj.b'.*ones(1,Nexx)));
                      zeros(2*obj.Nex,nx),obj.lsy.B];
                SMz = [zeros(nx,obj.Nex);obj.lsy.SM];
                Sz = obj.lsy.SFF;
                syz = LinearSystem_1O(Az,Bz,SMz,Sz);
            elseif strcmp(syx.name,'NLS1O')
                nx = syx.nx;
                Nexx = syx.nex;
                nz = nx+2*obj.Nex;
                if Nexx~=obj.Nex
                    error('System input does not match filter output!');
                end
                Az = [syx.A,zeros(nx,2*obj.Nex);
                      zeros(2*obj.Nex,nx),eye(2*obj.Nex)];
                Bz = [syx.B,-syx.SM.*(ones(nx,1)*(obj.a'.*ones(1,Nexx))),-syx.SM.*(ones(nx,1)*(obj.b'.*ones(1,Nexx)));
                      zeros(2*obj.Nex,nx),obj.lsy.B];
                SMz = [syx.SM,zeros(nx,obj.Nex);zeros(obj.lsy.Ndof,obj.Nex),obj.lsy.SM];
                Sz = [zeros(obj.Nex,2*obj.Nex);zeros(obj.Nex),obj.lsy.SFF];
                mufexz = [syx.fex.mux;zeros(obj.Nex,1)];
                fnlz = syx.fnl.expand(nz);
                fexz = GaussianProcess(mufexz,Sz);
                dim_infoz = [syx.dim_info;3*ones(2*obj.Nex,1)];
                syz = NonLinearSystem_1O(Az,Bz,fnlz,SMz,fexz,dim_infoz);
            else
                error('unknown system type');
            end
        end
    end
end