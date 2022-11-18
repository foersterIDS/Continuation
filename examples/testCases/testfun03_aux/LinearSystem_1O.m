classdef LinearSystem_1O < handle
    
    properties
        %
        % A*dxdt+B*x=SM*F
        %
        name = 'LS1O';
        A % Einheitsmassenmatrix
        invA % Inverse der Einheitsmassenmatrix
        B % Zustandsraummatrix
        SM % Sortiermatrix
        SFF % Spektraldichtematrix der Anregung (Funktion von Om)
        KXX % Kovarianzmatrix
        KdXdX % Kovarianzmatrix der 1. Ableitung
        RXX % Korrelationskoeffizientenmatrix
        Ndof % Anzahl der Dimensionen
        Nex % Anzahl der verschiedenen Anregungen
        EW % Eigenwertmatrix
        ew % Eigenwertvektor
        EVr % rechte Eigenvektormatrix
        invEVr % Inverse der rechten Eigenvektormatrix
        EVl % linke Eigenvektormatrix
        invEVrinvASM
        isGWN
        G
        D
    end
    
    methods
        function obj = LinearSystem_1O(A,B,SM,SFF,invA)
            % A*dxdt+B*x=SM*F
            % SM: Sortiermatrix
            % SFF: Spektraldichte (Anregung)
            obj.SM = SM;
            obj.SFF = SFF;
            [obj.Ndof,obj.Nex] = size(SM);
            obj.A = A;
            if nargin<=4
                obj.invA = A\eye(size(A));
            else
                obj.invA = invA;
            end
            obj.B = B;
            if isa(SFF,'function_handle')
                obj.isGWN = false;
            else
                obj.isGWN = true;
                obj.G = -obj.invA*obj.B;
                obj.D = conj(obj.invA*obj.SM)*obj.SFF*(obj.invA*obj.SM).';
            end
        end
        
        function updateB(obj,Bnew)
            if prod(size(Bnew)==(obj.Ndof*[1,1]))
                obj.B = Bnew;
                if obj.isGWN
                    obj.G = -obj.invA*obj.B;
                end
                obj.EW = [];
                obj.ew = [];
                obj.EVr = [];
                obj.EVl = [];
                obj.invEVr = [];
                obj.invEVrinvASM = [];
            else
                error('Bnew has the wrong size');
            end
        end
        
        function calcEP(obj)
            if isempty(obj.EW)
                [obj.EVr,obj.EW,obj.EVl] = eig(obj.invA*obj.B);
                obj.invEVr = obj.EVr\eye(size(obj.EVr));
                obj.ew = zeros(size(obj.A(:,1)));
                for h=1:obj.Ndof
                    obj.ew(h) = obj.EW(h,h);
                end
                if sum(real(obj.ew)>0)<obj.Ndof
                    warning('possible rigid body motion (root with negative real part)');
                end
                obj.invEVrinvASM = obj.invEVr*obj.invA*obj.SM;
            end
        end
        
        function H = getH(obj,omega)
            obj.calcEP();
            % Übertragungsfunktion
            o = length(omega);
            H = zeros(obj.Ndof,obj.Nex,o);
            for k=1:o
%                 Hk = ((1i*omega(k)*obj.A+obj.B)\eye(obj.Ndof))*obj.SM; % normal
                Hk = obj.EVr*diag(1./(1i*omega(k)+obj.ew))*obj.invEVrinvASM; % modal (noch schneller)
                H(:,:,k) = Hk;
            end
        end
        
        function SXX = getSXX(obj,omega)
            obj.calcEP();
            % Spektraldichtematrix der Systemantwort
            o = length(omega);
            SXX = zeros(obj.Ndof,obj.Ndof,o);
            for k=1:o
%                 Hk = ((1i*omega(k)*obj.A+obj.B)\eye(obj.Ndof))*obj.SM; % normal
                Hk = obj.EVr*diag(1./(1i*omega(k)+obj.ew))*obj.invEVrinvASM; % modal (noch schneller)
                if obj.isGWN
                    SXX(:,:,k) = conj(Hk)*obj.SFF*Hk.';
                else
                    SXX(:,:,k) = conj(Hk)*obj.SFF(omega(k))*Hk.';
                end
            end
        end
        
        function KXX = getKXX(obj)
            % Kovarianzmatrix der Systemantwort
            if obj.isGWN
                KXX = lyap(obj.G,obj.D);
            else
                KXX = (1/(2*pi))*integral(@(Om) obj.getSXX(Om),-inf,inf,'ArrayValued',true);
            end
            errimag = sqrt(sum(imag(KXX(:)).^2));
%             if errimag>10^-10
%                 warning('Imaginaeranteil in KXX!');
%             end
%             KXX = real(KXX);
        end
        
        function KdXdX = getKdXdX(obj)
            % Kovarianzmatrix der Systemantwort
            KdXdX = ((1/(2*pi))*integral(@(Om) Om.^2.*obj.getSXX(Om),-inf,inf,'ArrayValued',true));
            errimag = sqrt(sum(imag(KdXdX(:)).^2));
            if errimag>10^-10
                warning('Imaginaeranteil in KdXdX!');
            end
%             KdXdX = real(KdXdX);
        end
        
        function RXX = getRXX(obj,varargin)
            % Korrelationskoeffizientenmatrix der Systemantwort
            if nargin==1
                KXX = obj.getKXX();
            else
                KXX = varargin{1};
            end
            RXX = getCorrCoef(KXX);
            RXX = real(RXX);
        end
        
        function EW = getEW(obj)
            obj.calcEP();
            EW = obj.EW;
        end
        
        function [EVr,EVl] = getEV(obj)
            obj.calcEP();
            EVr = obj.EVr;
            EVl = obj.EVl;
        end
    end
end