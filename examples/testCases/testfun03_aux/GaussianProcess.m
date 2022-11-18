classdef GaussianProcess < RandomProcess
    
    properties
        Sxx
    end
    
    methods
        function obj = GaussianProcess(mu,S0)
            name = 'gp';
            obj = obj@RandomProcess(mu,name);
            obj.Sxx = S0; % Spektrale Leistungsdichte des Mittelwertfreien Prozesses
        end
        
        function Sxx = getPSD(obj)
            Sxx = obj.Sxx;
        end
    end
end