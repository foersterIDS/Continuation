classdef RandomProcess < handle
    % Superclass of random processes
    %   GaussianProcess, FilteredProcess
    
    properties
        mux
        n
        name
    end
    
    methods
        function obj = RandomProcess(mu,name)
            obj.mux = mu; % Mittelwert
            obj.n = length(mu);
            obj.name = name;
        end
        
        function mux = getMean(obj)
            mux = obj.mux;
        end
        
        function n = length(obj)
            n = obj.n;
        end
    end
end