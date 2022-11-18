classdef NL_Element
    %% Super class of nonlinear elements
    
    properties
        ne; % number of equations
        nv; % number of variables
    end
    
    methods
        function obj = NL_Element(ne,nv)
            obj.ne = ne;
            obj.nv = nv;
        end
    end
end

