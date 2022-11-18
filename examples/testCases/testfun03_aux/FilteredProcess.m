classdef FilteredProcess < RandomProcess
    
    properties
        fi
    end
    
    methods
        function obj = FilteredProcess(mu,fi)
            name = 'fp';
            obj = obj@RandomProcess(mu,name);
            obj.fi = fi; % Filter
        end
        
        function fi = getFilter(obj)
            fi = obj.fi;
        end
        
        function nls = applyExcitation(obj,nls)
            nlsh = nls;
            nls = obj.fi.merge(nls);
            nls.nlsh = nlsh;
        end
    end
end