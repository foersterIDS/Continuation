classdef RnRotation < handle
    
    properties
        vx
        n
        TR
    end
    
    methods
        function obj = RnRotation(vx)
            vx = vx(:);
            n = length(vx);
            obj.vx = vx(:);
            obj.n = n;
            ex1 = [1;zeros(n-1,1)];
            ex2 = [0;1;zeros(n-2,1)];
            alpha = acos((vx'*ex1)/(sqrt(vx'*vx)*sqrt(ex1'*ex1)));
            
            g1 = vx/sqrt(vx'*vx);
            
            if alpha==0
                g2 = ex2;
            else
                a1 = vx(1)/(vx'*vx);
                g2 = a1*vx-ex1;
                g2 = g2/sqrt(g2'*g2);
            end
            
            V = g1*g1'+g2*g2';
            W = g1*g2'-g2*g1';
            I = eye(n);
            
            obj.TR = I+(cos(-alpha)-1)*V-sin(-alpha)*W;
        end
        
        function TR = getTR(obj,vx)
            if length(vx)==obj.n && prod(vx==obj.vx)
                TR = obj.TR;
            else
                vx = vx(:);
                n = length(vx);
                obj.vx = vx(:);
                obj.n = n;
                ex1 = [1;zeros(n-1,1)];
                ex2 = [0;1;zeros(n-2,1)];
                alpha = acos((vx'*ex1)/(sqrt(vx'*vx)*sqrt(ex1'*ex1)));
                
                g1 = vx/sqrt(vx'*vx);
                
                if alpha==0
                    g2 = ex2;
                else
                    a1 = vx(1)/(vx'*vx);
                    g2 = a1*vx-ex1;
                    g2 = g2/sqrt(g2'*g2);
                end
                
                V = g1*g1'+g2*g2';
                W = g1*g2'-g2*g1';
                I = eye(n);
                
                TR = I+(cos(-alpha)-1)*V-sin(-alpha)*W;
                obj.TR = TR;
            end
        end
    end
end

