%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for Playoperator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [yk] = fun_SuperpositionsOperator(xk , r_S)

    if (r_S > 0)
        yk = max([xk - r_S, 0]);
    elseif (r_S < 0)
        yk = min([xk - r_S, 0]);
    else
        yk = xk;
    end
    
end