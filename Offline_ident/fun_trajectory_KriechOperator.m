%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for Trajectory of KriechOperator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = fun_trajectory_KriechOperator(x, y_K0, a_K, r_K, T_s)

% Trajektorie des Kriechoperators

y = [length(x),1];

    for k=1:length(x)
            if (k == 1)
                %STARTWERTE VON X=??????
                y(k) = fun_KriechOperator(x(k), y_K0, r_K, a_K, T_s);
            else  
                y(k) = fun_KriechOperator(x(k-1), y(k-1), r_K, a_K, T_s);
            end
    end
end


