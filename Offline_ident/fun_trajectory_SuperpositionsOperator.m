%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for Trajectory of SuperpositionsOperator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = fun_trajectory_SuperpositionsOperator(x, r_S)

% Trajektorie des Superpositionsoperators

y = [length(x),1];

    for k=1:length(x)
            y(k) = fun_SuperpositionsOperator(x(k), r_S);
    end

end