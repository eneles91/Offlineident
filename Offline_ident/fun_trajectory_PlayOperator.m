%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for Trajectory of Playoperator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y] = fun_trajectory_PlayOperator(x, r_H, y0)

% Trajektorie des Playoperators

y = [length(x),1];

    for k=1:length(x)
        if (k==1)
            y(k) = fun_PlayOperator(x(k), y0, r_H);
        else
            y(k) = fun_PlayOperator(x(k), y(k-1), r_H); 
        end    
    end

end