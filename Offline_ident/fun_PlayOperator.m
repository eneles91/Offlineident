%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for Playoperator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [yk] = fun_PlayOperator(xk , yk_1, r_H)

yk = max ([(xk - r_H), min([(xk + r_H), yk_1])]);

end