%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for Creepoperator %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [yk] = fun_KriechOperator(xk_1, yk_1, r_K, a_K, T_s)

yk = yk_1 + (1 - exp(-a_K*T_s)) * fun_PlayOperator((xk_1 - yk_1),0, r_K);
%yk = yk_1 + (1 - exp(-a_K*T_s)) * max([(xk_1-yk_1-r_K), min([(xk_1 - yk_1 + r_K), 0])]);
end