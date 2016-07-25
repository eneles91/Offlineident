%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for the inverse of r_H                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r_H_inv] = fun_r_H_inv(r_H, w_H)
    r_H_inv = ones(length(r_H),1);

    for i=1:length(r_H)
        sum1 = 0;
        for j=1:i
            sum1 = sum1 + w_H(j)*(r_H(i)-r_H(j));
        end
        r_H_inv(i) = sum1;
    end
end