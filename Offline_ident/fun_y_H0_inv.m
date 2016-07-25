%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for the inverse of y_H0                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y_H0_inv] = fun_y_H0_inv(w_H, y_H0)
    y_H0_inv = zeros(length(y_H0), 1);

    for i=1:length(w_H)
        sum1 = 0;
        sum2 = 0;
        for j=1:i
           sum1 = sum1 + w_H(j)*y_H0(i);
        end

        for j=(i+1):length(w_H)
            sum2 = sum2 + w_H(j)*y_H0(j);
        end
        y_H0_inv(i) = sum1+sum2;
    end
end