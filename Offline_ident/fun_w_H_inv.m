%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for the inverse of w_H                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w_H_inv] = fun_w_H_inv(w_H)
    
    w_H_inv = ones(length(w_H), 1);

    w_H_inv(1) = 1/w_H(1);        
    for i=2:length(w_H)
        sum1 = 0;
        sum2 = 0;
        zaehler = w_H(i);
        
        for j=2:i
            sum1 = sum1 + w_H(j);
        end
        nenner1 = (w_H(1) + sum1);
        
        for j=2:i-1
            sum2 = sum2 + w_H(j);
        end
        nenner2 = (w_H(1) + sum2);
        
        w_H_inv(i,1) = -zaehler/(nenner1 * nenner2);
    end
end