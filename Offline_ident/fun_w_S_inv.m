%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for the inverse of w_S                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w_S_inv] = fun_w_S_inv(w_S, L)

    w_S_inv = ones(length(w_S), 1);
    
    w_S0 = w_S(L+1);
    
    %Berechnung von w_S0_inv
    w_S_inv(L+1) = 1/w_S0; 
    
    % analog zu: i= -L ... -1
    for i=1:L
        sum1 = 0;
        sum2 = 0;
        zaehler = w_S(i);
        
        parfor j=i:L
            sum1 = sum1 + w_S(j);
        end
        nenner1 = (w_S0 + sum1);
        
        parfor j=(i+1):L
            sum2 = sum2 + w_S(j);
        end
        nenner2 = (w_S0 + sum2);
        
        w_S_inv(i,1) = -zaehler/(nenner1 * nenner2);
    end
    
    % analog zu: i= 0 ... L
    for i=(L+2):((2*L)+1)
        sum1 = 0;
        sum2 = 0;
        zaehler = w_S(i);
        
        parfor j=(L+2):i
            sum1 = sum1 + w_S(j);
        end
        nenner1 = (w_S0 + sum1);
        
        parfor j=(L+2):(i-1)
            sum2 = sum2 + w_S(j);
        end
        nenner2 = (w_S0 + sum2);
        
        w_S_inv(i,1) = -zaehler/(nenner1 * nenner2);
    end
end