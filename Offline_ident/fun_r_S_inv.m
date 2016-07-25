function [r_S_invf] = fun_r_S_inv(r_S, w_S, L)
    
    r_S_invf = ones(2*L+1, 1);
    
    %Berechnung von r_S0_inv
    r_S_invf(L+1) = 0;
    
    % analog zu: i= -L ... -1
    for i=1:L
        sum1 = 0;
        parfor j=i:(L+1)
            sum1 = sum1 + w_S(j) * (r_S(i) - r_S(j));
        end
        r_S_invf(i,1) =r_S(i)+ sum1;
    end
    
    % analog zu: i= 1 ... L
    for i=(L+2):((2*L)+1)
        sum1 = 0;
        parfor j=(L+1):i
            sum1 = sum1 + w_S(j) * (r_S(i) - r_S(j));
        end
        r_S_invf(i,1) =r_S(i)+ sum1;
    end
    

end