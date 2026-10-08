function impulseResponse = GiannakisFormula(c3x,q)
    %This function calculates the impulse response of the system based on
    %the Giannakis formula using the order q of the MA process that the
    %signal is given by and the third order cumulants of the same signal.
    
    L3 = (length(c3x)-1)/2;
    
    impulseResponse = zeros(1,q+1);
    for i=0:q
        impulseResponse(i+1) = c3x(L3+1+q,L3+1+i)/c3x(L3+1+q,L3+1);
    end
end
        
        