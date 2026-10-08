function shifted_array = arrayShift(array,tau)
    %This function implements the shiftings of an array 'tau' possitions to
    %the left or to the right regarding the polarity of the value (if it is
    %possitive or negative).
    
    if tau == 0
        shifted_array = array;
    elseif tau>0
        shifted_array = circshift(array,-abs(tau));
        shifted_array(length(array)-abs(tau):length(array)) = 0;
    else
        shifted_array = circshift(array,abs(tau));
        shifted_array(1:abs(tau)) = 0;
    end
end