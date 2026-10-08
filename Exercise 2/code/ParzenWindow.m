function parzenwindow = ParzenWindow(shiftings)
    %This function creates a parzen window taking into accound the number
    %of the maximum shifitngs.
    
    L3 = (length(shiftings)-1)/2;
    
    parzenwindow = NaN([length(shiftings) length(shiftings)]);
    
    for i=1:length(shiftings)
        for j=1:length(shiftings)
            parzenwindow(i,j) = dtaf(shiftings(i),L3)*dtaf(shiftings(j),L3)*dtaf(shiftings(j)-shiftings(i),L3);
        end
    end

    function d_taf = dtaf(lag,L3)
        if abs(lag)<= L3/2
            d_taf = 1-6*(abs(lag)/L3)^2+6*(abs(lag)/L3)^3;
        elseif abs(lag)>=L3/2 && abs(lag)<=L3
            d_taf = 2*(1-abs(lag)/L3)^3;
        elseif abs(lag)>L3
            d_taf = 0;
        end
    end
end