function rectwindow = RectangularWindow(shiftings,height)
    %This function creates a rectangular window taking into account the
    %number of maximum shiftings. The 'height' parameter defines the height
    %of the window. If not defined takes automatically the value one.
    
    if nargin == 1
        height = 1;
    end
    
    rectwindow = height*ones([length(shiftings) length(shiftings)]);
end