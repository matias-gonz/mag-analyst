function MH = mvshanhpos(Hanh,Mq,Moffset,Hoffset)
%MVSHANPOS Local function that calculates the Manh[Hanh] distribution for positive Hanh and Manh (i.e., 1st quadrant).
%The function assumes the loop is centered, so it only consideres the positive H interval.
    
    % Hanh is the array of Hanh
    % Mq is the array of Manh(Hanh)
    % Moffset and Hoffset are the vertical and horizontal offsets, respectively, to add to the hysteresis loop so that M(H)=-M(-H)
    
    MqShifted = (Mq + Moffset);
    HanhShifted = @(Hoffset)(Hanh + Hoffset);
    
    tmpHanhShifted = HanhShifted(Hoffset);   % temporary variable to store the function of shifted Hanh at Hoffset input. 
                                             % Although Octave does handle a syntaxis such as HanhShifted(Hoffset)(1:end_Hanhneg_shifted)
                                             % -having two sets of parentheses immediately following one another-, Matlab doesn't.

    % First, find where we need to "fold" Manh(H) distribution (ideally, at Hanh=0, but that value may not be present in the Hanh(M) table)
    % We'll do it both for the original Hanh(M) distribution, as for the shifted distribution, just to compare them both in a plot
    % (we actually only need it for the shifted distribution)

    [Hanhzero_shifted, Izero_shifted] = min(abs(tmpHanhShifted));   % Find |Hanh_shifted| minimum (i.e., where Hanh is zero o changes sign) 
                                                                    % as well as the row index of data in which it appears
    
    if tmpHanhShifted(Izero_shifted) == 0;
	    start_Hanhpos_shifted = Izero_shifted;
    elseif tmpHanhShifted(Izero_shifted) > 0;
	    start_Hanhpos_shifted = Izero_shifted;
    else tmpHanhShifted(Izero_shifted) < 0;
	    start_Hanhpos_shifted = Izero_shifted + 1;
    end

    absHanhpos_shifted = tmpHanhShifted(start_Hanhpos_shifted:end);
    absManhpos_shifted = abs(MqShifted(start_Hanhpos_shifted:end));         % array of M(Hanh) with Hanhshifted>0

    % We need to remove duplicates using 'unique' function. C = unique(A) returns the same data as in A, but with no repetitions. C is in sorted order.
    [ManhUnique, index] = unique(absManhpos_shifted);
    HanhUnique = absHanhpos_shifted(index);
    
    MH = [transpose(HanhUnique),transpose(ManhUnique)];
    if HanhUnique(1) > 0
        MH = [[0,0]; MH];
    end
end