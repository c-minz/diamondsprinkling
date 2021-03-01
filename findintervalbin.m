function [ bins, bincount, min ] = findintervalbin( values, max, ...
    bincount, min )
% FINDINTERVALBIN returns the bin indices for the vector VALUES such that
% the lowest bin index is the underflow and the highest bin index is the
% overflow. 
% The range from MIN (excluded) to MAX (included) is split into BINCOUNT
% many equivalent bins. Underflows (including exact values of MIN) are put 
% in the first bin (underflow bin), overflows in the last bin (overflow 
% bin).
% The first inputs are mandatory. BINCOUNT is by default 102, MIN is by
% default 0.
    
    if nargin < 3
        bincount = 102;
    end
    if nargin < 4
        min = 0;
    end
    underflow = 1;
    overflow = bincount;
    bins = ceil( ( bincount - 2 ) * double( values - min ) ...
                                  / ( max - min ) ) + 1;
    bins( bins < 1 ) = underflow;
    bins( bins > overflow ) = overflow;
end

