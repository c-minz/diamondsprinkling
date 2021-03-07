function intervalincrement = convert2interval( values, max, ...
    bincount, min )
% CONVERT2INTERVAL converts the vector VALUES into a length of VALUES 
% times ( INTERVALS + 2 ) matrix. 
% The range from MIN (excluded) to MAX (included) is split into BINCOUNT
% many equivalent bins. Underflows (including exact values of MIN) are put 
% in the first bin (underflow bin), overflows in the last bin (overflow 
% bin).
% The first inputs are mandatory. Default for BINCOUNT and MIN are set by 
% FINDINTERVALBIN.
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.
    
    if nargin < 3
        [ bins, bincount ] = findintervalbin( values, max );
    elseif nargin < 4
        bins = findintervalbin( values, max, bincount );
    else
        bins = findintervalbin( values, max, bincount, min );
    end
    intervalincrement = zeros( length( values ), bincount );
    for i = 1 : length( values )
        intervalincrement( i, bins( i ) ) = 1;
    end
end

