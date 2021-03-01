function g = metric( d, spacetime )
%METRIC returns the metric for SPACETIME (Default: Minkowski) in D
% dimensions.
    
    if nargin < 2
        spacetime = 'Minkowski';
    end
    if strcmp( spacetime, 'Minkowski' )
        g = metric_flat( d );
    else % not supported spacetime
        g = eye( d );
    end
end

