function chains = causet_get_chains( C, degree )
%CAUSET_GET_CHAINS the numbers of total ordered k-chains in the 
% (sub)-causet with causal matrix C up to k = DEGREE.
% 
% Arguments:
% C                   upper triangular (logical) causals matrix.
%
% Optional arguments:
% DEGREE              maximal chain to be returned (length of return
%                     vector). (Default: 2)
% 
% Returns:
% CHAINS              row vector of length DEGREE with counts of total
%                     ordered k tuples (k-chains).
    
    %% set default values:
    if nargin < 2 || degree < 1
        degree = 2;
    end
    %% compute link matrix powers and chain counts:
    chains = zeros( 1, degree );
    chains( 1 ) = length( C );
    C = double( C );
    Cpower = C;
    for i = 2 : ( degree - 1 )
        chains( i ) = sum( sum( Cpower > 0, 1 ), 2 );
        Cpower = Cpower * C;
    end
    chains( degree ) = sum( sum( Cpower > 0, 1 ), 2 );
end

