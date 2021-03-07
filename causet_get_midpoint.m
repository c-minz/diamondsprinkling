function [ midpoints, mincard ] = causet_get_midpoint( C )
%CAUSET_GET_MIDPOINT finds the central-most event(s) of the causet.
% 
% Arguments:
% C                   logical upper triangular causal matrix.
%
% Returns:
% MIDPOINTS           row vector of events, which maximize the product of
%                     the event counts in their past and future lightcone.
% MINCARD             minimum of the past or future volumes for each event.
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.
    
    pcard = sum( C, 1 );
    fcard = transpose( sum( C, 2 ) );
    [ m, midpoints ] = max( pcard .* fcard ); %#ok<ASGLU>
    mincard = min( pcard, fcard );
    mincard = mincard( midpoints );
end

