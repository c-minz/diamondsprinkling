function C = causet_edit_relate( coordinates, metric )
%CAUSET_EDIT_RELATE is a quick version of CAUSET_EDIT_LINK and only 
% returns the causal matrix.
% 
% Arguments:
% COORDINATES         positions of the elements.
% METRIC              matrix for measurement.
% 
% Returns:
% C                   upper triangular (logical) causal matrix.
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.
    
    N = size( coordinates, 1 );
    C = false( N );
    for j = 2 : N
        for i = 1 : ( j - 1 )
            dcoordinates = coordinates( j, : ) - coordinates( i, : );
            causaldistanceIJ = dcoordinates * metric * transpose( dcoordinates );
            if causaldistanceIJ >= 0
                C( i, j ) = true;
            end
        end
    end
end

