function preftimes = causet_get_preftimes( L )
%CAUSET_GET_PREFTIMES calculates the matrix of preferred future/past from the 
% causet link matrix. 
% 
% Arguments: 
% L                   logical link matrix.
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.

    preftimes = L^2 >= 2; % preferred future/past in two link distance
end
