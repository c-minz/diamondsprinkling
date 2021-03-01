function preftimes = causet_get_preftimes( L )
%CAUSET_GET_PREFTIMES calculates the matrix of preferred future/past from the 
% causet link matrix. 
% 
% Arguments: 
% L                   logical link matrix.

    preftimes = L^2 >= 2; % preferred future/past in two link distance
end
