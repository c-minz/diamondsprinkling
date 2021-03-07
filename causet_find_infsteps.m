function infevents = causet_find_infsteps( L, steps )
%CAUSET_FIND_INFSTEPS returns the indices of the k-step future/past infinity by
% applying a find to the function CAUSET_SELECT_INF.
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.
    
    infevents = find( ~isnan( causet_select_infsteps( L, steps ) ) );
end

