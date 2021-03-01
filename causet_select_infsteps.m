function infsteps = causet_select_infsteps( L, steps )
%CAUSET_SELECT_INFSTEPS returns a N-vector holding the number of steps of 
% the future/past infinity up to STEPS for each element in the causet. 
% Unidentified elements are returned as NaN.
% Use Inf and -Inf to address 0 step future/past infinity and values 
% STEPS > 0 for future, STEPS < 0 for past.
    
    N = size( L, 1 );
    infsteps = NaN * zeros( 1, N ); % pre-allocate memory
    if steps >= 0 % find future infinity steps:
        direction = 1;
        infsteps( sum( L, 2 ) == 0 ) = 0;
        linkcounts = sum( L, 2 );
    else % find past infinity steps:
        direction = -1;
        steps = -steps;
        infsteps( sum( L, 1 ) == 0 ) = 0;
        linkcounts = sum( L, 1 );
    end
    if isinf( steps )
        steps = 0; % +/-Inf represents +/-0
    end
    for s = 1 : steps
        select = ~isnan( infsteps );
        if direction > 0
            isinfstep = transpose( sum( L( :, select ), 2 ) == linkcounts ) ...
                & ~select;
        else
            isinfstep = ( sum( L( select, : ), 1 ) == linkcounts ) ...
                & ~select;
        end
        infsteps( isinfstep ) = s;
    end
end

