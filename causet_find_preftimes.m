function preftime = causet_find_preftimes( perimetrals, internals, criterion, minprm )
%CAUSET_FIND_PREFTIMES gives the indices of the preferred time
% (past/future) under a CRITERION (Default: 1) from the vector of 
% PURITYCOUNTS with respective IMPURITYCOUNTS.
% 
% Arguments:
% PERIMETRALS         column vector of rank 2 path counts.
% INTERNALS           column vector of number of elements in paths with 
%                     rank > 2.
% 
% Optional arguments:
% CRITERION           one of the following integers (Default: 2):
%  1:                 largest diamonds.
%  2:                 smallest diamonds.
%  3:                 largest pure diamonds.
%  4:                 diamonds with firstly most 2-step paths, and 
%                     secondly least impurities.
%  5:                 diamonds with firstly most 2-step paths, and 
%                     secondly most impurities.
%  6:                 searches for all singletons and chooses the one with
%                     firstly most 2-step paths, and secondly least
%                     impurities. If there is no singelton set, it returns
%                     criterion 5.
%  7:                 diamond improving the second criterion as follows. 
%                     Let i be the number of impurities and p the number 
%                     of pure 2-step paths in the diamond according to 
%                     criterion 2. This criterion allows for at most 
%                     max(0,i-k) impurities in an (i+p-k)-diamond. For 
%                     m \in [1,i+p], it returns the smallest m-diamond 
%                     with this internals restriction and if there is only 
%                     one with size m. If it does not find a unique 
%                     diamond smaller than the (i+p)-diamond, this 
%                     criterion is identical to criterion two.
% MINPRM              minimal number of rank 2 paths. (Default: 1).
% 
% Returns:
% PREFTIME            vector of indexes of preferred time (past/future). 
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.
    
    %% set default values:
    if nargin < 3 || criterion < 1 || criterion > 6
        criterion = 6;
    end
    if nargin < 4 || minprm < 1
        minprm = 1;
    end
    %% pick indexes of preferred time by criterion:
    %  include +-Inf in min/max functions to avoid empty returns from empty
    %  vectors
    if criterion == 1
        dsizes = perimetrals + internals;
        preftime = find( perimetrals >= minprm & ...
            dsizes == max( [ dsizes( perimetrals >= minprm ), -Inf ] ) );
        return
    elseif criterion == 2
        dsizes = perimetrals + internals;
        preftime = find( perimetrals >= minprm & ...
            dsizes == min( [ dsizes( perimetrals >= minprm ), Inf ] ) );
        return
    elseif criterion == 3
        sel_onlypure = ( internals == 0 ) & ( perimetrals >= minprm );
        preftime = find( ...
            sel_onlypure ...
          & perimetrals == max( [ perimetrals( sel_onlypure ), -Inf ] ) );
        return
    elseif criterion == 4
        sel_maxprm = ( perimetrals == max( perimetrals ) ) ...
            & ( perimetrals >= minprm );
        preftime = find( ...
            sel_maxprm ...
          & internals == max( [ internals( sel_maxprm ), -Inf ] ) );
      return
    elseif criterion == 6
        sel_prm = find( perimetrals >= minprm );
        [ i_all, p_all ] = find( sparse( internals( sel_prm ) + 1, ...
            perimetrals( sel_prm ), ones( 1, length( sel_prm ) ) ) == 1 );
        i = min( i_all );
        if ~isempty( i )
            p = max( p_all( i_all == i ) );
            preftime = find( ( perimetrals == p ) & ...
                ( internals == i - 1 ) );
            return
        end
    elseif criterion == 7
        % find internals maximum, which shall decrease with decresing
        % rank 2 path count:
        maxprm = max( perimetrals );
        maxitn = min( internals( perimetrals == maxprm ) );
        if isempty( maxitn )
            maxitn = 0;
        end
        % decrease internals maximum with rank 2 path count:
        preftime = [];
        for i = minprm : 1 : maxprm
            sel_prm = ( perimetrals == i );
            thismaxitn = max( 0, maxitn - maxprm + i );
            preftime = find( ...
                sel_prm ...
              & internals == min( [ internals( sel_prm ), thismaxitn ] ) );
            if length( preftime ) == 1
                break
            end
        end
        return
    end
    % criterion == 5 or (criterion == 6 and no singeltons)
    sel_maxprm = ( perimetrals == max( perimetrals ) ) ...
        & ( perimetrals >= minprm );
    preftime = find( ...
        sel_maxprm ...
      & internals == min( [ internals( sel_maxprm ), Inf ] ) );
end

