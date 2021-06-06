function results = causet_get_statistics( C, L, maxsizes, ...
    flags, events, criterias, coordinates, g )
%CAUSET_GET_STATISTICS counts the link relations in the link matrix L 
% between all the element pairs, which have a diamond arrangement, starting 
% with 1 intermediate element and going up to MAXSIZES. 
% 
% Arguments:
% C                   logical upper triangular causal matrix.
% L                   logical upper triangular link matrix.
% MAXSIZES            3 element vector: 
%  1: MAXPURITIES     positive number for the maximal number of
%                     center elements in a diamond.
%  2: MAXIMPURITY     positive number for the maximal number of
%                     impurities in a diamond. Set to 0 if diamond-link
%                     impurity should not matter for the counting.
%  3: MAXPREFFUTURES  positive integer for the maximal number of
%                     preferred futures in the statistics per event.
%  4: MAXPROPERTIME   maximal proper time. Only necessary if coordinates
%                     are specified.
%  5: PT_BINCOUNT     number of bins for the proper time statistics, 
%                     excluding overflow bins. Only necessary if 
%                     coordinates are specified.
%  6: MAXHYPERBCOORD  maximal hyperboloic coordinate value. Only necessary 
%                     if coordinates are specified.
%  7: HB_BINCOUNT     number of bins for the unit hyperboloid statistics, 
%                     excluding overflow bins. Only necessary if 
%                     coordinates are specified.
%  8: HBSPEED_BINCOUNT number of bins for the unit hyperboloid speed 
%                     statistics, excluding overflow bins. Only necessary 
%                     if coordinates are specified.
% 
% Optional arguments:
% FLAGS               string of return fields and flags separated by 
%                     spaces. The sizes of the returns contain an
%                     overflow row (+ 1).
% -allfields          include all of the following fields (Default).
% -fields             followed by the list of result fields to be included:
% .diamonds           diamond links 
%                     int64 matrix, [ MAXIMPURITY + 1, MAXPURITY + 1 ]
% .simplices          counters of sub-causets that are resemble simplices. 
%                     [ 1+0, 1+1, 1+2, 1+3 ] dim. counter 
%                     (See also CAUSET_GET_SIMPLICES) 
%                     int64 matrix, [ 4, 1 ]
% .chains             k-chain counts for k in [ 0, 4 ] (the zeroth entry
%                     holds the average counter) for spacetime dimension 
%                     estimators.
%                     int64 matrix, [ 5, 1 ]
% .diamondtimes       sum of all diamond link proper time separations
%                     matrix, [ MAXIMPURITY + 1, MAXPURITY + 1 ]
% .propertimes        proper time separation distribution
%                     int64 column vector, [ PT_BINCOUNT, 1 ]
% .hyperbdistribution radial unit hyperboloid distribution
%                     int64 column vector, [ HB_BINCOUNT, 1 ]
% .hyperbspeeddistribution radial unit hyperboloid speed distribution
%                     int64 column vector, [ HBSPEED_BINCOUNT, 1 ]
% .preffutures        analysis of preferred future criteria
%                     structure with the automatically included fields:
% (..counts)          counts of elements with preferred future.
%                     (Automatically included with .preffutures)
%                     int64 tensor, [ MAXPREFFUTURES + 1, 
%                                     criterias( 1 ), criterias( 2 ) ]
% ..diamonds          diamond links to preferred future
%                     int64 tensor, [ MAXIMPURITY + 1, MAXPURITY + 1, 
%                                     criterias( 1 ), criterias( 2 ) ]
% ..chains            k-chain counts for k in [ 0, 4 ] - the zeroth entry
%                     holds the average counter
%                     (for spacetime estimators)
%                     matrix, [ 5, criterias( 1 ), criterias( 2 ) ]
% ..propertimes       proper time separation distribution to preferred
%                     future
%                     int64 row vector, [ PT_BINCOUNT, 
%                                         criterias( 1 ), criterias( 2 ) ]
% ..unithyperboloid   projected coordinates on the unit hyperboloid for 
%                     each preferred future, determined by the coordinates 
%                     of the preferred future relative to the event. 
%                     cell array, [ criterias( 1 ), criterias( 2 ) ]
%                       each element: matrix, 
%                         [ total number of preferred futures, dim - 1 ]
% ..hyperbdistribution radial unit hyperboloid distributions
%                     int64 tensor, [ HB_BINCOUNT, 
%                                     criterias( 1 ), criterias( 2 ) ]
% ..hyperbspeeddistribution radial unit hyperboloid speed distributions
%                     int64 tensor, [ HBSPEED_BINCOUNT, 
%                                     criterias( 1 ), criterias( 2 ) ]
% -set                analyse all diamond links (Default).
% -chain              only analyse diamond links which contain three 
%                     totally ordered events from the list EVENTS.
% EVENTS              vector of causet events to be considered in the
%                     counting. Default: 1 : N (all elements).
%                     Here it is possible to specify the list of an event 
%                     chain with the flag '-chain' to focus on diamond 
%                     links along the chain only.
% CRITERIAS           two element vector. Use [] to exclude preferred 
%                     future counting.
%                     1. criteria count (Default: 5)
%                     2. maximum of the minimal 2-step count per diamond
%                     (Default: 3).
% COORDINATES         positions of the elements. Use [] (default) to skip
%                     calculation of proper time separations.
% G                   metric for the computation of proper times. Has to be
%                     specified when COORDINATES is set. Use a function
%                     pointer for curved spacetimes. The function has to 
%                     compute the square of the geodesic length between a 
%                     pair of coordinates, arguments A and B, which are d
%                     dimensional row vectors.
% 
% Returns: 
% RESULTS             structure with the fields specified above.
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.
    
    %% set start up parameter:
    if nargin < 4 || isempty( flags )
        flags = '-allfields';
    end
    [ is, proc, get ] = causet_get_statistics_flags( flags, ...
        ( nargin >= 7 ) && ~isempty( coordinates ), ...
        ~isempty( criterias ) );
    maxpurity = maxsizes( 1 );
    maximpurity = maxsizes( 2 );
    maxpreffutures = maxsizes( 3 );
    if proc.coordinates
        dim = size( coordinates, 2 );
        maxpropertime = maxsizes( 4 );
        pt_bincount = maxsizes( 5 );
        maxhyperbcoord = maxsizes( 6 );
        hb_bincount = maxsizes( 7 );
        hbspeed_bincount = maxsizes( 8 );
    end
    N = size( L, 1 );
    if nargin < 5
        events = 1 : N;
    end
    if nargin < 6
        criterias = [];
    elseif ~isempty( criterias ) && length( criterias ) < 2
        criterias = [ criterias, 1 ];
    end
    %% allocate memory for results fields:
    results = struct();
    if get.diamonds
        results.diamonds = int64( zeros( maximpurity + 1, maxpurity + 1 ) );
    end
    if get.simplices
        results.simplices = int64( zeros( 4, 1 ) );
        C_t = transpose( C );
        L_t = transpose( L );
    end
    if get.diamondtimes
        results.diamondtimes = zeros( maximpurity + 1, maxpurity + 1 );
    end
    if get.chains
        results.chains = int64( zeros( 5, 1 ) );
    end
    if get.propertimes
        results.propertimes = int64( zeros( pt_bincount, 1 ) );
    end
    if get.hyperbdistribution
        results.hyperbdistribution = int64( zeros( hb_bincount, 1 ) );
        results.hyperbspeeddistribution = ...
            int64( zeros( hbspeed_bincount, 1 ) );
    end
    if get.preffutures
        criteria_count = criterias( 1 );
        criteria_minsize_max = criterias( 2 );
        results.preffutures.counts = ...
            int64( zeros( maxpreffutures + 1, ...
                          criteria_count, criteria_minsize_max ) );
        if get.preffutures_diamonds
            results.preffutures.diamonds = ...
                int64( zeros( maximpurity + 1, maxpurity + 1, ...
                              criteria_count, criteria_minsize_max ) );
        end
        if get.preffutures_chains
            results.preffutures.chains = ...
                int64( zeros( 5, criteria_count, criteria_minsize_max ) );
        end
        if get.preffutures_propertimes
            results.preffutures.propertimes = ...
                int64( zeros( pt_bincount, ...
                              criteria_count, criteria_minsize_max ) );
        end
        if get.preffutures_unithyperboloid
            results.preffutures.unithyperboloid = ...
                cell( criteria_count, criteria_minsize_max );
            results.preffutures.unithyperboloid_rowcount = ...
                int64( zeros( criteria_count, criteria_minsize_max ) );
        end
        if get.preffutures_hyperbdistribution
            results.preffutures.hyperbdistribution = ...
                int64( zeros( hb_bincount, ...
                              criteria_count, criteria_minsize_max ) );
            results.preffutures.hyperbspeeddistribution = ...
                int64( zeros( hbspeed_bincount, ...
                              criteria_count, criteria_minsize_max ) );
        end
    end
    %% statistics for selected events:
    if is.eventchain
        chain_sel = false( 1, N );
        chain_sel( events ) = true;
    end
    for i = events
        %% record simplex dimension:
        if get.simplices
            stdim = causet_get_simplices( C, L, { i } ) + 1;
            if stdim < 4
                stdim = max( stdim, ...
                    causet_get_simplices( C_t, L_t, { i } ) + 1 );
            end
            results.simplices( stdim ) = ...
                results.simplices( stdim ) + 1;
            if get.onlysimplices
                continue
            end
        end
        %% get diamond link properties:
        rank1_sel = L( i, : ); % links from bottoms to centers
        if is.eventchain % only single diamond top along chain
            rank2_sel = L( chain_sel & rank1_sel, : );
            puritycounts = chain_sel .* sum( rank2_sel, 1 );
        else
            rank2_sel = L( rank1_sel, : ); % links from centers to tops
            puritycounts = sum( rank2_sel, 1 ); % numbers of pure links
        end
        impuritycounts = zeros( 1, N ); % pre-allocate
        %% step through all diamond link top events:
        future_idx = find( puritycounts ); % tops
        for j = future_idx
            % causal from i-th AND causal to j-th:
            center_sel = C( i, 1:j ) & transpose( C( 1:j, j ) );
            % size of diamond link:
            puritycount = puritycounts( j );
            impuritycount = sum( center_sel ) - puritycount;
            impuritycounts( j ) = impuritycount;
            %% record diamond link type:
            if get.diamonds || get.diamondtimes
                % limit index range:
                if impuritycount > maximpurity % overflow
                    impurity_idx = maximpurity + 1;
                else
                    impurity_idx = impuritycount + 1;
                end
                if puritycount > maxpurity % overflow
                    purity_idx = maxpurity + 1;
                else
                    purity_idx = puritycount;
                end
            end
            if get.diamonds
                results.diamonds( impurity_idx, purity_idx ) = ...
                    results.diamonds( impurity_idx, purity_idx ) + 1;
            end
            %% record diamond link proper time seperation:
            if proc.coordinates
                if ismatrix( g )
                    dcoordinates = coordinates( j, : ) - coordinates( i, : );
                    propertime = sqrt( dcoordinates * g * transpose( dcoordinates ) );
                    unithyperbcoords = dcoordinates( 2 : dim ) / propertime;
                else
                    propertime = sqrt( feval( g, coordinates( i, : ), coordinates( j, : ) ) );
                    unithyperbcoords = zeros( 1, dim - 1 );
                end
            end
            if get.diamondtimes
                results.diamondtimes( impurity_idx, purity_idx ) = ...
                    results.diamondtimes( impurity_idx, purity_idx ) + propertime;
            end
            if get.propertimes
                idx = findintervalbin( propertime, maxpropertime, pt_bincount );
                results.propertimes( idx ) = ...
                    results.propertimes( idx ) + 1;
            end
            if get.hyperbdistribution
                [ hb_rad, hb_speed ] = hyperb_rescaled( unithyperbcoords );
                idx = findintervalbin( hb_rad, maxhyperbcoord, hb_bincount );
                results.hyperbdistribution( idx ) = ...
                    results.hyperbdistribution( idx ) + 1;
                idx = findintervalbin( hb_speed, 1, hbspeed_bincount );
                results.hyperbspeeddistribution( idx ) = ...
                    results.hyperbspeeddistribution( idx ) + 1;
            end
            %% record chains for dimension estimators of diamond link:
            if get.chains
                diamond_sel = center_sel;
                diamond_sel( i ) = true;
                diamond_sel( j ) = true;
                results.chains( 1:5 ) = results.chains( 1:5 ) + ...
                    int64( [ 1; transpose( causet_get_chains( C( diamond_sel, diamond_sel ), 4 ) ) ] );
            end
        end % j = future_idx
        %% statistics for preferred future criteria:
        if ~get.preffutures
            continue
        end
        for c1 = 1 : criteria_count
            for c2 = 1 : criteria_minsize_max
                preffutures_idx = ...
                    causet_find_preftimes( puritycounts, impuritycounts, c1, c2 );
                preftimecount_idx = length( preffutures_idx ) + 1;
                if preftimecount_idx > maxpreffutures
                    results.preffutures.counts( maxpreffutures + 1, c1, c2 ) = ...
                        results.preffutures.counts( maxpreffutures + 1, c1, c2 ) + 1;
                else
                    results.preffutures.counts( preftimecount_idx, c1, c2 ) = ...
                        results.preffutures.counts( preftimecount_idx, c1, c2 ) + 1;
                end
                for j = preffutures_idx
                    % causal from i-th AND causal to j-th:
                    center_sel = C( i, 1:j ) & transpose( C( 1:j, j ) );
                    % size of diamond link:
                    puritycount = puritycounts( j );
                    impuritycount = impuritycounts( j );
                    %% record diamond link type:
                    if get.preffutures_diamonds
                        % limit index range:
                        if impuritycount > maximpurity % overflow
                            impurity_idx = maximpurity + 1;
                        else
                            impurity_idx = impuritycount + 1;
                        end
                        if puritycount > maxpurity % overflow
                            purity_idx = maxpurity + 1;
                        else
                            purity_idx = puritycount;
                        end
                        results.preffutures.diamonds( impurity_idx, purity_idx, c1, c2 ) = ...
                            results.preffutures.diamonds( impurity_idx, purity_idx, c1, c2 ) + 1;
                    end
                    %% record diamond link proper time seperation:
                    if proc.coordinates
                        if ismatrix( g )
                            dcoordinates = coordinates( j, : ) - coordinates( i, : );
                            propertime = ...
                                sqrt( dcoordinates * g * transpose( dcoordinates ) );
                            unithyperbcoords = dcoordinates( 2 : dim ) / propertime;
                        else
                            propertime = ...
                                sqrt( feval( g, coordinates( i, : ), coordinates( j, : ) ) );
                            unithyperbcoords = zeros( 1, dim - 1 );
                        end
                    end
                    if get.preffutures_propertimes
                        idx = findintervalbin( propertime, ...
                            maxpropertime, pt_bincount );
                        results.preffutures.propertimes( idx, c1, c2 ) = ...
                            results.preffutures.propertimes( idx, c1, c2 ) + 1;
                    end
                    if get.preffutures_unithyperboloid
                        unithyperboloid_appendrow( c1, c2, unithyperbcoords );
                    end
                    if get.preffutures_hyperbdistribution
                        [ hb_rad, hb_speed ] = ...
                            hyperb_rescaled( unithyperbcoords );
                        idx = findintervalbin( hb_rad, ...
                            maxhyperbcoord, hb_bincount );
                        results.preffutures.hyperbdistribution( idx, c1, c2 ) = ...
                            results.preffutures.hyperbdistribution( idx, c1, c2 ) + 1;
                        idx = findintervalbin( hb_speed, 1, hbspeed_bincount );
                        results.preffutures.hyperbspeeddistribution( idx, c1, c2 ) = ...
                            results.preffutures.hyperbspeeddistribution( idx, c1, c2 ) + 1;
                    end
                    %% record spacetime dimension of diamond link:
                    if get.preffutures_chains
                        diamond_sel = center_sel;
                        diamond_sel( i ) = true;
                        diamond_sel( j ) = true;
                        results.preffutures.chains( 1:5, c1, c2 ) = ...
                            results.preffutures.chains( 1:5, c1, c2 ) + ...
                            int64( [ 1; transpose( causet_get_chains( C( diamond_sel, diamond_sel ), 4 ) ) ] );
                    end
                end % j = preffutures_idx
            end % c2 = 1 : criteria_minsize_max
        end % c1 = 1 : criteria_count
    end % i = events
    %% remove allocated excess memory:
    if get.preffutures
        if get.preffutures_unithyperboloid
            for c1 = 1 : criteria_count
                for c2 = 1 : criteria_minsize_max
                    results.preffutures.unithyperboloid{ c1, c2 } = ...
                        results.preffutures.unithyperboloid{ c1, c2 }( ...
                            1 : results.preffutures.unithyperboloid_rowcount( c1, c2 ), : );
                end
            end
            results.preffutures = rmfield( results.preffutures, ...
                'unithyperboloid_rowcount' );
        end
    end
    
    %% add row, but allocate memory in chunks:
    function unithyperboloid_appendrow( c1, c2, row )
        r = results.preffutures.unithyperboloid_rowcount( c1, c2 ) + 1;
        if r == 1
            results.preffutures.unithyperboloid{ c1, c2 } = row;
        else
            if r > size( results.preffutures.unithyperboloid{ c1, c2 }, 1 )
                % allocate memory in chunks of 1000 rows:
                results.preffutures.unithyperboloid{ c1, c2 }( r + 1000, 1 ) = int64( 0 );
            end
            results.preffutures.unithyperboloid{ c1, c2 }( r, : ) = row;
        end
        results.preffutures.unithyperboloid_rowcount( c1, c2 ) = r;
    end
end
