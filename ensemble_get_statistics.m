function results = ensemble_get_statistics( N, d, shape, runs, ...
    flags, maxsizes, shapeparam, spacetime )
%ENSEMBLE_GET_STATISTICS repeats the CAUSET_GET_STATISTICS functions 
% for the SPACETIME SHAPE and averages over RUNS repetitions.
% 
% Arguments:
% N                   number of events in the causets if int32 value. If
%                     double value it is the lambda parameter for the
%                     Poisson distribution.
% D                   number of dimensions for the SPACETIME SHAPE.
% SHAPE               name of shape to be generated.
% RUNS                ensemble size.
% 
% Optional arguments:
% FLAGS               (see CAUSET_GET_STATISTICS)
% -showprogress       prints a progress percentage bar in the standard 
%                     output.
% -fields             followed by the list of result fields to be included:
%                     see structure RESULTS below.
% MAXSIZES            3 element vector: 
%  1: MAXPURITY       positive integer for the maximal number of
%                     center elements in a diamond. 
%                     (Detault: min( 2000, 60 * 5^( d - 2 ) ))
%  2: MAXIMPURITY     positive integer for the maximal number of
%                     impurities in a diamond. Set to 0 if diamond-link
%                     impurity should not matter for the counting.
%                     (Default: min( 2000, 4 * maxpurity ))
%  3: MAXPREFFUTURES  positive integer for the maximal number of
%                     preferred futures in the statistics per event. 
%                     (Default: 200)
%  4: MAXDIMENSION    maximal dimension. (Default: 8)
%  5: BINCOUNT        number of histogram bins for the dimension spectra. 
%                     (Default: 202)
% SHAPEPARAM          shape parameter, see CAUSET_NEW_SPRINKLE.
% SPACETIME           name of spacetime.
% 
% Returns: 
% RESULTS             structure with the following fields:
% .eventcounts        number of events in the respective subvolumes.
%                     int64 vector, [ 1, subvolumecount ]
% .distribution       if N is given as double value, i.e. as parameter for
%                     the Poisson distribution, this records the
%                     distribution and its properties:
% ..min, ..max        lower and upper bound of the distribution range
%                     (without under-/overflow bins).
% ..binsize           size of each bin.
% ..bincount          number of bins (including under-/overflow bins).
% ..bins              distribution counters.
%                     int64 vector, [ 1, ..bincount ]
% ..Ns                sum (int64 value) of all causet cardinalities. It 
%                     gives the expectation value of the distribution when 
%                     dividing it by RUNS.
% .comptime           computation time in seconds.
% .futureinfinities   [ imax, subvolumes ] matrix counting the number of
%                     elements in the i-step (i = 1 to 6) future infinity 
%                     falling inside the subvolumes, 1 (top shell), 2
%                     (middle shell), 3 (past bulk). 
% .volume             compact spacetime volume.
% .chains             k-chain counters for k in [ 1, 4 ] (for computation 
%                     of spacetime estimators). 
%                     int64 matrix, [ 3, 4 ]
%                     1st dimension index: 
%                     1: causal interval between events 1 and N
%                     2: larger cau. int. between midpoint and event 1 / N
%                     3: smaller cau. int. between midpoint and event 1 / N
% .maxpropertime      maximal proper time from MAXSIZES( 5 ).
% .maxdimension       maximal dimension from MAXSIZES( 6 ).
% .reducevolumeparams volume reduction parameters for the sprinkling
%                     routine (if used in the simulations).
% .dimestimators      dimension distributions, replaces the .CHAINS results
%                     from CAUSET_GET_STATISTICS.
%                     int64 tensor, [ 3, BINCOUNT, VOLUMESPLITS + 1, 3 ]
%                     first dimension index: 
%                     1: Myrheim-Meyer estimator of flat dimension
%                     2: extended estimator of curved dimension
%                     3: midpoint estimator
%                     fourth dimension index: 
%                     1: entire causet
%                     2: geodesic from event 1 to N by longest link chain
%                     3: geodesic from event 1 to N by maximal volume
% .preffutures
% ..dimestimators     dimension distributions, replaces the ..CHAINS 
%                     results from CAUSET_GET_STATISTICS.
%                     int64 tensor, [ 2, BINCOUNT, criterias( 1 ), 
%                                     criterias( 2 ), VOLUMESPLITS + 1 ]
%                     first dimension index: 
%                     1: Myrheim-Meyer estimator of flat dimension
%                     2: extended estimator of curved dimension
% .(...)              All other results fields are as in 
%                     CAUSET_GET_STATISTICS. Each field has an additional
%                     dimension with index up tp VOLUMESPLITS + 1. The
%                     statistics of the entire causet are stored under
%                     index = 1. Higher indexes correspond to a shell by 
%                     shell reduced volume. 
    
    %% set start up parameter:
    if nargin < 5
        flags = '';
    end
    [ is, proc, get ] = causet_get_statistics_flags( flags, true, true ); %#ok<ASGLU>
    if proc.showprogress
        runjobs_progress_initialize();
        runjobs_progress_update();
    end
    if nargin < 6 || isempty( maxsizes ) || maxsizes( 1 ) < 1
        maxpurity = min( 2000, ceil( 60 * 5^( d - 2 ) ) );
    else
        maxpurity = maxsizes( 1 );
    end
    if nargin >= 6 && length( maxsizes ) >= 2
        maximpurity = maxsizes( 2 );
    else
        maximpurity = min( 2000, 4 * maxpurity );
    end
    if nargin >= 6 && length( maxsizes ) >= 3
        maxpreffutures = maxsizes( 3 );
    else
        maxpreffutures = 200;
    end
    if nargin >= 6 && length( maxsizes ) >= 4
        maxdimension = maxsizes( 4 );
    else
        maxdimension = 8;
    end
    if nargin >= 6 && length( maxsizes ) >= 5
        bincount = maxsizes( 5 );
    else
        bincount = 202;
    end
    if nargin < 8
        spacetime = 'Minkowski';
    end
    propertime_binsize = 0.05; % in [sprinkling units]
    hyperbbincount = 102; % for radial unit hyperboloid distribution
    hyperbbincount_speed = 102; % for radial speed unit hyperboloid distr.
    results.maxhyperbcoord = 20.0;
    results.maxhyperbspeed = 1.0;
    startTimeVal = tic();
    g = metric( d, spacetime );
    %% allocate memory for results fields:
    subvolumecount = 1;
    if strcmp( shape, 'bicone' )
        if isfield( proc, 'reduceto' )
            volumesplits = -1; % single reduced volume
            reducevolumeparams = [ proc.reduceto, -1, proc.volalign ];
        else
            volumesplits = 5; % up to this many volume splits for bicone
            subvolumecount = volumesplits + 1;
            reducevolumeparams = [ 2^(-d/4), subvolumecount, proc.volalign ];
        end
        results.reducevolumeparams = reducevolumeparams;
    else
        volumesplits = 0;
        reducevolumeparams = [ 0, 1 ];
    end
    results.eventcounts = int64( zeros( 1, subvolumecount ) );
    if get.futureinfinities
        futureinfinitymax = 5; % number of k-step future infinities to compute
        results.futureinfinities = ...
            int64( zeros( futureinfinitymax + 1, subvolumecount ) );
    end
    if get.chains
        results.chains = ...
            int64( zeros( 3, 5 ) );
        flags = replace( flags, ' .chains', '' );
    end
    if get.dimestimators
        results.dimestimators = ...
            int64( zeros( 3, bincount, subvolumecount, proc.arrangementcount ) );
        % chain statistics are needed from CAUSET_GET_STATISTICS:
        flags = replace( flags, ' .dimestimators', ' .chains' );
    end
    if get.diamonds
        results.diamonds = ...
            int64( zeros( maximpurity + 1, maxpurity + 1, subvolumecount, proc.arrangementcount ) );
    end
    if get.simplices
        results.simplices = ...
            int64( zeros( 4, bincount, subvolumecount, proc.arrangementcount ) );
    end
    if get.diamondtimes
        results.diamondtimes = ...
            zeros( maximpurity + 1, maxpurity + 1, subvolumecount );
    end
    if get.preffutures
        criterias = [ 6, 1 ]; % check A criteria, each up to B sub-cases: [ A, B ]
        results.preffutures.counts = ...
            int64( zeros( maxpreffutures + 1, criterias( 1 ), criterias( 2 ), subvolumecount ) );
        if get.preffutures_diamonds
            results.preffutures.diamonds = ...
                int64( zeros( maximpurity + 1, maxpurity + 1, criterias( 1 ), criterias( 2 ), subvolumecount ) );
        end
        if get.preffutures_dimestimators
            results.preffutures.dimestimators = ...
                int64( zeros( 2, bincount, criterias( 1 ), criterias( 2 ), subvolumecount ) );
            flags = replace( flags, ' .preffutures.dimestimators', ' .preffutures.chains' );
        end
    else
        criterias = [];
    end
    %% repeat sprinkling process:
    resetN = isfloat( N );
    if ~resetN
        N = double( N );
    else
        results.distribution.binsize = max( 1, floor( sqrt( N ) / 4 ) );
        results.distribution.min = ...
            max( 0, floor( N - 30.5 * results.distribution.binsize ) ) - 0.5;
        results.distribution.bincount = 2 * ceil( ( round( N ) ...
            - results.distribution.min ) / results.distribution.binsize ) + 2;
        results.distribution.max = results.distribution.min ...
            + results.distribution.binsize * ( results.distribution.bincount - 2 );
        results.distribution.bins = int64( zeros( 1, results.distribution.bincount ) );
        results.distribution.Ns = int64( 0 );
    end
    lambda = N;
    for i = 1 : runs
        if resetN
            N = poissrnd( lambda );
            results.distribution.Ns = ...
                results.distribution.Ns + int64( N );
            results.distribution.bins = results.distribution.bins ...
                + int64( convert2interval( N, results.distribution.max, ...
                    results.distribution.bincount, results.distribution.min ) );
        end
        events = 1 : N;
        noevents_sel = false( N, 1 );
        chains = int64( zeros( 5, subvolumecount, proc.arrangementcount ) );
        simplices = int64( zeros( 4, subvolumecount, proc.arrangementcount ) );
        if get.preffutures
            preffutures_chains = ...
                int64( zeros( 5, criterias( 1 ), criterias( 2 ), subvolumecount ) );
        end
        %% sprinkle causet, possibly with further input parameters:
        if nargin < 7 % use default shape parameters:
            [ coordinates, cranges, volume, events_sel ] = ...
                causet_new_sprinkle( N, d, shape, reducevolumeparams );
        else
            [ coordinates, cranges, volume, events_sel ] = ...
                causet_new_sprinkle( N, d, shape, reducevolumeparams, ...
                'global', shapeparam, spacetime );
        end
        %% use sprinkle coordinates for proper times and unit hyperboloid:
        if proc.coordinates
            if ( i == 1 )
                results.unitfactor = ( lambda / volume )^( 1 / d );
                results.maxpropertime = results.unitfactor ...
                    * ( cranges( 2, 1 ) - cranges( 1, 1 ) );
                propertimebincount = min( 2000, ... % limit to 2002 bins
                    ceil( results.maxpropertime / propertime_binsize ) ) + 2;
                maxsizes = [ maxpurity, maximpurity, maxpreffutures, ...
                    results.maxpropertime, propertimebincount, ...
                    results.maxhyperbcoord, hyperbbincount, ...
                    hyperbbincount_speed ];
                if get.propertimes
                    results.propertimes = int64( zeros( propertimebincount, ...
                        subvolumecount, proc.arrangementcount ) );
                end
                if get.hyperbdistribution
                    results.hyperbdistribution = int64( zeros( hyperbbincount, ...
                        subvolumecount, proc.arrangementcount ) );
                    results.hyperbspeeddistribution = ...
                        int64( zeros( hyperbbincount_speed, ...
                        subvolumecount, proc.arrangementcount ) );
                end
                if get.preffutures_propertimes
                    results.preffutures.propertimes = ...
                        int64( zeros( propertimebincount, ...
                        criterias( 1 ), criterias( 2 ), subvolumecount ) );
                end
                if get.preffutures_unithyperboloid
                    results.preffutures.unithyperboloid = ...
                        cell( criterias( 1 ), criterias( 2 ), subvolumecount );
                end
                if get.preffutures_hyperbdistribution
                    results.preffutures.hyperbdistribution = ...
                        int64( zeros( hyperbbincount, ...
                        criterias( 1 ), criterias( 2 ), subvolumecount ) );
                    results.preffutures.hyperbspeeddistribution = ...
                        int64( zeros( hyperbbincount_speed, ...
                        criterias( 1 ), criterias( 2 ), subvolumecount ) );
                end
            end
            % convert to sprinkling units:
            coordinates = coordinates * results.unitfactor;
        end
        %% compute causality and identify future infinity:
        C = causet_edit_relate( coordinates, g );
        L = causet_get_links( C );
        if get.futureinfinities
            infsteps = causet_select_infsteps( L, futureinfinitymax );
        end
        %% set up geodesic chains:
        if ~proc.geodesics
            geodesics = {};
        else
            if proc.linkgeodesic
                [ geodesics1, cardinalities, gevents1 ] = ...
                    causet_find_linkgeodesics( C, L, 1, N ); %#ok<ASGLU>
            end
            if proc.volumegeodesic
                [ geodesics2, cardinalities, gevents2, midpointdim ] = ...
                    causet_find_volumegeodesics( C, L, 1, N ); %#ok<ASGLU>
            end
            if proc.linkgeodesic && proc.volumegeodesic
                geodesics = { geodesics1, geodesics2 };
            elseif proc.linkgeodesic
                geodesics = { geodesics1, {} };
            else
                geodesics = { {}, geodesics2 };
            end
        end
        %% start with past bulk and continue shell by shell to entire set:
        thiseventcounts = zeros( 1, subvolumecount );
        for k = 1 : subvolumecount
            this_events_sel = events_sel( :, k );
            thiseventcounts( k ) = sum( this_events_sel );
            if k <= volumesplits % unselect next shells:
                this_events_sel = this_events_sel ...
                    & ~events_sel( :, k + 1 );
            end
            %% get diamond links statistics for past bulk and future shells:
            thisresults = causet_get_statistics( C, L, maxsizes, ...
                sprintf( '-set %s', replace( flags, ' .chains(geodesics)', '' ) ), ...
                events( this_events_sel ), criterias, coordinates, g );
            if get.preffutures_unithyperboloid
                for c1 = 1 : criterias( 1 )
                    for c2 = 1 : criterias( 2 )
                        results.preffutures.unithyperboloid{ c1, c2, k } = [ ...
                            results.preffutures.unithyperboloid{ c1, c2, k }; 
                            thisresults.preffutures.unithyperboloid{ c1, c2 } ];
                    end
                end
            end
            %% future infinity sizes for past bulk and future shells:
            if get.futureinfinities
                thisfutureinfinities = zeros( futureinfinitymax + 1, 1 );
                for j = 0 : 1 : futureinfinitymax
                    thisfutureinfinities( j + 1, 1 ) = ...
                        sum( ( infsteps <= j ) & transpose( this_events_sel ) );
                end
                thisfutureinfinities = int64( thisfutureinfinities );
            end
            %% add diamond links statistics:
            for j = 1 : k
                if get.dimestimators && ~get.dimestimators_onlyforgeodesics
                    chains( :, j, 1 ) = chains( :, j, 1 ) ...
                      + thisresults.chains;
                end
                if get.simplices
                    simplices( :, j, 1 ) = simplices( :, j, 1 ) ...
                        + thisresults.simplices;
                end
                if get.futureinfinities
                    results.futureinfinities( :, j ) = ...
                        results.futureinfinities( :, j ) ...
                      + thisfutureinfinities;
                end
                if get.diamonds
                    results.diamonds( :, :, j, 1 ) = ...
                        results.diamonds( :, :, j, 1 ) ...
                      + thisresults.diamonds;
                end
                if get.diamondtimes
                    results.diamondtimes( :, :, j ) = ...
                        results.diamondtimes( :, :, j ) ...
                      + thisresults.diamondtimes;
                end
                if get.propertimes
                    results.propertimes( :, j, 1 ) = ...
                        results.propertimes( :, j, 1 ) ...
                      + thisresults.propertimes;
                end
                if get.hyperbdistribution
                    results.hyperbdistribution( :, j, 1 ) = ...
                        results.hyperbdistribution( :, j, 1 ) ...
                      + thisresults.hyperbdistribution;
                    results.hyperbspeeddistribution( :, j, 1 ) = ...
                        results.hyperbspeeddistribution( :, j, 1 ) ...
                      + thisresults.hyperbspeeddistribution;
                end
                if get.preffutures
                    results.preffutures.counts( :, :, :, j ) = ...
                        results.preffutures.counts( :, :, :, j ) ...
                      + thisresults.preffutures.counts;
                    if get.preffutures_dimestimators
                        preffutures_chains( :, :, :, j ) = ...
                            preffutures_chains( :, :, :, j ) ...
                          + thisresults.preffutures.chains;
                    end
                    if get.preffutures_diamonds
                        results.preffutures.diamonds( :, :, :, :, j ) = ...
                            results.preffutures.diamonds( :, :, :, :, j ) ...
                          + thisresults.preffutures.diamonds;
                    end
                    if get.preffutures_propertimes
                        results.preffutures.propertimes( :, :, :, j ) = ...
                            results.preffutures.propertimes( :, :, :, j ) ...
                          + thisresults.preffutures.propertimes;
                    end
                    if get.preffutures_hyperbdistribution
                        results.preffutures.hyperbdistribution( :, :, :, j ) = ...
                            results.preffutures.hyperbdistribution( :, :, :, j ) ...
                          + thisresults.preffutures.hyperbdistribution;
                        results.preffutures.hyperbspeeddistribution( :, :, :, j ) = ...
                            results.preffutures.hyperbspeeddistribution( :, :, :, j ) ...
                          + thisresults.preffutures.hyperbspeeddistribution;
                    end
                end
            end
            %% get geodesics statistics:
            if proc.geodesics
                for m = 2 : proc.arrangementcount
                    for i_g = 1 : length(geodesics{ m - 1 })
                        this_chain_sel = noevents_sel;
                        this_chain_sel( geodesics{ m - 1 }{ i_g } ) = true;
                        thisresults = causet_get_statistics( C, L, maxsizes, ...
                            sprintf( '-chain %s', flags ), ...
                            events( this_events_sel & this_chain_sel ), ...
                            [], coordinates, g );
                        %% add diamond links statistics:
                        for j = 1 : k
                            if get.dimestimators
                                chains( :, j, m ) = chains( :, j, m ) ...
                                  + thisresults.chains;
                            end
                            if get.simplices
                                simplices( :, j, m ) = simplices( :, j, m ) ...
                                    + thisresults.simplices;
                            end
                            if get.diamonds
                                results.diamonds( :, :, j, m ) = ...
                                    results.diamonds( :, :, j, m ) ...
                                  + thisresults.diamonds;
                            end
                            if get.propertimes
                                results.propertimes( :, j, m ) = ...
                                    results.propertimes( :, j, m ) ...
                                  + thisresults.propertimes;
                            end
                            if get.hyperbdistribution
                                results.hyperbdistribution( :, j, m ) = ...
                                    results.hyperbdistribution( :, j, m ) ...
                                  + thisresults.hyperbdistribution;
                                results.hyperbspeeddistribution( :, j, m ) = ...
                                    results.hyperbspeeddistribution( :, j, m ) ...
                                  + thisresults.hyperbspeeddistribution;
                            end
                        end
                    end
                end
            end
        end % k = 1 : subvolumecount
        results.eventcounts = results.eventcounts ...
            + int64( thiseventcounts );
        %% add dimension estimators to statistics:
        for k = 1 : subvolumecount
            if get.dimestimators
                mstart = 1;
                if get.dimestimators_onlyforgeodesics
                    mstart = 2;
                end
                for m = mstart : proc.arrangementcount
                    thismidpointdim = 0;
                    thischains = double( chains( 2:5, k, m ) ) ...
                      / max( 1, double( chains( 1, k, m ) ) );
                    if m == 3
                        thismidpointdim = midpointdim;
                    end
                    results.dimestimators( :, :, k, m ) = ...
                        results.dimestimators( :, :, k, m ) ...
                      + int64( convert2interval( ...
                        [ causet_get_MMdim( thischains( 1:2 ) ), ...
                          causet_get_MMdim( thischains( 1:4 ) ), ...
                          thismidpointdim ], maxdimension, bincount ) );
                end
            end
            if get.simplices
                for m = 1 : proc.arrangementcount
                    results.simplices( :, :, k, m ) = ...
                        results.simplices( :, :, k, m ) ...
                      + int64( convert2interval( ...
                            double( simplices( :, k, m ) ) ...
                          / thiseventcounts( k ), 1, bincount ) );
                end
            end
            if get.preffutures_dimestimators
                for c1 = 1 : criterias( 1 )
                    for c2 = 1 : criterias( 2 )
                        thischains = double( preffutures_chains( 2:5, c1, c2, k ) ) ...
                          / max( 1, double( preffutures_chains( 1, c1, c2, k ) ) );
                        results.preffutures.dimestimators( :, :, c1, c2, k ) = ...
                            results.preffutures.dimestimators( :, :, c1, c2, k ) ...
                          + int64( convert2interval( ...
                            [ causet_get_MMdim( thischains( 1:2 ) ), ...
                              causet_get_MMdim( thischains( 1:4 ) ) ], ...
                            maxdimension, bincount ) );
                    end
                end
            end
        end
        %% add chain counts for global spacetime estimation to statistics:
        if get.chains
            midpoints = causet_get_midpoint( C );
            tot_sel = causet_find_Aset( C, 1, N );
            bot_sel = causet_find_Aset( C, 1, midpoints( 1 ) );
            top_sel = causet_find_Aset( C, midpoints( 1 ), N );
            tot_chains = int64( causet_get_chains( C( tot_sel, tot_sel ), 4 ) );
            bot_chains = int64( causet_get_chains( C( bot_sel, bot_sel ), 4 ) );
            top_chains = int64( causet_get_chains( C( top_sel, top_sel ), 4 ) );
            results.chains( 1, : ) = results.chains( 1, : ) ...
              + int64( [ 1, tot_chains ] );
            if bot_chains( 1 ) > top_chains( 1 )
                results.chains( 2, : ) = results.chains( 2, : ) ...
                  + int64( [ 1, bot_chains ] );
                results.chains( 3, : ) = results.chains( 3, : ) ...
                  + int64( [ 1, top_chains ] );
            else
                results.chains( 2, : ) = results.chains( 2, : ) ...
                  + int64( [ 1, top_chains ] );
                results.chains( 3, : ) = results.chains( 3, : ) ...
                  + int64( [ 1, bot_chains ] );
            end
        end
        if proc.showprogress
            runjobs_progress_update( i / runs );
        end
    end % i = 1 : runs
    if proc.showprogress
        runjobs_progress_finalize();
    end
    %% record further info:
    results.maxdimension = maxdimension;
    results.comptime = toc( startTimeVal );
end

