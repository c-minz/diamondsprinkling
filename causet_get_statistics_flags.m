function [ is, proc, get ] = causet_get_statistics_flags( flagsstring, hascoordinates, hascriterias )
%CAUSET_GET_STATISTICS_FLAGS splits the flags for the causet statistics 
% function in a structure of booleans.
    
    is.eventchain = ~isempty( strfind( flagsstring, '-chain' ) ); %#ok<*STREMP>
    is.eventset = ~is.eventchain;
    
    proc.coordinates = hascoordinates;
    proc.criterias = hascriterias;
    if ~isempty( strfind( flagsstring, '-volalign=future' ) )
        proc.volalign = 1;
    elseif ~isempty( strfind( flagsstring, '-volalign=center' ) )
        proc.volalign = 0;
    else
        proc.volalign = -1;
    end
    key = '-reduceto=';
    if ~isempty( strfind( flagsstring, key ) )
        flagsstring = [ flagsstring, ' ' ];
        i = strfind( flagsstring, key ) + length( key );
        j_sel = strfind( flagsstring, ' ' );
        j_sel = j_sel( j_sel >= i );
        j = j_sel( 1 ) - 1;
        proc.reduceto = str2double( flagsstring( i : j ) );
    end
    key = '-geoendtime=';
    if ~isempty( strfind( flagsstring, key ) )
        flagsstring = [ flagsstring, ' ' ];
        i = strfind( flagsstring, key ) + length( key );
        j_sel = strfind( flagsstring, ' ' );
        j_sel = j_sel( j_sel >= i );
        j = j_sel( 1 ) - 1;
        proc.geoendtime = str2double( flagsstring( i : j ) );
    end
    proc.showprogress = ~isempty( strfind( flagsstring, '-showprogress' ) );
    proc.linkgeodesic = ~isempty( strfind( flagsstring, '-linkgeodesic' ) ) ...
        || ~isempty( strfind( flagsstring, '-geodesics' ) );
    proc.volumegeodesic = ~isempty( strfind( flagsstring, '-volumegeodesic' ) ) ...
        || ~isempty( strfind( flagsstring, '-geodesics' ) );
    proc.geodesics = proc.linkgeodesic || proc.volumegeodesic;
    if proc.volumegeodesic
        proc.arrangementcount = 5; % 4: all, 5: single volume geodesic
    elseif proc.linkgeodesic
        proc.arrangementcount = 3; % 2: all, 3: single link geodesic
    else
        proc.arrangementcount = 1; % 1: entire causet
    end
    
    get.allfields = isempty( strfind( flagsstring, '-fields' ) );
    get.futureinfinities = get.allfields ...
        || ~isempty( strfind( flagsstring, '.futureinfinities' ) );
    get.dimestimators = get.allfields ...
        || ~isempty( strfind( flagsstring, '.dimestimators' ) );
    get.dimestimators_onlyforgeodesics = get.allfields ...
        || ~isempty( strfind( flagsstring, '.dimestimators(geodesics)' ) );
    get.diamonds = get.allfields ...
        || ~isempty( strfind( flagsstring, '.diamonds' ) );
    get.simplices = get.allfields ...
        || ~isempty( strfind( flagsstring, '.simplices' ) );
    get.chains = get.allfields ...
        || ~isempty( strfind( flagsstring, '.chains' ) );
    get.diamondtimes = hascoordinates ...
        && ( ( get.allfields && is.eventset ) ...
           || ( ~get.allfields && ~isempty( strfind( flagsstring, '.diamondtimes' ) ) ) );
    get.propertimes = hascoordinates ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.propertimes' ) ) );
    get.hyperbdistribution = hascoordinates ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.hyperbdistribution' ) ) );
    get.preffutures = hascriterias ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.preffutures' ) ) );
    get.preffutures_dimestimators = get.preffutures ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.preffutures.dimestimators' ) ) ...
           || ~isempty( strfind( flagsstring, '.preffutures.*' ) ) );
    get.preffutures_diamonds = get.preffutures ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.preffutures.diamonds' ) ) ...
           || ~isempty( strfind( flagsstring, '.preffutures.*' ) ) );
    get.preffutures_propertimes = get.preffutures ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.preffutures.propertimes' ) ) ...
           || ~isempty( strfind( flagsstring, '.preffutures.*' ) ) );
    get.preffutures_unithyperboloid = get.preffutures ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.preffutures.unithyperboloid' ) ) ...
           || ~isempty( strfind( flagsstring, '.preffutures.*' ) ) );
    get.preffutures_hyperbdistribution = get.preffutures ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.preffutures.hyperbdistribution' ) ) ...
           || ~isempty( strfind( flagsstring, '.preffutures.*' ) ) );
    %% getter for CAUSET_GET_STATISTICS only:
    get.onlysimplices = get.simplices ...
        && ~( get.diamonds || get.chains || get.diamondtimes ...
             || get.propertimes || get.preffutures );
    get.preffutures_chains = get.preffutures ...
        && ( get.allfields ...
           || ~isempty( strfind( flagsstring, '.preffutures.chains' ) ) ...
           || ~isempty( strfind( flagsstring, '.preffutures.*' ) ) );
end
