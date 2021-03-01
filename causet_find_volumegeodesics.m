function [ geodesics, cardinalities, allevents, midpointdim ] = ...
    causet_find_volumegeodesics( C, L, pastevent, futureevent )
%CAUSET_FIND_VOLUMEGEODESICS finds all geodesics between the events with 
% indexes PASTEVENT and FUTUREVENT. At each point the geodesic maximizes 
% the volume between the point and FUTUREEVENT and, secondly, minimizes 
% the volume between the point and PASTEVENT.
% 
% Arguments:
% C                   logical upper triangular causal matrix.
% L                   logical upper triangular link matrix.
% PASTEVENT           index of the event where the geodesics start.
% FUTUREEVENT         index of the event where the geodesics end.
% 
% Returns:
% GEODESICS           cell column vector with geodesics as a row vectors.
%                     The geodesics are sorted descending by their 
%                     cardinality, so that the first geodesic is the 
%                     longest. If the events are causally disconnected then
%                     the cell vector is empty.
% CARDINALITIES       is a column vector holding the length of each
%                     geodesic.
% ALLEVENTS           is a vector of the indexes of all events for any
%                     bifurcation.
% MIDPOINTDIM         dimension estimate using the midpoint approach, 
%                     d = log2( N / max( N_small ) ), where N is the
%                     cardinality of the Alexandrov set between PASTEVENT
%                     and FUTUREEVENT and N_small is the cardinality of the
%                     smaller Alexandrov set between the midpoint and
%                     PASTEVENT or FUTUREEVENT.
    
    %% analyse elements causally between past- and future-event:
    Aset = causet_find_Aset( C, pastevent, futureevent );
    Aset_card = length( Aset );
    all_sel = false( 1, Aset_card );
    all_sel( Aset == pastevent ) = true;
    all_sel( Aset == futureevent ) = true;
    geodesics = cell( 0, 1 );
    cardinalities = zeros( Aset_card, 1 );
    N_small = 0;
    if Aset_card == 0
        bifurcations = 0;
    else
        %% pre-allocate memory:
        b = 1;
        bifurcations = 1;
        geodesics{ 1 } = zeros( 1, Aset_card );
        event = pastevent;
        card = 1;
        %% step through (bifurcations) of geodesics:
        geodesics{ 1 }( 1 ) = pastevent;
        cardinalities( 1 ) = 1;
        while b <= bifurcations
            nextevents = find( L( event, : ) & transpose( C( :, futureevent ) ) );
            if isempty( nextevents )
                %% add last event:
                card = card + 1;
                geodesics{ b }( card ) = futureevent;
                cardinalities( b ) = card;
                %% continue with next bifurcation:
                b = b + 1;
                if b <= bifurcations
                    event = geodesics{ b }( cardinalities( b ) );
                    card = cardinalities( b );
                end
                continue
            end
            nextevents_count = length( nextevents );
            nextevents_pvolumes = zeros( 1, nextevents_count );
            nextevents_fvolumes = zeros( 1, nextevents_count );
            for n = 1 : nextevents_count
                nextevent = nextevents( n );
                nextevents_pvolumes( n ) = ...
                    sum( C( pastevent, : ) & transpose( C( :, nextevent ) ) ) + 2;
                nextevents_fvolumes( n ) = ...
                    sum( C( nextevent, : ) & transpose( C( :, futureevent ) ) ) + 2;
            end
            %% use volume maximization to the future, minimizing to the past:
            nextevents_sel = nextevents_fvolumes == max( nextevents_fvolumes );
            nextevents_sel = nextevents_sel ...
                & ( nextevents_pvolumes == min( nextevents_pvolumes( nextevents_sel ) ) );
            nextevents = nextevents( nextevents_sel );
            nextevents_count = length( nextevents );
            N_small = max( N_small, ...
                max( min( nextevents_pvolumes( nextevents_sel ), ...
                          nextevents_fvolumes( nextevents_sel ) ) ) );
            %% store bifurcations for later loop runs:
            if nextevents_count > 1
                bifurcations_new = bifurcations + nextevents_count - 1;
                m = 1;
                for b2 = ( bifurcations + 1 ) : 1 : bifurcations_new
                    m = m + 1;
                    card2 = cardinalities( b ) + 1;
                    geodesics{ b2 } = geodesics{ b };
                    geodesics{ b2 }( card2 ) = nextevents( m );
                    cardinalities( b2 ) = card2;
                end
                bifurcations = bifurcations_new;
            end
            %% add next element (of current bifurcation):
            event = nextevents( 1 );
            card = card + 1;
            geodesics{ b }( card ) = event;
            cardinalities( b ) = card;
            all_sel( Aset == event ) = true;
        end
    end
    %% remove pre-allocation of each geodesic:
    cardinalities = cardinalities( 1 : bifurcations );
    for b = 1 : bifurcations
        geodesic = geodesics{ b };
        geodesics{ b } = geodesic( 1 : cardinalities( b ) );
    end
    %% sort geodesics by cardinalities:
    if bifurcations > 1
        [ cardinalities, indexes ] = sort( cardinalities, 'descend' );
        geodesics = geodesics( indexes );
    end
    allevents = Aset( all_sel );
    %% get midpoint dimension estimate:
    if N_small > 0
        midpointdim = log2( Aset_card / N_small );
    else
        midpointdim = 0;
    end
end
