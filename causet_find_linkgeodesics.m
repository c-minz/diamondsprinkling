function [ geodesics, cardinality, allevents ] = ...
    causet_find_linkgeodesics( C, L, pastevent, futureevent )
%CAUSET_FIND_LINKGEODESICS finds all geodesics between the events with 
% indexes PASTEVENT and FUTUREVENT, which maximize the link count. 
% 
% Arguments:
% C                   logical upper triangular causal matrix.
% L                   logical upper triangular link matrix.
% PASTEVENT           index of the event where the geodesics start.
% FUTUREEVENT         index of the event where the geodesics end.
% 
% Returns:
% GEODESICS           cell vector of geodesics, each as a row 
%                     vector of event indexes. 
% CARDINALITY         is a vector holding the length of each geodesic.
% ALLEVENTS           is a vector of the indexes of all events for any
%                     bifurcation.
    
    %% use Alexandrov set between past- and future-event:
    Aset = causet_find_Aset( C, pastevent, futureevent );
    Aset_card = length( Aset );
    if Aset_card == 0
        geodesics = cell( 0, 1 );
        cardinality = 0;
        allevents = [];
    elseif Aset_card == 1
        geodesics = cell( 1, 1 );
        geodesics{ 1 } = pastevent;
        cardinality = 1;
        allevents = pastevent;
    else
        Aset_L = L( Aset, Aset );
        all_sel = false( 1, Aset_card );
        all_sel( Aset == pastevent ) = true;
        all_sel( Aset == futureevent ) = true;
        %% compute maximal path-lengths from each event to the top:
        maxpathlengths = zeros( 1, Aset_card );
        maxpathlengths( Aset == futureevent ) = 1;
        for i = Aset_card : -1 : 2
            past_sel = Aset_L( :, i );
            maxpathlengths( past_sel ) = ...
                max( maxpathlengths( past_sel ), maxpathlengths( i ) + 1 );
        end
        %% get path along maximal path-length:
        bifurcations = 1;
        geodesics = cell( 1, 1 );
        cardinality = max( maxpathlengths );
        geodesics{ 1 } = zeros( 1, cardinality );
        geodesics{ 1 }( 1 ) = 1;
        for i = 1 : ( cardinality - 2 )
            for b = 1 : bifurcations
                prevevent = geodesics{ b }( i );
                this_sel = Aset_L( prevevent, : );
                this_sel = this_sel ...
                  & ( maxpathlengths == max( maxpathlengths( this_sel ) ) );
                future_sel = logical( sum( Aset_L( this_sel, : ), 1 ) );
                future_sel = future_sel ...
                  & ( maxpathlengths == max( maxpathlengths( future_sel ) ) );
                all_sel = all_sel | this_sel | future_sel;
                nowbifurcate = false;
                for j = find( this_sel )
                    j_future = find( Aset_L( j, : ) & future_sel );
                    if ~isempty( j_future )
                        if nowbifurcate
                            bifurcations = bifurcations + 1;
                            geodesics{ bifurcations } = geodesics{ b };
                            geodesics{ bifurcations }( i + 1 ) = j;
                        else
                            geodesics{ b }( i + 1 ) = j;
                            nowbifurcate = true;
                        end
                        future_sel( j_future ) = false;
                    end
                end
            end
        end
        %% translate to indexes of entire causet:
        for b = 1 : bifurcations
            for i = 1 : cardinality - 1
                geodesics{ b }( i ) = Aset( geodesics{ b }( i ) );
            end
            geodesics{ b }( cardinality ) = futureevent;
        end
        allevents = Aset( all_sel );
    end
end
