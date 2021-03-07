function [ spacedim, simplexcount, events ] = causet_get_simplices( C, L, preselected, ...
    maxsimplex )
%CAUSET_GET_SIMPLICES identifies the space dimension (up to 3) and returns 
% the events that form the respective simplex(es). 
% 
% Arguments:
% C                   logical upper triangular causal matrix.
% L                   logical upper triangular link matrix.
% PRESELECTED         [ 1, n ] cell array for the event indexes that are in
%                     the simplex(es). The first cell is a vector of 0-face 
%                     events (vertices), the second cell is a vector 
%                     1-faces (edges), and so on. 
%                     Use { x } to get the largest simplex at event x. Then
%                     the events A1, A2, A3 are srictly increasing. 
% 
% Optional arguments:
% MAXSIMPLEX          largest simplex to search. Use 
%                     n == causet_get_simplices( C, L, { x }, n )
%                     to quickly check if at least one n-simplex is present
%                     at event x.
%                     Default: 3 (supported maximum)
% 
% Returns:
% SPACEDIM            index of the largest simplex found at POSITION and
%                     also the number of faces n.
% SIMPLEXCOUNT        number m of n-simplices.
% EVENTS              [ m, n ] cell matrix for the events in the
%                     simplex(es). There are m rows for each n-simplex.
%                     The first column holds the 0-face events, the second
%                     column the 1-face events, etc.
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.
    
    %% initialize:
    if nargin < 4
        maxsimplex = 3;
    end
    spacedim = 0; % initialize as time-line
    presel_facecount = length( preselected );
    presel_0face_count = length( preselected{ 1 } );
    if presel_facecount >= 2
        presel_1face_count = length( preselected{ 2 } );
    else
        presel_1face_count = 0;
    end
    if presel_facecount >= 3
        presel_2face_count = length( preselected{ 3 } );
    else
        presel_2face_count = 0;
    end
    break_simplexsearch = false;
    sel_all = true( 1, size( L, 1 ) );
    %% select event A1:
    newnexttoA1 = sel_all; % stores previous A1 indexes
    A1 = preselected{ 1 }( 1 );
    newnexttoA1( A1 ) = false; % further ignore event
    %% select events B1:
    newnexttoB1 = sel_all;
    afterA1 = L( A1, : );
    if ( presel_1face_count >= 1 ) ...
      && afterA1( preselected{ 2 }( 1 ) )
        set_B1 = preselected{ 2 }( 1 );
    else
        set_B1 = find( afterA1 );
    end
    for B1 = set_B1
    newnexttoB1( B1 ) = false;
    %% select events A2:
    newnexttoA2 = newnexttoA1;
    beforeB1 = transpose( L( :, B1 ) );
    sel_A2 = newnexttoA1 & beforeB1;
    if ( presel_0face_count >= 2 ) ...
      && sel_A2( preselected{ 1 }( 2 ) )
        set_A2 = preselected{ 1 }( 2 );
    else
        set_A2 = find( sel_A2 );
    end
    for A2 = set_A2
    newnexttoA2( A2 ) = false;
    %% at least 1-simplex:
    if spacedim < 1
        events = cell( 1, 1 );
        spacedim = 1;
        simplexcount = 0;
    end
    if maxsimplex > 1
    %% select events B2:
    newnexttoB2 = newnexttoB1;
    ninfutrofA2 = ~C( A2, : );
    ninpastofB1 = ~transpose( C( :, B1 ) );
    septoA2B1 = ninfutrofA2 & ninpastofB1;
    sel_B2 = newnexttoB1 & afterA1 & septoA2B1;
    if ( presel_1face_count >= 2 ) ...
      && sel_B2( preselected{ 2 }( 2 ) )
        set_B2 = preselected{ 2 }( 2 );
    else
        set_B2 = find( sel_B2 );
    end
    for B2 = set_B2
    newnexttoB2( B2 ) = false;
    %% select events A3:
    newnexttoA3 = newnexttoA2;
    beforeB2 = transpose( L( :, B2 ) );
    sel_A3 = newnexttoA2 & beforeB2 & septoA2B1;
    if ( presel_0face_count >= 3 ) ...
      && sel_A3( preselected{ 1 }( 3 ) )
        set_A3 = preselected{ 1 }( 3 );
    else
        set_A3 = find( sel_A3 );
    end
    for A3 = set_A3
    newnexttoA3( A3 ) = false;
    %% select events B3:
    newnexttoB3 = newnexttoB2;
    septoA1 = ~( C( A1, : ) | transpose( C( :, A1 ) ) );
    afterA2 = L( A2, : );
    afterA3 = L( A3, : );
    sel_B3 = newnexttoB2 & afterA2 & afterA3 & septoA1;
    if ( presel_1face_count >= 3 ) ...
      && sel_B3( preselected{ 2 }( 3 ) )
        set_B3 = preselected{ 2 }( 3 );
    else
        set_B3 = find( sel_B3 );
    end
    for B3 = set_B3
    newnexttoB3( B3 ) = false;
    %% at least 2-simplex:
    if spacedim < 2
        events = cell( 1, 2 );
        spacedim = 2;
        simplexcount = 0;
    end
    if maxsimplex > 2
    %% select events B4:
    newnexttoB4 = newnexttoB3;
    ninfutrofA3 = ~C( A3, : );
    ninpastofB3 = ~transpose( C( :, B3 ) );
    septoA3B3 = ninfutrofA3 & ninpastofB3;
    sel_B4 = newnexttoB3 & afterA1 & septoA2B1 & septoA3B3;
    if ( presel_1face_count >= 4 ) ...
      && sel_B4( preselected{ 2 }( 4 ) )
        set_B4 = preselected{ 2 }( 4 );
    else
        set_B4 = find( sel_B4 );
    end
    for B4 = set_B4
    newnexttoB4( B4 ) = false;
    %% select events B5:
    newnexttoB5 = newnexttoB4;
    ninpastofB2 = ~transpose( C( :, B2 ) );
    ninpastofB4 = ~transpose( C( :, B4 ) );
    septoA1B2B4 = septoA1 & ninpastofB2 & ninpastofB4;
    sel_B5 = newnexttoB4 & afterA2 & septoA1B2B4 & septoA3B3;
    if ( presel_1face_count >= 5 ) ...
      && sel_B5( preselected{ 2 }( 5 ) )
        set_B5 = preselected{ 2 }( 5 );
    else
        set_B5 = find( sel_B5 );
    end
    for B5 = set_B5
    newnexttoB5( B5 ) = false;
    %% select events A4:
    newnexttoA4 = newnexttoA3;
    ninpastofB1 = ~transpose( C( :, B1 ) );
    septoA3B2B3 = septoA3B3 & ninpastofB2;
    beforeB4 = transpose( L( :, B4 ) );
    beforeB5 = transpose( L( :, B5 ) );
    sel_A4 = newnexttoA3 & beforeB4 & beforeB5 & septoA3B2B3 & ninpastofB1;
    if ( presel_0face_count >= 4 ) ...
      && sel_A4( preselected{ 1 }( 4 ) )
        set_A4 = preselected{ 1 }( 4 );
    else
        set_A4 = find( sel_A4 );
    end
    for A4 = set_A4
    newnexttoA4( A4 ) = false;
    %% select events B6:
    afterA4 = L( A4, : );
    septoA1B1B2B4 = septoA1B2B4 & ninpastofB1;
    sel_B6 = newnexttoB5 & afterA3 & afterA4 & septoA1B1B2B4;
    if ( presel_1face_count >= 6 ) ...
      && sel_B6( preselected{ 2 }( 6 ) )
        set_B6 = preselected{ 2 }( 6 );
    else
        set_B6 = find( sel_B6 );
    end
    for B6 = set_B6
    %% select events C1:
    newnexttoC1 = newnexttoA4;
    beforeB3 = transpose( L( :, B3 ) );
    ninfutrofA4 = ~C( A4, : );
    ninpastofB5 = ~transpose( C( :, B5 ) );
    ninpastofB6 = ~transpose( C( :, B6 ) );
    sel_C1 = newnexttoA4 & beforeB1 & beforeB2 & beforeB3 ...
            & ninfutrofA4 & ninpastofB4 & ninpastofB5 & ninpastofB6;
    if ( presel_2face_count >= 1 ) ...
      && sel_C1( preselected{ 3 }( 1 ) )
        set_C1 = preselected{ 3 }( 1 );
    else
        set_C1 = find( sel_C1 );
    end
    for C1 = set_C1
    newnexttoC1( C1 ) = false;
    %% select events C2:
    newnexttoC2 = newnexttoC1;
    sel_C2 = newnexttoC1 & beforeB1 & beforeB4 & beforeB5 ...
            & septoA3B2B3 & ninpastofB6;
    if ( presel_2face_count >= 2 ) ...
      && sel_C2( preselected{ 3 }( 2 ) )
        set_C2 = preselected{ 3 }( 2 );
    else
        set_C2 = find( sel_C2 );
    end
    for C2 = set_C2
    newnexttoC2( C2 ) = false;
    %% select events C3:
    newnexttoC3 = newnexttoC2;
    beforeB6 = transpose( L( :, B6 ) );
    sel_C3 = newnexttoC2 & beforeB2 & beforeB4 & beforeB6 ...
            & septoA2B1 & ninpastofB3 & ninpastofB5;
    if ( presel_2face_count >= 3 ) ...
      && sel_C3( preselected{ 3 }( 2 ) )
        set_C3 = preselected{ 3 }( 3 );
    else
        set_C3 = find( sel_C3 );
    end
    for C3 = set_C3
    newnexttoC3( C3 ) = false;
    %% select events C4:
    sel_C4 = newnexttoC3 & beforeB3 & beforeB5 & beforeB6 & septoA1B1B2B4;
    if ( presel_2face_count >= 4 ) ...
      && sel_C4( preselected{ 3 }( 2 ) )
        set_C4 = preselected{ 3 }( 4 );
    else
        set_C4 = find( sel_C4 );
    end
    for C4 = set_C4
    %% at least 3-simplex:
    if spacedim < 3
        events = cell( 1, 3 );
        spacedim = 3;
        simplexcount = 0;
    end
    if maxsimplex > 3
    %% search for 4-simplices not implemented:
    warning( 'A search for 4-simplices is not implemented.' );
    end
    if spacedim == 3
        %% if 3-simplex is the largest:
        if nargout > 1
            simplexcount = simplexcount + 1;
            if ( nargout > 2 )
                events{ simplexcount, 1 } = [ A1 A2 A3 A4 ];
                events{ simplexcount, 2 } = [ B1 B2 B4 B3 B5 B6 ];
                events{ simplexcount, 3 } = [ C1 C2 C3 C4 ];
            end
        elseif maxsimplex == 3
            break_simplexsearch = true;
        end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    end % maxsimplex > 2
    if spacedim == 2
        %% if 2-simplex is the largest:
        if nargout > 1
            simplexcount = simplexcount + 1;
            if ( nargout > 2 )
                events{ simplexcount, 1 } = [ A1 A2 A3 ];
                events{ simplexcount, 2 } = [ B1 B2 B3 ];
            end
        elseif maxsimplex == 2
            break_simplexsearch = true;
        end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    end % maxsimplex > 1
    if spacedim == 1
        %% if 1-simplex is the largest:
        if nargout > 1
            simplexcount = simplexcount + 1;
            if ( nargout > 2 )
                events{ simplexcount, 1 } = [ A1 A2 ];
            end
        elseif maxsimplex == 1
            break_simplexsearch = true;
        end
    end
    if break_simplexsearch; break; end
    end
    if break_simplexsearch; break; end
    end
    %% if 0-simplex is the largest:
    if ( spacedim == 0 ) && ( nargout > 1 )
        simplexcount = 1;
        if ( nargout > 2 )
            events = { A1 };
        end
    end
end
