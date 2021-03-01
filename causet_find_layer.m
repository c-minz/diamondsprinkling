function events = causet_find_layer( C, event, layer )
%CAUSET_FIND_LAYER returns the indexes of the elements in the layer of 
% index LAYER for EVENT.
    
    N = size( C, 1 );
    events = zeros( 1, N );
    j = 0;
    if layer >= 0 % future (or present) layer:
        for i = event : N
            i_layer = length( causet_find_Aset( C, event, i ) ) - 1;
            if i_layer == layer
                j = j + 1;
                events( j ) = i;
            end
        end
    else % past layer:
        for i = 1 : ( event - 1 )
            i_layer = length( causet_find_Aset( C, i, event ) ) - 1;
            if i_layer == -layer
                j = j + 1;
                events( j ) = i;
            end
        end
    end
    events = sort( events( 1 : j ) );
end

