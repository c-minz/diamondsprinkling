function L = causet_get_links( C )
%CAUSET_GET_LINKS converts the causals matrix into a links matrix.
% 
% Arguments:
% C                   upper triangular (logical) causals matrix.
% 
% Returns:
% L                   upper triangular (logical) links matrix.
    
    if ~islogical( C )
        C = logical( C ); % only for backwards compatibility, 
        % causal and link matrixes are now always logical
    end
    N = size( C, 1 );
    L = false( N );
    for i = 1 : N
        %  Selector for the i-th future light cone:
        causalsel = C( i, : );
        %  Vector of the causal connection count TO each event in the 
        %  future light cone FROM any other event in the cone:
        connections = sum( C( causalsel, : ), 1 );
        %  If such a number is greater than 1, then the respective event 
        %  is not linked because there is a longer path to it. Linked are
        %  only those events that do not have any connection by other
        %  events in the future light cone:
        L( i, : ) = causalsel & ( connections == 0 );
    end
end

