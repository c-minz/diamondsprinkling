function events = causet_find_Aset( C, pastevent, futureevent )
%CAUSET_FIND_ASET returns the indexes of the elements in the Alexandrov set
% between the elements with indexes PASTEVENT and FUTUREEVENT. 
% 
% Arguments:
% C                   logical upper triangular causal matrix.
% PASTEVENT           index of the first event in the set.
% FUTUREEVENT         index of the last event in the set.
% 
% Returns:
% EVENTS              list of events, which are causally between PASTEVENT
%                     and FUTUREEVENT. If they are not causally related the
%                     return is an empty row vector. Otherwise the return
%                     also includes PASTEVENT as well as FUTUREEVENT and
%                     returns only a single index if both events are the
%                     same.
    
    if futureevent == pastevent % only one element
        events = pastevent;
    elseif ~C( pastevent, futureevent ) % no causal relation between the events
        events = zeros( 1, 0 );
    else % pastevent, all elements causally between, futureevent
        events = [ pastevent, ...
            find( C( pastevent, 1:futureevent ) ...
                & transpose( C( 1:futureevent, futureevent ) ) ), ...
            futureevent ];
    end
end

