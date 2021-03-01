function runjobs_progress_update( progress )
%RUNJOBS_PROGRESS_UPDATE Extends the print out of the progress.
% 
% Arguments:
% PROGRESS            progress percentage.
    
    persistent prgcounter
    if ( nargin < 1 ) || ( progress <= 0 )
        prgcounter = 0;
    else
        newprgcounter = floor( 20 * progress );
        fprintf( repmat( '#', 1, newprgcounter - prgcounter ) );
        prgcounter = newprgcounter;
    end
end

