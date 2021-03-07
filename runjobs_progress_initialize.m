function runjobs_progress_initialize()
%RUNJOBS_PROGRESS_INITIALIZE Prints the first lines of the progress bar to
% the std output.
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.

    fprintf( 'Progress\n' );
    fprintf( 'Percent: 0                   100\n' );
    fprintf( '         |____________________|\n' );
    fprintf( '         [' );
end

