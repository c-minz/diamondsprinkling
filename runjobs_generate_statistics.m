function runjobs_generate_statistics( task, N, dims, shapes, runs, flags, maxsizes )
%RUNJOBS_GENERATE_STATISTICS is a wrapper function for 
% ENSEMBLE_GET_STATISTICS as its n-th task of a task sequence.
% 
% The TASKINDEX is used to select one dimension and one shape from the 
% array DIMS and the cell list SHAPES by modulo and integer division, 
% respectively. 
% 
% Arguments:
% TASK                task index, positive integer. The generated data is 
%                     stored in a file with a name including the job index. 
%                     A string can be specified to add to the filename
%                     instead of the job index.
% N                   number of elements in the causet.
% DIMS                vector of spacetime dimensions.
% SHAPES              cell array of shape names or string with comma 
%                     separated list of names.
% RUNS                number of loop iterations.
% 
% Optional arguments:
% MAXSIZES            memory dimensions of the results, see
%                     ENSEMBLE_GET_STATISTICS.

    %% set start up parameter:
    isSingleJob = false;
    if ~isnumeric( task )
        task_num = str2double( task );
        if isnan( task_num )
            isSingleJob = true;
            jobname = task;
            task = 1;
        else
            task = task_num;
        end
    end
    if ~iscell( shapes )
        shapes = strsplit( shapes, ',' );
    end
    %% convert taskindex to job parameters:
    count_d = length( dims );
    count_shape = length( shapes );
    dshape = mod( task - 1, count_d * count_shape ) + 1;
    % job index:
    j = floor( ( task - 1 ) / ( count_d * count_shape ) ) + 1;
    % char sequence as shape name, selected from the shapes array:
    shape = char( shapes( floor( ( dshape - 1 ) / count_d ) + 1 ) );
    % spacetime dimension, selected from the dimensions array:
    d = dims( mod( dshape - 1, count_d ) + 1 );
    %% initialize pseudo-random generator:
    stream = RandStream( 'mlfg6331_64' );
    stream.Substream = task;
    RandStream.setGlobalStream( stream );
    %% generate causal sets and count diamonds:
    if nargin < 6
        results = ensemble_get_statistics( N, d, shape, runs );
    elseif nargin < 7
        results = ensemble_get_statistics( N, d, shape, runs, flags );
    else
        results = ensemble_get_statistics( N, d, shape, runs, flags, maxsizes );
    end
    results.runs = runs;
    results.N = N;
    results.d = d;
    results.shape = shape;
    results.comptimestr = char( seconds( results.comptime ), 'hh:mm:ss' );
    %% save to file:
    if isSingleJob
        filename = sprintf( 'data/N%dD%d%sJ%s.mat', N, d, shape, jobname );
    else
        filename = sprintf( 'data/N%dD%d%sJ%d.mat', N, d, shape, j );
    end
    save( filename, '-struct', 'results' );
    info = whos( 'results' );
    fprintf( 'Task %d has finished after %s. Results (%0.3f kiB) written to file ''%s''.\n', ...
        task, results.comptimestr, info.bytes / 1024, filename );
end

