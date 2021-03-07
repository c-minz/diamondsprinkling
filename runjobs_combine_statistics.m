function runjobs_combine_statistics( Ns, dims, shapes, jobs, extend, divisions, deljobs )
%RUNJOBS_COMBINE_STATISTICS loads single jobs from a job sequence with the
% same N's, dimensions and shapes and combines them to a single data file.
% All files are read from and written to a subdirectory 'data'.
% 
% Arguments:
% NS                  array of numbers of elements in the causet.
% DIMS                array of spacetime dimensions.
% SHAPES              string array of shape names.
% JOBS                array of job-indexes to be combined.
% 
% Optional arguments:
% EXTEND              boolean flag to extend existing files.
%                     true:  job file data is added to existing files
%                            (Default)
%                     false: combined job file data overwrites existing
%                            files
% DIVISIONS           create new job file by combining only a subset of all
%                     files; necessary for combining large data in 
%                     parallel. New job files get the name extension JAxB 
%                     where A stands for the new division order and 
%                     B stands for the new job number.
%                     A:           A is the new devision order (Default: 1)
%                     [ A, m, n ]: A is the new devision order (in 
%                     the first call it shall be 1), m is subfile index of 
%                     n new job files in this division.
%                     n many files with job indexes in JOBS (and in 
%                     division A-1 if A>0) are read.
% DELJOBS             boolean flag to delete job files. 
%                     false: job files remain (Default)
%                     true:  job files will be deleted after loading data
% 
% Copyright 2021, C. Minz. BSD 3-Clause License.

    %% set start up parameter:
    if nargin < 5
        extend = true;
    end
    if nargin < 6
        divisions = 1;
    end
    if nargin < 7
        deljobs = false;
    end
    if length( divisions ) == 1
        divisionorder = divisions;
        subfile = 0;
    else
        divisionorder = divisions( 1 );
        subfile = divisions( 2 );
        firstjobidx = length( jobs ) / divisions( 3 ) * ( subfile - 1 ) + 1;
        lastjobidx = length( jobs ) / divisions( 3 ) * subfile;
        jobs = jobs( firstjobidx : 1 : lastjobidx );
    end
    %% run through Ns, dimensions and shapes:
    for N = Ns
        for d = dims
            for shape = shapes
                fprintf( 'Combining files for the parameters N%dD%d%s:\n', ...
                    N, d, lower( char( shape ) ) );
                runjobs_progress_initialize();
                runjobs_progress_update();
                %% load previous master data:
                if ~subfile
                    masterfilename = sprintf( 'data/N%dD%d%s.mat', ...
                        N, d, lower( char( shape ) ) );
                else
                    masterfilename = sprintf( 'data/N%dD%d%sJ%dx%d.mat', ...
                        N, d, lower( char( shape ) ), divisionorder, subfile );
                end
                if ( exist( masterfilename, 'file' ) > 0 ) && extend
                    masterdata = load( masterfilename );
                else
                    masterdata = struct();
                end
                %% load job data file and add to master data:
                jobcounter = 0;
                jobsnumber = length( jobs );
                jobsexist = zeros( size( jobs ) );
                jobsincompatible = zeros( size( jobs ) );
                incompatiblefield = '';
                for j = jobs
                    jobcounter = jobcounter + 1;
                    if divisionorder > 1
                        filename = sprintf( 'data/N%dD%d%sJ%dx%d.mat', ...
                            N, d, lower( char( shape ) ), divisionorder - 1, j );
                    else
                        filename = sprintf( 'data/N%dD%d%sJ%d.mat', ...
                            N, d, lower( char( shape ) ), j );
                    end
                    jobsexist( jobcounter ) = exist( filename, 'file' );
                    if jobsexist( jobcounter )
                        jobdata = load( filename );
                        for fcell = transpose( fieldnames( jobdata ) )
                            fieldname = char( fcell );
                            jobvalue = jobdata.(fieldname);
                            if ~isfield( masterdata, fieldname )
                                % fields to copy:
                                if strcmp( fieldname, 'shape' )
                                    masterdata.(fieldname) = char( jobvalue );
                                else
                                    masterdata.(fieldname) = jobvalue;
                                end
                            elseif sum( strcmp( fieldname, { 'runs', ...
                                    'futureinfinities', 'eventcounts', ...
                                    'chains', 'simplices', 'dimestimators', ...
                                    'diamonds', 'diamondtimes', 'propertimes', ...
                                    'hyperbdistribution', 'comptime' } ) ) == 1
                                % fields to add up:
                                mastervalue = masterdata.(fieldname);
                                if size( jobvalue ) == size( mastervalue )
                                    mastervalue = mastervalue + jobvalue;
                                    masterdata.(fieldname) = mastervalue;
                                else
                                    jobsincompatible( jobcounter ) = 1;
                                    incompatiblefield = fieldname;
                                end
                            elseif sum( strcmp( fieldname, { ...
                                    'distribution', 'preffutures' } ) ) == 1
                                % field with subfields:
                                for subfcell = transpose( fieldnames( jobdata.(fieldname) ) )
                                    subfieldname = char( subfcell );
                                    jobvalue = jobdata.(fieldname).(subfieldname);
                                    if ~isfield( masterdata.(fieldname), subfieldname )
                                        masterdata.(fieldname).(subfieldname) = jobvalue;
                                    elseif sum( strcmp( subfieldname, ...
                                        { 'Ns', 'bins', 'counts', 'diamonds', ...
                                          'dimestimators', 'propertimes', ...
                                          'hyperbdistribution' } ) ) == 1
                                        % subfields to add up:
                                        mastervalue = masterdata.(fieldname).(subfieldname);
                                        if size( jobvalue ) == size( mastervalue )
                                            mastervalue = mastervalue + jobvalue;
                                        else
                                            jobsincompatible( jobcounter ) = 1;
                                            incompatiblefield = ...
                                                sprintf( '%s.%s', fieldname, subfieldname );
                                        end
                                        masterdata.(fieldname).(subfieldname) = mastervalue;
                                    elseif strcmp( subfieldname, 'unithyperboloid' )
                                        % subfields to extend:
                                        mastervalue = masterdata.(fieldname).(subfieldname);
                                        if size( jobvalue, 2 ) == size( mastervalue, 2 )
                                            mastervalue = [ mastervalue; jobvalue ]; %#ok<AGROW>
                                        else
                                            jobsincompatible( jobcounter ) = 1;
                                            incompatiblefield = ...
                                                sprintf( '%s.%s', fieldname, subfieldname );
                                        end
                                        masterdata.(fieldname).(subfieldname) = mastervalue;
                                    end
                                end
                            end
                        end
                        if deljobs
                            delete( filename );
                        end
                    end
                    runjobs_progress_update( jobcounter / jobsnumber );
                end
                runjobs_progress_finalize();
                %% if any jobs had been processed:
                if sum( jobsexist ) > 0
                    %% set default max. dimension:
                    if ~isfield( masterdata, 'maxdimension' )
                        masterdata.maxdimension = 8;
                    end
                    %% set default runs:
                    if ~isfield( masterdata, 'runs' )
                        masterdata.runs = 1;
                    end
                    %% save to combined file:
                    masterdata.comptimestr = ...
                        char( seconds( masterdata.comptime ), 'hh:mm:ss' ); 
                    save( masterfilename, '-struct', 'masterdata' );
                end
                %% show list of missing jobs:
                missingjobs = jobs( jobsexist == 0 );
                missingjobcount = length( missingjobs );
                if missingjobcount > 0
                    fprintf( 'N%dD%d%s had %d / %d missing job files: ', ...
                        N, d, lower( char( shape ) ), missingjobcount, length( jobs ) );
                    jobcounter = 0;
                    for j = missingjobs
                        jobcounter = jobcounter + 1;
                        if jobcounter < missingjobcount
                            fprintf( '%d, ', j );
                        else
                            fprintf( '%d.\n', j );
                        end
                    end
                else
                    fprintf( 'N%dD%d%s had %d job files.\n', ...
                        N, d, lower( char( shape ) ), length( jobs ) );
                end
                %% show list of jobs with incompatible fields:
                incompatiblejobs = jobs( jobsincompatible == 1 );
                incompatiblejobcount = length( incompatiblejobs );
                if incompatiblejobcount > 0
                    fprintf( 'N%dD%d%s had %d / %d job files with incompatible sizes for some fields (last error for .%s): ', ...
                        N, d, lower( char( shape ) ), incompatiblejobcount, length( jobs ), incompatiblefield );
                    jobcounter = 0;
                    for j = incompatiblejobcount
                        jobcounter = jobcounter + 1;
                        if jobcounter < incompatiblejobcount
                            fprintf( '%d, ', j );
                        else
                            fprintf( '%d.\n', j );
                        end
                    end
                end
            end
        end
    end
end

