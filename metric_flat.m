function g = metric_flat( d )
%METRIC returns the Minkowski metric in D dimensions as matrix. 
    
    metrictime = zeros( d );
    metrictime( 1, 1 ) = 1;
    g = 2 * metrictime - eye( d );
end

