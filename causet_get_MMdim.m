function MMdim = causet_get_MMdim( chains, newaccuracy )
%CAUSET_GET_MMDIM returns the Myrheim-Meyer estimator for the spacetime
% dimension from the given k-chain averages CHAINS (row vector). 
% 
% Arguments:
% CHAINS              2-element row vector to get original flat spacetime 
%                     dimension estimator. 
%                     Solution to equation (11) in [1]. 
%                     4-element row vector to get modified curved spacetime
%                     dimension estimator. 
%                     Solution to equation (59) in [1].
%
% Optional arguments:
% NEWACCURACY         sets the precision for the return value (number 
%                     of digits), the default is 3.
% 
% Returns:
% MMDIM               (lowest) minimum of the estimation, which is referred
%                     to as the Mzrheim-Meyer dimension.
% 
% References: 
% [1] author = {Roy, Mriganko and Sinha, Debdeep and Surya, Sumati},
%     journal = {Physical Review D},
%     number = {4},
%     pages = {044046},
%     title = {{The Discrete Geometry of a Small Causal Diamond}},
%     volume = {87},
%     year = {2013}
% 
    
    %% set default values:
    if nargin < 2
        newaccuracy = 3;
    end
    %% generate coefficients and keep in local memory:
    persistent accuracy dim coeffflat coeffcurved
    if isempty( dim ) || accuracy ~= newaccuracy 
        accuracy = newaccuracy;
        dimlen = ( 10 - 0 ) * 10^accuracy + 1;
        dim = linspace( 0, 10, dimlen );
        % coeffflat is 1 / chi2:
        coeffflat = gamma( dim + 1 ) .* gamma( dim / 2 ) ./ 4 ./ gamma( 1.5 * dim );
        coeffcurved = zeros( 4, dimlen );
        coeffcurved( 1, : ) = ( dim + 2 ) .* ( 2 * dim + 2 );
        coeffcurved( 2, : ) = -3 * ( 2 * dim + 2 ) .* ( 3 * dim + 2 ) ...
            ./ coeffflat.^2;
        coeffcurved( 3, : ) = 3 * ( 3 * dim + 2 ) .* ( 4 * dim + 2 ) ...
            ./ ( 1/12 * dim.^2 .* gamma( dim / 2 ) .* gamma( dim ).^3 ...
                ./ gamma( 1.5 .* dim ) ./ gamma( 2 * dim ) ).^( 4 / 3 );
        coeffcurved( 4, : ) = - ( 4 * dim + 2 ) .* ( 5 * dim + 2 ) ...
            ./ ( 1/32 * dim.^2 .* gamma( dim / 2 ) .* gamma( dim ).^4 ...
                ./ gamma( 2.5 .* dim ) ./ gamma( 2 * dim ) );
    end
    %% compute estimators for the dimension values:
    if length( chains ) == 2 % flat estimator
        estimators = coeffflat .* chains( 1 )^2 - chains( 2 );
    elseif length( chains ) == 4 % curved estimator
        % The coefficient curves are logarithmically increasing such that
        % a step by step addition is necesary.
        estimators = coeffcurved( 1, : ) .* chains( 1 )^4;
        estimators = estimators + coeffcurved( 2, : ) .* chains( 2 )^2;
        estimators = estimators + coeffcurved( 3, : ) .* chains( 3 )^( 4 / 3 );
        estimators = estimators + coeffcurved( 4, : ) .* chains( 4 );
    else
        estimators = 0;
    end
    %% return first result:
    estimators = abs( estimators ); % closest to zero
    MMdim = dim( estimators == min( estimators ) );
    MMdim = MMdim( 1 );
end

