% parameter sets to get geodesics with an expected length of 10 diamonds:
%                                         length  starts    ends
d = 2; N = 80; t = 0.035 * N^(1/d); %         11       8       8
d = 3; N = 440; t = 0.035 * N^(1/d); %        11       8       8
d = 4; N = 2420; t = 0.035 * N^(1/d); %       11       8       8
d = 5; N = 13310; t = 0.035 * N^(1/d); %      11       8       8
d = 6; N = 73205; t = 0.035 * N^(1/d);

% test simulations:
runs = 5;
average = 0;
fails = 0;
averageA = 0;
averageB = 0;
for i = 1 : runs
    [ coordinates, coordinateranges, volume, events, geoendevents ] = ...
        causet_new_sprinkle(N, d, 'bicone', [1, -1, 0, t]);
    C = causet_edit_relate(coordinates, metric_flat(d)); 
    L = causet_get_links(C);
    A = find( geoendevents(:, 1) );
    B = find( geoendevents(:, 2) );
    if isempty(A) || isempty(B)
        fails = fails + 1;
    else
        averageA = averageA + length(A);
        averageB = averageB + length(B);
        a = randi(length(A), 1);
        b = randi(length(B), 1);
        [ geodesics, cardinality ] = ...
            causet_find_linkgeodesics(C, L, A(a), B(b));
        average = average + cardinality;
    end
end
average = average / (runs - fails);
averageA = averageA / (runs - fails);
averageB = averageB / (runs - fails);

% plot:
scatter(coordinates(geoendevents(:, 1), 2), ...
        coordinates(geoendevents(:, 1), 1), 'ob');
hold on;
scatter(coordinates(geoendevents(:, 2), 2), ...
        coordinates(geoendevents(:, 2), 1), 'or');
scatter(coordinates(geodesics{1}, 2), ...
        coordinates(geodesics{1}, 1), '.k');
title(sprintf('length: %.2f, bottoms: %.2f, tops: %.2f, fails: %d / %d', ...
    average, averageA, averageB, fails, runs));
xlim([-1, 1]);
ylim([-1, 1]);
hold off;
