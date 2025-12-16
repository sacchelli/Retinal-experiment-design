clc

%rng(10)

d = 6;

points = randn(10000,d);

[minVals, indexMin] = min(points, [], 1);
[maxVals, indexMax] = max(points, [], 1);


hullPoints = [points(indexMax,:); points(indexMin,:)];


K = convhulln(hullPoints);

% Build Delaunay triangulation from hull points
DT = delaunayn(hullPoints);

tot=2000;
compteur=0;

for k =1:tot
% Example point to test
    pt = randn(1,d);
    %pt = [0 0 0 0 0 0];
    
    % Check if pt is inside the convex hull
    id = tsearchn(hullPoints, DT, pt);
    
    if ~isnan(id)
        compteur=compteur+1;
        %disp('Point is inside the hull')
    else
        %disp('Point is outside the hull')
    end

end
display(100*compteur/tot)
