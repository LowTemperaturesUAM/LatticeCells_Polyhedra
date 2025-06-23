function [V_voronoi,P,Vol] = plotVoronoiCell3D(points,center,ax)
arguments
    points (:,3) double
    center (:,3) double = [0 0 0];
    ax = [];
end

points = uniquetol(points,1e-10,'ByRows',true);
if isempty(ax) % In case there were not axes given
    figure;
    ax = axes;
    scatter3(ax,points(:,1),points(:,2),points(:,3),20, ...
        'or','filled','MarkerEdgeColor','k');
end
hold(ax,'on')

dt = delaunayTriangulation(points);
[V, Reg] = voronoiDiagram(dt);

% Initialize vertices, path array and volumes
V_voronoi = cell(size(center,1),1);
P = gobjects(size(V_voronoi));
Vol = nan(size(V_voronoi));

for n = 1:size(center,1)
% Obtain lattice point nearest to center
tid = nearestNeighbor(dt,center(n,1),center(n,2),center(n,3));
% Calculate vertices of voronoi cell around the center
V_voronoi{n} = uniquetol(V(Reg{tid},:),1e-6,'ByRows',true);
[F_voronoi,Vol(n)] = facesPatch3D(V_voronoi{n});

% Plot polyhedra on axes
    
    P(n) = patch(ax,'vertices',V_voronoi{n},'faces',F_voronoi, ...
        'faceAlpha',.2,'FaceColor',[.1 .5 .7] );
end
    daspect([1 1 1])
    view([10 2])
    xlabel('X')
    ylabel('Y')
    hold(ax,'off')

end