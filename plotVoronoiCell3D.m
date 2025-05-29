function [V_voronoi,P] = plotVoronoiCell3D(points,center,ax)
arguments
    points (:,3) double
    center (1,3) double = [0 0 0];
    ax = [];
end

points = uniquetol(points,1e-10,'ByRows',true);
dt = delaunayTriangulation(points);
[V, Reg] = voronoiDiagram(dt);
% Obtain lattice point nearest to center
tid = nearestNeighbor(dt,center(1),center(2),center(3));
% Calculate vertices of voronoi cell around the center
V_voronoi = uniquetol(V(Reg{tid},:),1e-6,'ByRows',true);
[F_voronoi,Vol] = facesPatch3D(V_voronoi);

if isempty(ax) % In case there were not axes given
    figure;
    ax = axes;
    scatter3(ax,points(:,1),points(:,2),points(:,3),20, ...
        'or','filled','MarkerEdgeColor','k');
end
    hold(ax,'on')
    P = patch(ax,'vertices',V_voronoi,'faces',F_voronoi, ...
        'faceAlpha',.2,'FaceColor',[.1 .5 .7] );

    daspect([1 1 1])
    view([10 2])
    xlabel('X')
    ylabel('Y')
    hold(ax,'off')

end