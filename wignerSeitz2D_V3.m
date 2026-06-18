function [P,Area] = wignerSeitz2D_V3(Bragg,center,opts)
arguments
    Bragg (:,2) double
    center (:,2) double = [0 0] % center point
    opts.Method {mustBeMember(opts.Method,{'voronoi','Intersection'})} = 'voronoi'
    opts.Output {mustBeMember(opts.Output,{'struct','verticesSorted'})} = 'verticesSorted'
end
% Given a set of 2D Bragg peaks, wignerSeitz2D_V3(Bragg,center) calculates
%  the vertices of a WignerSeitz (Brillouin in reciprocal space) cell around 
% center (or the origin if none is given), i.e. the area enclosed by these 
% Bragg peaks whose points lie closer to the origin than to any of the peaks.

% INPUT
% Bragg is a 2 column matrix with coordinates of the peaks
% center is a 2 column matrix, with as many rows as cells you want
% 'Method' = 'voronoi','Intersection' determines algorithm. By default, set
% to 'voronoi'
% 'Output' = 'verticesSorted','struct' determines type of output variable.
% By default, it is 'verticesSorted', where a matrix (cell if more than 1
% center) is given. 'Struct' would offer a struct, so that patch(P) already
% plots the Wigner-Seitz cell.

%Example: BrillouinPeaks = wignerSeitz2D_V3([0 1;1 0]).

% Zero = mean(Bragg,1); %Center of polygon as center of mass
% Bragg = Bragg - Zero;

for nCenter = 1:size(center,1)
    ctr = center(nCenter,:);

    % Sort peaks by angle
    a = atan2d(Bragg(:,2),Bragg(:,1));
    [~,idx] = sort(a);
    Bragg = Bragg(idx,:);

switch opts.Method
    case 'Intersection'
        Bragg = Bragg - ctr; % move points around zero

        normBragg = 0.5*vecnorm(Bragg,2,2).^2;
        % Number of peaks
        numPeaks = numel(normBragg);

        % Solve linear system for every adjacent pair of peaks
        Pn = zeros(size(Bragg));
    for i = 1:numPeaks
        j = mod(i,numPeaks)+1;
    
        b = [normBragg(i); normBragg(j)];
        A = [Bragg(i,1), Bragg(i,2);
            Bragg(j,1), Bragg(j,2)];
    
        Pn(i,:) = (A\b).';
    end
P{nCenter} = Pn + ctr; % return polygon to original position

    case 'voronoi'
        % Triangulate input points
            dt = delaunayTriangulation(Bragg);
            [verts,region] = voronoiDiagram(dt);
            % Obtain lattice point nearest to input center, in case it doeasn't match
            tid = nearestNeighbor(dt,ctr(1),ctr(2));

            % Calculate vertices of Voronoi cell around center
            P{nCenter} = uniquetol(verts(region{tid},:),1e-10,'ByRows',true);
            % Sort vertices by angle in 2D
            a = atan2d(P{nCenter}(:,2),P{nCenter}(:,1));
            [~,idx] = sort(a);
            P{nCenter} = P{nCenter}(idx,:);
end
% Compute cell area defined by P
    [~,Area(nCenter)] = convexHull(delaunayTriangulation(P{nCenter}));
end

if contains(opts.Output,'struct')
    C = tab10(nCenter);
    for n = 1:nCenter
        outStruct(n).vertices = P{n};
        outStruct(n).faces = 1:size(P{n},1); % face connectivity after sorting
        outStruct(n).faceColor = 'none'; % assign different colors

        outStruct(n).UserData.Area = Area(n);% save the cell area
    end
    % Replace first output by struct
    P = outStruct;

% simplify result if only 1 center
elseif nCenter == 1
    P = P{1};

end

end