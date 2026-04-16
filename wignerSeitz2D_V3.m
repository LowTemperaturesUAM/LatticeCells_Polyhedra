function [P] = wignerSeitz2D_V3(Bragg,center,opts)
arguments
    Bragg (:,2) double
    center (:,2) double = [0 0];
    opts.Method {mustBeMember(opts.Method,{'voronoi','Intersection'})} = 'voronoi'
    opts.Output {mustBeMember(opts.Output,{'struct','verticesSorted'})} = 'verticesSorted'

end

% Given a set of 2D Bragg peaks, wignerSeitz2D:V3(Bragg,center) calculates
%  the vertices of a WignerSeitz (Brillouin in reciprocal space) cell around 
% center (or the origin if none is given), i.e. the area enclosed by these 
% Bragg peaks whose points lie closer to the origin than to any of the peaks.

% INPUT
% Bragg is a 2 column matrix with coordinates of the peaks
% center is a 2 column matrix, with as many rows as cells you want
% 'Method' = 'voronoi','Intersection' determines algorithm. By default, set
% to 'voronoi'
% 'Output' = 'verticesSorted','struct' determines type of output variable.
% By default, it is 'vertocesSorted', where a matrix (cell if more than 1
% center) is given. 'Struct' would offer a struct, so that patch(P) already
% plots the Wigner-Seitz cell.

%Example: BrillouinPeaks = wignerSeitz_V3([0 1;1 0]).

switch opts.Method
    case 'Intersection' % Linear system -------------------------
        Zero = mean(Bragg,1); %Center of polygon as center of mass
        Bragg = Bragg - Zero;
        % Npeaks = size(Bragg,1);

        % Sort peaks by angle
        a = atan2d(Bragg(:,2),Bragg(:,1));
        [~,idx] = sort(a);
        Bragg = Bragg(idx,:);

        normBragg = 0.5*vecnorm(Bragg,2,2).^2;
        % Number of peaks
        numPeaks = numel(normBragg);

        % Solve linear system for every adjacent pair of peaks
        P = zeros(size(Bragg));
        for i = 1:numPeaks
            j = mod(i,numPeaks)+1;

            b = [normBragg(i); normBragg(j)];
            A = [Bragg(i,1), Bragg(i,2);
                Bragg(j,1), Bragg(j,2)];

            P(i,:) = (A\b).';
        end
        % Return polygon to original position
        P = P + Zero;

    case 'voronoi' % Voronoi teselation--------------------------------
        % Initialize variables
        P = cell(size(center,1),1);

        % get output for every center provided
        for nCenter = 1:size(center,1)
            ctr = center(nCenter,:);

            % Triangulate input points
            dt = delaunayTriangulation(Bragg);
            [verts,region] = voronoiDiagram(dt);

            % Obtain lattice point nearest to input center, in case it doeasn't match
            tid = nearestNeighbor(dt,ctr(1),ctr(2));

            % Calculate vertices of Voronoi cell around center
            tempVertices = ...
                uniquetol(verts(region{tid},:),1e-10,'ByRows',true);

            %Sort vertices by angle, so that connectivity is index vector
            a = atan2d(tempVertices(:,2),tempVertices(:,1));
            [~,idx] = sort(a);
            tempVertices = tempVertices(idx,:);

            P{nCenter} = tempVertices;
        end

        if contains(opts.Output,'struct')
            % C = tab10(nCenter);
            for n = 1:nCenter
                outStruct(n).vertices = P{n};
                outStruct(n).faces = 1:length(P{n});
                outStruct(n).faceColor = 'none'; % assign different colors
            end
            % Replace first output by struct
            P = outStruct;
        end

        % if there is only 1 center, use matrix instead of cells
        if nCenter == 1
            P = cell2mat(P);
        end

end
end