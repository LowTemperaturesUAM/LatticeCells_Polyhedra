function [Vertices,Faces,Vol] = wignerSeitz3D(Points,center,opts)
arguments
    Points (:,:) double
    center (:,:) double = [0 0 0]; % center point
    opts.Method {mustBeMember(opts.Method,{'voronoi','planeIntersection'})} = 'voronoi'
    opts.Output {mustBeMember(opts.Output,{'struct','verticesAndFaces'})} = 'verticesAndFaces'
end
% Given a set of points in 3D space around a center points, gives the vertices
% of the Wigner-Seitz Cell. Second output is a Face matrix, indicating 
% which vertex is in which face, used to plot with patch.
% EXAMPLE: [Vertices,Faces] = wignerSeitz3D([0 0 1;0 1 0;1 0 0;
% 0 0 -1;0 -1 0;-1 0 0], [0 0 0])

% remove duplicates
Points = uniquetol(Points,1e-10,"ByRows",true);
% Initialize vertices, faces and volumes
Vertices = cell(size(center,1),1);
Faces = cell(size(center,1),1);
Vol = nan(size(center,1),1);

flag2D = false;

% If there are many lattice points, it is faster to triangulate. If 2D,
% avoid Intersection
if length(Points) > 50 || size(Points,2)<3
    opts.Method = 'voronoi';
    warning('Switching calculation method to Voronoi teselation')
end

for nCenter = 1:size(center,1)
    ctr = center(nCenter,:);

    switch opts.Method
        case 'planeIntersection'
            sol = [];
            vertexPlanes = [];
            % Every combination of 3 points
            combin = nchoosek(1:size(Points,1),3);
            % Calculate every coefficient for the plane equation: n(x-m) = 0
            n = Points-ctr; % Normal vector
            m = (Points+ctr)/2; % Midpoint
            d = -sum(n.*m,2); % independent factor - position of plane

            for i = 1:size(combin,1) % loop over all 3point combinations
                P = Points(combin(i,:),:);
                N = P-ctr;
                M = (P+ctr)/2;
                D = -sum(N.*M,2);

                %MatrixCoefs = n; if determinant is non-zero, there is a single
                %solution
                if abs(det(N)) > 1e-6
                    solved = N\(-D);
                    sol = [sol; solved'];
                    vertexPlanes = [vertexPlanes; combin(i,:)];
                else
                end
            end
            % Leave only those that satisfy ALL inequations
            valCond = any(sol * n' + d' >eps,2);
            sol(valCond,:) = [];
            vertexPlanes(valCond,:) = [];

            Vertices{nCenter} = uniquetol(sol,1e-10,'ByRows',true);

        case 'voronoi' % Use Voronoi teselation from triangulation
            % IF POINTS ARE COPLANAR (null space of 1D)
            % tempPoints = Points(vecnorm(Points,2,2)~=0,:);
            relatPoints = Points-Points(1,:);

            if size(null(relatPoints,1e-4),2) == 1
                flag2D = true;
                % Project points to 2D beforehand, get the base coordinate,
                % and return to 3D after calculation
                tempPoints = Points; % save 3D for later
                warning("Coplanar points. Switched to 2D for Voronoi calulations")

                % Ensure that vectors for calculating normal are not
                % parallel
                idPoint = 2;
                while size(null(relatPoints([2 idPoint],:),1e-4),2) ~= 1
                    idPoint = idPoint + 1;
                end
              % normal vector
                n = cross(relatPoints(2,:),relatPoints(idPoint,:));
                n = n/vecnorm(n);
                disp("Normal vector: "+mat2str(round(n,3)))

                [~,Points,newBase] = projectToPlane(Points,n,Points(1,:));
            end
            % Triangulate input points
            dt = delaunayTriangulation(Points);

            [verts,region] = voronoiDiagram(dt);
            % Obtain lattice point nearest to input center, in case it doeasn't match
            if ~flag2D
                tid = nearestNeighbor(dt,ctr(1),ctr(2),ctr(3));
            else
                ctr = ctr*newBase.vectors';
                tid = nearestNeighbor(dt,ctr(1),ctr(2));
            end

            % Calculate vertices of Voronoi cell around center
            Vertices{nCenter} = uniquetol(verts(region{tid},:),1e-10, ...
                'ByRows',true);

            % Return to 3D
            if size(null(relatPoints,1e-4),2) == 1
                % Sort by angle
                new0 = mean(Vertices{nCenter}); % Center points around mean
                a = atan2d(Vertices{nCenter}(:,2)-new0(2), ...
                    Vertices{nCenter}(:,1)-new0(1));
                [~,idx] = sort(a);
                Vertices{nCenter} = Vertices{nCenter}(idx,:);

                Vertices{nCenter} = Vertices{nCenter}*newBase.vectors...
                    + newBase.offset;
                Faces{nCenter} = 1:size(Vertices{nCenter},1);
            end
    end
    % Compute Connectivity for faces and volume of cell
    if ~flag2D
        [Faces{nCenter}, Vol(nCenter)] = facesPatch3D(Vertices{nCenter});
    end

end

if nCenter == 1 % Remove cell structure if there is only 1
    Vertices = Vertices{1};
    Faces = Faces{1};
end

    if contains(opts.Output,'struct')
        C = tab10(nCenter);
        for n = 1:nCenter
            outStruct(n).vertices = Vertices{n};
            outStruct(n).faces = Faces{n};
            outStruct(n).faceColor = 'none'; % assign different colors

            outStruct(n).UserData.rAtom = center(n,:);% save the atom position
        end
        % Replace first output by struct, and the second by the volumes
        Vertices = outStruct;
        Faces = Vol;
        Vol = [];
    end

end


