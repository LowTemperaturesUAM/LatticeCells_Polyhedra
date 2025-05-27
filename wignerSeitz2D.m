function [P] = wignerSeitz2D(Bragg)
arguments
    Bragg (:,2) double
end

% Given a set of Bragg peaks, wignerSeitz2D(Bragg) calculates the vertices of a WignerSeitz
% (Brillouin in reciprocal space) cell around the origin, i.e. the area 
% enclosed by these Bragg peaks whose points lie closer to the origin than
% to any of the peaks.

% INPUT
% Bragg is a 2 column matrix with coordinates of the peaks

%% Rewrite using voronoi

Zero = mean(Bragg,1); %Center of polygon as center of mass
% Zero = [0 0];
Bragg = Bragg - Zero;
Npeaks = size(Bragg,1);

% Sort peaks by angle
a = atan2d(Bragg(:,2),Bragg(:,1));
 [~,idx] = sort(a);
 Bragg = Bragg(idx,:);

% Normal vectors to All peaks
v = cell2mat(rowfun(@null,table(Bragg),"OutputFormat","cell")')';
M = Bragg/2; % Middle point between Origin and Peak

% Solve linear system for every adjacent pair
P = zeros(size(Bragg));
for i = 1:Npeaks
    k = mod(i,Npeaks)+1;

    A = [v(i,2), -v(i,1) ;  v(k,2), -v(k,1)];
    b = [v(i,2)*M(i,1)-v(i,1)*M(i,2); v(k,2)*M(k,1)-v(k,1)*M(k,2)];

    P(i,:) = (A\b)';
end

% Return polygon to original position
P = P + Zero;

end