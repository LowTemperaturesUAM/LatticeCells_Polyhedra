function CrystalPolys(crystal, C, axPoly)

%This fucntion CrystalPolys calculates and paints the polyhedras of the
%atoms in a cell and does nothing more than that
n = [0 0 1 ] % Miller indices
 

% START ---------------------------------------------------------------
% Check if there is incorrect number of colors given
% If more colors are given, the next one is used for Void polyhedra
C = [0.2 0.2 0.2];
numC = numel(crystal.Data.atom_site_fract_x);
if numC > size(C,1)
    C = tab10(numC+1);
    voidColor = C(end,:);
elseif numC < size(C,1)
    voidColor = C(end,:);
else
    voidColor = [.4 .7 .8];
end
C = C(1:numC,:);
clearvars numC
mvId = [0 0 0;n(1:3);1 0 0;0 1 0];



% Create lattice vectors
LatR = latticeVectors(crystal,"real","conventional");
% create convenional cell patch
[cellP] = cellPatch(crystal);

% Load 1st cell
[~,R0] = breakSymm(crystal,1,'IncludeBoundary',true,'NormalizedPosition',false);
% Decide atoms I want to paint
Rpaint = R0;
Rpaint = cellfun(@(x) x.*[-1 -1 -1]+sum(LatR),R0,'UniformOutput',false);
%Rpaint = cellfun(@(x) x.*[1 1 1]+sum(LatR),R0,'UniformOutput',false);
R0=Rpaint;
% ALL atoms, for Voronoi. Use 3 cells to each side
Rall = cellfun(@(x) replicateCell(x,1,LatR),R0, ...
    'UniformOutput',false);

% Create polyhedra of R0 atoms
Fpoly = cellfun(@(R0) wignerSeitz3D(vertcat(Rall{:}),R0, ...
    "Method","voronoi","Output","struct"),R0,'UniformOutput',false);
for i = 1:numel(R0) % assign colors to the atoms
    tempColor = repmat({C(i,:)},size(Fpoly{i}));
    [Fpoly{i}.faceColor] = tempColor{:};
end

Fpoly = [Fpoly{:}];

% PLOT--------------------------------------------------------------------
% figure();


cla(axPoly)
hold(axPoly,'on')

set(axPoly,'DataAspectRatio',[1 1 1],'view',[47 17])
xlabel('X (Å)'); ylabel('Y (Å)');zlabel('Z (Å)')
% zlim([0 crystal.lengths(3)*1.1])

% Plot cell countours
arrayfun(@(x) patch(axPoly,x,'facecolor','none','Linewidth',3),cellP);


% plot atoms
for i = 1:length(Rpaint)
    scatter3(axPoly,Rpaint{i}(:,1), Rpaint{i}(:,2), Rpaint{i}(:,3),100, ...
        C(i,:).*ones(size(Rpaint{i},1),1),'filled','MarkerEdgeColor','k');
end
% plot polyhedra
% arrayfun(@(x) patch(axPoly,x,'facealpha',.3),Spoly);
arrayfun(@(x) patch(axPoly,x,'facealpha',.3),Fpoly);

% Plot CMs
% if ~isempty(ogCM)&& ~isempty(mvdCM)
%     scatter3(axPoly,ogCM(:,1), ogCM(:,2), ogCM(:,3),180, ...
%         [0 .8 0].*ones(size(ogCM,1),1),'filled','^', ...
%         'MarkerEdgeColor','k');
%     scatter3(axPoly,mvdCM(:,1), mvdCM(:,2), mvdCM(:,3),180, ...
%         [0 .8 .8].*ones(size(ogCM,1),1),'filled','v','MarkerEdgeColor','k');
% end

end