function plotLattices(latticeData,axReal,axRecip)
arguments
    latticeData
    axReal
    axRecip
end

% Given a struct with lattice data and 2 axes, plots Real and reciprocal
% lattices with their cells.

% ADD OPTIONAL INPUT FOR CELL, VECTORS AND UNITS


%---------------------------------------------------------------------
% get data from struct
fullR = latticeData.Real;
fullK = latticeData.Reciprocal;
Pr = latticeData.realCellVertices;
Pk = latticeData.reciprCellVertices;
Fr = latticeData.realCellFaces;
Fk = latticeData.reciprCellFaces;

if nargin<2
    figure;
    tiledlayout(1,2)
    axReal = nexttile;
    axRecip = nexttile;
end

% Clear and plot real lattice
vectorR = latticeData.realLatticeVector;
cla(axReal)
hold(axReal,'on')
% box(axReal,'on')
% Plot lattice
scatter3(axReal,fullR(:,1),fullR(:,2),fullR(:,3),50,'r','filled');
% Plot cell with conventional vectors
patch(axReal,'Vertices',Pr,'faces',Fr,'FaceColor','none','Linewidth',2)
% xlim(axReal,[-0.5 1.5].*max(Pr(:,1)))
% ylim(axReal,[-.5 1.5].*max(Pr(:,2)))
% zlim(axReal,[-.5 1.5].*max(Pr(:,3)))
xlim(axReal,[-1 2].*max([vectorR(:,1);abs(Pr(:,1))]))
ylim(axReal,[-1 2].*max(vectorR(:,2)))
zlim(axReal,[-0 2].*max(vectorR(:,3)))
% view(axReal,[110 10])
xlabel(axReal,'X')
ylabel(axReal,'Y')
zlabel(axReal,'Z')
hold(axReal,'off')

% Clear and plot reciprocal lattice
vectorK = latticeData.reciprLatticeVector;
cla(axRecip)
hold(axRecip,'on')
% box(axRecip,'on')
scatter3(axRecip,fullK(:,1),fullK(:,2),fullK(:,3),50,'b','filled');
patch(axRecip,'Vertices',Pk,'faces',Fk,'FaceColor','none','Linewidth',2)
% xlim(axRecip,[-1 1].*max(Pk(:,1)))
% ylim(axRecip,[-1 1].*max(Pk(:,2)))
% zlim(axRecip,[-1 1].*max(Pk(:,3)))
xlim(axRecip,[-1 1].*max(vectorK(:,1)))
ylim(axRecip,[-1 1].*max(vectorK(:,2)))
zlim(axRecip,[-1 1].*max(vectorK(:,3)))
% view(axRecip,[110 10])
xlabel(axRecip,'k_x')
ylabel(axRecip,'k_y')
zlabel(axRecip,'k_z')
hold(axRecip,'off')

set([axReal,axRecip],'DataAspectRatio',[1 1 1],'View',[110,10])
end