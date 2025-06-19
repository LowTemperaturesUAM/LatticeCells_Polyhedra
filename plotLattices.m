function plotLattices(latticeData,axReal,axRecip)

% Given a struct with lattice data and 2 axes, plots Real and reciprocal
% lattices with their cells.

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
cla(axReal)
hold(axReal,'on')
scatter3(axReal,fullR(:,1),fullR(:,2),fullR(:,3),50,'r','filled');
patch(axReal,'Vertices',Pr,'faces',Fr,'FaceColor','none','Linewidth',2)
daspect([1 1 1])
xlim([-0.5 1.5].*max(Pr(:,1)))
ylim([-.5 1.5].*max(Pr(:,2)))
zlim([-.5 1.5].*max(Pr(:,3)))
% xlim([-1 2].*max(R(:,1)))
% ylim([-1 2].*max(R(:,2)))
% zlim([-1 2].*max(R(:,3)))
view([110 10])
hold(axReal,'off')

% Clear and plot real lattice
cla(axRecip)
hold(axRecip,'on')
scatter3(axRecip,fullK(:,1),fullK(:,2),fullK(:,3),50,'b','filled');
patch(axRecip,'Vertices',Pk,'faces',Fk,'FaceColor','none','Linewidth',2)
daspect([1 1 1])
xlim([-1 1].*max(Pk(:,1)))
ylim([-1 1].*max(Pk(:,2)))
zlim([-1 1].*max(Pk(:,3)))
view([110 10])
hold(axRecip,'off')


end