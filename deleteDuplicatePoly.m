function newFPoly = deleteDuplicatePoly(newFPoly)
% Removes duplicate polyhedrons from a patch array
% EXAMPLE: newFPoly = deleteDuplicatePoly(newFPoly)
%
%

if isfield(newFPoly(1),"faceColor")
 C = unique(vertcat(newFPoly.faceColor),'rows','stable');
else
   C = [.2 .2 .2];
   warning('It will be considered as all polyhedra of same color')
end

% Select only unique patches from atomic subset
% save index to delete 
tempPoly = [];
for jj = 1:size(C,1)
    idDelete = [];
    if isfield(newFPoly(1),"faceColor")
    subPoly = newFPoly(all(vertcat(newFPoly.faceColor)==C(jj,:),2));
    else
        subPoly = newFPoly;
    end

% Loop over all the subset to detect duplicates avoid repeating indices 
% or combinations (upper triangular matrix of combinations)
    for subId = 1:numel(subPoly)
        for subId2 = subId+1:numel(subPoly)
           if all(sortrows(round(subPoly(subId).vertices,5),[3 2 1])- ...
           sortrows(round(subPoly(subId2).vertices,5),[3 2 1])==0,'all')
                idDelete(end+1) = subId;
            end
        end
    end
    subPoly(idDelete) = [];
    tempPoly = [tempPoly, subPoly];
end
newFPoly= tempPoly;
end