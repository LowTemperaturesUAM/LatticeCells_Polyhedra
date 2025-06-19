function cifData = readcif(path)

% Reads the .cif file in path and sves it as a struct. 

% Import function from
% Sandor Toth (2025). Class handling .cif formatted files 
% (https://www.mathworks.com/matlabcentral/fileexchange/43266-class-handling-cif-formatted-files),
%  MATLAB Central File Exchange

% use function to load .cif as atrut array
cifT = importcif([],path);

cifData = struct();
% Reshape it to a 1x1 struct with field names
fieldNames = {cifT.name};
fieldNames = cellfun(@(x) replace(x,{'-','.'},'_'),fieldNames, ...
    'UniformOutput',false);

for i = 1:numel(fieldNames)
    cifData.(fieldNames{i}) = cifT(i).val;
end
end