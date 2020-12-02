function [fixDOFs] = BoundaryCoundition(model,S)

switch S
    case 'C'    % Cantilever
        fixDOFs = 1:3;
    case 'B'    % Simply supported
        fixDOFs = [1 2 model.nDOF-2 model.nDOF-1];
    case 'DC'  % Doubly clamped
        fixDOFs = [1:3 (model.nDOF-2):model.nDOF];
end
        