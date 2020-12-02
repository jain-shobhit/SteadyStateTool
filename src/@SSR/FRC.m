function [om, N] = FRC(O,T_range,ncontsteps,varargin)
% function that computes the FRC of S for a given range

if nargin == 3
    type = 'original';
else
    type = varargin{1};
end

switch type
    case 'original'
        [om,cont_sol,~] = O.sequential_continuation(T_range,'ncontsteps',ncontsteps);
    case 'reformulated'
        [om,cont_sol,~] = O.sequential_continuation_reformulated(T_range,'ncontsteps',ncontsteps);
end

N = nan(length(om),1);
for j = 1:length(om)
    if ~isempty(cont_sol{j})
        N(j) = O.dt*norm(cont_sol{j},'fro');
    end
end    
end