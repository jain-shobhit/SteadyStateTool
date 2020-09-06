function [om, N] = FRC(O,T_range,ncontsteps)
% function that computes the FRC of S for a given range

[om,cont_sol,~] = O.sequential_continuation(T_range,'ncontsteps',ncontsteps);

N = nan(length(om),1);
for j = 1:length(om)
    if ~isempty(cont_sol{j})
        N(j) = O.dt*norm(cont_sol{j},'fro');
    end
end    
end