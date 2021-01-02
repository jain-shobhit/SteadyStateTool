function [om, N] = FRC(O,T_range,ncontsteps)
% function that computes the FRC of S for a given range

[om,cont_sol,~] = O.sequential_continuation(T_range,'ncontsteps',ncontsteps);

N = nan(length(om),1);
for j = 1:length(om)
    if ~isempty(cont_sol{j})
        O.T = 2*pi/om(j);
        N(j) = sqrt(O.dt)*norm(cont_sol{j}(1:end/2,:),'fro')/O.T;
    end
end    
end