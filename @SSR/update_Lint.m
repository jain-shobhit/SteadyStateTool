function update_Lint(O)
% this is a function that does the "pre integration" of convolution
% integrals described in the Picard Iteration section
% currently at one point / integral
Lint = zeros(O.n,length(O.t),length(O.t));
for i = 1:length(O.t)
        Lint(:,i,:) = O.dt * O.Lvec(:,ceil(end/2)-i+1:ceil(end/2)-i+length(O.t));
end
Lint(:,1,:) = Lint(:,1,:)/2;
Lint(:,end,:) = Lint(:,end,:)/2;
O.Lint = Lint;
O.isupdated.Lint = true;
end