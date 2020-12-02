function [OMEGA, SOL, PICARD] = sequential_continuation(O,range,varargin)
% this method performs sequential continuation of periodic
% response assuming the time-period (T) of the solution to be
% the continuation parameter. Warning: this scheme cannot
% capture fold points in the frequency response curves.

%% Input parsing
defaultInit = [];
% these inputs are only needed in case of quasiperiodic
% solutions
defaultOmega0 = O.Omega;
defaultparamidx = 1;
ncont = 50; % number of continuation steps between each FREQ

p = inputParser;
addOptional(p,'Omega0',defaultOmega0)
addOptional(p,'init',defaultInit)
addOptional(p,'idx',defaultparamidx)
addOptional(p,'ncontsteps',ncont)
parse(p,varargin{:});
x0 = p.Results.init;
paramidx = p.Results.idx;
Omega0 = p.Results.Omega0;
n_cont_steps = p.Results.ncontsteps;

%% check inputs
if ~isrow(range)
    range = range.';
end
switch O.type
    case 'p'
        Omega_range = 2*pi./range;
    case 'qp'
        Omega_range = range;
end
% check which eigenfrequencies are in this range
check = ((Omega_range(1) - O.omega).*(Omega_range(2) - O.omega)) <0;
idx = find(check); % indices of natural frequencies in the range
if length(idx)>1
    mid_freq = 0.5* (O.omega(idx(1:end-1)) +...
        O.omega(idx(2:end)) );
else
    mid_freq = [];
end
FREQ = sort([Omega_range, O.omega(idx).', mid_freq.']);
n_intervals = length(FREQ) - 1;
SOL = cell(n_intervals*n_cont_steps + 1,1);
OMEGA = zeros(n_intervals*n_cont_steps + 1,1);
PICARD = zeros(n_intervals*n_cont_steps + 1,1);
for j = 1:n_intervals
    interval = FREQ(j:j+1);
    d_Omega = (interval(2)-interval(1))/n_cont_steps;
    if mod(j,2)
        % go left to right in odd-numbered intervals
        Omega_array = interval(1):d_Omega:interval(2);
    else
        % go right to left in even-numbered intervals
        Omega_array = interval(2):-d_Omega:interval(1);
        Omega_array = Omega_array(1:end-1);
    end
    n_Omega = length(Omega_array);
    sol = cell(1,n_Omega);
    pic = zeros(1,n_Omega);
    try_Picard = true;
    for k = 1:n_Omega
        switch O.type
            case 'p'
                O.T = 2*pi./Omega_array(k);
            case 'qp'
                Omega0(paramidx) = Omega_array(k);
                O.Omega = Omega0;
        end
        if isempty(x0)
            x0 = LinearResponse(O);
        end
        if try_Picard
            try
                [x0, xd0] = Picard(O,'init',x0);
            catch
                disp('Switching to Newton--Raphson iteration')
                try_Picard = false;
                [x0, xd0] = NewtonRaphson(O,'init',x0);
            end
        else
            [x0, xd0] = NewtonRaphson(O,'init',x0);
        end
        sol{k} = [x0;xd0];
        pic(k) = try_Picard;
    end
    if mod(j,2)
        % go left to right in odd-numbered intervals
        SOL((j-1)*n_cont_steps + (1:n_Omega)) = sol;
        OMEGA((j-1)*n_cont_steps + (1:n_Omega)) = Omega_array;
        PICARD((j-1)*n_cont_steps + (1:n_Omega)) = pic;
    else
        % go right to left in even-numbered intervals
        OMEGA((j-1)*n_cont_steps + 1 + (1:n_Omega)) = Omega_array(end:-1:1);
        SOL((j-1)*n_cont_steps + 1 + (1:n_Omega)) = sol(end:-1:1);
        PICARD((j-1)*n_cont_steps + 1+ (1:n_Omega)) = pic(end:-1:1);
    end
end
end
