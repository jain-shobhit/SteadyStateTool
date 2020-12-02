function class = copyfromfull(class,full)
% copies properties of the second entry to the first entry
if ~isempty(full.f)
    class.f = full.f;
end
if ~isempty(full.T)
    class.T = full.T;
    class.n_steps = full.n_steps;
end
if ~isempty(full.order)
    class.order = full.order;
end
    class.S = full.S;
    class.DS = full.DS;
end