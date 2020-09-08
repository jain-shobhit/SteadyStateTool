function modes = modeselect(full,param,X) 
% automated nonlinear mode selection procedure

try
    modes = param.initialmodes; 
catch    
    % initial set based on spectral quotient analysis
    realpart = min(abs(real(-full.zeta.*full.omega0 + sqrt(full.zeta.^2-1).*full.omega0)),...
        abs(real(-full.zeta.*full.omega0 - sqrt(full.zeta.^2-1).*full.omega0)));
    [rpsorted,ind] = sort(realpart);
    ratio = zeros(1,length(rpsorted)-1);
    for i = 1:length(rpsorted)-1
        ratio(i) = rpsorted(i+1)/rpsorted(i);
    end
    cut = param.nmodes;
    while cut > ceil(param.nmodes/3)    % limit the number of modes selected via s.q.a.
        [~,cut] = max(ratio);
        if cut ~= 1
        ratio = ratio(1:cut-1);
        end
    end
    modes = ind(1:cut).';

    % refinement based on linear response to forcing (if given)
    if ~isempty(full.f)
        [x_linf,~] = full.LinearResponse();
        z = full.U.' * full.M * x_linf;
        b = size(z);
        NORM = zeros(1,b(1));
        for j = 1:b(1)
            NORM(j) = norm(z(j,:));
        end
        while length(modes) < ceil(param.nmodes*2/3)    % limit the number of modes selected via linear response   
            if sum(NORM(modes)) > 0.9*sum(NORM)
                break
            end
            modesC = setdiff(1:full.n,modes);
            [~,ind2] = max(NORM(modesC));
            modes = [modes ind2];
        end
        percent = 100*sum(NORM(modes))/sum(NORM);
        disp(['The initial mode selection recovers ' num2str(percent) '% of the linear response.'])
    else
        disp('Since no forcing was provided for the SSR class, it is not guaranteed that the mode selection recovers the linear response.')
    end
end
disp('Initial linear selection:')
disp(modes)

% nonlinear selection
itnumber = 0;
while length(modes) < param.nmodes
    red = SSR(full.M,full.C,full.K,modes);
    red = copyfromfull(red,full);       % copies general props
    switch param.example                % obtaining the nonlinear coefficients is problem dependent
        case 'simple'                   % heavy on memory
            R = S11reduct(red,X);
        case 'plate'                    % memory preserving approach
            V = zeros(param.fulldof,red.n);
            VC = zeros(param.fulldof,length(red.M)-red.n);
            V(param.freedofs,:) = red.U;
            VC(param.freedofs,:) = red.U2;
            [K2,~] = AssembleTensors(param.model,V,VC); % computes quadratic coeffs R_i
            R = cell(1,size(red.U2,2));
            for i = 1:size(red.U2,2)
                R{i} = K2(i,:,:);
            end
    end
    W = red.SSM2(R);                    % computing SSM
    A = [zeros(red.n),eye(red.n);-red.K1,-red.C1];
    [~ ,scaldir] = scalg(A,W);          % computing directional curvatures
    [scnorm,scindex] = sort(abs(scaldir));
    Jmodes = setdiff(1:full.n,modes);
    scindex = fliplr(scindex);
    crit = sum(scnorm) * (1-param.curvtolerance);
    Z = 0;
    numrec = 0;
    while Z < crit 
        Z = Z + max(scnorm);
        scnorm = setdiff(scnorm,max(scnorm));
        numrec = numrec + 1;
    end
    recomm = Jmodes(scindex(1:numrec)); % recommnended modes
    check = param.nmodes - (length(modes) + length(recomm)); % check if it overshoots the number of allowed modes
    if check < 0
        modes = sort([modes recomm(1:(param.nmodes - length(modes)))]);
    else
        modes = sort([modes recomm]);
    end
    disp(modes)
    itnumber = itnumber + 1;
    if ~param.type
        break
    end
end
disp(['The nonlinear mode selection procedure was performed ' num2str(itnumber) ' time(s).'])
end