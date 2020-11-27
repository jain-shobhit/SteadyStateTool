function X = AssembleQuadCoeff(model,Misc)
ndof=model(1).ndof;
nelem=model(1).nelem;
% K = tensor(@zeros,[ndof ndof ndof]);
X = cell(1,length(Misc.freedofs));
for j = 1:length(Misc.freedofs)
    X{j} = zeros(ndof);
end


for i = 1:nelem
    disp([num2str(i/nelem*100) '% \t completed'])
    [K2el, K3el]=Allman_nonlinear_stiffness(model(i).BL,model(i).Kxx,model(i).Kyy,model(i).Kxy,model(i).Kt,model(i).Am,model(i).area);
    [K2el, ~] = RotateK2K3(K2el,K3el, model(i).T);
    index = model(i).index;
    for k = 1:length(index)
        j = index(k);
        XX = X{j};
        XX(index,index) = XX(index,index) + K2el(k,:,:);
        X{j} = XX;
    end
end
end