function K2 = AssembleTensors(model,V,VC)
% assembles nonlinear coefficients R_i

M = size(V,2);
N = size(VC,2);
nelem = model(1).nelem;
K2 = tensor(@zeros,[N M M]);


for i=1:nelem
    disp([num2str(i/nelem*100) '% \t completed'])
    
    % element level coefficients / replace these with your specific FE
    % model
    [K2el, K3el] = Allman_nonlinear_stiffness(model(i).BL,model(i).Kxx,model(i).Kyy,model(i).Kxy,model(i).Kt,model(i).Am,model(i).area);
    [K2el, ~] = RotateK2K3(K2el,K3el, model(i).T);
    
    index = model(i).index;   
    Vel = tensor(V(index,:),[18,M]);
    VelC = tensor(VC(index,:),[18,N]);

    % coefficients in the enslaved eqs (needed for SSM computation)
    K2 = K2 + ttt(ttt(ttt(VelC,K2el,1,1),Vel,2,1),Vel,2,1); % quadratic         
end  

