function eta = convolution_x(O,phi)
% convolution with L to obtain velocities
eta = zeros(size(phi));

switch O.sys_order
    case 'first'
        for j = 1:O.N    
            eta(j,:) = conv(O.weights.* phi(j,:),O.Gvec(j,:),'same');
        end        
        
    case 'second'
        for j = 1:O.n
            eta(j,:) = conv(O.weights.* phi(j,:),O.Lvec(j,:),'same');
        end

end