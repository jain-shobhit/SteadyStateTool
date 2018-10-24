function eta = convolution_x(O,phi)
% convolution with L to obtain velocities
eta = zeros(size(phi));
for j = 1:O.n
    eta(j,:) = conv(O.weights.* phi(j,:),O.Lvec(j,:),'same');
end
end