function eta = convolution_dT(O,phi)
eta = zeros(size(phi));
for j = 1:size(eta,1)
    eta(j,:) = conv(O.weights.* phi(j,:),O.DLvec(j,:),'same');
end
end