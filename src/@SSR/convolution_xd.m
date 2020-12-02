function etad = convolution_xd(O,phi)
% convolution with J to obtain velocities
etad = zeros(size(phi));
for j = 1:O.n
    etad(j,:) = conv(O.weights.* phi(j,:),O.Jvec(j,:),'same');
end
end