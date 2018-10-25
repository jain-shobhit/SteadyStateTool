function w = NewtonCotes(order)
% returns the weights corresponding to the appropriate order for performing
% Newton-Cotes numerical integration on a uniformly spaced grid
switch order
    case 1
        w = [1 1]/2;
    case 2
        w = [1 4 1]/3;
    case 3
        w = [1 3 3 1]*(3/8);
    case 4
        w = [7 32 12 32 7]*(2/45);
end
end