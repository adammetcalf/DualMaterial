% Function to find the closest grid point
function [idx, idy, idz] = findClosestGridPoint(x, y, z, p)
    [~, idx] = min(abs(x(1,:,1) - p(1)));
    [~, idy] = min(abs(y(:,1,1) - p(2)));
    [~, idz] = min(abs(z(1,1,:) - p(3)));
end