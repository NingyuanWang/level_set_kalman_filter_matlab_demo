function [dxMdt] = advdiffvel_mc_avg(t,xM,K,advvel,alpha)
%Finds advection diffusion velocity given the center x, a square root
% of the covariance matrix M, dynamics advvel, diffusion
% K, and time t (not used). xM is given by a n by n+1 matrix
% This version uses the average velocity on both sides.
n = size(K,1);%xM is of size n*(n+1)
dMdt_adv_p = zeros(n);
dMdt_adv_n = zeros(n);
for j = 1:n
    x_p = alpha * xM(:,j+1) + xM(:,1);
    x_n = -alpha * xM(:,j+1) + xM(:,1);
    dxdt_adv_p = advvel(t,x_p);
    dMdt_adv_p(:,j) = dxdt_adv_p;
    dxdt_adv_n = advvel(t,x_n);
    dMdt_adv_n(:,j) = dxdt_adv_n;
end
dMdt_diff = diffvel(xM(:,2:(n+1)),K);
dxdt_center = 0.5 * mean(dMdt_adv_n + dMdt_adv_p,2);
dMdt = (dMdt_adv_p - dMdt_adv_n) ./ (2.0 * alpha) + dMdt_diff;
dxMdt = cat(2,dxdt_center,dMdt);
end
