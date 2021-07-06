function [dMdt] = diffvel(M,K)
%Find the diffusion velocity given a square root of 
% the covariance matrix M, and diffusion coefficient matrix K
%Check: Should K be square root?
n = size(M,1);
% Sigma = M*M';
for j = 1:n
    %As derived originally, less stable.
    % v = M(:,j);
    % dvdt = K * (Sigma\v);
    %New method, more stable when M is near singular
    %b = zeros(n,1);
    %b(j) = 1;
    %dvdt_new = K * (M' \ b);
    %dMdt(:,j) = dvdt_new; 
    
end
dMdt = 0.5 * K / M';
end

