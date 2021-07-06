function [xS_initial] = get_initial_condition()
%Generates initial condition used for the coordinated turns model
%   w0: rotation speed in degrees per second
x0_initial = [1;0;0];
Sigma0 = diag([1e-2;1e-2;3e-2].^2);
U0 = chol(Sigma0);
xS_initial = mean_covariance_sqrt_cls(x0_initial,U0);
end

