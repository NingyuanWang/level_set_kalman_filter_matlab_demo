function [xS_initial] = get_initial_condition(w0_degree,varargin)
%Generates initial condition used for the coordinated turns model
%   w0: rotation speed in degrees per second
w0_sigma_degree=5.72958;
if nargin >= 2
    w0_sigma_degree = varargin{1};
end
x0_initial = [1000;0;2650;150;200;0;w0_degree*pi/180];
Sigma0 = diag([1e1;1;1e1;1;1e1;1;w0_sigma_degree*pi/180].^2);
U0 = chol(Sigma0);
xS_initial = mean_covariance_sqrt_cls(x0_initial,U0);
end

