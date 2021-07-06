function [xspan] = euler_maruyama(tspan,x0,K,advvel_fn,step_subdivision)
%Euler Maruyama method for update of state variable x0 at time points tspan
%   K: diffusion matrix, N: time subdivision countadvvel_fn: ODE of dynamic
%   returns x1: state variable at tspan
n = numel(x0);%dimension of space
N = numel(tspan)-1;
xspan = zeros(n,N);
xspan(:,1) = x0(:);
dt_span = tspan(2:N+1)-tspan(1:N);
for stepcount = 1:N
    t = tspan(stepcount);
    dt = dt_span(stepcount) / step_subdivision;
    dK = 1.0 * K *dt;
    x_cur = xspan(:,stepcount);
    for j=1:step_subdivision
        R = mvnrnd(zeros(n,1),dK)';
        advvel = advvel_fn(t,x_cur);
        x_cur = x_cur + dt * advvel + R;
    end
    xspan(:,stepcount+1) = x_cur;
end
end

