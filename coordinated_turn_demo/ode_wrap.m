function [varargout] = ode_wrap(odesolver_fun,odefun,tspan,y0,dim_state,use_adaptive_solver)
%ode_wrap a workaround to use standard matlab ode solver for matrix valued
%state variable and derivatives, which only takes vector-valued states.
% [varargout] = [y_end,t_int,y_int]
%y0 and y_end is of mean_covariance_sqrt_cls class
%The function gives optional output of intermediate timesteps and
%intermediate state at these steps t_int, y_int. 
%NOTE: performance is not considered, and there is probably a lot of
%overhead cost.
    y0_l = y0.data(:);
    if use_adaptive_solver
        [t_int,y_ll] = odesolver_fun(@(t,y)odefun_l(t,y,odefun,dim_state),tspan,y0_l);
    else
        y_ll = odesolver_fun(@(t,y)odefun_l(t,y,odefun,dim_state),tspan,y0_l);
        t_int = tspan;
    end
    y_l = y_ll(end,:);
    y_end = reshape(y_l,dim_state,dim_state+1);
    nout = max(nargout,1);
    varargout{1} = mean_covariance_sqrt_cls(y_end);
    if nout >= 2
        varargout{2} = t_int;
    end
    if nout >= 3
        nn = numel(t_int);
        varargout{3} = reshape(y_ll',[dim_state,dim_state+1,nn]);
    end

end
function [dydt_l] = odefun_l (t,y_l,odefun,dim_state)
    y = reshape(y_l,dim_state,dim_state+1);
    dydt = odefun(t,y);
    dydt_l = dydt(:);
end