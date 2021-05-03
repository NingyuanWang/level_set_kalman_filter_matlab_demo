function [varargout] = level_set_filter(xS_guess,measurements,continuous_discrete_model,instance_parameter)
%Finds the predicted state variables and its predicted variance at tspan,
%given initial state guess and measurements at the times by level set method. 
%WARNING: ode_wrap is MORE object-oriented than previous implementations
%and cannot be substituted.
%	[varargout] = [X_corrected, var_corrected, t_int, X_int, var_int]
%	only X_corrected is given by default.
%	the mean, variance at measurements. States at intermediate steps t_int. 
%	initial_state_guess: mean_covariance_cls type
%	measurements: matrix containing all the measurements made, including an
%	initial measurement at start time that is not used.
%	continuous_discrete_model: of continuous_discrete_model_cls type
%	instance_parameter: of instance_parameter_cls type
    xS = xS_guess;
    X_corrected = zeros(continuous_discrete_model.dim_state,instance_parameter.observation_count+1);
    var_corrected = zeros(continuous_discrete_model.dim_state,instance_parameter.observation_count+1);
    t_int_full = [0];
    X_int_full = xS_guess.mean;
    var_int_full = diag(xS_guess.covariance);
    X_corrected(:,1) = xS_guess.mean;
    var_corrected(:,1) = diag(xS_guess.covariance);
    X_predicted = X_corrected;
    var_predicted = var_corrected;
    t = 0;
    dt = instance_parameter.T_simulation / instance_parameter.observation_count;
    for n = 1:instance_parameter.observation_count
        %Prediction part
        if nargout <=2
            if instance_parameter.use_adaptive_solver
                [xS] = ode_wrap(@ode113,@(t,x)advdiffvel_mc_avg(t,x,continuous_discrete_model.K_process,@continuous_discrete_model.dxdt,instance_parameter.alpha),linspace(t,t+dt,instance_parameter.ode_subdivision+1),xS,continuous_discrete_model.dim_state,instance_parameter.use_adaptive_solver);
            else
                [xS] = ode_wrap(@ode4,@(t,x)advdiffvel_mc_avg(t,x,continuous_discrete_model.K_process,@continuous_discrete_model.dxdt,instance_parameter.alpha),linspace(t,t+dt,instance_parameter.ode_subdivision+1),xS,continuous_discrete_model.dim_state,instance_parameter.use_adaptive_solver);
            end
        else
            if instance_parameter.use_adaptive_solver
                [xS,t_int,X_int] = ode_wrap(@ode113,@(t,x)advdiffvel_mc_avg(t,x,continuous_discrete_model.K_process,@continuous_discrete_model.dxdt,instance_parameter.alpha),linspace(t,t+dt,instance_parameter.ode_subdivision+1),xS,continuous_discrete_model.dim_state,instance_parameter.use_adaptive_solver);
            else
                [xS,t_int,X_int] = ode_wrap(@ode4,@(t,x)advdiffvel_mc_avg(t,x,continuous_discrete_model.K_process,@continuous_discrete_model.dxdt,instance_parameter.alpha),linspace(t,t+dt,instance_parameter.ode_subdivision+1),xS,continuous_discrete_model.dim_state,instance_parameter.use_adaptive_solver);
            end
            t_int_full = cat(1,t_int_full,t_int(2:end));
            X_int_full = cat(2,X_int_full,squeeze(X_int(:,1,2:end)));
            if nargout == 5
                var_int = zeros(continuous_discrete_model.dim_state,numel(t_int)-1);
                for j = 1:numel(t_int) - 1
                    S_sqrt = X_int(:,2:end,j+1);
                    var_int(:,j) = diag(S_sqrt*S_sqrt');
                end
                var_int_full = cat(2,var_int_full,var_int);
            end
        end
        X_predicted(:,n+1) = xS.mean;
        if nargout >=2
            var_predicted(:,n+1) = diag(xS.covariance);
        end
        %Update part (Rewritten to use square root update form.
        xS = square_root_cubature_measurement_update(xS,measurements(:,n+1),@continuous_discrete_model.ob_fun,continuous_discrete_model.K_observe);
        %Data_output
        X_corrected(:,n+1) = xS.mean;
        if nargout >=2
            var_corrected(:,n+1) = diag(xS.covariance);
        end
        if nargout >= 3
            if instance_parameter.return_uncorrected_prediction
                X_int_full(:,end) = X_corrected(:,n+1);
                var_int_full(:,end) = var_corrected(:,n+1);
            else
                X_int_full(:,end) = X_predicted(:,n+1);
                var_int_full(:,end) = var_predicted(:,n+1);
            end
        end
        t = t + dt;
    end
    %End cycle, parse output:
    if instance_parameter.return_uncorrected_prediction
        X_output = X_predicted;
        var_output = var_predicted;
    else
        X_output = X_corrected;
        var_output = var_corrected;
    end
    varargout{1} = X_output; 
    if nargout>= 2
        varargout{2} = var_output;
    end
    if nargout >= 3
        varargout{3} = t_int_full;
    end
    if nargout >= 4
        varargout{4} = X_int_full;
    end
    if nargout >= 5
        varargout{5} = var_int_full;
    end
end

