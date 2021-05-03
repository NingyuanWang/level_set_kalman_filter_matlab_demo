function [xspan,yspan] = gen_traj_and_meas(true_initial_state,continuous_discrete_model,instance_parameter)
%Generate a trajectory and measurements given initial state
    N = instance_parameter.observation_count;
    tspan = linspace(0,instance_parameter.T_simulation,N+1);
    xspan = euler_maruyama(tspan,true_initial_state,continuous_discrete_model.K_process,continuous_discrete_model.dxdt,ceil(4096/instance_parameter.ode_subdivision));
    yspan = zeros(continuous_discrete_model.dim_ob,N+1);
    for n = 1:N+1
        v = mvnrnd(zeros(continuous_discrete_model.dim_ob,1),continuous_discrete_model.K_observe)';
        yspan(:,n) = continuous_discrete_model.ob_fun(xspan(:,n)) + v;
    end
end

