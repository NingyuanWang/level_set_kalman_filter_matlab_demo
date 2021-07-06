function [x_output,var_output,RMSE] = run_test_case(model,instance_parameter,initial_state,xspan_true,measurements,filter_method)
%Runs the test case using the given parameters and evaluates performance
%indicarors.
%   model: custom class, generate using get_test_model
%   instance_parameter: custom class, get_instance_parameter
%   initial_condition: custom class, get_initial_condition
%   xspan_true:true trajectory generated using gen_traj_and_meas
%   measurements: measurements generated using gen_traj_and_meas
%   filter_method: discrete_CKF, level_set_filter, etc. as a function
%   handle
%   Outputs: 
%   x_output: full predicted (TODO: before/after correction) state at
%   measurement times.
%   var_output: variance vector at measurement times
%   RMSE: mean root squared error averaged over all measurement times.
[x_output,var_output] = filter_method(initial_state,measurements,model,instance_parameter);
N = instance_parameter.observation_count;
x_diff = x_output(:,2:N+1) - xspan_true(:,2:N+1);
RMSE = sqrt(sum(x_diff.^2,2)) / sqrt(N);
end

