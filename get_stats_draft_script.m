%Setting constants:
sampling_interval = 2;
m_list = [64,32,16,8,4,2];
w0_degree = 3;
N_repeat = 64;
filter_method = @continuous_dicrete_cubature_kalman_filter;
%Setup problem:
xS0 = get_initial_condition(w0_degree);
x0 = xS0.mean;
test_model = get_test_model();
traj_generating_instance_parameter = get_instance_parameter(2,w0_degree,sampling_interval);
if exist('xtrue_list','var') && numel(xtrue_list) == test_model.dim_state*(traj_generating_instance_parameter.observation_count+1)*N_repeat
    fprintf('Using existing simulated measurements...\n');
else
    xtrue_list = zeros(test_model.dim_state,traj_generating_instance_parameter.observation_count+1,N_repeat);
    ymeasure_list = zeros(test_model.dim_ob,traj_generating_instance_parameter.observation_count+1,N_repeat);
    fprintf('Generating simulated measurements...\n')
    sc = parallel.pool.Constant(RandStream('Threefry'));%Fixes random seed.
    parfor n = 1:N_repeat
        %Fixing random seed:
        stream = sc.Value;   % Extract the stream from the Constant
        stream.Substream = n;
        [xspan,yspan] = gen_traj_and_meas(x0,test_model,traj_generating_instance_parameter);    
        xtrue_list(:,:,n) = xspan;
        ymeasure_list(:,:,n) = yspan;
    end
end    
%Containers for results:
RMSE_vm = zeros(test_model.dim_state,numel(m_list));
fprintf('Starting tracking filter..\n')
for j = 1:numel(m_list)
    m = m_list(j);
    fprintf('m = %g\n',m)
    instance_parameter = get_instance_parameter(m,w0_degree,sampling_interval);
    RMSE_list = zeros(test_model.dim_state,N_repeat);
    parfor n = 1:N_repeat
        [~,~,RMSE] = run_test_case(test_model,instance_parameter,xS0,xtrue_list(:,:,n),ymeasure_list(:,:,n),filter_method);
        RMSE_list(:,n) = RMSE;
    end
    RMSE_vm(:,j) = sqrt(sum(RMSE_list.^2,2)) / sqrt(N_repeat);
end