%This script compares the performance between CD-CKF and level set filter
%for the coordinated turn problem. 

%Setting constants:
sampling_interval = 4;
m_list = [64,32,16,8,4,2,1];
w0_degree = 12;
sample_count = 96;
%Setup problem:
xS0 = get_initial_condition(w0_degree);
x0 = xS0.mean;
test_model = get_test_model();
traj_generating_instance_parameter = get_instance_parameter(2,w0_degree,sampling_interval);
%Use existing test case or generate test case if non exist:
if exist('xtrue_list','var') && numel(xtrue_list) == test_model.dim_state*(traj_generating_instance_parameter.observation_count+1)*sample_count
    fprintf('Using existing simulated measurements...\n');
else
    fprintf('Generating simulated measurements...\n')
    [xtrue_list,ymeasure_list] = get_traj_and_measure_high_level(x0,test_model,traj_generating_instance_parameter,sample_count);
end
fprintf('Starting tracking filters..\n')
RMSE_vm_CDCKF = find_RMSE_fun(test_model,w0_degree,sampling_interval,@continuous_discrete_cubature_kalman_filter,xS0,m_list,xtrue_list,ymeasure_list,false);
RMSE_vm_LSKF_fixed = find_RMSE_fun(test_model,w0_degree,sampling_interval,@level_set_filter,xS0,m_list,xtrue_list,ymeasure_list,false);
RMSE_vm_LSKF_adaptive = find_RMSE_fun(test_model,w0_degree,sampling_interval,@level_set_filter,xS0,m_list,xtrue_list,ymeasure_list,true);
figure(1)
clf
plot(log2(m_list),sqrt(RMSE_vm_CDCKF(1,:).^2 + RMSE_vm_CDCKF(3,:).^2 + RMSE_vm_CDCKF(5,:).^2))
hold on
plot(log2(m_list),sqrt(RMSE_vm_LSKF_fixed(1,:).^2 + RMSE_vm_LSKF_fixed(3,:).^2 + RMSE_vm_LSKF_fixed(5,:).^2))
plot(log2(m_list),sqrt(RMSE_vm_LSKF_adaptive(1,:).^2 + RMSE_vm_LSKF_adaptive(3,:).^2 + RMSE_vm_LSKF_adaptive(5,:).^2))
axis([0 7 0 100])
legend('CDCKF','LSKF-fixed RK4','LSKF-adaptive')
title('RMSE in position')
xlabel('log2(m)')
ylabel('RMSE (meter)')
figure(2)
clf
plot(log2(m_list),sqrt(RMSE_vm_CDCKF(2,:).^2 + RMSE_vm_CDCKF(4,:).^2 + RMSE_vm_CDCKF(6,:).^2))
hold on
plot(log2(m_list),sqrt(RMSE_vm_LSKF_fixed(2,:).^2 + RMSE_vm_LSKF_fixed(4,:).^2 + RMSE_vm_LSKF_fixed(6,:).^2))
plot(log2(m_list),sqrt(RMSE_vm_LSKF_adaptive(2,:).^2 + RMSE_vm_LSKF_adaptive(4,:).^2 + RMSE_vm_LSKF_adaptive(6,:).^2))
axis([0 7 0 100])
legend('CDCKF','LSKF-fixed RK4','LSKF-adaptive')
title('RMSE in velocity')
xlabel('log2(m)')
ylabel('RMSE (meter/second)')
figure(3)
clf
plot(log2(m_list),RMSE_vm_CDCKF(7,:))
hold on
plot(log2(m_list),RMSE_vm_LSKF_fixed(7,:))
plot(log2(m_list),RMSE_vm_LSKF_adaptive(7,:))
axis([0 7 0 0.1])
legend('CDCKF','LSKF-fixed RK4','LSKF-adaptive')
title('RMSE in turn rate')
xlabel('log2(m)')
ylabel('RMSE (degree/second)')
function [xtrue_list,ymeasure_list] = get_traj_and_measure_high_level(x0,test_model,traj_generating_instance_parameter,sample_count)
    xtrue_list = zeros(test_model.dim_state,traj_generating_instance_parameter.observation_count+1,sample_count);
    ymeasure_list = zeros(test_model.dim_ob,traj_generating_instance_parameter.observation_count+1,sample_count);
    sc = parallel.pool.Constant(RandStream('Threefry','Seed',0));%Fixes random seed.
    parfor n = 1:sample_count
        %Fixing random seed:
        stream = sc.Value;   % Extract the stream from the Constant
        stream.Substream = n;
        [xspan,yspan] = gen_traj_and_meas(x0,test_model,traj_generating_instance_parameter);    
        xtrue_list(:,:,n) = xspan;
        ymeasure_list(:,:,n) = yspan;
    end
end
function [RMSE_vm] = find_RMSE_fun(test_model,w0,sampling_interval,filter_method,xS0,m_list,xtrue_list,ymeasure_list,use_adaptive_solver)
%Finds the RMSE for filter_method given measurments, true location and
%initial condition.
%   Detailed explanation goes here
RMSE_vm = zeros(test_model.dim_state,numel(m_list));
sample_count = size(xtrue_list,3);
for j = 1:numel(m_list)
    m = m_list(j);
    fprintf('m = %g\n',m)
    instance_parameter = get_instance_parameter(m,w0,sampling_interval,1.0,use_adaptive_solver);
    RMSE_list = zeros(test_model.dim_state,sample_count);
    parfor n = 1:sample_count
        [~,~,RMSE] = run_test_case(test_model,instance_parameter,xS0,xtrue_list(:,:,n),ymeasure_list(:,:,n),filter_method);
        RMSE_list(:,n) = RMSE;
    end
    RMSE_vm(:,j) = sqrt(sum(RMSE_list.^2,2)) / sqrt(sample_count);
end
end

