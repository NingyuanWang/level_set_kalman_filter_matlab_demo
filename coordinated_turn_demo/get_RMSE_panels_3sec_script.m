%This script compares the performance between CD-CKF and level set filter
%for the coordinated turn problem. 

%Setting constants:
sampling_interval = 3;
m_list = [64,32,16,8,4,2,1];
w0_degree_list = [6,12,24];
divergence_bound = 500;
%Set grid: 
figure(1)
clf
t = tiledlayout(numel(w0_degree_list),4,'TileSpacing','compact');

sample_count = 100;
%Iterate through w0_degree_list:
for j = 1:numel(w0_degree_list)
    w0_degree = w0_degree_list(j);
    %Setup problem:
    xS0 = get_initial_condition(w0_degree);
    x0 = xS0.mean;
    test_model = get_test_model();
    traj_generating_instance_parameter = get_instance_parameter(2,w0_degree,sampling_interval);
    %generate test case:
    fprintf('Generating simulated measurements...\n')
    [xtrue_list,ymeasure_list] = get_traj_and_measure_high_level(x0,test_model,traj_generating_instance_parameter,sample_count);
    fprintf('Starting tracking filters..\n')
    [RMSE_vm_CDCKF,divergence_count_CDCKF] = find_RMSE_fun(test_model,w0_degree,sampling_interval,@continuous_discrete_cubature_kalman_filter,xS0,m_list,xtrue_list,ymeasure_list,false,divergence_bound);
    [RMSE_vm_LSKF_fixed,divergence_count_LSKF_fixed] = find_RMSE_fun(test_model,w0_degree,sampling_interval,@level_set_filter,xS0,m_list,xtrue_list,ymeasure_list,false,divergence_bound);
    [RMSE_vm_LSKF_adaptive,divergence_count_LSKF_adaptive] = find_RMSE_fun(test_model,w0_degree,sampling_interval,@level_set_filter,xS0,m_list,xtrue_list,ymeasure_list,true,divergence_bound);
    nexttile(4*(j-1)+1)
    
    plot(log2(m_list),sqrt(RMSE_vm_CDCKF(1,:).^2 + RMSE_vm_CDCKF(3,:).^2 + RMSE_vm_CDCKF(5,:).^2),'-o')
    hold on
    plot(log2(m_list),sqrt(RMSE_vm_LSKF_fixed(1,:).^2 + RMSE_vm_LSKF_fixed(3,:).^2 + RMSE_vm_LSKF_fixed(5,:).^2),'--x')
    plot(log2(m_list),sqrt(RMSE_vm_LSKF_adaptive(1,:).^2 + RMSE_vm_LSKF_adaptive(3,:).^2 + RMSE_vm_LSKF_adaptive(5,:).^2),'-.+')
    axis([0 7 20 45])
    %legend('CDCKF','LSKF-fixed RK4','LSKF-adaptive')
    title('RMSE in position')
    xlabel('log2(m)')
    ylabel('RMSE (meter)')
    nexttile(4*(j-1)+2)
    
    plot(log2(m_list),sqrt(RMSE_vm_CDCKF(2,:).^2 + RMSE_vm_CDCKF(4,:).^2 + RMSE_vm_CDCKF(6,:).^2),'-o')
    hold on
    plot(log2(m_list),sqrt(RMSE_vm_LSKF_fixed(2,:).^2 + RMSE_vm_LSKF_fixed(4,:).^2 + RMSE_vm_LSKF_fixed(6,:).^2),'--x')
    plot(log2(m_list),sqrt(RMSE_vm_LSKF_adaptive(2,:).^2 + RMSE_vm_LSKF_adaptive(4,:).^2 + RMSE_vm_LSKF_adaptive(6,:).^2),'-.+')
    axis([0 7 6.5 9.5])
    %legend('CDCKF','LSKF-fixed RK4','LSKF-adaptive')
    title('RMSE in velocity')
    xlabel('log2(m)')
    ylabel('RMSE (meter/second)')
    nexttile(4*(j-1)+3)
    
    plot(log2(m_list),RMSE_vm_CDCKF(7,:),'-o')
    hold on
    plot(log2(m_list),RMSE_vm_LSKF_fixed(7,:),'--x')
    plot(log2(m_list),RMSE_vm_LSKF_adaptive(7,:),'-.+')
    axis([0 7 0.014 0.018])
    %legend('CDCKF','LSKF-fixed RK4','LSKF-adaptive')
    title('RMSE in turn rate')
    xlabel('log2(m)')
    ylabel('RMSE (rad/second)')
    %Divergence count
    nexttile(4*(j-1)+4)
    plot(log2(m_list),divergence_count_CDCKF,'-o')
    hold on
    plot(log2(m_list),divergence_count_LSKF_fixed,'--x')
    plot(log2(m_list),divergence_count_LSKF_adaptive,'-.+')
    axis([0 7 -5 100])
    if j == 1
        legend('CDCKF','LSKF-fixed RK4','LSKF-adaptive')
    end
    title('Count of divergent results')
    xlabel('log2(m)')
    ylabel('Count')
    drawnow
end
%Reformatting titles, legends, etc.
%ylabel(t,'Root Mean Square Error')
titleletters = 'a':'z';
for j = 1:4*numel(w0_degree_list)
    nexttile(j)
    title(sprintf('(%s)',titleletters(j)))
end
column_title_list = {'RMSE in Position','RMSE in velocity','RMSE in turn rate','Divergent results count'};
for j = 1:4
    nexttile(j)
    title({column_title_list{j},sprintf('(%s)',titleletters(j))})
end
for j = 1:numel(w0_degree_list)
    nexttile(4*(j-1)+1)
    ylabel({sprintf('\\omega_0 = %i',w0_degree_list(j));'meter'})
    nexttile(4*(j-1)+2)
    ylabel('meter/second')
    nexttile(4*(j-1)+3)
    ylabel('rad/second')
    nexttile(4*(j-1)+4)
    ylabel('Count')
end
exportgraphics(t,'RMSE_3_sec.png','Resolution',450)
savefig('RMSE_3_sec.fig')

function [xtrue_list,ymeasure_list] = get_traj_and_measure_high_level(x0,test_model,traj_generating_instance_parameter,sample_count)
    xtrue_list = zeros(test_model.dim_state,traj_generating_instance_parameter.observation_count+1,sample_count);
    ymeasure_list = zeros(test_model.dim_ob,traj_generating_instance_parameter.observation_count+1,sample_count);
    sc = parallel.pool.Constant(RandStream('Threefry','Seed',180));%Fixes random seed.
    parfor n = 1:sample_count
        %Fixing random seed:
        stream = sc.Value;   % Extract the stream from the Constant
        stream.Substream = n;
        [xspan,yspan] = gen_traj_and_meas(x0,test_model,traj_generating_instance_parameter);    
        xtrue_list(:,:,n) = xspan;
        ymeasure_list(:,:,n) = yspan;
    end
end
function [RMSE_vm,divergence_count_vm] = find_RMSE_fun(test_model,w0,sampling_interval,filter_method,xS0,m_list,xtrue_list,ymeasure_list,use_adaptive_solver,divergence_bound)
%Finds the RMSE for filter_method given measurments, true location and
%initial condition.
%   Detailed explanation goes here
RMSE_vm = zeros(test_model.dim_state,numel(m_list));
divergence_count_vm = zeros(1,numel(m_list));
sample_count = size(xtrue_list,3);

for j = 1:numel(m_list)
    m = m_list(j);
    fprintf('m = %g\n',m)
    instance_parameter = get_instance_parameter(m,w0,sampling_interval,1.0,use_adaptive_solver);
    RMSE_list = zeros(test_model.dim_state,sample_count);
    divergence_list = zeros(1,sample_count);
    parfor n = 1:sample_count
        [~,~,RMSE] = run_test_case(test_model,instance_parameter,xS0,xtrue_list(:,:,n),ymeasure_list(:,:,n),filter_method);
        RMSE_list(:,n) = RMSE;
        if sqrt(sum(RMSE.^2,'all')) > divergence_bound
            divergence_list(n) = 1;
        end
    end
    divergence_count = sum(divergence_list);
    divergence_count_vm(j) = divergence_count;
    RMSE_vm(:,j) = sqrt(sum(RMSE_list.^2.*(1-divergence_list),2)) / sqrt(sample_count-divergence_count + 1e-6);
end
end

