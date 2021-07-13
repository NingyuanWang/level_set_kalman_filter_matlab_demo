%This script compares the performance between CD-CKF and level set filter
%for the coordinated turn problem. 

%Setting constants:
column_count = 3;
sampling_interval_list = [1,2,3,4,5,6,7];
m_CDCKF = 64;
m_LSKF = 1;
w0_degree_list = [6,12,24];
divergence_bound = 500;
%Set grid: 
figure(1)
clf
t = tiledlayout(numel(w0_degree_list),column_count,'TileSpacing','compact');

sample_count = 100;
%Iterate through w0_degree_list:
for j = 1:numel(w0_degree_list)
    w0_degree = w0_degree_list(j);
    xS0 = get_initial_condition(w0_degree,0.2);
    x0 = xS0.mean;
    test_model = get_test_model();
    RMSE_CDCKF_vT = zeros(test_model.dim_state,numel(sampling_interval_list));
    covar_CDCKF_vT = zeros(test_model.dim_state,numel(sampling_interval_list));
    divcount_CDCKF_vT = zeros(1,numel(sampling_interval_list));
    RMSE_LSKF_vT = zeros(test_model.dim_state,numel(sampling_interval_list));
    covar_LSKF_vT = zeros(test_model.dim_state,numel(sampling_interval_list));
    divcount_LSKF_vT = zeros(1,numel(sampling_interval_list));
    for k = 1:numel(sampling_interval_list)      
        %Setup problem:
        sampling_interval = sampling_interval_list(k);
        traj_generating_instance_parameter = get_instance_parameter(2,w0_degree,sampling_interval);
        %generate test case:
        fprintf('Generating simulated measurements...\n')
        [xtrue_list,ymeasure_list] = get_traj_and_measure_high_level(xS0,test_model,traj_generating_instance_parameter,sample_count);
        fprintf('Starting tracking filters..\n')
        xS0_guess = mean_covariance_sqrt_cls(mvnrnd(xS0.mean,xS0.covariance,1),xS0.c_sqrt);
        [RMSE_CDCKF,covar_CDCKF,divcount_CDCKF] = find_RMSE_covar_fun(test_model,w0_degree,sampling_interval,@continuous_discrete_cubature_kalman_filter,xS0_guess,m_CDCKF,xtrue_list,ymeasure_list,false,divergence_bound);
        RMSE_CDCKF_vT(:,k) = RMSE_CDCKF;
        covar_CDCKF_vT(:,k) = covar_CDCKF;
        divcount_CDCKF_vT(k) = divcount_CDCKF;
        [RMSE_LSKF_adaptive,covar_LSKF_adaptive,divcount_LSKF] = find_RMSE_covar_fun(test_model,w0_degree,sampling_interval,@level_set_filter,xS0_guess,m_LSKF,xtrue_list,ymeasure_list,true,divergence_bound);
        RMSE_LSKF_vT(:,k) = RMSE_LSKF_adaptive;
        covar_LSKF_vT(:,k) = covar_LSKF_adaptive;
        divcount_LSKF_vT(k) = divcount_LSKF;
    end
    nexttile(column_count*(j-1)+1);
    hold on
    plot(sampling_interval_list,norm_in_position(RMSE_CDCKF_vT),'o-')
    plot(sampling_interval_list,norm_in_position(RMSE_LSKF_vT),'-.+')
    %legend({'CDCKF-64','LSKF-adaptive'},'Location','Northwest')
    %plot(sampling_interval_list,norm_in_position(sqrt(covar_CDCKF_vT)),'x--')
    %plot(sampling_interval_list,norm_in_position(sqrt(covar_LSKF_vT)),'d:')
    %legend('RMSE: CDCKF-64','RMSE: LSKF adaptive','stdev: CDCKF-64','stdev: LSKF adaptive')
    xlabel('sampling interval (s)')
    ylabel('meters')
    nexttile(column_count*(j-1)+2);
    hold on
    plot(sampling_interval_list,norm_in_velocity(RMSE_CDCKF_vT),'o-')
    plot(sampling_interval_list,norm_in_velocity(RMSE_LSKF_vT),'-.+')
    %legend({'CDCKF-64','LSKF-adaptive'},'Location','Northwest')
    %plot(sampling_interval_list,norm_in_velocity(sqrt(covar_CDCKF_vT)),'x--')
    %plot(sampling_interval_list,norm_in_velocity(sqrt(covar_LSKF_vT)),'d:')
    %legend('RMSE: CDCKF-64','RMSE: LSKF adaptive','stdev: CDCKF-64','stdev: LSKF adaptive')
    xlabel('sampling interval (s)')
    ylabel('meters/second')
    nexttile(column_count*(j-1)+3);
    hold on
    plot(sampling_interval_list,norm_in_turnrate(RMSE_CDCKF_vT)*180/pi,'o-')
    plot(sampling_interval_list,norm_in_turnrate(RMSE_LSKF_vT)*180/pi,'-.+')
    %legend({'CDCKF-64','LSKF-adaptive'},'Location','Northwest')
    %plot(sampling_interval_list,norm_in_turnrate(sqrt(covar_CDCKF_vT))*180/pi,'x--')
    %plot(sampling_interval_list,norm_in_turnrate(sqrt(covar_LSKF_vT))*180/pi,'d:')
    %legend('RMSE: CDCKF-64','RMSE: LSKF adaptive','stdev: CDCKF-64','stdev: LSKF adaptive')
    xlabel('sampling interval (s)')
    ylabel('degrees/second')
    %Count of divergent results is not drawn since it is constantly 0
    %nexttile(column_count*(j-1)+4);
    %hold on
    %plot(sampling_interval_list,divcount_CDCKF_vT,'o-')
    %plot(sampling_interval_list,divcount_LSKF_vT,'-.+')
    %legend({'CDCKF-64','LSKF-adaptive'},'Location','Northwest')
    %xlabel('sampling interval (s)')
    %ylabel('count')
    drawnow
end
%Reformatting titles, legends, etc.
%ylabel(t,'Root Mean Square Error')
titleletters = 'a':'z';
for j = 1:column_count*numel(w0_degree_list)
    nexttile(j)
    title(sprintf('(%s)',titleletters(j)))
end
column_title_list = {'RMSE in Position','RMSE in velocity','RMSE in turn rate'};
for j = 1:column_count
    nexttile(j)
    title({column_title_list{j},sprintf('(%s)',titleletters(j))})
end
for j = 1:numel(w0_degree_list)
    nexttile(column_count*(j-1)+1)
    ylabel({sprintf('\\omega_0 = %i',w0_degree_list(j));'meters'})
end
nexttile(column_count)
%axis([1 8 0.14 0.2])
legend({'CDCKF-64','LSKF-adaptive'},'Location','Northeast')
exportgraphics(t,'RMSE_var_secs_n.png','Resolution',450)
savefig('RMSE_var_secs_n.fig')
function [val_in_position] = norm_in_position(values)
    val_in_position = sqrt(values(1,:).^2+values(3,:).^2+values(5,:).^2);
end
function [val_in_velocity] = norm_in_velocity(values)
    val_in_velocity = sqrt(values(2,:).^2+values(4,:).^2+values(6,:).^2);
end
function [val_in_turnrate] = norm_in_turnrate(values)
    val_in_turnrate = abs(values(7,:));
end
function [xtrue_list,ymeasure_list] = get_traj_and_measure_high_level(xS0,test_model,traj_generating_instance_parameter,sample_count)
    xtrue_list = zeros(test_model.dim_state,traj_generating_instance_parameter.observation_count+1,sample_count);
    ymeasure_list = zeros(test_model.dim_ob,traj_generating_instance_parameter.observation_count+1,sample_count);
    sc = parallel.pool.Constant(RandStream('Threefry','Seed',0));%Fixes random seed.
    parfor n = 1:sample_count
        %Fixing random seed:
        stream = sc.Value;   % Extract the stream from the Constant
        stream.Substream = n;
        x0 = xS0.mean;
        %x0 = mvnrnd(xS0.mean,xS0.covariance,1);
        [xspan,yspan] = gen_traj_and_meas(x0,test_model,traj_generating_instance_parameter);    
        xtrue_list(:,:,n) = xspan;
        ymeasure_list(:,:,n) = yspan;
    end
end
function [RMSE,covar,divergence_count] = find_RMSE_covar_fun(test_model,w0,sampling_interval,filter_method,xS0_guess,m,xtrue_list,ymeasure_list,use_adaptive_solver,divergence_bound)
%Finds the RMSE for filter_method given measurments, true location and
%initial condition.
%   Detailed explanation goes here
sample_count = size(xtrue_list,3);
instance_parameter = get_instance_parameter(m,w0,sampling_interval,1.0,use_adaptive_solver);
RMSE_list = zeros(test_model.dim_state,sample_count);
var_list = zeros(test_model.dim_state,sample_count);
divergence_list = zeros(1,sample_count);
parfor n = 1:sample_count
    [~,var_at_t,RMSE_single] = run_test_case(test_model,instance_parameter,xS0_guess,xtrue_list(:,:,n),ymeasure_list(:,:,n),filter_method);
    var_single = mean(var_at_t,2);
    var_list(:,n) = var_single;
    if sqrt(sum(RMSE_single.^2,'all')) > divergence_bound
        divergence_list(n) = 1;
    end
    RMSE_list(:,n) = RMSE_single;
end
divergence_count = sum(divergence_list);
covar = mean(var_list,2);
RMSE = sqrt(sum(RMSE_list.^2.*(1-divergence_list),2)) / sqrt(sample_count-divergence_count + 1e-6);
end

