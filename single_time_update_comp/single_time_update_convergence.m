%% Compares convergence of single time-update step
% CDCKF (IT-1.5 with Skk1 computed once)
% CDCKF-ts (IT-1.5 with Skk1 computed every substep)
% LSKF
sigma_1 = 0.01;
sigma_2 = 0.02;
dt = 0.2;
t0=0;
xS0 = get_initial_condition();
test_model = get_test_model(sigma_1,sigma_2);
%Get "accurate" result:
acc_instance_parameter = get_instance_parameter(1,dt,1.0,true);
[xS_end] = ode_wrap(@ode113,@(t,x)advdiffvel_mc_avg(t,x,test_model.K_process,@test_model.dxdt,acc_instance_parameter.alpha),linspace(t0,t0+dt,acc_instance_parameter.ode_subdivision+1),xS0,test_model.dim_state,acc_instance_parameter.use_adaptive_solver);
%Check for each case:
m_list = [32,16,8,4,2,1];
m_count = numel(m_list);
x_E_CDCKF = zeros(m_count,1);
S_E_CDCKF = zeros(m_count,1);
x_E_CDCKF_ts = x_E_CDCKF;
S_E_CDCKF_ts = S_E_CDCKF;
x_E_LSKF = x_E_CDCKF;
S_E_LSKF = S_E_CDCKF;
x_E_LSKF_RK4 = x_E_CDCKF;
S_E_LSKF_RK4 = S_E_CDCKF;
for j = 1:numel(m_list)
    m = m_list(j);
    instance_parameter = get_instance_parameter(m,dt,1.0,false);
    %CDCKF:
    [xS_CDCKF] = simple_harmonic_cubature_predict_wrap(test_model,linspace(t0,t0+dt,instance_parameter.ode_subdivision+1),xS0);
    %CDCKF-ts:
    xS_CDCKF_ts = xS0;
    for l = 1:instance_parameter.ode_subdivision
        ts = t0 + (l-1)*dt/instance_parameter.ode_subdivision;
        te = t0 + l*dt/instance_parameter.ode_subdivision;
        [xS_CDCKF_ts] = simple_harmonic_cubature_predict_wrap(test_model,[ts,te],xS_CDCKF_ts);
    end
    %LSKF-RK2:
    [xS_LSKF] = ode_wrap(@ode2,@(t,x)advdiffvel_mc_avg(t,x,test_model.K_process,@test_model.dxdt,instance_parameter.alpha),linspace(t0,t0+dt,instance_parameter.ode_subdivision+1),xS0,test_model.dim_state,instance_parameter.use_adaptive_solver);
    %LSKF-RK4:
    [xS_LSKF_RK4] = ode_wrap(@ode4,@(t,x)advdiffvel_mc_avg(t,x,test_model.K_process,@test_model.dxdt,instance_parameter.alpha),linspace(t0,t0+dt,instance_parameter.ode_subdivision+1),xS0,test_model.dim_state,instance_parameter.use_adaptive_solver);
    %Find error:
    [x_E_ss,S_E_ss] = find_error(xS_end,xS_CDCKF);
    x_E_CDCKF(j) = x_E_ss;
    S_E_CDCKF(j) = S_E_ss;
    [x_E_ss,S_E_ss] = find_error(xS_end,xS_CDCKF_ts);
    x_E_CDCKF_ts(j) = x_E_ss;
    S_E_CDCKF_ts(j) = S_E_ss;
    [x_E_ss,S_E_ss] = find_error(xS_end,xS_LSKF);
    x_E_LSKF(j) = x_E_ss;
    S_E_LSKF(j) = S_E_ss;
    [x_E_ss,S_E_ss] = find_error(xS_end,xS_LSKF_RK4);
    x_E_LSKF_RK4(j) = x_E_ss;
    S_E_LSKF_RK4(j) = S_E_ss;
end
clf
t = tiledlayout(1,2,'TileSpacing','compact');
nexttile(1)
semilogy(log2(m_list),x_E_CDCKF,'-o')
hold on
semilogy(log2(m_list),x_E_CDCKF_ts,'--x')
semilogy(log2(m_list),x_E_LSKF,'-.+')
semilogy(log2(m_list(3:end)),x_E_LSKF_RK4(3:end),':*')
xlabel('log2(m)')
ylabel('error (L^2 norm)')
legend('CDCKF','IT 1.5','LSKF-RK2','LSKF-RK4')
title('error in mean')
nexttile(2)
semilogy(log2(m_list),S_E_CDCKF,'-o')
hold on
semilogy(log2(m_list),S_E_CDCKF_ts,'--x')
semilogy(log2(m_list),S_E_LSKF,'-.+')
semilogy(log2(m_list(3:end)),S_E_LSKF_RK4(3:end),':*')
xlabel('log2(m)')
ylabel('error (Frobenius norm)')
title('error in covariance matrix')
savefig('linear_single_step_conv.fig')
exportgraphics(t,'linear_single_step_conv.png','Resolution',450)
function [x_E,S_E] = find_error(x_S_true,x_S_compute)
    x_E_v = x_S_true.mean - x_S_compute.mean;
    x_E = sqrt(sum(x_E_v.^2,'all'));
    S_E_m = x_S_true.covariance - x_S_compute.covariance;
    S_E = sqrt(sum(S_E_m.^2,'all'));
end