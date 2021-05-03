classdef continuous_discrete_model_cls
    % Description of a continuous-discrete model with its state space,
    %dynamics described by a convection-diffusion equation (dxdt and
    %continuous process noise), measurement model, and measurement noise
    properties
        dim_state %dimension of state space
        dxdt      %ODE defining velocity field
        K_process %continuous process noise
        dim_ob    %dimension of observation space
        ob_fun    %observation model function
        K_observe %observation noise
    end
    methods
        function obj = continuous_discrete_model_cls(Dim_state,Dxdt,KK_process,Dim_ob,Ob_fun,KK_observe)
            obj.dim_state = Dim_state;
            obj.dxdt = Dxdt;
            obj.K_process = KK_process;
            obj.dim_ob = Dim_ob;
            obj.ob_fun = Ob_fun;
            obj.K_observe = KK_observe;
        end
    end
end