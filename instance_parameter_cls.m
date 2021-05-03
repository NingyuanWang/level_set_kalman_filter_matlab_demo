classdef instance_parameter_cls
    %Collects parameters associated with a particular instance of
    %simulation. 
    properties
        T_simulation
        observation_count
        euler_subdivision
        ode_subdivision
        alpha
        use_adaptive_solver
        return_uncorrected_prediction
    end
    methods
        function obj = instance_parameter_cls(TT_simulation, Observation_count,Euler_subdivision,Ode4_subdivision,varargin)
            obj.T_simulation = TT_simulation;
            obj.observation_count = Observation_count;
            obj.euler_subdivision = Euler_subdivision;
            obj.ode_subdivision = Ode4_subdivision;
            if nargin >=5
                obj.alpha = varargin{1};
            else
                obj.alpha = 1e-2;
            end
            if nargin >= 6
                obj.use_adaptive_solver = varargin{2};
            else
                obj.use_adaptive_solver = true;
            end
            if nargin >= 7
                obj.return_uncorrected_prediction = varargin{3};
            else
                obj.return_uncorrected_prediction = false;
            end
        end
    end
end