function [instance_parameter] = get_instance_parameter(m,w0_degree,measurement_interval,varargin)
%Generates the instance parameter used for the coordinated turns test case
%   m: subdivision between measurements
%   w0: rotation rate (affects Total simulation time T)
%   measurement_interval: <-
%   varargin = alpha,use_adaptive_solver,return_uncorrected_prediction. Same
%   meaning as in the instance_parameter_cls
rotation_angle = 1800;
if (nargin >= 7)
    T = varargin{4};
else
    T = max(rotation_angle / w0_degree,120);
end
measurement_count = ceil(T/measurement_interval);
if (nargin >= 4)
    alpha = varargin{1};
else
    alpha = 1.0;
end
if (nargin >= 5)
    use_adaptive_solver = varargin{2};
else
    use_adaptive_solver = false;
end
if (nargin >= 6)
    return_uncorrected_prediction = varargin{3};
else
    return_uncorrected_prediction = false;
end

instance_parameter = instance_parameter_cls(T,measurement_count,4096,m,alpha,use_adaptive_solver,return_uncorrected_prediction);

end

