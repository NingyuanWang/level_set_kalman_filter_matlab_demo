function [test_model] = get_test_model(sigma_1,sigma_2)
%sigma_1: process noise in location and velocity
%sigma_2: process noise in acceleration
%Constants: 
sigma_x = 0.04;
sigma_v = 0.15;
%w0 = 3;%Testcase: 3,4.5,6,

c23 = sqrt(sigma_2.^2 - sigma_1.^2);
K_process = diag([sigma_1,sigma_1,sigma_2].^2);
K_observe = diag([sigma_x,sigma_v].^2);

test_model = continuous_discrete_model_cls(3,@oscillator_dynamic,K_process,2,@oscillator_observe,K_observe);
end
function [y] = oscillator_observe(x)
y = x(1:2);
end

function [advvel] = oscillator_dynamic(t,x)
%Coordinated turns radar tracking problem as given in CD-CKF paper
%  state variable x = (x,vx,y,vy,z,vz,w): position, velocity, "turn rate"
advvel = 0*x;
advvel(1) = x(2);
advvel(2) = x(3);
advvel(3) = -x(1);
end

