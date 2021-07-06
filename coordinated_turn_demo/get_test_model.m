function [test_model] = get_test_model()
%Constants: 
sigma_1 = sqrt(0.2);
sigma_2 = 7e-4;
sigma_r = 50;
sigma_theta = 0.1/180*pi;
sigma_phi = 0.1/180*pi;
%w0 = 3;%Testcase: 3,4.5,6,


K_process = diag([0,sigma_1,0,sigma_1,0,sigma_1,sigma_2].^2);
K_observe = diag([sigma_r,sigma_theta,sigma_phi].^2);

test_model = continuous_discrete_model_cls(7,@coordinated_turns_dynamic,K_process,3,@coordinated_turns_observe,K_observe);
end
function [y] = coordinated_turns_observe(x)
%Coordinated turns radar tracking problem as given in CD-CKF paper
%  state variable x = (x,vx,y,vy,z,vz,w): position, velocity, "turn rate" 
%  observation variable y = (r,theta,phi)
s = [1500;10;0];
%Location of radar station. Chosen to be the same as in the example code
xs = [x(1);x(3);x(5)] - s;
y = zeros(3,1);
y(1) = sqrt(xs(1)^2+xs(2)^2+xs(3)^2);
y(2) = atan2(xs(2),xs(1));
y(3) = atan2(xs(3),sqrt(xs(1)^2+xs(2)^2));
end

function [advvel] = coordinated_turns_dynamic(t,x)
%Coordinated turns radar tracking problem as given in CD-CKF paper
%  state variable x = (x,vx,y,vy,z,vz,w): position, velocity, "turn rate"
advvel = 0*x;
advvel(1) = x(2);
advvel(3) = x(4);
advvel(5) = x(6);
advvel(2) = -x(4)*x(7);
advvel(4) = x(2)*x(7);
end

