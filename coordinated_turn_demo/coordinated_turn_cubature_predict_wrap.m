function [xS_predict] = coordinated_turn_cubature_predict_wrap(continuous_discrete_model,timesteps,xS0)
%Gives the cubature updates of the state variable for the coordinated turn 
%problem
% [varargout] = [xS_predict,t_int,y_int]
%xS0 and xS_predict is of mean_covariance_sqrt_cls class
%The function gives optional output of intermediate timesteps and
%intermediate state at these steps t_int, y_int. 

%Step 0a: Extract information from model to eliminate global variables:
%WARNING: The implementation is specific to the problem.
nx = 7;%dimension of state space
nPts = 2 * nx;
sig1 = sqrt(continuous_discrete_model.K_process(2,2));%sigma1 in process noise
sig2 = sqrt(continuous_discrete_model.K_process(7,7));%sigma2 in process noise
T = timesteps(end) - timesteps(1);%Interval between measurements
[L,D] = ldl(continuous_discrete_model.K_process);
Qsqrt = L * sqrt(D);%Square root pf process noise
QPtArray = sqrt(nx) * cat(2,eye(nx),-eye(nx));%Quadrature point definition

%Step 0b: reformulate input parameters:
xkk = xS0.mean;
Skk = xS0.c_sqrt;
%Copying predict line by line:
%Xi = repmat(xkk,1,nPts) + Skk*QPtArray;
Xi = repmat(xkk,1,nPts) + Skk*QPtArray;
%Xi = StateEq_filt(Xi,M);
%WARNING: the implementation is problem specific
for j = 1:numel(timesteps)-1
    deltat = timesteps(j+1)-timesteps(j);
    Xi = StateEq_filt_singlestep(Xi,deltat);
end
%xkk1 = sum(Xi,2)/nPts; 
xkk1 = sum(Xi,2)/nPts; 
%X = (Xi-repmat(xkk1,1,nPts))/sqrt(nPts);
X = (Xi-repmat(xkk1,1,nPts))/sqrt(nPts);
%WARNING: Entries related to Li and L are problem-specific
%L2 = sig1*[1 0 0 xkk(7) 0 0 0]';
L2 = sig1*[1 0 0 xkk(7) 0 0 0]';
%L4 = sig1*[0 -xkk(7) 1 0 0 0 0]';
L4 = sig1*[0 -xkk(7) 1 0 0 0 0]';
%L6 = [0 0 0  0 sig1 0 0]';
L6 = [0 0 0  0 sig1 0 0]';
%L7 = sig2*[0 -xkk(4) 0 xkk(2) 0 0 0]';
L7 = sig2*[0 -xkk(4) 0 xkk(2) 0 0 0]';
%L = [zeros(7,1) L2 zeros(7,1) L4 zeros(7,1) L6 L7];
L = [zeros(7,1) L2 zeros(7,1) L4 zeros(7,1) L6 L7];
%[foo,Skk1] = qr([ X sqrt(T)*(Qsqrt+0.5*T*L)  sqrt(T/3)*0.5*T*L]',0);
[~,Skk1] = qr([ X sqrt(T)*(Qsqrt+0.5*T*L)  sqrt(T/3)*0.5*T*L]',0);
%Skk1 = Skk1'; 
Skk1 = Skk1'; 
%Reformulating output:
xS_predict = mean_covariance_sqrt_cls(xkk1,Skk1);
end
function Xi_next = StateEq_filt_singlestep(Xi,deltat)
Xi_next = Xi;
sigma_point_count = size(Xi,2);
    for i = 1:sigma_point_count
    x = Xi(:,i);
        f = [x(2)  -x(7)*x(4)  x(4) x(7)*x(2) x(6) 0 0]';
        L0 = x(7)*[-x(4) -x(7)*x(2) x(2) -x(7)*x(4) 0 0 0]';
        x = x+ deltat*f + 0.5*deltat^2*L0;
        Xi_next(:,i) = x;
    end
end
