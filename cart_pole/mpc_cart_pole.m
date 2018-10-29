clear;
clc;
close all;

 % Time initialization
 horizon = 300;
 dt = 0.01;
 T = 16; % shorter DDP time horizon
 
 % Trajectory initialization
 xo = zeros(4,1); % state: [x; x_dot; theta; theta_dot]
 x_traj = zeros(4,horizon); % trajectory
 
 % Final configuration
 p_target(1,1) = 10; 
 p_target(2,1) = 0;
 p_target(3,1) = pi;
 p_target(4,1) = 0;
 
% Weight in Final State:
Q_f(1,1) = 1;
Q_f(2,2) = 150;
Q_f(3,3) = 400;
Q_f(4,4) = 50;

% Weight in the Control:
R = 2 * eye(1,1);

% Learning Rate:
gamma = 0.5;

% Number of iterations
num_iter = 100;

% Changing the nominal model from the real one
uncert = 1.1;
 
 for i = 1:horizon-1
    
     [u] = fnDDP_cart_pole(xo,p_target,Q_f,R,T,dt,gamma,num_iter,1);
     [x_traj_next] = fnsimulate_2(xo,u(1),2,dt,0,uncert);
     x_traj(:,i+1) = x_traj_next(:,end);
     xo = x_traj_next(:,end);
     
     fprintf('iLQG Horizon %d,  Current State %d %d,   Control applied %d \n',i, xo(3), xo(4), u(1));
 end 



%% PLOTTING

% Setting up a time vector to span the horizon to plot against
time(1)=0;
for i= 2:horizon
time(i) =time(i-1) + dt;  
end

figure(2)
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 12)

subplot(3,2,1)
hold on
plot(time,x_traj(1,:),'b', 'linewidth', 1.5) 
plot(time,p_target(1,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('X')
xlabel('Time in sec')
hold off
grid

subplot(3,2,2)
hold on
plot(time,x_traj(2,:),'b', 'linewidth', 1.5);
plot(time,p_target(2,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('X dot')
xlabel('Time in sec')
hold off
grid

subplot(3,2,3)
hold on
plot(time,x_traj(3,:),'b', 'linewidth', 1.5);
plot(time,p_target(3,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('Theta')
xlabel('Time in sec')
hold off
grid

subplot(3,2,4)
hold on
plot(time,x_traj(4,:), 'b', 'linewidth', 1.5)
plot(time,p_target(4,1)*ones(1,horizon),'r', 'linewidth', 1.5)
title('Theta dot')
xlabel('Time in sec')
hold off
grid

